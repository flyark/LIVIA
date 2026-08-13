// LIVIA CORS proxy — Cloudflare Worker
//
// Usage from client:
//   fetch(`https://livia-proxy.<account>.workers.dev/?url=${encodeURIComponent(targetUrl)}`)
//
// Security model:
//   - Only allowed origins receive an Access-Control-Allow-Origin echo
//   - Only allow-listed target hosts can be proxied (no open relay)
//   - GET only

const ALLOWED_ORIGINS = [
  'https://flyark.github.io',
  'http://localhost:8765',
  'http://localhost:8000',
  'http://127.0.0.1:8765',
  'http://127.0.0.1:8000',
];

const ALLOWED_TARGET_HOSTS = [
  'www.flyrnai.org',
  'flyrnai.org',
  // AFDB already sends CORS headers, so direct fetch works without the proxy.
  // Listed here as a resilience fallback only.
  'alphafold.ebi.ac.uk',
  // OSF-hosted prediction bundles (FlyPredictome cluster structures) — OSF sends no CORS header, so LIVIA
  // (universal.html?data=<osf_url>) must fetch them through this proxy. fetch() follows the OSF 302 to the file.
  'osf.io',
  'files.osf.io',
  'files.us.osf.io',
];

function pickAllowOrigin(reqOrigin) {
  return ALLOWED_ORIGINS.includes(reqOrigin) ? reqOrigin : ALLOWED_ORIGINS[0];
}

// Row cap — a safety valve so a pathological page can't run the parser unbounded.
const FP_SUMMARY_MAX_ROWS = 8000;

// The #famdb_summary_results tbody has a fixed 56-column layout (decoded from the rendered rows). We map
// by this list rather than the page's grouped <thead> (which has a 2-row header that doesn't align 1:1)
// or its per-cell debug <script> markers (present on some genes, absent on others — unreliable). If
// FlyPredictome changes the summary columns, update this list; `theadCols` is returned for cross-checking.
const FP_SUMMARY_COLS = [
  'partner_display', 'CompositeScore',
  'iLIS_best', 'iLIS_avg', 'iLISA_best', 'iLISA_avg', 'ipSAE_best', 'ipSAE_avg', 'actifpTM_best', 'actifpTM_avg', 'ipTM_best', 'ipTM_avg',
  'Length_ij', 'fPPI', 'fGI', 'iPPI', 'iGI',
  'id', 'ranks_counts', 'taxon1', 'taxon2', 'fbgn1', 'fbgn2', 'ncbi_gene_id1', 'ncbi_gene_id2', 'gene_symbol1', 'gene_symbol2', 'Protein_Len_A', 'Protein_Len_B',
  'LIS_avg', 'LIS_best', 'LIA_avg', 'LIA_best', 'LIR_avg', 'LIR_best', 'LIpLDDT_avg', 'LIpLDDT_best',
  'cLIS_avg', 'cLIS_best', 'cLIA_avg', 'cLIA_best', 'cLIR_avg', 'cLIR_best', 'cLIpLDDT_avg', 'cLIpLDDT_best',
  'pLDDT_avg', 'pLDDT_best', 'pTM_avg', 'pTM_best', 'iLIA_best', 'iLIA_avg', 'Confidence_avg', 'Confidence_best',
  'directory_name', 'Protein_1_size', 'Protein_2_size',
];

// Fast partner-list endpoint (?fpSummary=<FBgn>).
//
// FlyPredictome's famdb_summary page can be ~80 MB (the partner table is server-rendered into the HTML,
// bloated by a debug <script> per cell; the real data is only ~350 KB). Pulling that into the browser
// per gene-search is a non-starter. Here we fetch it SERVER-SIDE and run it through HTMLRewriter, which
// parses the response as a stream — we never hold the 80 MB in memory, only the extracted cells. The
// result (compact JSON: {cols, rows}) is cached on Cloudflare, so the 80 MB is paid at most once per gene.
//
// Caveat: HTMLRewriter still has to CHUNK through 80 MB, which costs CPU. This needs Workers Paid with a
// raised CPU limit (wrangler.toml `limits.cpu_ms`); on the free 10 ms tier it will time out. On failure
// we return an error status and LIVIA falls back to the "search unavailable" note.
async function handleFpSummary(fbgn, allowOrigin, ctx) {
  const jsonHeaders = (extra) => ({ ...corsHeaders(allowOrigin), 'Content-Type': 'application/json; charset=utf-8', ...(extra || {}) });
  if (!/^FBgn\d+$/i.test(fbgn)) {
    return new Response(JSON.stringify({ error: 'fpSummary expects an FBgn id, e.g. FBgn0010379' }), { status: 400, headers: jsonHeaders() });
  }
  fbgn = fbgn.replace(/^fbgn/i, 'FBgn');   // canonical case — FlyBase ids are case-sensitive, flyrnai 500s on FBGN…

  // Cache lookup (Cloudflare edge cache, keyed by fbgn). Re-echo CORS for the requesting origin.
  const cache = caches.default;
  const cacheKey = new Request('https://livia-proxy.internal/fpSummary/' + fbgn);
  const hit = await cache.match(cacheKey);
  if (hit) {
    const h = new Headers(hit.headers);
    for (const [k, v] of Object.entries(corsHeaders(allowOrigin))) h.set(k, v);
    h.set('X-FP-Cache', 'HIT');
    return new Response(hit.body, { status: hit.status, headers: h });
  }

  const target = 'https://www.flyrnai.org/tools/fly_predictome/web/famdb_summary/' + fbgn + '/';
  let upstream;
  try {
    // Match the plain ?url= proxy exactly — flyrnai 500s on `Accept: text/html`, and cacheEverything lets
    // Cloudflare cache the big upstream fetch too.
    upstream = await fetch(target, {
      headers: { 'User-Agent': 'Mozilla/5.0 (compatible; LIVIA-proxy/1.0; +https://flyark.github.io/LIVIA)', 'Accept': '*/*' },
      cf: { cacheTtl: 3600, cacheEverything: true },
    });
  } catch (e) {
    return new Response(JSON.stringify({ error: 'upstream fetch failed', detail: String(e && e.message || e) }), { status: 502, headers: jsonHeaders() });
  }
  if (!upstream.ok) {
    return new Response(JSON.stringify({ error: 'upstream returned ' + upstream.status }), { status: 502, headers: jsonHeaders() });
  }

  // Streaming extraction of #famdb_summary_results: the thead names (for cross-checking) and the tbody
  // rows (cell text). State is shared across handlers; HTMLRewriter fires them in document order.
  const theadCols = [], rows = [];
  let curRow = null, cell = null, truncated = false;
  const rewriter = new HTMLRewriter()
    .on('#famdb_summary_results thead th', {
      element(el) { cell = ''; el.onEndTag(() => { theadCols.push((cell || '').replace(/\s+/g, ' ').trim()); cell = null; }); },
      text(t) { if (cell !== null) cell += t.text; },
    })
    .on('#famdb_summary_results tbody tr', {
      element(el) {
        if (rows.length >= FP_SUMMARY_MAX_ROWS) { truncated = true; curRow = null; return; }
        curRow = [];
        el.onEndTag(() => { if (curRow && curRow.length) rows.push(curRow); curRow = null; });
      },
    })
    .on('#famdb_summary_results tbody td', {
      element(el) { if (!curRow) return; cell = ''; el.onEndTag(() => { if (curRow) curRow.push((cell || '').replace(/\s+/g, ' ').trim()); cell = null; }); },
      text(t) { if (cell !== null) cell += t.text; },
    });

  try {
    const transformed = rewriter.transform(upstream);
    // Drive the stream to completion but DISCARD the (80 MB) transformed output — reading chunk-by-chunk
    // keeps peak memory to the extracted arrays, not the whole page.
    const reader = transformed.body.getReader();
    for (;;) { const { done } = await reader.read(); if (done) break; }
  } catch (e) {
    return new Response(JSON.stringify({ error: 'extraction failed (likely the Worker CPU/memory limit on a large page)', detail: String(e && e.message || e), rowsSoFar: rows.length }), { status: 507, headers: jsonHeaders() });
  }

  // Name the columns with the known fixed layout when a row has the expected width; otherwise fall back
  // to the raw thead so the caller can at least see what changed.
  const rowWidth = rows[0] ? rows[0].length : 0;
  const cols = (rowWidth === FP_SUMMARY_COLS.length) ? FP_SUMMARY_COLS : theadCols;
  const body = JSON.stringify({ fbgn, cols, colsMatched: rowWidth === FP_SUMMARY_COLS.length, rowWidth, theadCols, rows, count: rows.length, truncated });
  const resp = new Response(body, { status: 200, headers: jsonHeaders({ 'Cache-Control': 'public, max-age=86400', 'X-FP-Cache': 'MISS' }) });
  if (ctx && ctx.waitUntil) ctx.waitUntil(cache.put(cacheKey, resp.clone()));
  return resp;
}

function corsHeaders(allowOrigin) {
  return {
    'Access-Control-Allow-Origin': allowOrigin,
    'Access-Control-Allow-Methods': 'GET, OPTIONS',
    'Access-Control-Allow-Headers': 'Content-Type',
    'Access-Control-Max-Age': '86400',
    'Vary': 'Origin',
  };
}

export default {
  async fetch(request, env, ctx) {
    const reqOrigin = request.headers.get('Origin') || '';
    const allowOrigin = pickAllowOrigin(reqOrigin);

    if (request.method === 'OPTIONS') {
      return new Response(null, { status: 204, headers: corsHeaders(allowOrigin) });
    }

    if (request.method !== 'GET') {
      return new Response('Method not allowed', { status: 405, headers: corsHeaders(allowOrigin) });
    }

    const reqUrl = new URL(request.url);

    // Fast FlyPredictome partner-list endpoint (server-side extract + cache; see handleFpSummary).
    const fpSummaryFbgn = reqUrl.searchParams.get('fpSummary');
    if (fpSummaryFbgn) return handleFpSummary(fpSummaryFbgn, allowOrigin, ctx);

    const target = reqUrl.searchParams.get('url');
    if (!target) {
      return new Response('LIVIA CORS proxy. Usage: ?url=<encoded target URL>', {
        status: 400, headers: { ...corsHeaders(allowOrigin), 'Content-Type': 'text/plain' },
      });
    }

    let targetUrl;
    try { targetUrl = new URL(target); }
    catch { return new Response('Invalid target URL', { status: 400, headers: corsHeaders(allowOrigin) }); }

    if (targetUrl.protocol !== 'https:' && targetUrl.protocol !== 'http:') {
      return new Response('Only http/https targets allowed', { status: 400, headers: corsHeaders(allowOrigin) });
    }
    if (!ALLOWED_TARGET_HOSTS.includes(targetUrl.hostname)) {
      return new Response(`Target host not allowed: ${targetUrl.hostname}`, { status: 403, headers: corsHeaders(allowOrigin) });
    }

    try {
      const upstream = await fetch(targetUrl.toString(), {
        method: 'GET',
        headers: {
          'User-Agent': 'Mozilla/5.0 (compatible; LIVIA-proxy/1.0; +https://flyark.github.io/LIVIA)',
          'Accept': '*/*',
        },
        // Cloudflare cache for hot assets (PAE PNG, PDB CIF) — 1 hour
        cf: { cacheTtl: 3600, cacheEverything: true },
      });

      const headers = new Headers(corsHeaders(allowOrigin));
      const ct = upstream.headers.get('Content-Type');
      if (ct) headers.set('Content-Type', ct);
      const cl = upstream.headers.get('Content-Length');
      if (cl) headers.set('Content-Length', cl);
      headers.set('Cache-Control', 'public, max-age=3600');

      return new Response(upstream.body, {
        status: upstream.status,
        statusText: upstream.statusText,
        headers,
      });
    } catch (e) {
      return new Response(`Upstream fetch failed: ${e && e.message ? e.message : 'unknown'}`, {
        status: 502, headers: corsHeaders(allowOrigin),
      });
    }
  },
};
