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
  async fetch(request) {
    const reqOrigin = request.headers.get('Origin') || '';
    const allowOrigin = pickAllowOrigin(reqOrigin);

    if (request.method === 'OPTIONS') {
      return new Response(null, { status: 204, headers: corsHeaders(allowOrigin) });
    }

    if (request.method !== 'GET') {
      return new Response('Method not allowed', { status: 405, headers: corsHeaders(allowOrigin) });
    }

    const reqUrl = new URL(request.url);
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
