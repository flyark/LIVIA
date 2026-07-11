/*
 * cLIP data pipeline — convert-lis +
 * gene_search (search / orient / filter_and_deduplicate) + build_contact_map.
 * Turns raw lis.py output (or pre-converted) rows into per-partner contact
 * fingerprints for clustering.
 */
(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) module.exports = factory();
  else root.CLIPPipeline = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // ---- minimal RFC4180 CSV parser (handles quoted fields w/ commas & newlines) ----
  function parseCSV(text) {
    const rows = [];
    let field = '', row = [], inQ = false;
    for (let i = 0; i < text.length; i++) {
      const c = text[i];
      if (inQ) {
        if (c === '"') { if (text[i + 1] === '"') { field += '"'; i++; } else inQ = false; }
        else field += c;
      } else if (c === '"') inQ = true;
      else if (c === ',') { row.push(field); field = ''; }
      else if (c === '\n') { row.push(field); rows.push(row); field = ''; row = []; }
      else if (c === '\r') { /* skip */ }
      else field += c;
    }
    if (field.length || row.length) { row.push(field); rows.push(row); }
    if (!rows.length) return { header: [], rows: [] };
    const header = rows[0];
    const out = [];
    for (let r = 1; r < rows.length; r++) {
      if (rows[r].length === 1 && rows[r][0] === '') continue; // blank line
      const o = {};
      for (let c = 0; c < header.length; c++) o[header[c]] = rows[r][c];
      out.push(o);
    }
    return { header, rows: out };
  }

  // ---- convert-lis: raw lis.py → cLIP column names ----
  const COLUMN_MAP = {
    rank: 'Rank', len_i: 'Protein_Len_A', len_j: 'Protein_Len_B',
    cLIR_indices_i: 'cLIR_indice_A', cLIR_indices_j: 'cLIR_indice_B',
    LIR_indices_i: 'LIR_indice_A', LIR_indices_j: 'LIR_indice_B',
    LIR_i: 'LIR_A', LIR_j: 'LIR_B', cLIR_i: 'cLIR_A', cLIR_j: 'cLIR_B',
    pLDDT_i: 'pLDDT_A', pLDDT_j: 'pLDDT_B',
    LIpLDDT_i: 'LIpLDDT_A', LIpLDDT_j: 'LIpLDDT_B',
    cLIpLDDT_i: 'cLIpLDDT_A', cLIpLDDT_j: 'cLIpLDDT_B',
  };
  function parseProteinNames(name, sep) {
    const seps = [sep, '___', ' & ', '_vs_', ' vs ', '_VS_', '--', '__'];
    for (const s of seps) if (s && name.includes(s)) { const p = name.split(s); return [p[0].trim(), p.slice(1).join(s).trim()]; }
    return [name, name];
  }
  const DROP_COLS = ['structure_file', 'pae_plot', 'deeploc_1', 'deeploc_2', 'iLIS_x_iLIA'];   // columns cLIP never renders — dropped on parse to cut memory on large multi-file loads
  function convertRows(rows, sep) {
    sep = sep || '_vs_';
    if (!rows.length) return rows;
    const isRaw = ('name' in rows[0]) && ('cLIR_indices_i' in rows[0] || 'cLIR_indices_j' in rows[0]);
    if (!isRaw) { for (const r of rows) for (const d of DROP_COLS) delete r[d]; return rows; } // already converted — just trim
    return rows.map((r) => {
      const o = {};
      for (const k in r){ const nk = COLUMN_MAP[k] || k; if (DROP_COLS.indexOf(nk) >= 0) continue; const v = r[k]; o[nk] = (v !== '' && v != null && !isNaN(v)) ? +v : v; }   // numeric strings -> Numbers (~8 B vs ~40 B) — cLIR_indice '[...]' / names stay strings (isNaN)
      const [s1, s2] = parseProteinNames(String(r.name), sep);
      o.Symbol_1 = s1; o.Symbol_2 = s2;
      return o;
    });
  }

  // ---- gene search / orient ----
  function searchGene(rows, gene) {
    const g = String(gene);
    return rows.filter((r) => String(r.Symbol_1) === g || String(r.Symbol_2) === g);
  }
  const SWAP = [['Symbol_1', 'Symbol_2'], ['cLIR_indice_A', 'cLIR_indice_B'],
    ['LIR_indice_A', 'LIR_indice_B'], ['Protein_Len_A', 'Protein_Len_B']];
  function orient(rows, gene) {
    const g = String(gene);
    return rows.map((r) => {
      if (String(r.Symbol_2) !== g) return r;
      const o = Object.assign({}, r);
      for (const [a, b] of SWAP) if (a in o || b in o) { const t = o[a]; o[a] = o[b]; o[b] = t; }
      return o;
    });
  }

  const num = (x) => { const v = parseFloat(x); return isNaN(v) ? null : v; };
  const DEFAULT_EXCL = ['phospho', 'mutant'];

  // ---- filter_and_deduplicate ----
  function filterAndDedup(rows, opts) {
    opts = opts || {};
    const ilisCut = opts.ilisCutoff == null ? 0.223 : opts.ilisCutoff;
    const excl = opts.exclusionKeywords || DEFAULT_EXCL;
    const sortBy = opts.sortBy || 'iLIS';
    const topN = opts.topN == null ? 50 : opts.topN;

    let r = rows.map((x) => Object.assign({}, x));
    for (const x of r) { for (const c of ['iLIS', 'iLIA', 'iLISA', 'ipTM', 'Rank']) if (c in x) x['_' + c] = num(x[c]); }
    r.forEach((x) => { x._iLISxiLIA = x._iLISA != null ? x._iLISA : (x._iLIS || 0); });

    // iLIS cutoff
    r = r.filter((x) => x._iLIS != null && x._iLIS >= ilisCut);
    // exclusion keywords on Symbol_2 (case-insensitive)
    const el = excl.map((k) => k.toLowerCase());
    r = r.filter((x) => { const s = String(x.Symbol_2).toLowerCase(); return !el.some((k) => s.includes(k)); });
    if (!r.length) return r;

    // best Rank-1 per (Symbol_1,Symbol_2): keep pairs that have a rank-1 row
    const key = (x) => x.Symbol_1 + '' + x.Symbol_2;
    const hasRank1 = new Set();
    for (const x of r) if (x._Rank === 1) hasRank1.add(key(x));
    if (hasRank1.size) r = r.filter((x) => hasRank1.has(key(x)));
    if (!r.length) return r;

    // stable sort by iLIS desc, then dedup by (Symbol_1,Symbol_2,Rank) keep first
    const sv = (x) => (x['_' + sortBy] != null ? x['_' + sortBy] : (x._iLIS || 0));
    r = r.map((x, i) => [x, i]).sort((A, B) => (sv(B[0]) - sv(A[0])) || (A[1] - B[1])).map((p) => p[0]);
    const seen = new Set();
    const dedup = [];
    for (const x of r) { const k = key(x) + '' + x._Rank; if (!seen.has(k)) { seen.add(k); dedup.push(x); } }
    r = dedup;

    if (topN != null) r = r.slice(0, topN);
    return r;
  }

  // ---- parse residue-index string → int list (matches parse_indice_string) ----
  function parseIndice(s) {
    if (s == null) return [];
    s = String(s).trim();
    if (!s) return [];
    const cleaned = s.replace(/^\[|\]$/g, '');            // strip brackets
    const out = [];
    for (const tok of cleaned.split(/[,\s]+/)) {           // comma- or space-separated
      if (!tok) continue;
      const rng = tok.match(/^(\d+)\s*-\s*(\d+)$/);        // range "253-258" -> 253..258 (lis.py compresses runs)
      if (rng) { const a = +rng[1], b = +rng[2]; for (let v = a; v <= b; v++) out.push(v); }
      else { const v = parseFloat(tok); if (!isNaN(v)) out.push(Math.trunc(v)); }
    }
    return out;
  }

  // ---- build fingerprints (drop rows with no contacts) ----
  function buildFingerprints(rows, indiceCol, proteinLen) {
    indiceCol = indiceCol || 'cLIR_indice_A';
    const parsed = rows.map((r) => parseIndice(r[indiceCol]));
    let maxPos = proteinLen || 0;
    if (!maxPos) for (const idx of parsed) for (const v of idx) if (v > maxPos) maxPos = v;
    const fingerprints = [], keptRows = [];
    for (let i = 0; i < rows.length; i++) {
      const valid = parsed[i].filter((v) => v >= 1 && v <= maxPos);
      if (valid.length) { fingerprints.push(valid); keptRows.push(rows[i]); }
    }
    return { fingerprints, keptRows, proteinLen: maxPos };
  }

  // ---- residue list → compact "[a,b-c]" string (inverse of parseIndice) ----
  function serialiseIndice(nums) {
    const s = [...new Set(nums.filter((n) => n >= 1))].sort((a, b) => a - b);
    if (!s.length) return '[]';
    const out = []; let a = s[0], p = s[0];
    for (let i = 1; i < s.length; i++) { if (s[i] === p + 1) p = s[i]; else { out.push(a === p ? '' + a : a + '-' + p); a = p = s[i]; } }
    out.push(a === p ? '' + a : a + '-' + p);
    return '[' + out.join(',') + ']';
  }

  // ---- normalize mixed full/partial constructs ----
  // A screen may fold the same protein at different lengths across predictions (e.g. a receptor
  // as a full-length chain for some partners, an ectodomain/internal fragment for others). lis.py
  // residue indices are relative to whatever construct was folded, so those rows sit on different
  // coordinate frames and clustering/heatmap/sequence-viewer mis-register. Using the paired FASTA
  // (bySymLen: every construct of a symbol, by length) + alignMap(pred,target)->{map}, remap each
  // row's residue indices onto that symbol's LONGEST local construct so every row shares one frame.
  // Rows whose construct sequence isn't in the FASTA are left as-is (and flagged mixed). Mutates rows.
  function normalizeConstructs(rows, gene, parse, alignMapFn) {
    if (!parse || !parse.bySymLen || !alignMapFn) return { mixed: false, remapped: 0 };
    const refOf = (sym) => (parse.map && parse.map[sym]) || null;
    const cache = new Map();                              // (sym|len) → residue map (pred→ref) or null
    const mapFor = (sym, len) => {
      const k = sym + '|' + len; if (cache.has(k)) return cache.get(k);
      const ref = refOf(sym), byL = parse.bySymLen[sym]; let m = null;
      if (ref && byL && byL[len] && byL[len] !== ref) { const a = alignMapFn(byL[len], ref); m = (a && a.map) || null; }
      cache.set(k, m); return m;
    };
    let mixed = false, remapped = 0;
    for (const sym in parse.bySymLen) if (Object.keys(parse.bySymLen[sym]).length > 1) { mixed = true; break; }
    const side = (row, idxCols, lenCol, sym) => {
      const len = parseInt(row[lenCol], 10), ref = refOf(sym);
      if (!len || !ref || len === ref.length) return;    // no ref, or already full-length
      const m = mapFor(sym, len); if (!m) return;        // construct sequence not in FASTA → leave as-is
      for (const col of idxCols) { const out = []; for (const v of parseIndice(row[col])) { const t = (v >= 1 && v <= m.length) ? m[v - 1] : null; if (t) out.push(t); } row[col] = serialiseIndice(out); }
      row[lenCol] = ref.length; remapped++;
    };
    for (const row of rows) {
      side(row, ['cLIR_indice_A', 'LIR_indice_A'], 'Protein_Len_A', row.Symbol_1);
      side(row, ['cLIR_indice_B', 'LIR_indice_B'], 'Protein_Len_B', row.Symbol_2);
    }
    return { mixed, remapped };
  }

  // Best prediction per partner (max iLIS, any rank) — no cutoff/exclusion/topN.
  // Used as the scatter background so ALL predictions show, incl. below-cutoff ones.
  function bestPerPartner(rows) {
    const best = {};
    for (const r of rows) {
      const il = num(r.iLIS);
      if (il == null) continue;
      const s2 = String(r.Symbol_2);
      if (!best[s2] || il > best[s2]._il) best[s2] = Object.assign({}, r, { Symbol_2: s2, _il: il });  // keep full row → all metrics available on the scatter
    }
    return Object.values(best).map((o) => { const c = Object.assign({}, o); delete c._il; return c; });
  }

  // ---- pipeline on already-parsed+converted rows (fast path — parse once, reuse) ----
  function runPipelineRows(rows, gene, opts) {
    opts = opts || {};
    const hits = orient(searchGene(rows, gene), gene);
    const allPoints = bestPerPartner(hits);               // background: every partner
    const filt = filterAndDedup(hits, opts);              // (clones rows — safe to mutate below)
    const norm = normalizeConstructs(filt, gene, opts.fastaParse, opts.alignMap);   // unify mixed full/partial constructs onto one frame
    let plen = 0;
    for (const x of filt) { const v = num(x.Protein_Len_A); if (v && v > plen) plen = v; }
    const fp = buildFingerprints(filt, 'cLIR_indice_A', plen);
    return { rows: fp.keptRows, fingerprints: fp.fingerprints, proteinLen: fp.proteinLen, allPoints, allPredictions: hits, constructNorm: norm };
  }

  // ---- full pipeline: raw CSV text(s) + gene → {rows, fingerprints, proteinLen, allPoints} ----
  function runPipeline(csvTexts, gene, opts) {
    opts = opts || {};
    let all = [];
    for (const t of [].concat(csvTexts)) all = all.concat(convertRows(parseCSV(t).rows, opts.separator));
    return runPipelineRows(all, gene, opts);
  }

  return { parseCSV, convertRows, searchGene, orient, filterAndDedup, parseIndice, serialiseIndice, normalizeConstructs, buildFingerprints, bestPerPartner, runPipeline, runPipelineRows };
});
