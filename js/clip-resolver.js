/*
 * cLIP structure resolver — bait identity → UniProt accession → AlphaFold DB CIF.
 * Two paths, tried in order:
 *   1. exact sequence (FASTA)  → SwissProt CRC64 → UniParc checksum search → accession
 *   2. gene symbol + organism + length → UniProtKB search → accession (length-disambiguated)
 * Then accession → AFDB api/prediction → current cifUrl. All client-side (CORS-open).
 */
(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) module.exports = factory();
  else root.CLIPResolver = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';
  const _fetch = (typeof fetch !== 'undefined') ? fetch.bind(typeof self !== 'undefined' ? self : globalThis) : null;

  // ---- SwissProt CRC64 (matches UniProt's sequence.crc64 exactly) ----
  const POLY = 0xd800000000000000n;
  const TBL = Array.from({ length: 256 }, (_, i) => { let c = BigInt(i); for (let j = 0; j < 8; j++) c = (c & 1n) ? (c >> 1n) ^ POLY : (c >> 1n); return c; });
  function crc64(s) { let c = 0n; for (let i = 0; i < s.length; i++) c = TBL[Number((c ^ BigInt(s.charCodeAt(i))) & 0xffn)] ^ (c >> 8n); return c.toString(16).toUpperCase().padStart(16, '0'); }

  const SPECIES = { fly: '7227', human: '9606', mouse: '10090', yeast: '559292', worm: '6239', zebrafish: '7955', arabidopsis: '3702', ecoli: '83333' };

  // ---- FASTA → {symbol: sequence}. Handles ColabFold pair style (>p1___p2 / seqA:seqB)
  //      and plain monomer FASTA (>symbol / seq). ----
  function parseFastaToSeqMap(text, sep) {
    sep = sep || '___';
    const recs = []; let h = null, s = [];
    for (const line of text.split(/\r?\n/)) {
      if (line[0] === '>') { if (h !== null) recs.push({ h, s: s.join('') }); h = line.slice(1).trim(); s = []; }
      else if (line.trim()) s.push(line.trim());
    }
    if (h !== null) recs.push({ h, s: s.join('') });
    const clean = (x) => (x || '').replace(/[^A-Za-z]/g, '').toUpperCase();
    const map = {}, warn = [], allSeqs = [];
    const put = (sym, seq) => {
      sym = (sym || '').trim(); seq = clean(seq);
      if (!sym || !seq) return;
      if (map[sym] && map[sym] !== seq) warn.push(sym);      // inconsistent sequence across pairs
      if (!map[sym]) map[sym] = seq;
    };
    for (const r of recs) {
      const parts = r.h.split(sep).map((x) => x.trim()).filter(Boolean);
      const chains = r.s.split(':').map(clean).filter(Boolean);
      for (const c of chains) allSeqs.push(c);               // index every chain (length fallback)
      // confident symbol→seq only when header parts line up 1:1 with chains
      // (monomer ">sym", or named pair ">sym1___sym2" with seqA:seqB). Arbitrary headers
      // like ">prediction1" don't map by name — they're resolved by length below.
      if (parts.length === chains.length && parts.length >= 1) parts.forEach((sym, i) => put(sym, chains[i]));
    }
    const uniq = [...new Set(allSeqs)], byLen = {};
    for (const q of uniq) (byLen[q.length] = byLen[q.length] || []).push(q);
    return { map, byLen, seqs: uniq, warn: [...new Set(warn)] };
  }

  // Pick the bait's sequence from a parsed FASTA: by symbol/pair-name, else by unique length
  // (Protein_Len from the CSV) — handles arbitrary headers like ">prediction1".
  function baitSequence(parse, gene, len) {
    if (!parse || !parse.map) return null;
    const m = parse.map, g = String(gene);
    if (m[g]) return m[g];
    const gl = g.toLowerCase();
    for (const key in m) {                                   // case-insensitive + first token
      const k = key.toLowerCase();                           // (">Sym", ">sym desc", "sp|ACC|SYM_ORG")
      if (k === gl || k.split(/[\s|]/)[0] === gl) return m[key];
    }
    const byL = parse.byLen && parse.byLen[len];             // fall back to unique length (arbitrary headers)
    if (byL && byL.length === 1) return byL[0];
    return null;                                             // ambiguous / not present
  }

  const _sleep = (ms) => new Promise((r) => setTimeout(r, ms));
  async function _json(url, tries) {                 // retry on 429/5xx/network (transient rate-limits)
    tries = tries || 3;
    for (let i = 0; i < tries; i++) {
      try {
        const r = await _fetch(url);
        if (r.ok) return r.json();
        if (r.status !== 429 && r.status < 500) throw new Error(url + ' → ' + r.status);
      } catch (e) { if (i === tries - 1) throw e; }
      await _sleep(400 * (i + 1));
    }
  }

  async function _bySequence(seq, orgId, expectLen) {
    const crc = crc64(seq);
    const d = await _json(`https://rest.uniprot.org/uniparc/search?query=checksum:${crc}&fields=upi&format=json`);
    if (!d.results || !d.results.length) return null;
    const e = await _json(`https://rest.uniprot.org/uniparc/${d.results[0].uniParcId}?format=json`);
    const xrefs = (e.uniParcCrossReferences || []).filter((x) => String(x.database).includes('UniProtKB'));
    if (!xrefs.length) return null;
    // prefer: active · canonical (no -isoform) · organism match · Swiss-Prot
    const score = (x) => (x.active ? 8 : 0) + (!/-\d+$/.test(x.id) ? 4 : 0) + (orgId && String(x.organism && x.organism.taxonId) === String(orgId) ? 2 : 0) + (/Swiss-Prot/i.test(x.database) ? 1 : 0);
    const sorted = xrefs.slice().sort((a, b) => score(b) - score(a));
    return { by: 'sequence (exact)', candidates: sorted.map((x) => ({ acc: x.id.split('-')[0], length: seq.length, organism: x.organism && x.organism.scientificName, gene: null })) };
  }

  // Tier 2 (fast, ~1s): FlyBase / construct translations are usually the UniProt canonical
  // +/- a few terminal residues, so the full-length checksum misses. Checksum trimmed
  // prefixes/suffixes and match them all in ONE batched UniParc OR-query.
  async function _byTrimmedSequence(seq, orgId) {
    const crcs = new Set();
    for (let k = 0; k <= 25; k++) { const c = seq.slice(0, seq.length - k); if (c.length >= 50) crcs.add(crc64(c)); }
    for (let k = 1; k <= 8; k++) { const n = seq.slice(k); if (n.length >= 50) crcs.add(crc64(n)); }
    const q = '(' + [...crcs].map((c) => 'checksum:' + c).join(' OR ') + ')';
    let d; try { d = await _json('https://rest.uniprot.org/uniparc/search?query=' + encodeURIComponent(q) + '&fields=upi&format=json'); } catch (e) { return null; }
    if (!d.results || !d.results.length) return null;
    const e = await _json('https://rest.uniprot.org/uniparc/' + d.results[0].uniParcId + '?format=json');
    const xrefs = (e.uniParcCrossReferences || []).filter((x) => String(x.database).includes('UniProtKB'));
    if (!xrefs.length) return null;
    const score = (x) => (x.active ? 8 : 0) + (!/-\d+$/.test(x.id) ? 4 : 0) + (orgId && String(x.organism && x.organism.taxonId) === String(orgId) ? 2 : 0) + (/Swiss-Prot/i.test(x.database) ? 1 : 0);
    const sorted = xrefs.slice().sort((a, b) => score(b) - score(a));
    return { by: 'sequence (trimmed ends)', candidates: sorted.map((x) => ({ acc: x.id.split('-')[0], length: seq.length, organism: x.organism && x.organism.scientificName, gene: null })) };
  }

  async function _bySymbol(gene, orgId, expectLen) {
    if (!gene || gene === 'undefined' || gene === 'null') return null;   // seq-only resolves pass no gene; text search on "undefined" returns spurious hits
    // Try, in order: UniProt entry name (e.g. NAME_SPECIES) → accession →
    // accession embedded in entry name → gene name. Many datasets label proteins by entry
    // name, which gene: search misses.
    const acc0 = String(gene).split('_')[0];
    const ACC_RE = /^([OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9]([A-Z][A-Z0-9]{2}[0-9]){1,2})$/;
    const org = orgId ? ` AND organism_id:${orgId}` : '';
    const queries = [`id:${gene}`];
    if (ACC_RE.test(gene)) queries.push(`accession:${gene}`);
    if (ACC_RE.test(acc0) && acc0 !== gene) queries.push(`accession:${acc0}`);
    queries.push(`gene:${gene}${org}`);
    queries.push(`xref:${gene}${org}`); queries.push(`${gene}${org}`);   // cross-reference IDs: FlyBase FBgn, WormBase WBGene, Ensembl, …
    for (const q of queries) {
      let d;
      try { d = await _json(`https://rest.uniprot.org/uniprotkb/search?query=${encodeURIComponent(q)}&fields=accession,length,reviewed,organism_name,gene_names&format=json&size=25`); } catch (e) { continue; }
      const cands = (d.results || []).map((x) => ({ acc: x.primaryAccession, length: x.sequence && x.sequence.length, reviewed: /Swiss-Prot/i.test(x.entryType || ''), organism: x.organism && x.organism.scientificName, gene: x.genes && x.genes[0] && x.genes[0].geneName && x.genes[0].geneName.value }));
      if (!cands.length) continue;
      const score = (c) => (expectLen && c.length === expectLen ? 8 : 0) + (c.reviewed ? 2 : 0) + (expectLen ? -Math.min(2, Math.abs((c.length || 0) - expectLen) / 50) : 0);
      const via = q.startsWith('id:') ? 'entry name' : q.startsWith('accession:') ? 'accession' : q.startsWith('xref:') ? 'cross-ref (FlyBase/etc.)' : q.startsWith('gene:') ? 'symbol + length' : 'text search';
      return { by: via, candidates: cands.slice().sort((a, b) => score(b) - score(a)) };
    }
    return null;
  }

  // Tier 2/3 fallback: EBI BLAST when the exact checksum + name lookups miss (isoform /
  // construct sequences that differ from the UniProt canonical). Submit -> poll -> result;
  // accept only a confident, well-covered top hit.
  async function _byBlast(seq, db) {
    const base = 'https://www.ebi.ac.uk/Tools/services/rest/ncbiblast';
    let jid;
    try {
      const run = await _fetch(base + '/run', { method: 'POST', body: new URLSearchParams({ email: 'livia@example.com', program: 'blastp', stype: 'protein', database: db, sequence: seq, scores: '5', alignments: '5' }) });
      if (!run.ok) return null;
      jid = (await run.text()).trim();
    } catch (e) { return null; }
    for (let i = 0; i < 45; i++) {                       // ~4s x 45 ceiling; typically 15-25s
      await new Promise((r) => setTimeout(r, 4000));
      let st; try { st = (await (await _fetch(base + '/status/' + jid)).text()).trim(); } catch (e) { continue; }
      if (st === 'FINISHED') break;
      if (st === 'ERROR' || st === 'FAILURE' || st === 'NOT_FOUND' || i === 44) return null;
    }
    let res; try { res = await (await _fetch(base + '/result/' + jid + '/json')).json(); } catch (e) { return null; }
    const h = (res.hits || [])[0], hsp = h && h.hit_hsps && h.hit_hsps[0];
    if (!h || !hsp) return null;
    const idPct = +hsp.hsp_identity, aln = +hsp.hsp_align_len;
    if (!(idPct >= 90) || !(aln >= 0.7 * seq.length)) return null;   // confident + well-covered only
    return { by: 'BLAST (' + (/swissprot/i.test(db) ? 'Swiss-Prot' : 'UniProtKB') + ', ' + idPct.toFixed(0) + '% id)', candidates: [{ acc: h.hit_acc, length: null, organism: h.hit_os, gene: null }] };
  }

  async function _afdbCif(acc) {
    try { const d = await _json(`https://alphafold.ebi.ac.uk/api/prediction/${acc}`); if (Array.isArray(d) && d.length) return d[0].cifUrl || null; } catch (e) {}
    return null;
  }

  // Partial-protein construct named BASE_start_end (+ optional _MUT…), e.g. GENE_1_300
  // or GENE_301_600_A45G (point mutation). Validated by (end-start+1) === fragment length. The AFDB
  // structure is the full BASE protein; cLIR coords get shifted by (start-1) to full numbering.
  function parseFragmentRange(gene, length) {
    const m = String(gene).match(/^(.+?)_(\d+)_(\d+)(?:_(.+))?$/);
    if (!m) return null;
    const start = +m[2], end = +m[3];
    if (!(start >= 1 && end > start)) return null;
    if (length != null && (end - start + 1) !== length) return null;   // range must equal the fragment length
    return { base: m[1], start, end, mutations: m[4] ? m[4].split('_').filter(Boolean) : [] };
  }

  // Resolve one bait → {accession, cifUrl, matchedBy, …, lengthOk, fragment} | null
  async function resolveStructure({ gene, length, seq, species, blast }) {
    if (!_fetch) throw new Error('fetch unavailable');
    const orgId = SPECIES[species] || species || '';
    const frag = parseFragmentRange(gene, length);      // partial protein? resolve the BASE, offset coords
    let found = null;
    // A fragment's sequence isn't a full UniParc entry → checksum won't match; resolve BASE by name.
    if (seq && seq.length && !frag) { try { found = await _bySequence(seq, orgId, length); } catch (e) {} }
    const _symName = frag ? frag.base : gene;
    if ((!found || !found.candidates.length) && _symName) { try { found = await _bySymbol(_symName, orgId, frag ? null : length); } catch (e) {} }
    if ((!found || !found.candidates.length) && seq && seq.length) { try { found = await _bySequence(seq, orgId, length); } catch (e) {} }
    if ((!found || !found.candidates.length) && seq && seq.length >= 60) { try { found = await _byTrimmedSequence(seq, orgId); } catch (e) {} }   // Tier 2: canonical +/- terminal residues
    if ((!found || !found.candidates.length) && blast && seq && seq.length >= 25) {   // Tier 2 (Swiss-Prot) -> Tier 3 (full UniProtKB) similarity BLAST
      try { found = await _byBlast(seq, 'uniprotkb_swissprot'); } catch (e) {}
      if (!found || !found.candidates.length) { try { found = await _byBlast(seq, 'uniprotkb'); } catch (e) {} }
    }
    if (!found || !found.candidates.length) return null;
    // Among tied candidates, prefer one that actually has an AFDB structure (e.g. two 112-aa
    // TrEMBL entries where only one is in AlphaFold DB).
    let chosen = null, cifUrl = null;
    for (const c of found.candidates.slice(0, 6)) { const cu = await _afdbCif(c.acc); if (cu) { chosen = c; cifUrl = cu; break; } }
    if (!chosen) chosen = found.candidates[0];
    const acc = chosen.acc;
    // canonical UniProt sequence == the AFDB structure's sequence/numbering; used to align
    // the predicted (isoform/construct) sequence and remap cLIR residue coordinates.
    let structSeq = null, structLen = null, geneName = chosen.gene || null, organism = chosen.organism || null;
    try {
      const sd = await _json(`https://rest.uniprot.org/uniprotkb/${acc}?fields=sequence,gene_names,organism_name&format=json`);
      if (sd && sd.sequence) { structSeq = sd.sequence.value; structLen = sd.sequence.length; }
      if (sd && sd.genes && sd.genes[0] && sd.genes[0].geneName) geneName = sd.genes[0].geneName.value;
      if (sd && sd.organism && sd.organism.scientificName) organism = sd.organism.scientificName;
    } catch (e) {}
    return {
      accession: acc, cifUrl, matchedBy: found.by,
      organism, gene: geneName,
      matchedLength: chosen.length, expectedLength: length,
      structureSequence: structSeq, structureLength: structLen,
      lengthOk: frag ? true : (structLen == null || length == null || structLen === length),
      fragment: frag ? { start: frag.start, end: frag.end, offset: frag.start - 1, mutations: frag.mutations } : null,
    };
  }

  // ---- align predicted sequence → structure (canonical) sequence; return per-residue map
  //      map[predIndex-1] = structure residue number (1-based) or null. Tiers: identical →
  //      substring offset → Needleman-Wunsch (isoforms with internal indels). ----
  function _identity(n) { const m = new Array(n); for (let i = 0; i < n; i++) m[i] = i + 1; return m; }
  function nwAlign(a, b) {
    const n = a.length, m = b.length, GAP = -1, MIS = -1, MAT = 2;
    if ((n + 1) * (m + 1) > 9e6) {                       // too big → k-mer anchored offset fallback
      const K = 15;
      for (let i = 0; i + K <= n; i += Math.max(1, (n / 50) | 0)) {
        const idx = b.indexOf(a.substr(i, K));
        if (idx >= 0) { const off = idx - i, map = new Array(n); let cov = 0; for (let p = 0; p < n; p++) { const q = p + off; map[p] = (q >= 0 && q < m) ? q + 1 : null; if (map[p] && a[p] === b[q]) cov++; } return { map, method: 'anchored', covered: cov }; }
      }
      return { map: new Array(n).fill(null), method: 'unaligned', covered: 0 };
    }
    const dir = new Uint8Array((n + 1) * (m + 1));
    let prev = new Int32Array(m + 1), cur = new Int32Array(m + 1);
    for (let j = 0; j <= m; j++) { prev[j] = j * GAP; dir[j] = 1; }
    for (let i = 1; i <= n; i++) {
      cur[0] = i * GAP; dir[i * (m + 1)] = 2;
      for (let j = 1; j <= m; j++) {
        const diag = prev[j - 1] + (a[i - 1] === b[j - 1] ? MAT : MIS), up = prev[j] + GAP, left = cur[j - 1] + GAP;
        let best = diag, d = 0;
        if (up > best) { best = up; d = 2; }
        if (left > best) { best = left; d = 1; }
        cur[j] = best; dir[i * (m + 1) + j] = d;
      }
      const t = prev; prev = cur; cur = t;
    }
    const map = new Array(n).fill(null); let cov = 0, i = n, j = m;
    while (i > 0 && j > 0) { const d = dir[i * (m + 1) + j]; if (d === 0) { map[i - 1] = j; if (a[i - 1] === b[j - 1]) cov++; i--; j--; } else if (d === 2) i--; else j--; }
    return { map, method: 'aligned', covered: cov };
  }
  function alignMap(pred, target) {
    if (!pred || !target) return null;
    if (pred === target) return { map: _identity(pred.length), method: 'identical', covered: pred.length };
    let idx = target.indexOf(pred);
    if (idx >= 0) { const map = new Array(pred.length); for (let i = 0; i < pred.length; i++) map[i] = idx + i + 1; return { map, method: idx ? `offset +${idx}` : 'identical', covered: pred.length }; }
    return nwAlign(pred, target);
  }

  // Fetch UniProt domain features for an accession → [{start, end, name}] (canonical coords)
  async function fetchDomains(acc) {
    try {
      const d = await _json(`https://rest.uniprot.org/uniprotkb/${acc}?fields=ft_domain,ft_dna_bind,ft_zn_fing&format=json`);
      const KEEP = { 'Domain': 1, 'DNA binding': 1, 'Zinc finger': 1 };   // homeobox/ZF TFs annotate their key domain as "DNA binding"/"Zinc finger", not "Domain"
      return (d.features || [])
        .filter((f) => KEEP[f.type] && f.location && f.location.start && f.location.end)
        .map((f) => ({ start: +f.location.start.value, end: +f.location.end.value, name: f.description || f.type }))
        .filter((f) => f.start && f.end);
    } catch (e) { return []; }
  }

  async function fetchAlphaMissense(acc) {
    try {
      const meta = await _json(`https://alphafold.ebi.ac.uk/api/prediction/${acc}`);
      const e = Array.isArray(meta) ? meta[0] : meta;
      const url = e && e.amAnnotationsUrl;
      if (!url) return null;
      const txt = await (await fetch(url)).text();
      const sum = {}, cnt = {}, lines = txt.split('\n');
      for (let i = 1; i < lines.length; i++){ const c = lines[i].split(','); if (c.length < 2) continue; const m = c[0].match(/^[A-Z](\d+)[A-Z]$/); const p = parseFloat(c[1]); if (!m || isNaN(p)) continue; const pos = +m[1]; sum[pos] = (sum[pos]||0) + p; cnt[pos] = (cnt[pos]||0) + 1; }
      const out = []; for (const pos in sum) out[+pos] = sum[pos] / cnt[pos];
      return out.length ? out : null;
    } catch (e) { return null; }
  }
  return { crc64, parseFastaToSeqMap, baitSequence, resolveStructure, parseFragmentRange, alignMap, nwAlign, fetchDomains, fetchAlphaMissense, SPECIES };
});
