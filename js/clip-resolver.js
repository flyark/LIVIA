/*
 * cLIP structure resolver — bait identity → UniProt accession → AlphaFold DB CIF.
 * Paths, tried in order:
 *   0. hosted sequence index (flyark.github.io/LIVIA-seqindex; tools/seqindex/build_seqindex.py):
 *      exact sequence → accession with no UniProt call; fragments / point mutants → seed placement
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
    const map = {}, warn = [], allSeqs = [], bySymLen = {};
    const put = (sym, seq) => {
      sym = (sym || '').trim(); seq = clean(seq);
      if (!sym || !seq) return;
      (bySymLen[sym] = bySymLen[sym] || {})[seq.length] = seq;   // every construct of this symbol, keyed by length (mixed full/partial folds)
      if (map[sym] && map[sym] !== seq) warn.push(sym);          // inconsistent sequence across pairs (e.g. a partial construct)
      if (!map[sym] || seq.length > map[sym].length) map[sym] = seq;   // keep the LONGEST (full-length) as the symbol's canonical local sequence
    };
    for (const r of recs) {
      const parts = r.h.split(sep === '___' ? /___|\s*&\s*|_vs_|\s+vs\s+|_VS_|--/ : sep).map((x) => x.trim()).filter(Boolean);   // handle old-style ' & ' and other pair separators
      const chains = r.s.split(':').map(clean).filter(Boolean);
      for (const c of chains) allSeqs.push(c);               // index every chain (length fallback)
      // confident symbol→seq only when header parts line up 1:1 with chains
      // (monomer ">sym", or named pair ">sym1___sym2" with seqA:seqB). Arbitrary headers
      // like ">prediction1" don't map by name — they're resolved by length below.
      if (parts.length === chains.length && parts.length >= 1) parts.forEach((sym, i) => put(sym, chains[i]));
    }
    const uniq = [...new Set(allSeqs)], byLen = {};
    for (const q of uniq) (byLen[q.length] = byLen[q.length] || []).push(q);
    return { map, byLen, bySymLen, seqs: uniq, warn: [...new Set(warn)] };
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
  // Serialize every UniProt/EBI call through one FIFO queue with a minimum gap, so a resolution's
  // many candidate queries — and any partner lookups queued behind them — go out one at a time
  // instead of as a burst. The query is enqueued before partner lookups (those wait until it
  // resolves), so "target first, the rest slowly" falls out of FIFO order. Gentler on the API and
  // avoids self-inflicted rate-limiting.
  const _GAP = 120; let _chain = Promise.resolve(), _lastAt = 0;
  function _throttle(fn) {
    const run = _chain.then(async () => {
      const wait = _GAP - (_now() - _lastAt);
      if (wait > 0) await _sleep(wait);
      _lastAt = _now();
      return fn();
    });
    _chain = run.then(() => {}, () => {});   // keep the queue alive across errors
    return run;
  }
  const _now = () => (typeof performance !== 'undefined' && performance.now) ? performance.now() : Date.now();
  async function _json(url, tries) {                 // retry on 429/5xx/network (transient rate-limits)
    tries = tries || 3;
    for (let i = 0; i < tries; i++) {
      try {
        const r = await _throttle(() => _fetch(url));
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

  // canonical sequence + primary gene + organism of an accession: one UniProt call, shared by the
  // seed check below and resolveStructure (a failed call is not remembered, so a retry can succeed)
  const _uniSeqMemo = new Map();
  function _uniSeq(acc) {
    if (!_uniSeqMemo.has(acc)) _uniSeqMemo.set(acc, _json(`https://rest.uniprot.org/uniprotkb/${acc}?fields=sequence,gene_names,organism_name&format=json`)
      .catch(() => { _uniSeqMemo.delete(acc); return null; }));
    return _uniSeqMemo.get(acc);
  }

  // ---- Hosted sequence index: plain-text shards on GitHub Pages, one small file per lookup, fetched
  //      outside the UniProt throttle (it is our own static host). File layout and hashes are defined
  //      in tools/seqindex/build_seqindex.py. Unreachable or no match → null/[], and the caller falls
  //      through to the UniProt path exactly as before. ----
  const _idxBase = () => (typeof self !== 'undefined' && self.LIVIA_SEQINDEX_BASE) || 'https://flyark.github.io/LIVIA-seqindex/';
  const _idxFiles = new Map();
  function _idxText(path) {
    if (!_idxFiles.has(path)) _idxFiles.set(path, _fetch(_idxBase() + path).then((r) => (r.ok ? r.text() : null), () => null));
    return _idxFiles.get(path);
  }
  const _cleanSeq = (s) => String(s || '').toUpperCase().replace(/[^A-Z]/g, '');
  const _rows = (text, key) => (text || '').split('\n').filter((l) => l.startsWith(key + '\t')).map((l) => l.split('\t'));
  async function _sha256hex(s) {
    const c = (typeof crypto !== 'undefined') && crypto.subtle; if (!c) return null;
    const b = new Uint8Array(await c.digest('SHA-256', new TextEncoder().encode(s)));
    let h = ''; for (const x of b) h += x.toString(16).padStart(2, '0'); return h;
  }
  // When identical sequences exist in several species and nothing else decides: model organisms first.
  const _ORG_ORDER = ['9606', '10090', '7227', '559292', '6239', '7955', '10116', '3702', '83333', '284812'];
  const _orgOrder = (tax) => { const i = _ORG_ORDER.indexOf(String(tax)); return i < 0 ? 99 : i; };
  // Same preference as the UniParc path (canonical > organism match > Swiss-Prot), then a fixed
  // organism order and the accession, so a tie is decided the same way every time.
  const _scoreBy = (orgId) => (x) => (!/-\d+$/.test(x.acc) ? 4 : 0) + (orgId && String(x.tax) === String(orgId) ? 2 : 0) + (x.reviewed ? 1 : 0);
  const _rankBy = (orgId) => {
    const s = _scoreBy(orgId);
    return (a, b) => s(b) - s(a) || _orgOrder(a.tax) - _orgOrder(b.tax) || (a.acc < b.acc ? -1 : a.acc > b.acc ? 1 : 0);
  };
  async function _orgName(tax) {
    const t = await _idxText('t/' + String(Number(tax) % 100).padStart(2, '0') + '.txt');
    const f = _rows(t, String(tax))[0]; return f ? f[1] : null;
  }
  async function _withOrganisms(list) {
    const names = new Map(await Promise.all([...new Set(list.map((x) => x.tax))].map(async (t) => [t, await _orgName(t)])));
    return list.map((x) => Object.assign({}, x, { organism: names.get(x.tax) || null }));
  }

  // Every index entry whose sequence is exactly `seq` → [{acc, gene, tax, reviewed, length}];
  // [] = not in the index; null = index unreachable.
  async function indexExact(seq) {
    seq = _cleanSeq(seq); if (!seq) return [];
    const h = await _sha256hex(seq); if (!h) return null;
    const t = await _idxText('x/' + h.slice(0, 3) + '.txt'); if (t == null) return null;
    return _rows(t, h.slice(3, 16)).filter((f) => +f[5] === seq.length)
      .map((f) => ({ acc: f[1], gene: f[2], tax: f[3], reviewed: f[4] === 's', length: +f[5], afdb: f[6] == null ? undefined : f[6] === '1' }));
  }
  async function _byIndexExact(seq, orgId) {
    const hits = await indexExact(seq); if (!hits || !hits.length) return null;
    hits.sort(_rankBy(orgId));
    const cands = await _withOrganisms(hits.slice(0, 12));
    return { by: 'sequence (exact, LIVIA index)', exact: true,
             candidates: cands.map((x) => ({ acc: x.acc.split('-')[0], isoform: /-\d+$/.test(x.acc) ? x.acc : null,
               length: x.length, organism: x.organism, gene: x.gene || null, tax: x.tax, reviewed: x.reviewed, afdb: x.afdb })) };
  }

  // Seeds: a 12-mer is one when the low 4 bits of H1 are zero (1 in 16 positions, chosen by content,
  // so a fragment picks the same seeds as its parent). FNV-1a + murmur3 fmix32, bit-identical to the build.
  const _FNV = 16777619, _B1 = 2166136261, _B2 = (2166136261 ^ 0x5BD1E995) >>> 0;
  function _kmerHash(s, i, basis) {
    let h = basis;
    for (let j = 0; j < 12; j++) { h ^= s.charCodeAt(i + j); h = Math.imul(h, _FNV); }
    h ^= h >>> 16; h = Math.imul(h, 0x85EBCA6B); h ^= h >>> 13; h = Math.imul(h, 0xC2B2AE35); h ^= h >>> 16;
    return h >>> 0;
  }
  // Place a fragment / point mutant / tagged construct: look up up to `max` of its seeds, spread along the
  // query (a mutation or a tag only spoils the seeds it overlaps), and vote per (protein, diagonal).
  // → [{acc, gene, tax, reviewed, length, votes, offset}] best first; [] no seed hit; null index unreachable.
  async function indexSeeds(seq, max) {
    seq = _cleanSeq(seq); max = max || 12;
    const all = [];
    for (let i = 0; i + 12 <= seq.length; i++) {
      const h1 = _kmerHash(seq, i, _B1);
      if ((h1 & 15) === 0) all.push({ file: ((h1 >>> 4) & 0x3FFF).toString(16).padStart(4, '0'), key: _kmerHash(seq, i, _B2).toString(16).padStart(8, '0'), q: i });
    }
    if (!all.length) return [];
    const pick = all.length <= max ? all : Array.from({ length: max }, (_, k) => all[Math.floor(k * all.length / max)]);
    const texts = await Promise.all(pick.map((s) => _idxText('s/' + s.file + '.txt')));
    if (texts.every((t) => t == null)) return null;
    const votes = new Map();
    pick.forEach((s, k) => { for (const f of _rows(texts[k], s.key)) { const id = f[1] + ':' + (parseInt(f[2], 16) - s.q); votes.set(id, (votes.get(id) || 0) + 1); } });
    const out = [];
    for (const [id, n] of [...votes.entries()].sort((a, b) => b[1] - a[1]).slice(0, 6)) {
      const [pidHex, off] = id.split(':'), pid = parseInt(pidHex, 16);
      const line = ((await _idxText('p/' + (pid >> 8).toString(16).padStart(3, '0') + '.txt')) || '').split('\n')[pid & 255];
      const f = line ? line.split('\t') : [];
      if (f.length >= 5) out.push({ acc: f[0], gene: f[1], tax: f[2], reviewed: f[3] === 's', length: +f[4], afdb: f[5] == null ? undefined : f[5] === '1', votes: n, offset: +off });
    }
    return out;
  }
  // Accept a seed placement only if the query really is that protein: along the voted diagonal, the
  // best-scoring stretch (+1 match, -3 mismatch; a tag or linker falls outside it) must cover half the
  // query at >= 90% identity. Checked against the parent's canonical UniProt sequence.
  async function _byIndexSeeds(seq, orgId) {
    const q = _cleanSeq(seq);
    const hits = await indexSeeds(q); if (!hits || !hits.length || hits[0].votes < 2) return null;   // one lone seed is not evidence
    const best = hits.filter((x) => x.votes === hits[0].votes).sort(_rankBy(orgId))[0];
    const sd = await _uniSeq(best.acc), p = sd && sd.sequence && sd.sequence.value;
    if (!p) return null;
    let run = 0, runM = 0, runN = 0, top = 0, topM = 0, topN = 0;
    for (let i = 0; i < q.length; i++) {
      const j = i + best.offset; if (j < 0 || j >= p.length) continue;
      const m = q[i] === p[j];
      run += m ? 1 : -3; runM += m ? 1 : 0; runN++;
      if (run <= 0) { run = 0; runM = 0; runN = 0; } else if (run > top) { top = run; topM = runM; topN = runN; }
    }
    if (!(topN >= 0.5 * q.length && topM >= 0.9 * topN)) return null;
    const [c] = await _withOrganisms([best]);
    return { by: 'sequence (seed match, LIVIA index)', seed: true, parentSeq: p,
             candidates: [{ acc: c.acc.split('-')[0], length: c.length, organism: c.organism, gene: c.gene || null, tax: c.tax, reviewed: c.reviewed, afdb: c.afdb }] };
  }
  // A complex is almost always one organism: the species the unambiguous chains resolve to settles the
  // chains whose exact sequence exists in several species. → { hint: taxid|null, ambiguous: Map(seq → bool) }
  async function organismHint(seqs) {
    const tally = new Map(), ambiguous = new Map();
    const all = await Promise.all(seqs.map((s) => indexExact(s).catch(() => null)));
    all.forEach((hits, i) => {
      const taxa = new Set((hits || []).map((h) => h.tax));
      ambiguous.set(seqs[i], taxa.size > 1);
      if (taxa.size === 1) { const t = [...taxa][0]; tally.set(t, (tally.get(t) || 0) + 1); }
    });
    let hint = null, n = 0; for (const [t, c] of tally) if (c > n) { hint = t; n = c; }
    return { hint, ambiguous };
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
  // needCif: false for a caller that never uses the AFDB model URL (universal); an index candidate
  // already known to have a model is then taken without asking AFDB.
  async function resolveStructure({ gene, length, seq, species, blast, needCif }) {
    if (!_fetch) throw new Error('fetch unavailable');
    const orgId = SPECIES[species] || species || '';
    const frag = parseFragmentRange(gene, length);      // partial protein? resolve the BASE, offset coords
    let found = null;
    // Hosted index first: an exact match needs no UniProt call at all.
    if (seq && seq.length && !frag) { try { found = await _byIndexExact(seq, orgId); } catch (e) {} }
    // A fragment's sequence isn't a full UniParc entry → checksum won't match; resolve BASE by name.
    if ((!found || !found.candidates.length) && seq && seq.length && !frag) { try { found = await _bySequence(seq, orgId, length); } catch (e) {} }
    const _symName = frag ? frag.base : gene;
    if ((!found || !found.candidates.length) && _symName) { try { found = await _bySymbol(_symName, orgId, frag ? null : length); } catch (e) {} }
    if ((!found || !found.candidates.length) && seq && seq.length && frag) { try { found = await _bySequence(seq, orgId, length); } catch (e) {} }   // (without a fragment name UniParc was already asked above)
    // Point mutant / tagged construct / unnamed fragment: place it by the index's sampled 12-mer seeds.
    if ((!found || !found.candidates.length) && seq && seq.length >= 20 && !frag) { try { found = await _byIndexSeeds(seq, orgId); } catch (e) {} }
    if ((!found || !found.candidates.length) && seq && seq.length >= 60) { try { found = await _byTrimmedSequence(seq, orgId); } catch (e) {} }   // Tier 2: canonical +/- terminal residues
    if ((!found || !found.candidates.length) && blast && seq && seq.length >= 25) {   // Tier 2 (Swiss-Prot) -> Tier 3 (full UniProtKB) similarity BLAST
      try { found = await _byBlast(seq, 'uniprotkb_swissprot'); } catch (e) {}
      if (!found || !found.candidates.length) { try { found = await _byBlast(seq, 'uniprotkb'); } catch (e) {} }
    }
    if (!found || !found.candidates.length) return null;
    // Among tied candidates, prefer one that actually has an AFDB structure (e.g. two 112-aa
    // TrEMBL entries where only one is in AlphaFold DB).
    // First candidate (in preference order) that has an AFDB model. Index candidates carry UniProt's
    // AlphaFoldDB cross-reference, so one known to have no model is skipped without a request (no 404).
    const pickModel = async (cands) => {
      for (const c of cands.slice(0, 12)) {
        if (c.afdb === false) continue;
        if (c.afdb === true && needCif === false) return [c, null];
        const cu = await _afdbCif(c.acc); if (cu) return [c, cu];
      }
      return [null, null];
    };
    let [chosen, cifUrl] = await pickModel(found.candidates);
    // A gene's reference-proteome entry can lack the AFDB model that an identical-sequence entry outside
    // the reference proteome has (fly Orc3: A1Z996 has none, Q7K2L1 has one). The exact index holds those
    // entries too, so a seed match looks its parent's own sequence up there.
    if (!chosen && found.seed && found.parentSeq) {
      const org = found.candidates[0].organism, sib = await _byIndexExact(found.parentSeq, orgId).catch(() => null);
      if (sib) { [chosen, cifUrl] = await pickModel(sib.candidates.filter((c) => !org || c.organism === org)); if (chosen) chosen = Object.assign({}, chosen, { viaSibling: true }); }
    }
    // Still none: UniParc sees every UniProt entry; take the first of the same organism with a model.
    if (!chosen && (found.exact || found.seed)) {
      try {
        const same = found.exact ? _cleanSeq(seq) : found.parentSeq, org = found.candidates[0].organism;
        const up = same ? await _bySequence(same, orgId, same.length) : null;
        for (const c of ((up && up.candidates) || []).filter((c) => !org || c.organism === org).slice(0, 6)) {
          const cu = await _afdbCif(c.acc); if (cu) { chosen = Object.assign({}, c, { viaUniParc: true }); cifUrl = cu; break; }
        }
      } catch (e) {}
    }
    if (!chosen) chosen = found.candidates[0];
    const acc = chosen.acc;
    // canonical UniProt sequence == the AFDB structure's sequence/numbering; used to align
    // the predicted (isoform/construct) sequence and remap cLIR residue coordinates.
    let structSeq = null, structLen = null, geneName = chosen.gene || null, organism = chosen.organism || null;
    if (found.exact && !chosen.isoform && !chosen.viaUniParc && organism) {   // the chain IS this entry's canonical sequence: nothing to fetch
      structSeq = _cleanSeq(seq); structLen = structSeq.length;
    } else if (chosen.viaSibling && organism) {                                 // same canonical sequence as the verified parent
      structSeq = found.parentSeq; structLen = structSeq.length;
    } else try {
      const sd = await _uniSeq(acc);
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

  // Pfam/InterPro domains (fallback when UniProt has no curated Domain features) → [{start, end, name}]
  async function fetchPfam(acc) {
    try {
      const data = await _json(`https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/${acc}?format=json`);
      const out = [];
      for (const entry of (data.results || [])) {
        const name = (entry.metadata && entry.metadata.name) || 'Pfam domain';
        for (const prot of (entry.proteins || []))
          for (const loc of (prot.entry_protein_locations || []))
            for (const frag of (loc.fragments || []))
              if (frag.start != null && frag.end != null) out.push({ name, start: +frag.start, end: +frag.end });
      }
      return out.sort((a, b) => a.start - b.start);
    } catch (e) { return []; }
  }

  // TED structural domains from the AlphaFold DB → [{name, start, end, segments:[{start,end}], cath, tedNo}]
  async function fetchTed(acc) {
    try {
      const data = await _json(`https://alphafold.ebi.ac.uk/api/domains/${acc}`);
      const out = [];
      for (const ann of (data.annotations || [])) {
        const segments = (ann.segments || []).map((s) => ({ start: +s.af_start, end: +s.af_end })).filter((s) => s.start && s.end);
        if (!segments.length) continue;
        const start = Math.min(...segments.map((s) => s.start)), end = Math.max(...segments.map((s) => s.end));
        const cath = ann.cath_label || '';
        out.push({ name: cath ? `TED ${ann.ted_domain_no} (${cath})` : `TED ${ann.ted_domain_no}`, start, end, segments, cath, tedNo: ann.ted_domain_no });
      }
      return out.sort((a, b) => a.start - b.start);
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
  return { crc64, parseFastaToSeqMap, baitSequence, resolveStructure, parseFragmentRange, alignMap, nwAlign, fetchDomains, fetchPfam, fetchTed, fetchAlphaMissense, SPECIES,
           indexExact, indexSeeds, organismHint };
});
