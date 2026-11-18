// LIVIA's build: one string every page shows at its foot and every generated script, report and settings file carries, so an
// export can be traced to the page that made it. Also exposes the read-only window.LIVIA inspection and Markdown summary API.
(function () {
  window.LIVIA_BUILD = '2026-11-18';

  const CUTOFFS = { fpr10: 0.223, fpr5: 0.339, fpr1: 0.551 };

  function runsOf(iter) {
    if (!iter) return [];
    const s = [...iter].map(Number).filter((x) => Number.isFinite(x) && x >= 1).sort((a, b) => a - b);
    const out = [];
    if (!s.length) return out;
    let a = s[0], p = s[0];
    for (let i = 1; i < s.length; i++) {
      if (s[i] === p) continue;
      if (s[i] === p + 1) p = s[i];
      else { out.push([a, p]); a = p = s[i]; }
    }
    out.push([a, p]);
    return out;
  }

  function stretchesOf(rs) {
    const out = [];
    for (let i = 0; i < rs.length; i++) {
      const a = rs[i][0], b = rs[i][1];
      const last = out[out.length - 1];
      if (last && a - last[1] <= 4) last[1] = b;
      else out.push([a, b]);
    }
    return out;
  }

  function fmtRuns(rs) {
    return rs.map((x) => (x[0] === x[1] ? String(x[0]) : x[0] + '-' + x[1])).join(';');
  }

  function expandIndexStr(str) {
    if (!str) return [];
    const out = [];
    const clean = String(str).replace(/[\[\]\x22'\s]/g, '');
    if (!clean) return out;
    const parts = clean.split(/[;,]/);
    for (let i = 0; i < parts.length; i++) {
      const t = parts[i];
      if (!t) continue;
      const ab = t.split('-').map(Number);
      if (ab.length === 2 && Number.isFinite(ab[0]) && Number.isFinite(ab[1]) && ab[1] >= ab[0]) {
        for (let r = ab[0]; r <= ab[1]; r++) out.push(r);
      } else if (Number.isFinite(ab[0]) && ab[0] >= 1) {
        out.push(ab[0]);
      }
    }
    return out;
  }

  function seqMapOf(seqInput) {
    const map = new Map();
    if (!seqInput) return map;
    if (typeof seqInput === 'string') {
      for (let i = 0; i < seqInput.length; i++) map.set(i + 1, seqInput[i]);
      return map;
    }
    if (Array.isArray(seqInput)) {
      for (let i = 0; i < seqInput.length; i++) {
        const it = seqInput[i];
        if (it && it.resnum != null && it.aa) map.set(Number(it.resnum), String(it.aa));
      }
    }
    return map;
  }

  function seqStringOf(seqInput) {
    if (!seqInput) return '';
    if (typeof seqInput === 'string') return seqInput;
    if (Array.isArray(seqInput)) {
      if (!seqInput.length || (seqInput[0] && seqInput[0].type === 'ion')) return '';
      return seqInput.map((x) => (x && x.aa && x.aa.length === 1 ? x.aa : '')).join('');
    }
    return '';
  }

  function seqSnippets(seqInput, rs, maxStretches) {
    const lim = maxStretches || 8;
    const map = seqMapOf(seqInput);
    if (!map.size || !rs.length) return '';
    const st = stretchesOf(rs);
    const parts = [];
    for (let i = 0; i < Math.min(lim, st.length); i++) {
      const a = st[i][0], b = st[i][1];
      let s = '';
      for (let r = a; r <= b; r++) s += map.get(r) || 'X';
      if (!s || /^X+$/.test(s)) continue;
      const disp = s.length > 28 ? s.slice(0, 12) + '...' + s.slice(-8) : s;
      parts.push((a === b ? String(a) : a + '-' + b) + ':' + disp);
    }
    if (!parts.length) return '';
    return parts.join('; ') + (st.length > lim ? '; +' + (st.length - lim) + ' more' : '');
  }

  function rNum(v, d) {
    const n = Number(v);
    return Number.isFinite(n) ? +n.toFixed(d == null ? 3 : d) : null;
  }

  function fNum(v, d) {
    const n = Number(v);
    return Number.isFinite(n) ? n.toFixed(d == null ? 3 : d) : '–';
  }

  function bandTag(ilis) {
    const v = Number(ilis);
    if (!Number.isFinite(v)) return '';
    if (v >= CUTOFFS.fpr1) return '1% FPR';
    if (v >= CUTOFFS.fpr5) return '5% FPR';
    if (v >= CUTOFFS.fpr10) return '10% FPR';
    return 'below 10% FPR';
  }

  function detectPage() {
    const p = (location.pathname.split('/').pop() || 'index.html').replace(/\.html$/i, '').toLowerCase();
    if (['universal', 'clip', 'flypredictome', 'ortholog_predictome', 'dimer', 'monomer', 'network', 'about', 'tutorials', 'bookmarklet', 'index'].includes(p)) {
      return p;
    }
    if (typeof uniState !== 'undefined') return 'universal';
    if (typeof CLIPPipeline !== 'undefined' && document.getElementById('ptable')) return 'clip';
    if (typeof dimerPairResult !== 'undefined') return 'dimer';
    if (typeof monPairResults !== 'undefined') return 'monomer';
    if (typeof RAWPAIRS !== 'undefined' && document.getElementById('dlmenu')) return 'network';
    if (typeof flyParsed !== 'undefined') return 'flypredictome';
    return p || 'index';
  }

  function pageStatus() {
    const page = detectPage();
    let ready = true;
    if (page === 'universal') ready = typeof uniState !== 'undefined' && !!(uniState.modelPairs && uniState.modelPairs.length && uniState.chainNames && uniState.chainNames.length);
    else if (page === 'clip') ready = typeof LAST !== 'undefined' && !!(LAST && LAST.pipe && LAST.cl);
    else if (page === 'flypredictome' || page === 'ortholog_predictome') ready = typeof flyParsed !== 'undefined' && !!(flyParsed && flyParsed.rows && flyParsed.rows.length);
    else if (page === 'dimer') ready = typeof dimerPairResult !== 'undefined' && !!dimerPairResult;
    else if (page === 'monomer') ready = typeof monPairResults !== 'undefined' && !!(monPairResults && monPairResults.length);
    else if (page === 'network') ready = typeof LAST !== 'undefined' && !!(LAST && LAST.nodes && LAST.edges);
    return { page, ready, title: document.title, build: window.LIVIA_BUILD };
  }

  function inspectUniversal() {
    if (typeof uniState === 'undefined' || !uniState.modelPairs || !uniState.modelPairs.length) {
      return { page: 'universal', ready: false };
    }
    const titleEl = document.getElementById('result-title');
    const title = (titleEl && titleEl.textContent ? titleEl.textContent : 'prediction').trim();
    const platform = typeof platformLabel === 'function' ? platformLabel(uniState.platform) : (uniState.platform || '');
    const paeCutEl = document.getElementById('pae-cutoff');
    const cbCutEl = document.getElementById('cb-cutoff');
    const paeCutoff = paeCutEl ? parseFloat(paeCutEl.value) || 12 : 12;
    const cbCutoff = cbCutEl ? parseFloat(cbCutEl.value) || 8 : 8;
    const seedEl = document.getElementById('seed-select');
    const selectedSeed = seedEl ? seedEl.value : '0';
    const structKeys = Object.keys(uniState.structureContents || {});
    const useIdx = selectedSeed === 'avg' ? structKeys[0] : selectedSeed;
    const structText = (uniState.structureContents && (uniState.structureContents[useIdx] || uniState.structureContents[structKeys[0]])) || '';
    let seqMap = new Map();
    if (structText) {
      if (uniState.structureFormat === 'pdb' && typeof parsePdbSequence === 'function') seqMap = parsePdbSequence(structText);
      else if (typeof parseCifSequence === 'function') seqMap = parseCifSequence(structText);
    }
    const md0 = (uniState.modelData && uniState.modelData[0]) || {};
    const cp0 = typeof computeChainPlddt === 'function' && md0.structText ? computeChainPlddt(md0.structText, md0.structFormat) : {};
    const chains = (uniState.chainNames || []).map((cn, idx) => {
      const gName = typeof chainGeneName === 'function' ? chainGeneName(cn) : '';
      const dInfo = (uniState.chainDomainInfo && uniState.chainDomainInfo[cn]) || null;
      const sq = seqMap.get(cn) || [];
      return {
        chain: cn,
        gene: gName || '',
        accession: dInfo && dInfo.accession ? dInfo.accession : '',
        length: (uniState.chainSizes && uniState.chainSizes[idx]) || sq.length || 0,
        type: (uniState.chainTypes && uniState.chainTypes[idx]) || (sq[0] && sq[0].type) || 'protein',
        meanPlddt: cp0[cn] != null ? rNum(cp0[cn], 1) : null,
        sequence: seqStringOf(sq)
      };
    });
    const serializePair = (v, label, ms) => {
      const lirI = runsOf(v.lirI), lirJ = runsOf(v.lirJ);
      const clirI = runsOf(v.clirI), clirJ = runsOf(v.clirJ);
      const sqI = seqMap.get(v.ci) || [], sqJ = seqMap.get(v.cj) || [];
      const gI = typeof chainGeneName === 'function' ? chainGeneName(v.ci) : '';
      const gJ = typeof chainGeneName === 'function' ? chainGeneName(v.cj) : '';
      return {
        model: label,
        chainI: v.ci,
        chainJ: v.cj,
        geneI: gI || v.ci,
        geneJ: gJ || v.cj,
        iLIS: rNum(v.iLIS, 3),
        band: bandTag(v.iLIS),
        iLISA: rNum(v.iLISA, 1),
        ipSAE: rNum(v.ipSAE, 3),
        actifpTM: rNum(v.actifpTM, 3),
        ipTM: rNum(v.ipTM != null ? v.ipTM : ms && ms.ipTM, 3),
        pTM: rNum(ms && ms.pTM, 3),
        LIS: rNum(v.LIS, 3),
        cLIS: rNum(v.cLIS, 3),
        pDockQ: rNum(v.pDockQ, 3),
        LIpDockQ: rNum(v.LIpDockQ, 3),
        pDockQ2: rNum(v.pDockQ2, 3),
        LIpDockQ2: rNum(v.LIpDockQ2, 3),
        lirCountI: v.lirI ? v.lirI.size : 0,
        lirCountJ: v.lirJ ? v.lirJ.size : 0,
        clirCountI: v.clirI ? v.clirI.size : 0,
        clirCountJ: v.clirJ ? v.clirJ.size : 0,
        lirRangesI: fmtRuns(lirI),
        lirRangesJ: fmtRuns(lirJ),
        clirRangesI: fmtRuns(clirI),
        clirRangesJ: fmtRuns(clirJ),
        clirStretchSeqI: seqSnippets(sqI, clirI, 8),
        clirStretchSeqJ: seqSnippets(sqJ, clirJ, 8)
      };
    };
    const models = [];
    for (let mi = 0; mi < uniState.modelPairs.length; mi++) {
      const md = (uniState.modelData && uniState.modelData[mi]) || {};
      const ms = md.scores || {};
      const sym = typeof symmetrizePairs === 'function' ? symmetrizePairs(uniState.modelPairs[mi]) : uniState.modelPairs[mi];
      const entries = Object.values(sym || {});
      for (let k = 0; k < entries.length; k++) {
        models.push(serializePair(entries[k], md.label || ('model_' + (mi + 1)), ms));
      }
    }
    const average = [];
    if (uniState.avgResults && uniState.modelPairs.length > 1) {
      const allMs = (uniState.modelData || []).map((md) => md.scores || {});
      const ptmArr = allMs.filter((s) => s.pTM != null);
      const avgPtm = ptmArr.length ? ptmArr.reduce((s, m) => s + m.pTM, 0) / ptmArr.length : null;
      const entries = Object.values(uniState.avgResults);
      for (let k = 0; k < entries.length; k++) {
        average.push(serializePair(entries[k], 'average', { pTM: avgPtm }));
      }
    }
    const activePair = uniState.currentPairResult ? serializePair(uniState.currentPairResult, selectedSeed === 'avg' ? 'average' : ((uniState.modelData[+selectedSeed] && uniState.modelData[+selectedSeed].label) || ('model_' + (+selectedSeed + 1))), {}) : null;
    return {
      page: 'universal',
      ready: true,
      title,
      platform,
      paeCutoff,
      cbCutoff,
      selectedModel: selectedSeed,
      chains,
      activePair,
      models,
      average
    };
  }

  function inspectClip() {
    if (typeof LAST === 'undefined' || !LAST || !LAST.pipe || !LAST.cl) {
      return { page: 'clip', ready: false };
    }
    const pipe = LAST.pipe, cl = LAST.cl;
    const plen = pipe.proteinLen || 0;
    const qSeq = (typeof seqForSymbol === 'function' ? seqForSymbol(LAST.gene, plen) : '') || LAST._uniprotSeq || '';
    const uniquePartners = new Set(pipe.rows.map((r) => r.Symbol_2)).size;
    const byC = new Map();
    for (let i = 0; i < pipe.rows.length; i++) {
      const c = cl.labels[i];
      let g = byC.get(c);
      if (!g) { g = { c, n: 0, hits: new Map(), best: new Map() }; byC.set(c, g); }
      g.n++;
      const fp = pipe.fingerprints[i] || [];
      for (let j = 0; j < fp.length; j++) g.hits.set(fp[j], (g.hits.get(fp[j]) || 0) + 1);
      const pt = pipe.rows[i].Symbol_2;
      const il = Number(pipe.rows[i].iLIS) || 0;
      if (pt && (!g.best.has(pt) || il > g.best.get(pt))) g.best.set(pt, il);
    }
    const clusters = [...byC.values()].sort((a, b) => b.best.size - a.best.size || b.n - a.n || a.c - b.c).map((g) => {
      const foot = [...g.hits.entries()].filter((kv) => kv[1] / g.n >= 0.3).map((kv) => kv[0]).sort((a, b) => a - b);
      const stretches = stretchesOf(runsOf(foot));
      return {
        cluster: 'C' + g.c,
        predictions: g.n,
        partnerCount: g.best.size,
        footprintStretches: stretches.map((ab) => (ab[0] === ab[1] ? String(ab[0]) : ab[0] + '-' + ab[1])),
        footprintSequence: seqSnippets(qSeq, stretches, 6),
        topPartners: [...g.best.entries()].sort((a, b) => b[1] - a[1]).slice(0, 15).map((kv) => ({ partner: kv[0], bestILIS: rNum(kv[1], 3) }))
      };
    });
    const bestByPartner = new Map();
    for (let i = 0; i < pipe.rows.length; i++) {
      const r = pipe.rows[i], pt = r.Symbol_2, il = Number(r.iLIS) || 0;
      const cur = bestByPartner.get(pt);
      if (!cur || il > cur.iLIS) {
        const qRuns = runsOf(pipe.fingerprints[i] || []);
        const pIdx = expandIndexStr(r.cLIR_indice_B || r.cLIR_indices_j || '');
        bestByPartner.set(pt, {
          partner: pt,
          cluster: 'C' + cl.labels[i],
          rank: r.Rank || r.rank || null,
          iLIS: rNum(il, 3),
          iLISA: rNum(r.iLISA, 1),
          ipTM: rNum(r.ipTM, 2),
          queryCLIRCount: (pipe.fingerprints[i] || []).length,
          queryCLIRRanges: fmtRuns(qRuns),
          partnerCLIRRanges: fmtRuns(runsOf(pIdx))
        });
      }
    }
    const partners = [...bestByPartner.values()].sort((a, b) => (b.iLIS || 0) - (a.iLIS || 0));
    const doms = (LAST._domains || []).map((d) => ({ name: d.name, start: d.start, end: d.end }));
    return {
      page: 'clip',
      ready: true,
      gene: LAST.gene,
      uniprotAccession: LAST._uniprotAcc || '',
      proteinLength: plen,
      sequence: qSeq,
      ilisCutoff: LAST.opts ? LAST.opts.ilisCutoff : CUTOFFS.fpr10,
      partnerCount: uniquePartners,
      predictionCount: pipe.fingerprints.length,
      clusterCount: cl.optimalK,
      domains: doms,
      clusters,
      partners
    };
  }

  function inspectPredictome(page) {
    if (typeof flyParsed === 'undefined' || !flyParsed || !flyParsed.rows || !flyParsed.rows.length) {
      return { page, ready: false };
    }
    const r0 = flyParsed.rows[0] || {};
    const sym1 = r0['Symbol 1'] || r0.chain_i || 'Protein A';
    const sym2 = r0['Symbol 2'] || r0.chain_j || 'Protein B';
    const models = flyParsed.rows.map((r, i) => ({
      rank: r.Rank || r.rank || (i + 1),
      model: r.Model || r.model || '',
      iLIS: rNum(r.iLIS, 3),
      band: bandTag(r.iLIS),
      iLISA: rNum(r.iLISA, 1),
      LIS: rNum(r.LIS, 3),
      cLIS: rNum(r.cLIS, 3),
      ipTM: rNum(r.ipTM, 3),
      pTM: rNum(r.pTM, 3),
      pLDDT: rNum(r.pLDDT, 1),
      lirRangesA: fmtRuns(runsOf(r._LIR_A)),
      lirRangesB: fmtRuns(runsOf(r._LIR_B)),
      clirRangesA: fmtRuns(runsOf(r._cLIR_A)),
      clirRangesB: fmtRuns(runsOf(r._cLIR_B))
    }));
    return { page, ready: true, proteinA: sym1, proteinB: sym2, models };
  }

  function inspectDimer() {
    if (typeof dimerPairResult === 'undefined' || !dimerPairResult) {
      return { page: 'dimer', ready: false };
    }
    const v = dimerPairResult;
    const sd = (typeof selectedDimer !== 'undefined' && selectedDimer) || {};
    const genes = sd.gene || [];
    const accs = sd.uniprotAccession || [];
    const protA = genes[0] || accs[0] || v.ci;
    const protB = genes[1] || genes[0] || accs[1] || accs[0] || v.cj;
    const seqMap = typeof dimerCifText !== 'undefined' && dimerCifText && typeof parseCifSequence === 'function' ? parseCifSequence(dimerCifText) : new Map();
    const sqA = seqMap.get(v.ci) || [], sqB = seqMap.get(v.cj) || [];
    const clirA = runsOf(v.clirI), clirB = runsOf(v.clirJ);
    return {
      page: 'dimer',
      ready: true,
      complexName: sd.complexName || sd.modelEntityId || (protA + ' x ' + protB),
      modelEntityId: sd.modelEntityId || '',
      proteinA: protA,
      proteinB: protB,
      accessionA: accs[0] || '',
      accessionB: accs[1] || accs[0] || '',
      lengthA: v.lenI || sqA.length || 0,
      lengthB: v.lenJ || sqB.length || 0,
      iLIS: rNum(v.iLIS, 3),
      band: bandTag(v.iLIS),
      iLISA: rNum(v.iLISA, 1),
      ipSAE: rNum(v.ipSAE, 3),
      actifpTM: rNum(v.actifpTM, 3),
      ipTM: rNum(v.ipTM, 3),
      LIS: rNum(v.LIS, 3),
      cLIS: rNum(v.cLIS, 3),
      pDockQ: rNum(v.pDockQ, 3),
      LIpDockQ: rNum(v.LIpDockQ, 3),
      pDockQ2: rNum(v.pDockQ2, 3),
      LIpDockQ2: rNum(v.LIpDockQ2, 3),
      lirRangesA: fmtRuns(runsOf(v.lirI)),
      lirRangesB: fmtRuns(runsOf(v.lirJ)),
      clirRangesA: fmtRuns(clirA),
      clirRangesB: fmtRuns(clirB),
      clirStretchSeqA: seqSnippets(sqA, clirA, 8),
      clirStretchSeqB: seqSnippets(sqB, clirB, 8)
    };
  }

  function inspectMonomer() {
    if (typeof monPairResults === 'undefined' || !monPairResults || !monPairResults.length) {
      return { page: 'monomer', ready: false };
    }
    const uid = typeof monUniprotId !== 'undefined' ? monUniprotId : '';
    const regs = (typeof monRegions !== 'undefined' && monRegions ? monRegions : []).map((r) => ({ name: r.name, start: r.start, end: r.end }));
    const pairs = monPairResults.map((v) => ({
      regionA: v.ri ? v.ri.name : '',
      regionB: v.rj ? v.rj.name : '',
      iLIS: rNum(v.iLIS, 3),
      iLISA: rNum(v.iLISA, 1),
      LIS: rNum(v.LIS, 3),
      cLIS: rNum(v.cLIS, 3),
      LIpDockQ: rNum(v.LIpDockQ, 3),
      LIpDockQ2: rNum(v.LIpDockQ2, 3),
      lirRangesA: fmtRuns(runsOf(v.lirI)),
      lirRangesB: fmtRuns(runsOf(v.lirJ)),
      clirRangesA: fmtRuns(runsOf(v.clirI)),
      clirRangesB: fmtRuns(runsOf(v.clirJ))
    }));
    return { page: 'monomer', ready: true, uniprotId: uid, regions: regs, pairs };
  }

  function inspectNetwork() {
    if (typeof LAST === 'undefined' || !LAST || !LAST.nodes || !LAST.edges) {
      return { page: 'network', ready: false };
    }
    const cutEl = document.getElementById('cutv');
    const cutoff = cutEl ? rNum(cutEl.value, 3) : CUTOFFS.fpr10;
    const degMap = new Map(LAST.nodes.map((s, i) => [s, LAST.degree[i]]));
    const cIdx = new Map(LAST.nodes.map((s, i) => [s, LAST.membership[i]]));
    const commObj = {};
    for (let i = 0; i < LAST.nodes.length; i++) {
      const c = LAST.membership[i], s = LAST.nodes[i];
      if (!commObj[c]) commObj[c] = { c, members: [], e: 0, sum: 0 };
      commObj[c].members.push(s);
    }
    for (let i = 0; i < LAST.edges.length; i++) {
      const e = LAST.edges[i];
      if (cIdx.get(e.a) === cIdx.get(e.b)) {
        const c = cIdx.get(e.a);
        commObj[c].e++;
        commObj[c].sum += e.iLIS;
      }
    }
    const labelNode = typeof nodeLabel === 'function' ? nodeLabel : (x) => x;
    const labelComm = typeof commLabel === 'function' ? commLabel : (c) => 'Community ' + (c + 1);
    const communities = Object.values(commObj).filter((m) => m.members.length >= 3).sort((a, b) => b.e - a.e || b.members.length - a.members.length).map((m) => {
      const sorted = m.members.slice().sort((a, b) => (degMap.get(b) || 0) - (degMap.get(a) || 0));
      return {
        id: m.c,
        label: labelComm(m.c),
        proteins: m.members.length,
        edges: m.e,
        meanILIS: m.e ? rNum(m.sum / m.e, 3) : 0,
        topMembers: sorted.slice(0, 15).map(labelNode)
      };
    });
    const topEdges = LAST.edges.slice().sort((a, b) => b.iLIS - a.iLIS).slice(0, 40).map((e) => ({
      a: labelNode(e.a),
      b: labelNode(e.b),
      iLIS: rNum(e.iLIS, 3)
    }));
    return {
      page: 'network',
      ready: true,
      ilisCutoff: cutoff,
      nodeCount: LAST.nodes.length,
      edgeCount: LAST.edges.length,
      communityCount: Math.max(0, ...LAST.membership) + 1,
      communities,
      topEdges
    };
  }

  function inspectCurrent() {
    const page = detectPage();
    if (page === 'universal') return inspectUniversal();
    if (page === 'clip') return inspectClip();
    if (page === 'flypredictome' || page === 'ortholog_predictome') return inspectPredictome(page);
    if (page === 'dimer') return inspectDimer();
    if (page === 'monomer') return inspectMonomer();
    if (page === 'network') return inspectNetwork();
    return { page, ready: true, title: document.title, build: window.LIVIA_BUILD };
  }

  function toMarkdownCurrent() {
    const d = inspectCurrent();
    if (!d.ready) return '# LIVIA (' + d.page + ')\n\nNo prediction is loaded yet.';
    const lines = [];

    if (d.page === 'universal') {
      lines.push('# LIVIA Prediction Summary: ' + d.title);
      lines.push('- **Platform**: ' + (d.platform || 'Structure + PAE') + ' · **Cutoffs**: PAE ≤ ' + d.paeCutoff + ' Å, Cβ ≤ ' + d.cbCutoff + ' Å');
      lines.push('- **iLIS thresholds**: ≥ 0.223 (10% FPR) · ≥ 0.339 (5% FPR) · ≥ 0.551 (1% FPR)');
      if (d.chains.length) {
        lines.push('\n## Chains');
        for (let i = 0; i < d.chains.length; i++) {
          const c = d.chains[i];
          const nm = c.gene && c.gene !== c.chain ? c.chain + ' (' + c.gene + (c.accession ? ', ' + c.accession : '') + ')' : c.chain + (c.accession ? ' (' + c.accession + ')' : '');
          lines.push('- **Chain ' + nm + '**: ' + c.length + ' ' + (c.type === 'dna' || c.type === 'rna' ? 'nt' : 'aa') + (c.meanPlddt != null ? ' · mean pLDDT ' + c.meanPlddt.toFixed(1) : ''));
        }
      }
      const rows = (d.average && d.average.length ? d.average.concat(d.models) : d.models);
      if (rows.length) {
        lines.push('\n## Pairwise Interaction Scores');
        lines.push('| Model | Pair | iLIS | Band | iLISA | ipSAE | actifpTM* | ipTM | LIS | cLIS | pDockQ2 | cLIR (i / j) | LIR (i / j) |');
        lines.push('|---|---|---:|---|---:|---:|---:|---:|---:|---:|---:|---|---|');
        for (let i = 0; i < rows.length; i++) {
          const r = rows[i];
          const pairName = (r.geneI || r.chainI) + ' (' + r.chainI + ') × ' + (r.geneJ || r.chainJ) + ' (' + r.chainJ + ')';
          lines.push('| ' + r.model + ' | ' + pairName + ' | ' + fNum(r.iLIS, 3) + ' | ' + (r.band || '–') + ' | ' + fNum(r.iLISA, 1) + ' | ' + fNum(r.ipSAE, 3) + ' | ' + fNum(r.actifpTM, 3) + ' | ' + fNum(r.ipTM, 3) + ' | ' + fNum(r.LIS, 3) + ' | ' + fNum(r.cLIS, 3) + ' | ' + fNum(r.pDockQ2, 3) + ' | ' + (r.clirRangesI || 'none') + ' / ' + (r.clirRangesJ || 'none') + ' | ' + (r.lirRangesI || 'none') + ' / ' + (r.lirRangesJ || 'none') + ' |');
        }
      }
      const focus = d.activePair || (d.average && d.average[0]) || d.models[0];
      if (focus) {
        lines.push('\n## Contact & Local Interaction Residues (' + focus.model + ': ' + focus.geneI + ' × ' + focus.geneJ + ')');
        lines.push('- **' + focus.geneI + ' (Chain ' + focus.chainI + ') cLIR** (' + focus.clirCountI + ' aa, PAE ≤ ' + d.paeCutoff + ' Å & Cβ ≤ ' + d.cbCutoff + ' Å): `' + (focus.clirRangesI || 'none') + '`' + (focus.clirStretchSeqI ? ' — `' + focus.clirStretchSeqI + '`' : ''));
        lines.push('- **' + focus.geneI + ' (Chain ' + focus.chainI + ') LIR** (' + focus.lirCountI + ' aa, PAE ≤ ' + d.paeCutoff + ' Å): `' + (focus.lirRangesI || 'none') + '`');
        lines.push('- **' + focus.geneJ + ' (Chain ' + focus.chainJ + ') cLIR** (' + focus.clirCountJ + ' aa, PAE ≤ ' + d.paeCutoff + ' Å & Cβ ≤ ' + d.cbCutoff + ' Å): `' + (focus.clirRangesJ || 'none') + '`' + (focus.clirStretchSeqJ ? ' — `' + focus.clirStretchSeqJ + '`' : ''));
        lines.push('- **' + focus.geneJ + ' (Chain ' + focus.chainJ + ') LIR** (' + focus.lirCountJ + ' aa, PAE ≤ ' + d.paeCutoff + ' Å): `' + (focus.lirRangesJ || 'none') + '`');
      }
      return lines.join('\n');
    }

    if (d.page === 'clip') {
      lines.push('# LIVIA cLIP Summary: ' + d.gene + (d.uniprotAccession ? ' (' + d.uniprotAccession + ')' : ''));
      lines.push('- **Query length**: ' + d.proteinLength + ' aa · **iLIS cutoff**: ' + fNum(d.ilisCutoff, 3) + ' · **Partners**: ' + d.partnerCount + ' (' + d.predictionCount + ' predictions across ' + d.clusterCount + ' clusters)');
      if (d.domains.length) {
        lines.push('- **Domains**: ' + d.domains.map((x) => x.name + ' (' + x.start + '-' + x.end + ')').join(', '));
      }
      if (d.clusters.length) {
        lines.push('\n## Predicted Binding-Site Clusters (footprint = ≥ 30% of cluster models)');
        for (let i = 0; i < d.clusters.length; i++) {
          const c = d.clusters[i];
          const fp = c.footprintStretches.length ? c.footprintStretches.join(', ') : 'none ≥ 30%';
          const tops = c.topPartners.slice(0, 8).map((p) => p.partner + ' (' + fNum(p.bestILIS, 3) + ')').join(', ');
          lines.push('- **' + c.cluster + '** (' + c.partnerCount + ' partners, ' + c.predictions + ' models): footprint `' + fp + '`' + (c.footprintSequence ? ' (`' + c.footprintSequence + '`)' : '') + (tops ? ' · top partners: ' + tops : ''));
        }
      }
      if (d.partners.length) {
        const cap = Math.min(30, d.partners.length);
        lines.push('\n## Top Partners (showing ' + cap + ' of ' + d.partners.length + ')');
        lines.push('| Partner | Cluster | Best iLIS | iLISA | ipTM | Query cLIR (' + d.gene + ') | Partner cLIR |');
        lines.push('|---|---|---:|---:|---:|---|---|');
        for (let i = 0; i < cap; i++) {
          const p = d.partners[i];
          lines.push('| ' + p.partner + ' | ' + p.cluster + ' | ' + fNum(p.iLIS, 3) + ' | ' + fNum(p.iLISA, 1) + ' | ' + fNum(p.ipTM, 2) + ' | ' + (p.queryCLIRRanges || 'none') + ' | ' + (p.partnerCLIRRanges || 'none') + ' |');
        }
      }
      return lines.join('\n');
    }

    if (d.page === 'flypredictome' || d.page === 'ortholog_predictome') {
      lines.push('# LIVIA Predictome Pair: ' + d.proteinA + ' × ' + d.proteinB);
      lines.push('| Rank | Model | iLIS | Band | iLISA | ipTM | pLDDT | ' + d.proteinA + ' cLIR | ' + d.proteinB + ' cLIR | ' + d.proteinA + ' LIR | ' + d.proteinB + ' LIR |');
      lines.push('|---|---|---:|---|---:|---:|---:|---|---|---|---|');
      for (let i = 0; i < d.models.length; i++) {
        const m = d.models[i];
        lines.push('| ' + m.rank + ' | ' + (m.model || '–') + ' | ' + fNum(m.iLIS, 3) + ' | ' + (m.band || '–') + ' | ' + fNum(m.iLISA, 1) + ' | ' + fNum(m.ipTM, 3) + ' | ' + fNum(m.pLDDT, 1) + ' | ' + (m.clirRangesA || 'none') + ' | ' + (m.clirRangesB || 'none') + ' | ' + (m.lirRangesA || 'none') + ' | ' + (m.lirRangesB || 'none') + ' |');
      }
      return lines.join('\n');
    }

    if (d.page === 'dimer') {
      lines.push('# LIVIA AFDB Dimer: ' + d.proteinA + ' × ' + d.proteinB + (d.modelEntityId ? ' (' + d.modelEntityId + ')' : ''));
      lines.push('- **iLIS**: ' + fNum(d.iLIS, 3) + ' (' + d.band + ') · **iLISA**: ' + fNum(d.iLISA, 1) + ' · **ipSAE**: ' + fNum(d.ipSAE, 3) + ' · **actifpTM***: ' + fNum(d.actifpTM, 3) + ' · **ipTM**: ' + fNum(d.ipTM, 3) + ' · **pDockQ2**: ' + fNum(d.pDockQ2, 3));
      lines.push('- **' + d.proteinA + '** (' + (d.accessionA || 'Chain A') + ', ' + d.lengthA + ' aa): cLIR `' + (d.clirRangesA || 'none') + '`' + (d.clirStretchSeqA ? ' (`' + d.clirStretchSeqA + '`)' : '') + ' · LIR `' + (d.lirRangesA || 'none') + '`');
      lines.push('- **' + d.proteinB + '** (' + (d.accessionB || 'Chain B') + ', ' + d.lengthB + ' aa): cLIR `' + (d.clirRangesB || 'none') + '`' + (d.clirStretchSeqB ? ' (`' + d.clirStretchSeqB + '`)' : '') + ' · LIR `' + (d.lirRangesB || 'none') + '`');
      return lines.join('\n');
    }

    if (d.page === 'monomer') {
      lines.push('# LIVIA Monomer Domain Analysis: ' + (d.uniprotId || 'Monomer'));
      if (d.regions.length) {
        lines.push('- **Defined regions**: ' + d.regions.map((r) => r.name + ' (' + r.start + '-' + r.end + ')').join(', '));
      }
      lines.push('\n| Region A | Region B | iLIS | iLISA | LIS | cLIS | LIpDockQ2 | Region A cLIR | Region B cLIR |');
      lines.push('|---|---|---:|---:|---:|---:|---:|---|---|');
      for (let i = 0; i < d.pairs.length; i++) {
        const p = d.pairs[i];
        lines.push('| ' + p.regionA + ' | ' + p.regionB + ' | ' + fNum(p.iLIS, 3) + ' | ' + fNum(p.iLISA, 1) + ' | ' + fNum(p.LIS, 3) + ' | ' + fNum(p.cLIS, 3) + ' | ' + fNum(p.LIpDockQ2, 3) + ' | ' + (p.clirRangesA || 'none') + ' | ' + (p.clirRangesB || 'none') + ' |');
      }
      return lines.join('\n');
    }

    if (d.page === 'network') {
      lines.push('# LIVIA Interaction Network');
      lines.push('- **iLIS cutoff**: ' + fNum(d.ilisCutoff, 3) + ' · **Proteins**: ' + d.nodeCount + ' · **Edges**: ' + d.edgeCount + ' · **Communities**: ' + d.communityCount);
      if (d.communities.length) {
        lines.push('\n## Modules (≥ 3 proteins)');
        for (let i = 0; i < d.communities.length; i++) {
          const m = d.communities[i];
          lines.push('- **' + m.label + '** (' + m.proteins + ' proteins, ' + m.edges + ' edges, mean iLIS ' + fNum(m.meanILIS, 3) + '): ' + m.topMembers.join(', '));
        }
      }
      if (d.topEdges.length) {
        lines.push('\n## Strongest Predicted Pairs (top ' + d.topEdges.length + ')');
        lines.push('| Protein A | Protein B | iLIS |');
        lines.push('|---|---|---:|');
        for (let i = 0; i < d.topEdges.length; i++) {
          const e = d.topEdges[i];
          lines.push('| ' + e.a + ' | ' + e.b + ' | ' + fNum(e.iLIS, 3) + ' |');
        }
      }
      return lines.join('\n');
    }

    return '# ' + (d.title || 'LIVIA');
  }

  function copyMarkdownWithFeedback(btn) {
    const md = toMarkdownCurrent();
    const prev = btn.textContent;
    const done = (ok) => {
      btn.textContent = ok ? 'Copied ✓' : 'Copy failed';
      setTimeout(() => { btn.textContent = prev; }, 1600);
    };
    if (navigator.clipboard && navigator.clipboard.writeText) {
      navigator.clipboard.writeText(md).then(() => done(true), () => done(false));
      return;
    }
    try {
      const ta = document.createElement('textarea');
      ta.value = md;
      ta.style.cssText = 'position:fixed;left:-9999px;top:0';
      document.body.appendChild(ta);
      ta.select();
      const ok = document.execCommand('copy');
      ta.remove();
      done(ok);
    } catch (_) {
      done(false);
    }
  }

  function mountCopyButton() {
    if (document.getElementById('livia-copy-md')) return;
    const page = detectPage();
    if (page === 'universal' || page === 'flypredictome' || page === 'ortholog_predictome') {
      const csvLink = document.querySelector('a[onclick^="downloadCSV"], a[onclick^="downloadFlyCSV"]');
      if (!csvLink || !csvLink.parentElement) return;
      const a = document.createElement('a');
      a.id = 'livia-copy-md';
      a.textContent = 'Copy Markdown';
      a.title = 'Copy a structured Markdown summary of scores and residue ranges for notes or LLMs';
      a.style.cssText = 'font-size:0.8rem; color:#2471A3; cursor:pointer; text-decoration:none; border-bottom:1px dashed #2471A3;';
      a.onclick = function () { copyMarkdownWithFeedback(a); };
      csvLink.insertAdjacentElement('afterend', a);
      return;
    }
    if (page === 'dimer' || page === 'monomer') {
      const csvBtn = document.querySelector('a.btn-download[onclick^="downloadCSV"]');
      if (!csvBtn || !csvBtn.parentElement) return;
      const a = document.createElement('a');
      a.id = 'livia-copy-md';
      a.className = 'btn-download cxc';
      a.textContent = 'Copy Markdown';
      a.title = 'Copy a structured Markdown summary of scores and residue ranges for notes or LLMs';
      a.style.cssText = 'background:#2471A3; border-color:#2471A3; font-size:0.82rem; cursor:pointer; margin-left:0.45rem;';
      a.onclick = function () { copyMarkdownWithFeedback(a); };
      csvBtn.insertAdjacentElement('afterend', a);
      return;
    }
    if (page === 'clip') {
      const expBtn = document.getElementById('export');
      if (expBtn && expBtn.parentElement) {
        const b = document.createElement('button');
        b.id = 'livia-copy-md';
        b.type = 'button';
        b.className = 'ghost';
        b.textContent = 'Copy Markdown';
        b.title = 'Copy a structured Markdown summary of clusters, footprints and top partners';
        b.style.marginLeft = 'auto';
        expBtn.style.marginLeft = '0';
        b.onclick = function () { copyMarkdownWithFeedback(b); };
        expBtn.insertAdjacentElement('beforebegin', b);
      }
      return;
    }
    if (page === 'network') {
      const menu = document.querySelector('#dlmenu .panel.menu');
      const expCsv = document.getElementById('expcsv');
      if (!menu || !expCsv) return;
      const b = document.createElement('button');
      b.id = 'livia-copy-md';
      b.type = 'button';
      b.className = expCsv.classList.contains('hidden') ? 'ghost hidden' : 'ghost';
      b.textContent = 'Copy Markdown';
      b.title = 'Copy a structured Markdown summary of modules and strongest pairs';
      b.onclick = function () { copyMarkdownWithFeedback(b); };
      expCsv.insertAdjacentElement('afterend', b);
      if (typeof MutationObserver === 'function') {
        new MutationObserver(function () {
          b.classList.toggle('hidden', expCsv.classList.contains('hidden'));
        }).observe(expCsv, { attributes: true, attributeFilter: ['class'] });
      }
    }
  }

  function stamp() {
    mountCopyButton();
    if (document.getElementById('livia-build')) return;
    const foot = document.querySelector('.footer');
    const d = document.createElement('div'); d.id = 'livia-build';
    if (foot) {
      d.style.cssText = 'margin-top:0.35rem;font-size:0.76rem;color:#6E7781';
      d.innerHTML = 'LIVIA build ' + window.LIVIA_BUILD + ' &middot; <a href="https://github.com/flyark/LIVIA/commits/main" target="_blank" rel="noopener" style="color:#2471A3">changes</a>';
      foot.appendChild(d);
      return;
    }
    d.style.cssText = 'max-width:1100px;margin:24px auto 18px;padding:0 20px;font:12px/1.5 system-ui,sans-serif;color:#5F6771;text-align:center';
    d.innerHTML = 'LIVIA build ' + window.LIVIA_BUILD + ' &middot; <a href="https://github.com/flyark/LIVIA/commits/main" target="_blank" rel="noopener" style="color:#2471A3">changes</a> &middot; Kim &amp; Perrimon (2026), <a href="https://doi.org/10.64898/2026.05.01.721633" target="_blank" rel="noopener" style="color:#2471A3">bioRxiv</a>';
    document.body.appendChild(d);
  }

  window.LIVIA = {
    version: window.LIVIA_BUILD,
    cutoffs: CUTOFFS,
    get page() { return detectPage(); },
    status: pageStatus,
    inspect: inspectCurrent,
    toMarkdown: toMarkdownCurrent,
    waitReady: function (timeoutMs) {
      const limit = timeoutMs == null ? 15000 : timeoutMs;
      return new Promise(function (resolve) {
        const t0 = Date.now();
        const tick = function () {
          if (pageStatus().ready || Date.now() - t0 >= limit) resolve(inspectCurrent());
          else setTimeout(tick, 80);
        };
        tick();
      });
    }
  };

  if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', stamp); else stamp();
})();
