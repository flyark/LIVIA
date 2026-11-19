/**
 * livia-viewer.js — 3D structure viewer utilities for LIVIA tool pages
 *
 * Provides:
 *   parseBfactorsPerResidue()  — extract per-residue B-factors from PDB text (CA atoms)
 *   plddtColor()               — map pLDDT B-factor value to AlphaFold confidence color
 *   buildMolstarPage()         — build Mol* iframe HTML page with MVS coloring
 *
 * Dependencies: none (self-contained)
 */

// ── Parse B-factors per residue (PDB or CIF, including HETATM for ions) ──
// Memoized by text and format: a page re-reads the same structure on every redraw (per chain, per model), so the parse runs
// once per structure and later calls get the same Map back. No caller writes to it. The cache holds the last few structures.
const _BF_CACHE = new Map(), _BF_CACHE_MAX = 24;   // text -> { format -> Map }; keyed by the text itself, whose hash the engine keeps
function parseBfactorsPerResidue(text, format) {
    if (!format) format = text.includes('_atom_site.') ? 'cif' : 'pdb';
    let by = _BF_CACHE.get(text);
    if (by && by[format]) return by[format];
    const m = _parseBfactorsPerResidue(text, format);
    if (!by) { if (_BF_CACHE.size >= _BF_CACHE_MAX) _BF_CACHE.delete(_BF_CACHE.keys().next().value); by = {}; _BF_CACHE.set(text, by); }
    by[format] = m;
    return m;
}
function _parseBfactorsPerResidue(text, format) {
    const m = new Map();
    if (format === 'pdb') {
        for (const line of text.split('\n')) {
            if (line.length < 66) continue;
            if (line.startsWith('ATOM') && line.substring(12, 16).trim() === 'CA') {
                const ch = line.substring(21, 22).trim() || 'A';
                const rn = parseInt(line.substring(22, 26).trim());
                const bf = parseFloat(line.substring(60, 66).trim());
                if (!isNaN(rn) && !isNaN(bf)) m.set(`${ch}:${rn}`, bf);
            } else if (line.startsWith('HETATM')) {
                const ch = line.substring(21, 22).trim() || 'A';
                const rn = parseInt(line.substring(22, 26).trim());
                const bf = parseFloat(line.substring(60, 66).trim());
                if (!isNaN(bf)) m.set(`${ch}:${isNaN(rn) ? 1 : rn}`, bf);
            }
        }
    } else {
        const lines = text.split('\n');
        let inA = false; const cols = []; let ix = null;
        for (const line of lines) {
            if (line.startsWith('_atom_site.')) { inA = true; cols.push(line.trim().split('.')[1]); ix = null; continue; }
            if (!inA || line.startsWith('#') || !line.trim()) continue;
            if (line.startsWith('loop_') || line.startsWith('_')) { inA = false; continue; }
            if (!ix) ix = { group: cols.indexOf('group_PDB'), atom: cols.indexOf('label_atom_id'), ch: cols.indexOf('label_asym_id'), bf: cols.indexOf('B_iso_or_equiv'), rn: cols.indexOf('label_seq_id') };
            // only CA atoms and HETATM records are kept, so other atom lines skip the split
            const het = line.trimStart().startsWith('HETATM');
            if (!het && !/(^|\s)CA(\s|$)/.test(line)) continue;
            const p = line.trim().split(/\s+/);
            if (p.length < cols.length) continue;
            const g = (i) => i >= 0 ? p[i] : '';
            const group = g(ix.group);
            const atom = g(ix.atom);
            const ch = g(ix.ch);
            const bf = parseFloat(g(ix.bf));
            if (isNaN(bf)) continue;
            if (group === 'ATOM' && atom === 'CA') {
                const rn = parseInt(g(ix.rn));
                if (!isNaN(rn)) m.set(`${ch}:${rn}`, bf);
            } else if (group === 'HETATM') {
                const rn = parseInt(g(ix.rn));
                m.set(`${ch}:${isNaN(rn) ? 1 : rn}`, bf);
            }
        }
    }
    // Auto-scale 0–1 pLDDT to 0–100 (e.g. ESMFold2 native PDB stores pLDDT on 0–1).
    // AlphaFold/ColabFold use 0–100; ESMFold uses 0–1. Threshold mx ≤ 1.0 catches the
    // 0–1 case without false-positives on legit low-confidence 0–100 structures.
    let mx = 0;
    for (const v of m.values()) if (v > mx) mx = v;
    if (mx > 0 && mx <= 1.0) {
        for (const [k, v] of m) m.set(k, v * 100);
    }
    return m;
}

// ── Map pLDDT B-factor to AlphaFold confidence color ──
function plddtColor(b) {
    if (b > 90) return '#0053D6';
    if (b > 70) return '#65CBF3';
    if (b > 50) return '#FFDB13';
    return '#FF7D45';
}

// Jmol/CPK element colors for atomic (ball-and-stick) representations of PTM/modified residues.
const ELEMENT_CPK = [['C', '#909090'], ['N', '#3050F8'], ['O', '#FF0D0D'], ['P', '#FF8000'], ['S', '#FFFF30'], ['H', '#FFFFFF']];
const V3D_HL_COLOR = '#E91E63';

// ── Extract 1-letter sequences and residue numbers per chain from PDB or mmCIF (memoized) ──
const _V3D_AA3_TO_1 = {
    ALA:'A', ARG:'R', ASN:'N', ASP:'D', CYS:'C', GLN:'Q', GLU:'E', GLY:'G', HIS:'H', ILE:'I',
    LEU:'L', LYS:'K', MET:'M', PHE:'F', PRO:'P', SER:'S', THR:'T', TRP:'W', TYR:'Y', VAL:'V',
    MSE:'M', SEP:'S', TPO:'T', PTR:'Y', CSO:'C', CSD:'C', CME:'C', PCA:'E', HYP:'P', MLY:'K',
    M3L:'K', ALY:'K', SEC:'U', PYL:'O', ASX:'B', GLX:'Z',
    DA:'A', DT:'T', DG:'G', DC:'C', A:'A', U:'U', G:'G', C:'C'
};
const _V3D_SEQ_CACHE = new Map(), _V3D_SEQ_CACHE_MAX = 24;
function parseChainSequences(text, fmt) {
    if (!text || text === _MOLSTAR_WARMUP_CIF) return { chains: [], byChain: new Map() };
    const f = fmt === 'mmcif' || fmt === 'cif' || (!fmt && text.includes('_atom_site.')) ? 'cif' : 'pdb';
    let by = _V3D_SEQ_CACHE.get(text);
    if (by && by[f]) return by[f];
    const res = _parseChainSequences(text, f);
    if (!by) { if (_V3D_SEQ_CACHE.size >= _V3D_SEQ_CACHE_MAX) _V3D_SEQ_CACHE.delete(_V3D_SEQ_CACHE.keys().next().value); by = {}; _V3D_SEQ_CACHE.set(text, by); }
    by[f] = res;
    return res;
}
function _parseChainSequences(text, format) {
    const raw = new Map();
    const add = (ch, rn, comp) => {
        const aa = _V3D_AA3_TO_1[comp];
        if (!aa || !Number.isFinite(rn)) return;
        let m = raw.get(ch);
        if (!m) { m = new Map(); raw.set(ch, m); }
        if (!m.has(rn)) m.set(rn, aa);
    };
    if (format === 'pdb') {
        for (const line of text.split('\n')) {
            if (line.length < 26) continue;
            if (!line.startsWith('ATOM') && !line.startsWith('HETATM')) continue;
            const atom = line.substring(12, 16).trim();
            if (atom !== 'CA' && atom !== "C4'" && atom !== 'P') continue;
            const comp = line.substring(17, 20).trim().toUpperCase();
            if (!_V3D_AA3_TO_1[comp]) continue;
            const ch = line.substring(21, 22).trim() || 'A';
            const rn = parseInt(line.substring(22, 26).trim(), 10);
            add(ch, rn, comp);
        }
    } else {
        const lines = text.split('\n');
        let inA = false; const cols = []; let ix = null;
        for (const line of lines) {
            if (line.startsWith('_atom_site.')) { inA = true; cols.push(line.trim().split('.')[1]); ix = null; continue; }
            if (!inA || line.startsWith('#') || !line.trim()) continue;
            if (line.startsWith('loop_') || line.startsWith('_')) { inA = false; continue; }
            if (!ix) ix = {
                atom: cols.indexOf('label_atom_id'),
                comp: cols.indexOf('label_comp_id'),
                ch: cols.indexOf('label_asym_id'),
                rn: cols.indexOf('label_seq_id')
            };
            if (!/(^|\s)(CA|C4'|P)(\s|$)/.test(line)) continue;
            const p = line.trim().split(/\s+/);
            if (p.length < cols.length) continue;
            const atom = ix.atom >= 0 ? p[ix.atom].replace(/^['"]|['"]$/g, '') : '';
            if (atom !== 'CA' && atom !== "C4'" && atom !== 'P') continue;
            const comp = (ix.comp >= 0 ? p[ix.comp] : '').toUpperCase();
            if (!_V3D_AA3_TO_1[comp]) continue;
            const ch = ix.ch >= 0 ? p[ix.ch] : 'A';
            const rn = parseInt(ix.rn >= 0 ? p[ix.rn] : '', 10);
            add(ch, rn, comp);
        }
    }
    const chains = [...raw.keys()];
    const byChain = new Map();
    for (const ch of chains) {
        const entries = [...raw.get(ch).entries()].sort((a, b) => a[0] - b[0]);
        const resNums = entries.map(e => e[0]);
        const seq = entries.map(e => e[1]).join('');
        byChain.set(ch, { seq, resNums, resSet: new Set(resNums) });
    }
    return { chains, byChain };
}

// ── Search loaded structure by amino acid motif (X/. wildcard) or residue/range spec ──
function searchStructure3d(seqData, rawQuery) {
    const q = String(rawQuery || '').trim();
    if (!q || !seqData || !seqData.chains.length) return { kind: 'empty', query: q, hits: [] };
    // 1. Try residue / range syntax: e.g. "345", "120-135", "A:345", "A:120-135, B:45"
    const parts = q.split(/[,;]+/).map(s => s.trim()).filter(Boolean);
    const resRe = /^(?:([A-Za-z0-9]{1,4})\s*[:\-]\s*)?(\d+)(?:\s*[-–—..]+\s*(\d+))?$/;
    if (parts.length && parts.every(p => resRe.test(p))) {
        const hits = [];
        for (const p of parts) {
            const m = resRe.exec(p);
            const chReq = m[1] ? m[1].toUpperCase() : '';
            let a = parseInt(m[2], 10), b = m[3] ? parseInt(m[3], 10) : a;
            if (a > b) { const t = a; a = b; b = t; }
            const targetChains = chReq
                ? seqData.chains.filter(c => c.toUpperCase() === chReq)
                : seqData.chains;
            for (const ch of targetChains) {
                const info = seqData.byChain.get(ch);
                if (!info) continue;
                const present = info.resNums.filter(r => r >= a && r <= b);
                if (!present.length) continue;
                const s = present[0], e = present[present.length - 1];
                const i0 = info.resNums.indexOf(s), i1 = info.resNums.indexOf(e);
                const sub = i0 >= 0 && i1 >= i0 ? info.seq.slice(i0, i1 + 1) : '';
                hits.push({ chain: ch, start: s, end: e, seq: sub, label: s === e ? `${ch}:${s}` : `${ch}:${s}\u2013${e}` });
            }
        }
        return { kind: 'residue', query: q, hits };
    }
    // 2. Sequence / motif search (strip FASTA headers, spaces, digits, dashes; X or . = wildcard)
    let clean = q.replace(/^>[^\n]*\n?/gm, '').replace(/[\s0-9-]+/g, '').toUpperCase();
    if (!clean) return { kind: 'sequence', query: q, clean: '', hits: [] };
    if (!/^[A-Z.*]+$/.test(clean)) return { kind: 'invalid', query: q, hits: [] };
    const pat = clean.replace(/[.*+?^${}()|[\]\\]/g, (c) => (c === '.' ? '.' : '\\' + c)).replace(/X/g, '.');
    const re = new RegExp(pat, 'g');
    const hits = [];
    for (const ch of seqData.chains) {
        const info = seqData.byChain.get(ch);
        if (!info || !info.seq) continue;
        let idx = 0;
        while (idx <= info.seq.length - clean.length && hits.length < 50) {
            re.lastIndex = idx;
            const m = re.exec(info.seq);
            if (!m) break;
            const i0 = m.index, i1 = m.index + m[0].length - 1;
            const s = info.resNums[i0], e = info.resNums[i1];
            hits.push({
                chain: ch, start: s, end: e, seq: m[0],
                label: s === e ? `${ch}:${s}` : `${ch}:${s}\u2013${e}`
            });
            idx = m.index + 1;
        }
    }
    return { kind: 'sequence', query: q, clean, hits };
}

// ── Merge vivid highlight spans into base colorComponents without cartoon z-fighting ──
function _subtractRanges(ranges, cuts) {
    let cur = (ranges || []).map(r => ({ start: +r.start, end: +r.end }));
    for (const c of cuts) {
        const next = [];
        for (const r of cur) {
            if (r.end < c.start || r.start > c.end) { next.push(r); continue; }
            if (r.start < c.start) next.push({ start: r.start, end: c.start - 1 });
            if (r.end > c.end) next.push({ start: c.end + 1, end: r.end });
        }
        cur = next;
    }
    return cur;
}
function _mergeSpansByChain(spans) {
    const by = new Map();
    for (const sp of (spans || [])) {
        if (!sp || !sp.chain || !(sp.start > 0)) continue;
        const ch = String(sp.chain), s = +sp.start, e = +(sp.end || sp.start);
        if (!by.has(ch)) by.set(ch, []);
        by.get(ch).push({ start: Math.min(s, e), end: Math.max(s, e) });
    }
    const out = new Map();
    for (const [ch, list] of by) {
        list.sort((a, b) => a.start - b.start || a.end - b.end);
        const merged = [];
        for (const r of list) {
            const last = merged[merged.length - 1];
            if (last && r.start <= last.end + 1) last.end = Math.max(last.end, r.end);
            else merged.push({ start: r.start, end: r.end });
        }
        out.set(ch, merged);
    }
    return out;
}
function _mergeHighlightComponents(baseComponents, spans) {
    if (!spans || !spans.length) return baseComponents || [];
    const byCh = _mergeSpansByChain(spans);
    if (!byCh.size) return baseComponents || [];
    const out = [];
    for (const comp of (baseComponents || [])) {
        if (comp.isIon || comp.stick || comp.isHighlightStick) { out.push(comp); continue; }
        const cuts = byCh.get(String(comp.chain));
        if (!cuts || !comp.ranges) { out.push(comp); continue; }
        const rem = _subtractRanges(comp.ranges, cuts);
        if (rem.length) out.push(Object.assign({}, comp, { ranges: rem }));
    }
    let totalRes = 0;
    for (const ranges of byCh.values()) for (const r of ranges) totalRes += (r.end - r.start + 1);
    for (const [ch, ranges] of byCh) {
        out.push({ chain: ch, ranges, color: V3D_HL_COLOR });
        if (totalRes <= 40) out.push({ chain: ch, ranges, color: V3D_HL_COLOR, isHighlightStick: true });
    }
    return out;
}
function _spansOutsideBase(baseComponents, spans) {
    if (!baseComponents || !baseComponents.length || !spans || !spans.length) return false;
    const drawn = new Map();
    for (const c of baseComponents) {
        if (c.isIon || c.stick || c.isHighlightStick) continue;
        const ch = String(c.chain);
        if (!drawn.has(ch)) drawn.set(ch, []);
        if (c.ranges) drawn.get(ch).push(...c.ranges);
    }
    if (!drawn.size) return false;
    for (const sp of spans) {
        const rs = drawn.get(String(sp.chain));
        if (!rs || !rs.length) return true;
        const rem = _subtractRanges([{ start: +sp.start, end: +(sp.end || sp.start) }], rs);
        if (rem.length) return true;
    }
    return false;
}

const _v3dState = new Map();
function _getV3dState(frameId) {
    const id = frameId || 'viewer3d-frame';
    let st = _v3dState.get(id);
    if (!st) {
        st = { structData: '', fmt: 'mmcif', baseComponents: [], highlightSpans: [], pendingFocus: null, query: '', hits: [], activeHitIdx: -1, autoCompleteNote: '' };
        _v3dState.set(id, st);
    }
    return st;
}

// ── Build the MVS structureChildren (per-component representations + colors) ──
function _buildMvsStructureChildren(colorComponents) {
    const structureChildren = [];
    for (const comp of colorComponents) {
        if (comp.isIon) {
            // Non-polymer chain (ion / ligand / glycan / nucleic) → one ball-and-stick component
            // scoped to THIS chain's label_asym_id, colored by its own chord/legend color.
            // (Scoping by chain avoids the global 'ion'/'ligand' selectors letting the last color win.)
            structureChildren.push({
                kind: 'component',
                params: { selector: { label_asym_id: comp.chain } },
                children: [
                    { kind: 'representation', params: { type: 'ball_and_stick' },
                      children: [{ kind: 'color', params: { color: comp.color } }] }
                ]
            });
            continue;
        }
        if (comp.isHighlightStick) {
            const selector = comp.ranges.map(r => ({
                label_asym_id: comp.chain,
                beg_label_seq_id: r.start,
                end_label_seq_id: r.end,
            }));
            structureChildren.push({
                kind: 'component',
                params: { selector: selector.length === 1 ? selector[0] : selector },
                children: [{
                    kind: 'representation',
                    params: { type: 'ball_and_stick' },
                    children: [{ kind: 'color', params: { color: comp.color } }]
                }]
            });
            continue;
        }
        if (comp.stick) {
            // PTM / modified residues: show the modified sidechain as ball-and-stick over the cartoon
            // (backbone stays cartoon) so the functional group stays bonded to the chain. baseColor =
            // the subunit color for the amino acid part; highlights = per-atom overrides for the
            // functional group (phosphate: P orange, O red). Falls back to CPK if neither is given.
            const oneSel = (a) => ({ label_asym_id: a.chain, label_seq_id: a.seq, label_atom_id: a.atom });
            const sel = comp.atoms.map(oneSel);
            let children;
            if (comp.baseColor || comp.highlights) {
                children = [];
                if (comp.baseColor) children.push({ kind: 'color', params: { color: comp.baseColor } });
                for (const h of (comp.highlights || [])) {
                    const hs = h.atoms.map(oneSel);
                    children.push({ kind: 'color', params: { selector: hs.length === 1 ? hs[0] : hs, color: h.color } });
                }
            } else {
                children = ELEMENT_CPK.map(([el, col]) => ({ kind: 'color', params: { selector: { type_symbol: el }, color: col } }));
            }
            const rep = { kind: 'representation', params: { type: 'ball_and_stick' }, children };
            structureChildren.push({ kind: 'component', params: { selector: sel.length === 1 ? sel[0] : sel }, children: [rep] });
            continue;
        }
        const selector = comp.ranges.map(r => ({
            label_asym_id: comp.chain,
            beg_label_seq_id: r.start,
            end_label_seq_id: r.end,
        }));
        structureChildren.push({
            kind: 'component',
            params: { selector: selector.length === 1 ? selector[0] : selector },
            children: [{
                kind: 'representation',
                params: { type: 'cartoon' },
                children: [{ kind: 'color', params: { color: comp.color } }]
            }]
        });
    }
    return structureChildren;
}

// ── Build full MVS JSON (with placeholder __STRUCT_BLOB_URL__ for the structure URL) ──
function buildMvsJson(colorComponents, fmt) {
    return {
        kind: 'single',
        root: {
            kind: 'root',
            children: [
                { kind: 'canvas', params: { background_color: 'white' } },
                {
                    kind: 'download',
                    params: { url: '__STRUCT_BLOB_URL__' },
                    children: [{
                        kind: 'parse',
                        params: { format: fmt },
                        children: [{
                            kind: 'structure',
                            params: { type: 'model' },
                            children: _buildMvsStructureChildren(colorComponents)
                        }]
                    }]
                }
            ]
        },
        metadata: { version: '1.6' }
    };
}

// ── Send new colors to a Mol* iframe: loadMvsData with the same cell names, so Mol* keeps the downloaded, parsed structure and
// updates only the representations whose colors or selections changed (measured: 14 of 15 cells reused on a recolor) ──
// Iframe must have been built with buildMolstarPage (which embeds the listener).
// fmt: 'pdb' or 'mmcif' (must match initial load).
// Returns true if message was sent; false if iframe isn't ready yet.
function applyColorsToMolstarFrame(frameId, colorComponents, fmt, focusSpans) {
    const id = frameId || 'viewer3d-frame';
    const frame = document.getElementById(id);
    const st = _getV3dState(id);
    st.baseComponents = Array.isArray(colorComponents) ? colorComponents : [];
    if (fmt) st.fmt = fmt;
    const focusTarget = focusSpans || st.pendingFocus || null;
    st.pendingFocus = null;
    if (!frame || !frame.contentWindow) return false;
    const merged = _mergeHighlightComponents(st.baseComponents, st.highlightSpans);
    const mvsJson = buildMvsJson(merged, st.fmt);
    frame.contentWindow.postMessage({
        type: 'updateColors',
        mvsStr: JSON.stringify(mvsJson),
        focusSpans: focusTarget || null,
    }, '*');
    return true;
}

// ── Swap the structure inside an already-booted Mol* iframe (warm-up pattern) ──
// The iframe must have been built with buildMolstarPage. Use this to replace the
// placeholder structure loaded at page-entry warmup with the real prediction structure
// once analysis finishes, without rebuilding the iframe (avoids re-downloading Mol* JS).
function swapMolstarStructure(frameId, structData, colorComponents, fmt, label) {   // label: a readable name for Mol*'s structure panel, in place of the blob URL
    const id = frameId || 'viewer3d-frame';
    const frame = document.getElementById(id);
    const st = _getV3dState(id);
    st.structData = structData || '';
    st.fmt = fmt || 'mmcif';
    st.baseComponents = Array.isArray(colorComponents) ? colorComponents : [];
    st.autoCompleteNote = '';
    if (st.query) {
        const res = searchStructure3d(parseChainSequences(st.structData, st.fmt), st.query);
        st.hits = res.hits;
        st.highlightSpans = res.hits.map(h => ({ chain: h.chain, start: h.start, end: h.end }));
        st.activeHitIdx = res.hits.length === 1 ? 0 : -1;
    } else {
        st.hits = [];
        st.highlightSpans = [];
        st.activeHitIdx = -1;
    }
    if (typeof _syncViewer3dFindUi === 'function') _syncViewer3dFindUi(id);
    if (!frame || !frame.contentWindow) return false;
    const merged = _mergeHighlightComponents(st.baseComponents, st.highlightSpans);
    const mvsJson = buildMvsJson(merged, st.fmt);
    frame.contentWindow.postMessage({
        type: 'loadStructure',
        structData: structData,
        mvsStr: JSON.stringify(mvsJson),
        fmt: st.fmt,
        label: label || '',
        focusSpans: st.highlightSpans.length ? st.highlightSpans : null,
    }, '*');
    return true;
}

// ── Warm-up Mol* in an iframe at page entry with a 1-atom placeholder structure.
// The 4.87 MB Mol* JS finishes downloading while the user is still uploading/
// analyzing; when real results arrive, swapMolstarStructure() updates the structure
// in-place via postMessage — no second iframe build, no second JS fetch. ──
const _MOLSTAR_WARMUP_CIF = `data_warmup
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
ATOM 1 C CA . ALA A 1 1 ? 0.0 0.0 0.0 1.0 50.0 1 A
`;
const _molstarWarmState = new Map();
const _molstarIdleTimers = new Map();
function _cancelMolstarIdle(frameId) {
    const t = _molstarIdleTimers.get(frameId);
    if (!t) return;
    _molstarIdleTimers.delete(frameId);
    if (t.idle && typeof cancelIdleCallback === 'function') cancelIdleCallback(t.id);
    else clearTimeout(t.id);
}
function _runWarmupMolstarNow(frameId, parentBaseUrl) {
    _molstarIdleTimers.delete(frameId);
    const frame = document.getElementById(frameId);
    if (!frame) return false;
    const prev = _molstarWarmState.get(frameId);
    if (prev && prev.warmedUp) return true;
    const curSrc = frame.getAttribute('src') || '';
    if (curSrc && curSrc !== 'about:blank') return false;
    try {
        const base = parentBaseUrl || (window.location.origin + window.location.pathname.replace(/[^/]*$/, ''));
        const page = buildMolstarPage(_MOLSTAR_WARMUP_CIF, 'mmcif', [], base);
        const blob = new Blob([page], { type: 'text/html' });
        if (prev && prev.blobUrl) { try { URL.revokeObjectURL(prev.blobUrl); } catch(_e) {} }
        const url = URL.createObjectURL(blob);
        frame.src = url;
        _molstarWarmState.set(frameId, { warmedUp: true, blobUrl: url });
        const onMsg = (ev) => {
            if (ev.source !== frame.contentWindow) return;
            if (ev.data && ev.data.type === 'molstarFailed') {
                _molstarWarmState.delete(frameId);
                window.removeEventListener('message', onMsg);
            } else if (ev.data && ev.data.type === 'molstarReady') {
                window.removeEventListener('message', onMsg);
            }
        };
        window.addEventListener('message', onMsg);
        return true;
    } catch (e) { console.warn('Mol* warm-up failed for', frameId, e); return false; }
}
function warmupMolstarFrame(frameId, parentBaseUrl) {
    const prev = _molstarWarmState.get(frameId);
    if (prev && prev.warmedUp) return true;
    if (_molstarIdleTimers.has(frameId)) return true;
    const schedule = () => {
        if (typeof requestIdleCallback === 'function') {
            const id = requestIdleCallback(() => _runWarmupMolstarNow(frameId, parentBaseUrl), { timeout: 1200 });
            _molstarIdleTimers.set(frameId, { idle: true, id });
        } else {
            const id = setTimeout(() => _runWarmupMolstarNow(frameId, parentBaseUrl), 550);
            _molstarIdleTimers.set(frameId, { idle: false, id });
        }
    };
    if (typeof document !== 'undefined' && document.readyState !== 'complete') {
        window.addEventListener('load', schedule, { once: true });
    } else {
        schedule();
    }
    return true;
}
function isMolstarFrameWarmedUp(frameId) {
    const s = _molstarWarmState.get(frameId);
    return !!(s && s.warmedUp);
}
let _hlRaf = 0, _hlNext = null, _hlLast = '';
function highlightMolstarResidue(frameId, chain, seqId, endSeqId) {
    const key = chain && seqId ? `${frameId}|${chain}|${seqId}|${endSeqId || seqId}` : `${frameId}|clear`;
    if (key === _hlLast) return;
    _hlNext = { frameId, chain: chain || '', seqId: seqId ? +seqId : 0, endSeqId: endSeqId ? +endSeqId : (seqId ? +seqId : 0), key };
    if (_hlRaf) return;
    _hlRaf = requestAnimationFrame(() => {
        _hlRaf = 0;
        const h = _hlNext; _hlNext = null;
        if (!h || h.key === _hlLast) return;
        const frame = document.getElementById(h.frameId);
        if (!frame || !frame.contentWindow) return;
        const rect = frame.getBoundingClientRect();
        if (rect.width === 0 || rect.height === 0) return;
        _hlLast = h.key;
        try {
            frame.contentWindow.postMessage({ type: 'highlightResidue', chain: h.chain, seqId: h.seqId, endSeqId: h.endSeqId }, '*');
        } catch (_e) {}
    });
}
function showMolstarFocusToast(frame, chain, seqId, endSeqId) {
    if (!frame || typeof document === 'undefined') return;
    const rect = frame.getBoundingClientRect();
    if (!rect.width || !rect.height) return;
    const visH = Math.max(0, Math.min(rect.bottom, window.innerHeight) - Math.max(rect.top, 0));
    if (visH >= Math.min(rect.height * 0.55, window.innerHeight * 0.45)) return;
    let toast = document.getElementById('livia-3d-focus-toast');
    if (!toast) {
        toast = document.createElement('div');
        toast.id = 'livia-3d-focus-toast';
        toast.setAttribute('role', 'status');
        toast.style.cssText = 'position:fixed; right:18px; bottom:18px; z-index:2500; display:none; align-items:center; gap:10px; background:rgba(23,38,58,0.95); color:#fff; padding:8px 12px; border-radius:8px; font-size:12.5px; box-shadow:0 6px 20px rgba(0,0,0,0.22); font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;';
        document.body.appendChild(toast);
    }
    const spanTxt = endSeqId && +endSeqId !== +seqId ? `${seqId}\u2013${endSeqId}` : `${seqId}`;
    const dir = rect.top > window.innerHeight / 2 ? '\u2193' : '\u2191';
    toast.innerHTML = '';
    const msg = document.createElement('span');
    msg.textContent = `Highlighted ${chain}:${spanTxt} in 3D`;
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.textContent = `Jump to 3D ${dir}`;
    btn.style.cssText = 'border:1px solid rgba(255,255,255,0.4); background:#2471A3; color:#fff; border-radius:5px; padding:3px 9px; font-size:12px; font-weight:600; cursor:pointer;';
    btn.onclick = () => {
        const card = frame.closest('.card') || frame;
        card.scrollIntoView({ block: 'center' });
        toast.style.display = 'none';
    };
    toast.append(msg, btn);
    toast.style.display = 'flex';
    clearTimeout(toast._tid);
    toast._tid = setTimeout(() => { toast.style.display = 'none'; }, 3500);
}
function highlightSpansInMolstar(frameId, spans, focusTarget, opts) {
    const id = frameId || 'viewer3d-frame';
    const frame = document.getElementById(id);
    if (!frame || !frame.contentWindow) return false;
    const st = _getV3dState(id);
    const clean = (spans || []).filter(s => s && s.chain && s.start > 0).map(s => ({ chain: String(s.chain), start: +s.start, end: +(s.end || s.start) }));
    const foc = (focusTarget && focusTarget.length ? focusTarget : clean).map(s => ({ chain: String(s.chain), start: +s.start, end: +(s.end || s.start) }));
    st.highlightSpans = clean;
    st.autoCompleteNote = '';
    if (clean.length && _spansOutsideBase(st.baseComponents, clean)) {
        const sc = document.querySelector('.show-complete-cb');
        if (sc && !sc.checked && typeof toggleShowComplete === 'function') {
            st.pendingFocus = foc;
            st.autoCompleteNote = 'Complete structure turned on to show highlighted residues';
            sc.checked = true;
            toggleShowComplete(sc);
            if (typeof _syncViewer3dFindUi === 'function') _syncViewer3dFindUi(id);
            return true;
        }
    }
    if (st.baseComponents && st.baseComponents.length) {
        applyColorsToMolstarFrame(id, st.baseComponents, st.fmt, foc.length ? foc : null);
    } else if (foc.length) {
        try { frame.contentWindow.postMessage({ type: 'focusSpans', spans: foc }, '*'); } catch (_e) {}
    }
    if (!opts || !opts.keepFindUi) {
        if (typeof _syncViewer3dFindUi === 'function') _syncViewer3dFindUi(id);
    }
    return true;
}
function clearMolstarHighlight(frameId) {
    const id = frameId || 'viewer3d-frame';
    const st = _getV3dState(id);
    st.highlightSpans = [];
    st.pendingFocus = null;
    st.query = '';
    st.hits = [];
    st.activeHitIdx = -1;
    st.autoCompleteNote = '';
    const inp = document.getElementById('v3d-find-in');
    if (inp && id === 'viewer3d-frame') inp.value = '';
    if (st.baseComponents && st.baseComponents.length) {
        applyColorsToMolstarFrame(id, st.baseComponents, st.fmt, null);
    }
    if (typeof _syncViewer3dFindUi === 'function') _syncViewer3dFindUi(id);
}
function resetMolstarCamera(frameId) {
    const frame = document.getElementById(frameId || 'viewer3d-frame');
    if (!frame || !frame.contentWindow) return;
    try { frame.contentWindow.postMessage({ type: 'resetCamera' }, '*'); } catch (_e) {}
}
function focusMolstarResidue(frameId, chain, seqId, endSeqId) {
    if (!chain || !seqId) return;
    const id = frameId || 'viewer3d-frame';
    const frame = document.getElementById(id);
    if (!frame || !frame.contentWindow) return;
    const ch = String(chain), s = +seqId, e = +(endSeqId || seqId);
    const st = _getV3dState(id);
    const spanTxt = e !== s ? `${s}\u2013${e}` : `${s}`;
    st.query = '';
    st.hits = [{ chain: ch, start: s, end: e, seq: '', label: `${ch}:${spanTxt}` }];
    st.activeHitIdx = 0;
    const inp = document.getElementById('v3d-find-in');
    if (inp && id === 'viewer3d-frame') inp.value = '';
    highlightSpansInMolstar(id, [{ chain: ch, start: s, end: e }], [{ chain: ch, start: s, end: e }]);
    showMolstarFocusToast(frame, ch, s, e);
}
function downloadMolstarPng(frameId, filename) {
    const frame = document.getElementById(frameId || 'viewer3d-frame');
    if (!frame || !frame.contentWindow) return;
    try {
        frame.contentWindow.postMessage({ type: 'downloadPng', filename: filename || 'livia_3d_structure.png' }, '*');
    } catch (_e) {}
}

// ── Build complete Mol* viewer HTML page with MVS coloring ──
// structData: raw structure text (PDB or mmCIF)
// fmt: 'pdb' or 'mmcif'
// colorComponents: array of { chain, ranges: [{start, end}], color: '#hex' }
// parentBaseUrl: absolute URL of the parent LIVIA page (used for self-hosted Mol* fallback)
const MOLSTAR_VERSION = '5.9.0';
function buildMolstarPage(structData, fmt, colorComponents, parentBaseUrl) {
    let comps = colorComponents || [];
    if (structData !== _MOLSTAR_WARMUP_CIF) {
        _cancelMolstarIdle('viewer3d-frame');
        const st = _getV3dState('viewer3d-frame');
        st.structData = structData || '';
        st.fmt = fmt || 'mmcif';
        st.baseComponents = Array.isArray(colorComponents) ? colorComponents : [];
        if (st.query) {
            const res = searchStructure3d(parseChainSequences(st.structData, st.fmt), st.query);
            st.hits = res.hits;
            st.highlightSpans = res.hits.map(h => ({ chain: h.chain, start: h.start, end: h.end }));
            st.activeHitIdx = res.hits.length === 1 ? 0 : -1;
        }
        comps = _mergeHighlightComponents(st.baseComponents, st.highlightSpans);
        if (typeof _syncViewer3dFindUi === 'function') _syncViewer3dFindUi('viewer3d-frame');
    }
    const mvsJson = buildMvsJson(comps, fmt);
    const selfHostBase = parentBaseUrl || '';

    return `<!DOCTYPE html>
<html><head>
<style>
#viewer1 { position:absolute; top:0; left:0; right:0; bottom:0; }
/* Mol*'s state-snapshot picker ("[1/1] <load time>" + play): every structure loads as one MVS snapshot, so the picker
   only ever offers that one entry, stamped with the load time. Trajectory controls in the same corner stay. */
.msp-state-snapshot-viewport-controls { display: none !important; }
.msp-viewport-controls { display: block !important; }
.msp-plugin .msp-viewport-controls-panel { max-width: calc(100vw - 50px); }
body:not(.livia-ctl) .msp-plugin .msp-layout-standard .msp-layout-right { display: none !important; }
body:not(.livia-ctl) .msp-plugin .msp-layout-standard .msp-layout-main { inset: 0 !important; }
@media (min-width: 821px) {
  body.livia-ctl .msp-plugin .msp-layout-standard .msp-layout-main { inset: 0 300px 0 0 !important; }
  body.livia-ctl .msp-plugin .msp-layout-standard .msp-layout-right { display: block !important; position: absolute !important; inset: 0 0 0 auto !important; width: 300px !important; height: auto !important; border-top: none !important; border-left: 1px solid #cec8ba !important; }
}
@media (min-width: 520px) and (max-width: 820px) {
  body.livia-ctl .msp-plugin .msp-layout-standard .msp-layout-main { inset: 0 280px 0 0 !important; }
  body.livia-ctl .msp-plugin .msp-layout-standard .msp-layout-right { display: block !important; position: absolute !important; inset: 0 0 0 auto !important; width: 280px !important; height: auto !important; border-top: none !important; border-left: 1px solid #cec8ba !important; }
}
@media (max-width: 519px) {
  body.livia-ctl .msp-plugin .msp-layout-standard .msp-layout-main { inset: 0 0 50% 0 !important; }
  body.livia-ctl .msp-plugin .msp-layout-standard .msp-layout-right { display: block !important; position: absolute !important; inset: 50% 0 0 0 !important; width: 100% !important; height: auto !important; border-left: none !important; border-top: 1px solid #cec8ba !important; }
}
#loading-overlay, #error-overlay {
  position: absolute; top: 0; left: 0; right: 0; bottom: 0;
  background: white; display: flex; flex-direction: column;
  align-items: center; justify-content: center; z-index: 1000;
  font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
}
#error-overlay { display: none; }
#loading-overlay .loading-text { font-size: 13px; color: #555; margin-bottom: 14px; }
#loading-overlay .loading-source { font-size: 11px; color: #707070; margin-top: 8px; }
#loading-overlay .spinner-bar { width: 220px; height: 3px; background: #eee; border-radius: 2px; overflow: hidden; position: relative; }
#loading-overlay .spinner-fill { position: absolute; width: 30%; height: 100%; background: #2471A3; animation: livia-slide 1.2s ease-in-out infinite; }
@keyframes livia-slide { 0% { left: -30%; } 100% { left: 100%; } }
#error-overlay .error-text { font-size: 13px; color: #c62828; margin-bottom: 12px; }
#error-overlay button { padding: 6px 16px; background: #2471A3; color: white; border: none; border-radius: 4px; cursor: pointer; font-size: 13px; }
#error-overlay button:hover { background: #1a5276; }
</style>
</head><body>
<div id="viewer1"></div>
<div id="loading-overlay">
  <div class="loading-text">Loading 3D viewer&hellip;</div>
  <div class="spinner-bar"><div class="spinner-fill"></div></div>
  <div class="loading-source" id="loading-source"></div>
</div>
<div id="error-overlay">
  <div class="error-text">Could not load 3D viewer.</div>
  <button onclick="parent.postMessage({type:'molstarRetry'},'*');">Retry</button>
</div>
<script>
var structData = ${JSON.stringify(structData)};
var fmt = "${fmt}";
var mvsTemplate = ${JSON.stringify(JSON.stringify(mvsJson))};
var SELF_HOST_BASE = ${JSON.stringify(selfHostBase)};
var MOLSTAR_VERSION = ${JSON.stringify(MOLSTAR_VERSION)};
var MOLSTAR_SOURCES = [
    {label:'jsdelivr', js:'https://cdn.jsdelivr.net/npm/molstar@'+MOLSTAR_VERSION+'/build/viewer/molstar.js', css:'https://cdn.jsdelivr.net/npm/molstar@'+MOLSTAR_VERSION+'/build/viewer/molstar.css'},
    {label:'unpkg',    js:'https://unpkg.com/molstar@'+MOLSTAR_VERSION+'/build/viewer/molstar.js',           css:'https://unpkg.com/molstar@'+MOLSTAR_VERSION+'/build/viewer/molstar.css'},
    {label:'self-host',js:SELF_HOST_BASE+'lib/molstar/molstar.js',                                            css:SELF_HOST_BASE+'lib/molstar/molstar.css'},
];
var _viewer = null;
var _structUrl = null;
var _ready = false;
var _pendingColorMvs = null;
var _pendingStructure = null;
// The structure arrives as a blob URL, which Mol*'s structure panel would show as its name: once a load has built the
// download cell, its displayed label is replaced with the name the page sent (e.g. "HGTX × Akt, rank 1"). Only the shown label
// changes (cell.obj.label): updating the transform's params instead would make Mol* rebuild the structure and every
// representation under it, a second full load on each structure or color change.
var _structLabel = '', _relabelT = 0, _lastWH = '', _ctlTouched = false;
function _syncLayoutCtl() {
    try {
        var pl = _viewer && _viewer.plugin;
        var on = !!(pl && pl.layout && pl.layout.state && pl.layout.state.showControls);
        var prev = document.body.classList.contains('livia-ctl');
        document.body.classList.toggle('livia-ctl', on);
        if (on !== prev) {
            requestAnimationFrame(function() {
                try { if (pl && pl.canvas3d) pl.canvas3d.handleResize(); } catch(_er) {}
            });
        }
    } catch(_el) {}
}
function _restoreCtl(keep) {
    try {
        var pl = _viewer && _viewer.plugin;
        if (pl && pl.layout && pl.layout.state && !!pl.layout.state.showControls !== !!keep) {
            pl.layout.setProps({ showControls: !!keep });
        }
    } catch(_ec) {}
    _syncLayoutCtl();
}
function _relabelSoon() {   // debounced: the selection cells have no label yet when they are created
    clearTimeout(_relabelT);
    _relabelT = setTimeout(function() {
        try {
            if (!_viewer || !_viewer.plugin) return;
            var pl = _viewer.plugin, st = pl.state.data, nth = {}, changed = [];
            st.cells.forEach(function(cell) {
                var t = cell.transform, pr = t && t.params, o = cell.obj, u = pr && pr.url && (typeof pr.url === 'string' ? pr.url : pr.url.url);
                if (!o) return;
                if (_structLabel && typeof u === 'string' && u.indexOf('blob:') === 0 && o.label !== _structLabel) { o.label = _structLabel; changed.push(cell); }
                // a highlighted residue set reads "Custom Selection: [{label_asym_id: …}]": named by its chain instead
                var m = typeof o.label === 'string' && o.label.indexOf('Custom Selection') === 0 && /label_asym_id: "([^"]+)"/.exec(o.label);
                if (m) { nth[m[1]] = (nth[m[1]] || 0) + 1; o.label = 'Chain ' + m[1] + ', set ' + nth[m[1]]; changed.push(cell); }
            });
            changed.forEach(function(cell) { try { st.events.cell.stateUpdated.next({ state: st, ref: cell.transform.ref, cell: cell }); } catch (_e) {} });   // the panel redraws the names
        } catch (_e) {}
    }, 300);
}

function _setLoadingSource(label) {
    var el = document.getElementById('loading-source');
    if (el) el.textContent = 'Source: ' + label;
}
function _hideLoading() {
    var l = document.getElementById('loading-overlay');
    if (l) l.style.display = 'none';
}
function _showError() {
    _hideLoading();
    var e = document.getElementById('error-overlay');
    if (e) e.style.display = 'flex';
    try { parent.postMessage({ type: 'molstarFailed' }, '*'); } catch(_) {}
}

function _loadMolstarLib(idx) {
    if (idx >= MOLSTAR_SOURCES.length) { _showError(); return; }
    var src = MOLSTAR_SOURCES[idx];
    _setLoadingSource(src.label);
    var link = document.createElement('link');
    link.rel = 'stylesheet'; link.type = 'text/css'; link.href = src.css;
    document.head.appendChild(link);
    var s = document.createElement('script');
    s.src = src.js;
    s.onload = function() {
        init().catch(function(e) { console.error('Mol* init error:', e); _showError(); });
    };
    s.onerror = function() {
        console.warn('Mol* source failed:', src.label, src.js);
        try { link.remove(); s.remove(); } catch(_) {}
        _loadMolstarLib(idx + 1);
    };
    document.head.appendChild(s);
}

// Apply Mol*'s "Illustrative" style preset (matches the UI's Apply Style → Illustrative button)
// Equivalent to: ignoreLight on components + outline + occlusion + shadow off
async function _applyIllustrativeStyle(viewer) {
    try {
        var plugin = viewer.plugin;
        // 1) ignoreLight on all structure components (gives the matte/flat color look)
        if (plugin.managers && plugin.managers.structure && plugin.managers.structure.component) {
            var compMgr = plugin.managers.structure.component;
            await compMgr.setOptions(Object.assign({}, compMgr.state.options, { ignoreLight: true }));
        }
        // 2) Postprocessing: outline ON, occlusion ON, shadow OFF
        // 3) cameraClipping.radius = 0 — keep the "Clipping" slider pinned at 0 so the
        //    whole scene is always visible (no surprise scene-cropping after camera resets
        //    or structure swaps).
        var c3d = plugin.canvas3d;
        if (c3d) {
            var pp = c3d.props.postprocessing || {};
            var clip = (c3d.props && c3d.props.cameraClipping) || {};
            c3d.setProps({
                cameraClipping: Object.assign({}, clip, { radius: 0 }),
                postprocessing: {
                    outline: {
                        name: 'on',
                        params: (pp.outline && pp.outline.name === 'on') ? pp.outline.params
                              : { scale: 1, color: 0x000000, threshold: 0.33, includeTransparent: true }
                    },
                    occlusion: {
                        name: 'on',
                        params: (pp.occlusion && pp.occlusion.name === 'on') ? pp.occlusion.params
                              : { multiScale: { name: 'off', params: {} }, radius: 5, bias: 0.8, blurKernelSize: 15, blurDepthBias: 0.5, samples: 32, resolutionScale: 1, color: 0x000000, transparentThreshold: 0.4 }
                    },
                    shadow: { name: 'off', params: {} }
                }
            });
        }
    } catch(e) { console.warn('Illustrative style error:', e); }
}

// Inject a camera node into an MVS document with the current viewer camera, so
// loadMvsData renders with the user's existing view rather than auto-fitting.
function _injectCurrentCameraIntoMvs(mvsStr) {
    try {
        var c3d = _viewer && _viewer.plugin && _viewer.plugin.canvas3d;
        var cs = c3d && c3d.camera && c3d.camera.state;
        if (!cs || !cs.position || !cs.target || !cs.up) return mvsStr;
        var camNode = {
            kind: 'camera',
            params: {
                target:   [Number(cs.target[0]),   Number(cs.target[1]),   Number(cs.target[2])],
                position: [Number(cs.position[0]), Number(cs.position[1]), Number(cs.position[2])],
                up:       [Number(cs.up[0]),       Number(cs.up[1]),       Number(cs.up[2])],
            }
        };
        var obj = JSON.parse(mvsStr);
        if (!obj.root || !Array.isArray(obj.root.children)) return mvsStr;
        // Strip any existing camera nodes, then insert ours immediately after canvas (or at front)
        obj.root.children = obj.root.children.filter(function(c) { return c.kind !== 'camera'; });
        var insertAt = 0;
        for (var i = 0; i < obj.root.children.length; i++) {
            if (obj.root.children[i].kind === 'canvas') { insertAt = i + 1; break; }
        }
        obj.root.children.splice(insertAt, 0, camNode);
        return JSON.stringify(obj);
    } catch (_e) { return mvsStr; }
}

var _camGen = 0;
function _pinNoClipping() {
    try {
        var c3d = _viewer && _viewer.plugin && _viewer.plugin.canvas3d;
        if (c3d && c3d.props) {
            var clip = c3d.props.cameraClipping || {};
            if (clip.radius !== 0) c3d.setProps({ cameraClipping: Object.assign({}, clip, { radius: 0 }) });
        }
    } catch (_ec) {}
}
function _computeFocusSphere(text, format, spans) {
    if (!text || !spans || !spans.length) return null;
    var byCh = {};
    for (var i = 0; i < spans.length; i++) {
        var sp = spans[i];
        if (!sp || !sp.chain || !(sp.start > 0)) continue;
        var ch = String(sp.chain);
        if (!byCh[ch]) byCh[ch] = [];
        byCh[ch].push([+sp.start, +(sp.end || sp.start)]);
    }
    var inSpan = function(ch, rn) {
        var list = byCh[ch];
        if (!list) return false;
        for (var j = 0; j < list.length; j++) if (rn >= list[j][0] && rn <= list[j][1]) return true;
        return false;
    };
    var pts = [];
    var isCif = format === 'mmcif' || format === 'cif' || text.indexOf('_atom_site.') >= 0;
    var lines = text.split('\\n');
    if (!isCif) {
        for (var k = 0; k < lines.length; k++) {
            var ln = lines[k];
            if (ln.length < 54) continue;
            if (ln.indexOf('ATOM') !== 0 && ln.indexOf('HETATM') !== 0) continue;
            var ch = ln.substring(21, 22).trim() || 'A';
            var rn = parseInt(ln.substring(22, 26).trim(), 10);
            if (!inSpan(ch, rn)) continue;
            var x = parseFloat(ln.substring(30, 38)), y = parseFloat(ln.substring(38, 46)), z = parseFloat(ln.substring(46, 54));
            if (isFinite(x) && isFinite(y) && isFinite(z)) pts.push([x, y, z]);
        }
    } else {
        var inA = false, cols = [], ix = null;
        for (var k2 = 0; k2 < lines.length; k2++) {
            var ln2 = lines[k2];
            if (ln2.indexOf('_atom_site.') === 0) { inA = true; cols.push(ln2.trim().split('.')[1]); ix = null; continue; }
            if (!inA || ln2.charAt(0) === '#' || !ln2.trim()) continue;
            if (ln2.indexOf('loop_') === 0 || ln2.charAt(0) === '_') { inA = false; continue; }
            if (!ix) ix = {
                ch: cols.indexOf('label_asym_id'), ach: cols.indexOf('auth_asym_id'),
                rn: cols.indexOf('label_seq_id'), arn: cols.indexOf('auth_seq_id'),
                x: cols.indexOf('Cartn_x'), y: cols.indexOf('Cartn_y'), z: cols.indexOf('Cartn_z')
            };
            var p = ln2.trim().split(/\\s+/);
            if (p.length < cols.length) continue;
            var ch2 = ix.ch >= 0 ? p[ix.ch] : (ix.ach >= 0 ? p[ix.ach] : 'A');
            var rn2 = parseInt(ix.rn >= 0 ? p[ix.rn] : (ix.arn >= 0 ? p[ix.arn] : ''), 10);
            if (!inSpan(ch2, rn2)) {
                var ach2 = ix.ach >= 0 ? p[ix.ach] : '', arn2 = parseInt(ix.arn >= 0 ? p[ix.arn] : '', 10);
                if (!ach2 || !inSpan(ach2, arn2)) continue;
            }
            var x2 = parseFloat(p[ix.x]), y2 = parseFloat(p[ix.y]), z2 = parseFloat(p[ix.z]);
            if (isFinite(x2) && isFinite(y2) && isFinite(z2)) pts.push([x2, y2, z2]);
        }
    }
    if (!pts.length) return null;
    var cx = 0, cy = 0, cz = 0;
    for (var m = 0; m < pts.length; m++) { cx += pts[m][0]; cy += pts[m][1]; cz += pts[m][2]; }
    cx /= pts.length; cy /= pts.length; cz /= pts.length;
    var r2 = 0;
    for (var n = 0; n < pts.length; n++) {
        var dx = pts[n][0] - cx, dy = pts[n][1] - cy, dz = pts[n][2] - cz;
        var d2 = dx * dx + dy * dy + dz * dz;
        if (d2 > r2) r2 = d2;
    }
    return { center: [cx, cy, cz], radius: Math.max(Math.sqrt(r2) + 10, 18) };
}
function _focusOnSpans(spans, dur) {
    if (!_ready || !_viewer || !spans || !spans.length) return;
    try {
        var pl = _viewer.plugin;
        if (pl && pl.managers && pl.managers.structure && pl.managers.structure.focus) {
            pl.managers.structure.focus.clear();
        }
        _pinNoClipping();
        var sph = _computeFocusSphere(structData, fmt, spans);
        if (sph && pl && pl.managers && pl.managers.camera && typeof pl.managers.camera.focusSphere === 'function') {
            pl.managers.camera.focusSphere(sph, { durationMs: dur == null ? 250 : dur });
            setTimeout(_pinNoClipping, 60);
            setTimeout(_pinNoClipping, 300);
        }
    } catch (_ef) {}
}
function _resetCameraNow() {
    if (!_ready || !_viewer) return;
    try {
        var pl = _viewer.plugin;
        if (pl && pl.managers && pl.managers.structure && pl.managers.structure.focus) pl.managers.structure.focus.clear();
        _pinNoClipping();
        var c3d = pl && pl.canvas3d;
        if (c3d && typeof c3d.requestCameraReset === 'function') c3d.requestCameraReset({ durationMs: 250 });
        else if (pl && pl.managers && pl.managers.camera && typeof pl.managers.camera.reset === 'function') pl.managers.camera.reset(undefined, 250);
    } catch (_er) {}
}

// Listen for parent postMessage: structure swap (warm-up pattern) + color update
window.addEventListener('message', function(ev) {
    if (!ev.data) return;
    if (ev.data.type === 'highlightResidue') {
        if (!_ready || !_viewer || typeof _viewer.structureInteractivity !== 'function') return;
        try {
            _viewer.structureInteractivity({ action: 'highlight' });
            if (ev.data.chain && ev.data.seqId > 0) {
                _viewer.structureInteractivity({
                    action: 'highlight',
                    elements: { label_asym_id: String(ev.data.chain), beg_label_seq_id: +ev.data.seqId, end_label_seq_id: +(ev.data.endSeqId || ev.data.seqId) }
                });
            }
        } catch (_eh) {}
        return;
    }
    if (ev.data.type === 'focusResidue') {
        if (!ev.data.chain || !(ev.data.seqId > 0)) return;
        _camGen++;
        _focusOnSpans([{ chain: String(ev.data.chain), start: +ev.data.seqId, end: +(ev.data.endSeqId || ev.data.seqId) }], 250);
        return;
    }
    if (ev.data.type === 'focusSpans') {
        _camGen++;
        _focusOnSpans(ev.data.spans, 250);
        return;
    }
    if (ev.data.type === 'resetCamera') {
        _camGen++;
        _resetCameraNow();
        return;
    }
    if (ev.data.type === 'downloadPng') {
        if (!_ready || !_viewer) return;
        (async function() {
            var fn = String(ev.data.filename || 'livia_3d_structure.png').replace(/[^\\w.-]+/g, '_');
            if (!/\\.png$/i.test(fn)) fn += '.png';
            try {
                var sh = _viewer.plugin && _viewer.plugin.helpers && _viewer.plugin.helpers.viewportScreenshot;
                if (sh) {
                    var cur = sh.behaviors.values.value;
                    sh.behaviors.values.next(Object.assign({}, cur, { transparent: true }));
                    if (typeof sh.download === 'function') { sh.download(fn); return; }
                    if (typeof sh.getImageDataUri === 'function') {
                        var uri = await sh.getImageDataUri();
                        var a0 = document.createElement('a');
                        a0.href = uri; a0.download = fn; a0.click();
                        return;
                    }
                }
                var cv = document.querySelector('canvas');
                if (cv) {
                    var a = document.createElement('a');
                    a.href = cv.toDataURL('image/png'); a.download = fn; a.click();
                }
            } catch (_ep) {}
        })();
        return;
    }
    if (ev.data.type === 'loadStructure' && ev.data.structData && ev.data.mvsStr) {
        if (!_ready || !_viewer) { _pendingStructure = ev.data; return; }
        _structLabel = ev.data.label || '';
        structData = ev.data.structData;
        if (ev.data.fmt) fmt = ev.data.fmt;
        var loadFocus = ev.data.focusSpans || null;
        var myGen = ++_camGen;
        try {
            if (_structUrl) { try { URL.revokeObjectURL(_structUrl); } catch(_) {} }
            var newBlob = new Blob([ev.data.structData], { type: 'text/plain' });
            _structUrl = URL.createObjectURL(newBlob);
            var newMvs = ev.data.mvsStr.replace('__STRUCT_BLOB_URL__', _structUrl);
            var keepCtl = !!(_viewer.plugin && _viewer.plugin.layout && _viewer.plugin.layout.state && _viewer.plugin.layout.state.showControls);
            var p = _viewer.loadMvsData(newMvs, 'mvsj');
            var refit = function() {
                if (myGen !== _camGen) return;
                _restoreCtl(keepCtl);
                _applyIllustrativeStyle(_viewer);
                var resetCamera = function() {
                    if (myGen !== _camGen) return;
                    try {
                        var c3d = _viewer && _viewer.plugin && _viewer.plugin.canvas3d;
                        var wh = window.innerWidth + 'x' + window.innerHeight;
                        if (c3d && typeof c3d.handleResize === 'function' && wh !== _lastWH) { c3d.handleResize(); _lastWH = wh; }
                        _pinNoClipping();
                        if (loadFocus && loadFocus.length) _focusOnSpans(loadFocus, 0);
                        else if (c3d && typeof c3d.requestCameraReset === 'function') c3d.requestCameraReset();
                        else if (_viewer && _viewer.plugin && _viewer.plugin.managers && _viewer.plugin.managers.camera && typeof _viewer.plugin.managers.camera.reset === 'function') _viewer.plugin.managers.camera.reset();
                    } catch(_e) {}
                };
                resetCamera();
                requestAnimationFrame(function() { resetCamera(); requestAnimationFrame(resetCamera); });
                setTimeout(resetCamera, 100);
                setTimeout(resetCamera, 300);
            };
            if (p && typeof p.then === 'function') p.then(refit, refit);
            else setTimeout(refit, 50);
        } catch (e) { console.warn('Structure swap failed:', e); }
        return;
    }
    if (ev.data.type !== 'updateColors' || !ev.data.mvsStr) return;
    if (!_ready || !_viewer || !_structUrl) {
        _pendingColorMvs = ev.data;
        return;
    }
    try {
        var updFocus = ev.data.focusSpans || null;
        var myGen2 = ++_camGen;
        var newMvs = ev.data.mvsStr.replace('__STRUCT_BLOB_URL__', _structUrl);
        newMvs = _injectCurrentCameraIntoMvs(newMvs);
        var savedCamera = null;
        if (!updFocus || !updFocus.length) {
            try {
                var c3d = _viewer.plugin && _viewer.plugin.canvas3d;
                if (c3d && c3d.camera && c3d.camera.state) {
                    savedCamera = JSON.parse(JSON.stringify(c3d.camera.state));
                }
            } catch(_e) {}
        }
        var restore = function() {
            if (myGen2 !== _camGen) return;
            try {
                _pinNoClipping();
                if (updFocus && updFocus.length) { _focusOnSpans(updFocus, 180); return; }
                var c3dr = _viewer && _viewer.plugin && _viewer.plugin.canvas3d;
                if (savedCamera && c3dr && c3dr.camera) c3dr.camera.setState(savedCamera, 0);
            } catch(_e2) {}
        };
        var keepCtl = !!(_viewer.plugin && _viewer.plugin.layout && _viewer.plugin.layout.state && _viewer.plugin.layout.state.showControls);
        var p = _viewer.loadMvsData(newMvs, 'mvsj');
        var reapply = function() {
            if (myGen2 !== _camGen) return;
            _restoreCtl(keepCtl);
            _applyIllustrativeStyle(_viewer);
            if (updFocus && updFocus.length) {
                _focusOnSpans(updFocus, 250);
                setTimeout(function() { if (myGen2 === _camGen) _focusOnSpans(updFocus, 0); }, 180);
                setTimeout(function() { if (myGen2 === _camGen) _pinNoClipping(); }, 350);
            } else {
                restore();
                requestAnimationFrame(function() { restore(); requestAnimationFrame(restore); });
                setTimeout(restore, 100);
                setTimeout(restore, 300);
            }
        };
        if (p && typeof p.then === 'function') p.then(reapply, reapply);
        else setTimeout(reapply, 50);
    } catch (e) { console.warn('Color update failed:', e); }
});

function _notifyReady() {
    _ready = true;
    try { parent.postMessage({ type: 'molstarReady' }, '*'); } catch(e) {}
    if (_pendingStructure && _viewer) {
        _structLabel = _pendingStructure.label || '';
        structData = _pendingStructure.structData;
        if (_pendingStructure.fmt) fmt = _pendingStructure.fmt;
        var psFocus = _pendingStructure.focusSpans || null;
        var psGen = ++_camGen;
        try {
            if (_structUrl) { try { URL.revokeObjectURL(_structUrl); } catch(_) {} }
            var psBlob = new Blob([_pendingStructure.structData], { type: 'text/plain' });
            _structUrl = URL.createObjectURL(psBlob);
            var psMvs = _pendingStructure.mvsStr.replace('__STRUCT_BLOB_URL__', _structUrl);
            var psKeepCtl = !!(_viewer.plugin && _viewer.plugin.layout && _viewer.plugin.layout.state && _viewer.plugin.layout.state.showControls);
            var psP = _viewer.loadMvsData(psMvs, 'mvsj');
            var psRefit = function() {
                if (psGen !== _camGen) return;
                _restoreCtl(psKeepCtl);
                _applyIllustrativeStyle(_viewer);
                var resetCamera = function() {
                    if (psGen !== _camGen) return;
                    try {
                        var c3d = _viewer && _viewer.plugin && _viewer.plugin.canvas3d;
                        var wh = window.innerWidth + 'x' + window.innerHeight;
                        if (c3d && typeof c3d.handleResize === 'function' && wh !== _lastWH) { c3d.handleResize(); _lastWH = wh; }
                        _pinNoClipping();
                        if (psFocus && psFocus.length) _focusOnSpans(psFocus, 0);
                        else if (c3d && typeof c3d.requestCameraReset === 'function') c3d.requestCameraReset();
                        else if (_viewer && _viewer.plugin && _viewer.plugin.managers && _viewer.plugin.managers.camera && typeof _viewer.plugin.managers.camera.reset === 'function') _viewer.plugin.managers.camera.reset();
                    } catch(_e) {}
                };
                resetCamera();
                requestAnimationFrame(function() { resetCamera(); requestAnimationFrame(resetCamera); });
                setTimeout(resetCamera, 100);
                setTimeout(resetCamera, 300);
            };
            if (psP && typeof psP.then === 'function') psP.then(psRefit, psRefit);
            else setTimeout(psRefit, 50);
        } catch (e) { console.warn('Pending structure load failed:', e); }
        _pendingStructure = null;
        _pendingColorMvs = null;
        return;
    }
    if (_pendingColorMvs && _viewer && _structUrl) {
        try {
            var pmStr = typeof _pendingColorMvs === 'string' ? _pendingColorMvs : _pendingColorMvs.mvsStr;
            var pmFocus = typeof _pendingColorMvs === 'object' ? _pendingColorMvs.focusSpans : null;
            var newMvs = pmStr.replace('__STRUCT_BLOB_URL__', _structUrl);
            var pcKeepCtl = !!(_viewer.plugin && _viewer.plugin.layout && _viewer.plugin.layout.state && _viewer.plugin.layout.state.showControls);
            var p = _viewer.loadMvsData(newMvs, 'mvsj');
            var reapply = function() { _restoreCtl(pcKeepCtl); _applyIllustrativeStyle(_viewer); if (pmFocus && pmFocus.length) _focusOnSpans(pmFocus, 200); };
            if (p && typeof p.then === 'function') p.then(reapply, reapply);
            else setTimeout(reapply, 50);
        } catch (e) { console.warn('Pending color update failed:', e); }
        _pendingColorMvs = null;
    }
}

async function init() {
    // Narrow-frame detection: during warm-up the iframe boots inside display:none so
    // window.innerWidth is 0; estimate the iframe's visible width from parent.innerWidth
    // minus page/card padding (~96 px). Below 820 px the right-hand controls panel starts
    // closed so the 3D structure fills the frame, while .msp-viewport-controls stays visible
    // in the top-right corner so the user can open Controls or Viewport Settings at any time.
    var pw = 0;
    try { pw = (parent && parent !== window && parent.innerWidth) || 0; } catch(_e) {}
    var iw = window.innerWidth;
    var effW = iw > 0 ? iw : (pw > 0 ? Math.max(0, pw - 96) : 1024);
    var startOpen = effW > 820;
    document.body.classList.toggle('livia-ctl', startOpen);
    var viewer = await molstar.Viewer.create('viewer1', {
        layoutIsExpanded: false,
        layoutShowControls: startOpen,
        layoutShowRemoteState: false,
        layoutShowSequence: false,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: false,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
    });
    _viewer = viewer;
    try {
        document.addEventListener('click', function(ev) {
            if (ev.target && ev.target.closest && ev.target.closest('.msp-viewport-controls')) _ctlTouched = true;
        }, true);
        viewer.plugin.layout.events.updated.subscribe(_syncLayoutCtl);
        _syncLayoutCtl();
    } catch (_e0) {}
    var _revealedOnce = iw > 0;
    window.addEventListener('resize', function() {
        var w = window.innerWidth;
        if (!_revealedOnce && w > 0) {
            _revealedOnce = true;
            if (!_ctlTouched && w <= 820 && viewer.plugin && viewer.plugin.layout && viewer.plugin.layout.state && viewer.plugin.layout.state.showControls) {
                try { viewer.plugin.layout.setProps({ showControls: false }); } catch (_er2) {}
            }
        }
        try {
            var c3d = viewer.plugin && viewer.plugin.canvas3d;
            if (c3d) c3d.handleResize();
        } catch (_er3) {}
    });
    try { viewer.plugin.state.data.events.cell.created.subscribe(_relabelSoon); } catch (_e) {}

    var structBlob = new Blob([structData], { type: 'text/plain' });
    var structUrl = URL.createObjectURL(structBlob);
    _structUrl = structUrl;
    var mvsStr = mvsTemplate.replace('__STRUCT_BLOB_URL__', structUrl);

    // Load structure, awaiting promise so reject (not just sync throw) triggers fallback.
    try {
        var p = viewer.loadMvsData(mvsStr, 'mvsj');
        if (p && typeof p.then === 'function') await p;
    } catch(e) {
        console.warn('MVS load failed, falling back to loadStructureFromData:', e);
        try {
            var p2 = viewer.loadStructureFromData(structData, fmt);
            if (p2 && typeof p2.then === 'function') await p2;
        } catch(e2) {
            console.error('Fallback structure load also failed:', e2);
            _showError();
            return;
        }
    }

    // Apply illustrative style after the actual structure load
    try { await _applyIllustrativeStyle(viewer); } catch(_) {}
    try {
        var sh = viewer.plugin.helpers.viewportScreenshot;
        if (sh) {
            var cur = sh.behaviors.values.value;
            sh.behaviors.values.next(Object.assign({}, cur, { transparent: true }));
        }
    } catch(e3) {}
    _hideLoading();
    _notifyReady();
}
_loadMolstarLib(0);
<\/script>
</body></html>`;
}

// ── Parent-side: handle retry requests from inside the iframe ──
// User clicks "Retry" in the iframe's error overlay → iframe posts 'molstarRetry'.
// We rebuild the iframe with the same blob URL so the fallback chain runs again.
if (typeof window !== 'undefined' && !window.__livia_molstar_retry_listener) {
    window.__livia_molstar_retry_listener = true;
    window.addEventListener('message', function(e) {
        if (!e.data || e.data.type !== 'molstarRetry') return;
        const frame = document.getElementById('viewer3d-frame');
        if (!frame || !frame.src) return;
        const parent = frame.parentNode;
        const newFrame = frame.cloneNode(false);
        newFrame.src = frame.src;
        parent.replaceChild(newFrame, frame);
    });
}

// ── 3D color bar: color and display controls beside the structure, plus a domain coloring for the 3D view ──
// Row 1 picks what to color by (Interface, pLDDT, Chain, Polymer, Domains); row 2 shows only that choice's options (the
// palettes, the domain source, the pLDDT scale); row 3 repeats the LIR display (fill gaps, min segment, complete structure,
// gray non-LIR). Every control drives its twin in the Visualization Scripts card (a click there and here run the same handler),
// so both places agree and the scripts follow. Domains color the 3D view only: each domain its own color (the community colors:
// Tableau 10, then Tableau 20's pale pairs), every other residue white, with a key under the viewer.
// Pages pass their components through applyViewer3dDomains(); the bar mounts itself on pages with the presets.
let viewer3dDomainMode = '';   // '' | 'uniprot' | 'ted'
const DOM3D_COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#bcbd22', '#17becf',
    '#aec7e8', '#ffbb78', '#98df8a', '#ff9896', '#c5b0d5', '#c49c94', '#f7b6d2', '#dbdb8d', '#9edae5'];
const DOM3D_NAMES = { uniprot: 'UniProt/Pfam domains', ted: 'TED domains' };

function mountViewer3dColorBar() {
    const frame = document.getElementById('viewer3d-frame');
    if (!frame || document.getElementById('viewer3d-colorbar')) return;
    const rows = [...document.querySelectorAll('.presets-row')];
    const palRow = rows.find(r => r.querySelector('[onclick^="applyPreset"]'));
    const modeRow = rows.find(r => r.querySelector('[onclick^="applyCxcPreset"]'));
    if (!palRow && !modeRow) return;
    if (!document.getElementById('v3d-colorbar-css')) {
        const st = document.createElement('style'); st.id = 'v3d-colorbar-css';
        st.textContent = `.v3d-bar { display:flex; flex-direction:column; gap:0.4rem; margin:0 0 0.55rem; font-size:0.82rem; color:#555; }
.v3d-row { display:flex; flex-wrap:wrap; align-items:center; gap:0.35rem 0.4rem; min-height:28px; }
.v3d-lab { font-size:0.7rem; font-weight:600; letter-spacing:0.06em; text-transform:uppercase; color:#7A8899; min-width:5.6rem; }
.v3d-seg { display:inline-flex; flex-wrap:wrap; border:1px solid #DCE3EB; border-radius:7px; overflow:hidden; }
.v3d-seg button { border:0; border-right:1px solid #DCE3EB; background:#fff; color:#4A596F; font-family:inherit; font-size:0.8rem; font-weight:600; padding:0.3rem 0.7rem; cursor:pointer; }
.v3d-seg button:last-child { border-right:0; }
.v3d-seg button.on { background:#1A5276; color:#fff; }
.v3d-chip { display:inline-flex; align-items:center; gap:0.35rem; padding:0.24rem 0.5rem; border:1px solid #DCE3EB; border-radius:6px; background:#fff; cursor:pointer; font:inherit; font-size:0.78rem; color:#17263A; }
.v3d-chip:hover, .v3d-seg button:not(.on):hover { border-color:#2471A3; color:#2471A3; }
.v3d-chip.active { border-color:#1A5276; box-shadow:0 0 0 1px #1A5276; background:#EEF4FA; }
.v3d-bar button:focus-visible, .v3d-bar input:focus-visible { outline:2px solid #2471A3; outline-offset:1px; }
.v3d-chip .preset-strip { width:34px; height:10px; }
.v3d-more { background:none; border:0; padding:0.2rem 0.3rem; color:#2471A3; cursor:pointer; font:inherit; font-size:0.78rem; text-decoration:underline; }
.v3d-note { color:#7A8899; font-size:0.78rem; }
.v3d-num, .v3d-apply { box-sizing:border-box; height:26px; margin:0; vertical-align:middle; border-radius:5px; font-family:inherit; font-size:0.8rem; line-height:1; }
.v3d-num { width:3.4rem; padding:0 0.3rem; border:1px solid #DCE3EB; text-align:center; }
.v3d-apply { display:inline-flex; align-items:center; border:1px solid #2C6E9F; background:#2C6E9F; color:#fff; padding:0 0.65rem; font-weight:600; cursor:pointer; }
.v3d-bar label.v3d-cb { display:inline-flex; align-items:center; gap:0.3rem; margin:0; font-size:0.8rem; cursor:pointer; user-select:none; }
.v3d-pipe { color:#ccc; }
.v3d-sep { width:1px; align-self:stretch; background:#DCE3EB; margin:0 0.2rem; }
.v3d-key { display:flex; flex-wrap:wrap; gap:0.25rem 0.9rem; align-items:center; margin-top:0.45rem; font-size:0.8rem; color:#444; }
.v3d-key i { display:inline-block; width:11px; height:11px; border-radius:2px; margin-right:0.3rem; vertical-align:-1px; border:1px solid rgba(0,0,0,0.25); }
.v3d-key .muted { color:#888; }
@media (max-width:600px) { .v3d-lab { min-width:0; width:100%; } .v3d-pipe, .v3d-sep { display:none; } }`;
        document.head.appendChild(st);
    }
    const bar = document.createElement('div'); bar.id = 'viewer3d-colorbar'; bar.className = 'v3d-bar';
    frame.parentElement.insertBefore(bar, frame);
    const key = document.createElement('div'); key.id = 'viewer3d-domkey'; key.className = 'v3d-key'; key.hidden = true;
    frame.insertAdjacentElement('afterend', key);
    const redraw = () => { if (typeof onColorChange === 'function') onColorChange(); };
    const esc = (t) => String(t).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
    // the palette row the page shows now (universal swaps in its multi-chain row for 3+ chains) and the mode chips by name
    const visPal = () => { const mc = document.getElementById('multichain-presets'); return mc && mc.style.display !== 'none' && mc.querySelector('.preset-chip') ? mc : palRow; };
    const modeChip = (m) => modeRow && modeRow.querySelector(`[onclick*="'${m}'"]`);
    const current = () => { if (viewer3dDomainMode) return 'domains';
        for (const m of ['plddt', 'bychain', 'bypolymer']) { const c = modeChip(m); if (c && c.classList.contains('active')) return m; }
        return 'interface'; };
    const MODES = [['interface', 'Interface'], ['plddt', 'pLDDT'], ['bychain', 'Chain'], ['bypolymer', 'Polymer'], ['domains', 'Domains']].filter(([m]) => m === 'interface' || m === 'domains' || modeChip(m));
    const gap = document.getElementById('gap-fill-input'), seg = document.getElementById('min-segment-input');
    // Shade: Gradient (LIR light, cLIR dark) or Solid (one color per chain) for any palette, so the bar needs no separate solid chips.
    // Two chains: Solid = each palette's cLIR color for both LIR and cLIR, read off the chip's strip (LIR A, cLIR A, cLIR B, LIR B);
    // the Scripts card's solid chip with those colors lights up when there is one. 3+ chains: the card's own Gradient/Solid chips.
    const mcRow = () => { const mc = document.getElementById('multichain-presets'); return mc && mc.style.display !== 'none' && mc.querySelector('.preset-chip') ? mc : null; };
    const hexOf = (c) => { const m = String(c).match(/\d+/g); return m && /rgb/.test(c) ? '#' + m.slice(0, 3).map((x) => (+x).toString(16).padStart(2, '0')).join('') : String(c).toLowerCase(); };
    const cols = (o) => [...o.querySelectorAll('.preset-strip span')].map((x) => hexOf(x.style.background));
    const val = (id) => { const e = document.getElementById(id); return e ? String(e.value).toLowerCase() : ''; };
    const pals = () => { const mc = mcRow(); return mc ? [...mc.querySelectorAll('.palette-chip')] : palRow ? [...palRow.querySelectorAll('.preset-chip')].slice(0, 6) : []; };
    // Four shades of a palette (2026-10-16, as the Atlas pair page): Gradient, Solid, Pale LIR (LIR = the cLIR color faded to 20%) and
    // Contacts only (LIR white), so the contact residues stand out where LIR covers most of a chain. Two chains: read off the four color
    // inputs; 3+ chains: the card's own shade chips (universal's multi-chain row, data-shade).
    const fade = (h) => '#' + [1, 3, 5].map((i) => Math.round(255 - (255 - parseInt(h.slice(i, i + 2), 16)) * 0.2).toString(16).padStart(2, '0')).join('');
    const SHADES = [['gradient', 'Gradient', 'LIR light, cLIR dark'], ['solid', 'Solid', 'one color per chain'], ['pale', 'Pale LIR', 'LIR faded close to white, cLIR in color'], ['clir', 'Contacts only', 'LIR white, only cLIR in color']];
    const chipShade = (x) => x.dataset.shade || 'gradient';   // universal's shade chips carry their shade
    const shadeNow = () => { const mc = mcRow(); if (mc) { const a = mc.querySelector('.mode-chip.active'); return a ? chipShade(a) : 'gradient'; }
        const la = val('color-lir-a'), ca = val('color-clir-a'), lb = val('color-lir-b'), cb = val('color-clir-b'); if (!ca || !cb) return 'gradient';
        if (la === ca && lb === cb) return 'solid'; if (la === '#ffffff' && lb === '#ffffff') return 'clir'; if (la === fade(ca) && lb === fade(cb)) return 'pale'; return 'gradient'; };
    const shadeCols = (sh, ca, cb) => (sh === 'solid' ? [ca, ca, cb, cb] : sh === 'pale' ? [fade(ca), ca, fade(cb), cb] : sh === 'clir' ? ['#ffffff', ca, '#ffffff', cb] : null);   // LIR A, cLIR A, LIR B, cLIR B
    const sameDark = (o) => { const c = cols(o); return c.length === 4 && c[1] === val('color-clir-a') && c[2] === val('color-clir-b'); };
    const setTwo = (la, ca, lb, cb) => {
        const twin = [...document.querySelectorAll('.presets-row .preset-chip')].find((x) => { const c = cols(x); return c.length === 4 && c[0] === la && c[1] === ca && c[2] === cb && c[3] === lb; });
        if (twin) { twin.click(); return; }
        document.querySelectorAll('.preset-chip.active').forEach((x) => x.classList.remove('active'));
        if (typeof applyPreset === 'function') applyPreset(la, ca, lb, cb, null);
    };
    const lighten = (h) => '#' + [1, 3, 5].map((i) => Math.round(parseInt(h.slice(i, i + 2), 16) * 0.45 + 255 * 0.55).toString(16).padStart(2, '0')).join('');
    const pickPal = (o) => { viewer3dDomainMode = ''; const sh = shadeNow();
        if (mcRow() || sh === 'gradient') { o.click(); return; }
        const c = cols(o); if (c.length === 4) setTwo(...shadeCols(sh, c[1], c[2])); else o.click(); };
    const setShade = (sh) => { const mc = mcRow();
        if (mc) { const o = [...mc.querySelectorAll('.mode-chip')].find((x) => chipShade(x) === sh); if (o) o.click(); return; }
        const ca = val('color-clir-a'), cb = val('color-clir-b'); if (!ca || !cb) return;
        if (sh !== 'gradient') { setTwo(...shadeCols(sh, ca, cb)); return; }
        const grad = [...document.querySelectorAll('.presets-row .preset-chip')].find((x) => { const c = cols(x); return c.length === 4 && c[0] !== c[1] && c[1] === ca && c[2] === cb; });
        if (grad) grad.click(); else setTwo(lighten(ca), ca, lighten(cb), cb); };
    if (palRow && !document.getElementById('presets-pale')) {   // the card carries the same shades: two preset rows under its solid presets, from its six gradient palettes
        const solidLab = [...document.querySelectorAll('label')].find((l) => /Presets\s*\u2014\s*Solid/i.test(l.textContent)), solidRow = solidLab && solidLab.nextElementSibling;
        if (solidRow && solidRow.classList.contains('presets-row')) {
            const grads = [...palRow.querySelectorAll('.preset-chip')].slice(0, 6).map((o) => ({ c: cols(o), t: (o.textContent || '').trim() })).filter((x) => x.c.length === 4);
            const span = (c) => `<span style="background:${c}${c === '#ffffff' ? ';box-shadow:inset 0 0 0 1px #C9D2DC' : ''}"></span>`;   // a white LIR outlined
            let after = solidRow;
            for (const [id, text, sh] of [['presets-pale', 'Presets \u2014 pale LIR (contacts stand out)', 'pale'], ['presets-contacts', 'Presets \u2014 contacts only (LIR white)', 'clir']]) {
                const lab = solidLab.cloneNode(false), row = document.createElement('div'); lab.textContent = text; row.className = solidRow.className; row.id = id;
                row.innerHTML = grads.map(({ c, t }) => { const [la, ca, lb, cb] = shadeCols(sh, c[1], c[2]);
                    return `<div class="preset-chip" onclick="applyPreset('${la}','${ca}','${lb}','${cb}',this)"><div class="preset-strip">${span(la)}${span(ca)}${span(cb)}${span(lb)}</div>${esc(t)}</div>`; }).join('');
                after.insertAdjacentElement('afterend', lab); lab.insertAdjacentElement('afterend', row); after = row;
            }
        }
    }
    function render() {
        const cur = current(), P = visPal();
        let h = `<div class="v3d-row"><span class="v3d-lab">Color by</span><span class="v3d-seg" role="group" aria-label="Color by">${MODES.map(([m, l]) => `<button type="button" data-mode="${m}" class="${cur === m ? 'on' : ''}" aria-pressed="${cur === m}">${l}</button>`).join('')}</span></div>`;
        if (cur === 'interface' && P) { const sh = shadeNow(), multi = !!mcRow();
            h += `<div class="v3d-row"><span class="v3d-lab">Palette</span><span class="v3d-seg" role="group" aria-label="Shade">${SHADES.map(([v, l, t]) => `<button type="button" data-shade="${v}" class="${sh === v ? 'on' : ''}" aria-pressed="${sh === v}" title="${t}">${l}</button>`).join('')}</span><span class="v3d-sep" aria-hidden="true"></span>`
              + pals().map((o, i) => { const s = o.querySelector('.preset-strip'), t = (o.textContent || '').trim(), c = cols(o);
                const on = multi || sh === 'gradient' ? o.classList.contains('active') : sameDark(o);
                const sc = !multi && sh !== 'gradient' && c.length === 4 ? shadeCols(sh, c[1], c[2]) : null, sp = (x) => `<span style="background:${x}${x === '#ffffff' ? ';box-shadow:inset 0 0 0 1px #C9D2DC' : ''}"></span>`;
                const strip = sc ? `<div class="preset-strip">${(sh === 'solid' ? [sc[1], sc[3]] : [sc[0], sc[1], sc[3], sc[2]]).map(sp).join('')}</div>` : s ? s.outerHTML : '';
                return `<button type="button" class="v3d-chip${on ? ' active' : ''}" data-pal="${i}" title="${esc(t)}">${strip}<span>${esc(t)}</span></button>`; }).join('')
              + `<button type="button" class="v3d-more" data-more="1" title="every palette and the custom colors, in Visualization Scripts">more ↓</button></div>`; }
        else if (cur === 'domains') h += `<div class="v3d-row"><span class="v3d-lab">Domains from</span><span class="v3d-seg" role="group" aria-label="Domains from">${['uniprot', 'ted'].map((m) => `<button type="button" data-dom="${m}" class="${viewer3dDomainMode === m ? 'on' : ''}" aria-pressed="${viewer3dDomainMode === m}">${m === 'ted' ? 'TED' : 'UniProt/Pfam'}</button>`).join('')}</span><span class="v3d-note">each domain its own color, the rest white · 3D view only</span></div>`;
        else if (cur === 'plddt') h += `<div class="v3d-row"><span class="v3d-lab">Scale</span><span class="v3d-note">AlphaFold confidence per residue: <b style="color:#0053D6">&gt;90</b> · <b style="color:#3BA6D9">70–90</b> · <b style="color:#C9A800">50–70</b> · <b style="color:#E8642E">≤50</b></span></div>`;
        else h += `<div class="v3d-row"><span class="v3d-lab">Colors</span><span class="v3d-note">one per ${cur === 'bychain' ? 'chain' : 'polymer'}</span></div>`;
        if (gap || seg || document.querySelector('.show-complete-cb')) {
            h += `<div class="v3d-row"><span class="v3d-lab">LIR display</span>`
              + (gap ? `<label class="v3d-cb" title="bridge breaks of up to this many residues, for a continuous cartoon">Fill gaps ≤ <input type="number" class="v3d-num" data-proxy="gap" min="0" max="200" value="${esc(gap.value)}" aria-label="Fill gaps up to this many residues"></label>` : '')
              + (seg ? `<span class="v3d-pipe">|</span><label class="v3d-cb" title="drop isolated LIR fragments shorter than this">Min segment ≥ <input type="number" class="v3d-num" data-proxy="seg" min="1" max="50" value="${esc(seg.value)}" aria-label="Minimum LIR segment length"></label>` : '')
              + (gap || seg ? `<button type="button" class="v3d-apply" data-apply="1">Apply</button><span class="v3d-pipe">|</span>` : '')
              + `<label class="v3d-cb" title="display all residues instead of LIR only"><input type="checkbox" class="show-complete-cb" onchange="toggleShowComplete(this)"${typeof showComplete !== 'undefined' && showComplete ? ' checked' : ''}> Complete structure</label>`
              + `<label class="v3d-cb" title="color non-interacting residues gray"><input type="checkbox" class="gray-nonlir-cb" onchange="toggleGrayNonLir(this)"${typeof grayNonLir !== 'undefined' && grayNonLir ? ' checked' : ''}> Gray non-LIR</label></div>`;
        }
        bar.innerHTML = h;
        bar.querySelectorAll('[data-mode]').forEach((b) => b.onclick = () => { const m = b.dataset.mode, was = viewer3dDomainMode;
            if (m === 'domains') { if (!viewer3dDomainMode) { viewer3dDomainMode = 'uniprot'; render(); redraw(); } return; }
            viewer3dDomainMode = '';
            if (m === 'interface') { const ps = P ? [...P.querySelectorAll('.preset-chip')] : []; const o = ps.find((x) => x.classList.contains('active')) || ps[0]; if (o) o.click(); else redraw(); }
            else { const o = modeChip(m); if (o) o.click(); }
            render(); if (was && m !== 'interface' && !modeChip(m)) redraw(); });
        bar.querySelectorAll('[data-dom]').forEach((b) => b.onclick = () => { viewer3dDomainMode = b.dataset.dom; render(); redraw(); });
        bar.querySelectorAll('[data-pal]').forEach((b) => b.onclick = () => { const o = pals()[+b.dataset.pal]; if (o) pickPal(o); render(); });
        bar.querySelectorAll('[data-shade]').forEach((b) => b.onclick = () => { setShade(b.dataset.shade); render(); });
        const more = bar.querySelector('[data-more]'); if (more) more.onclick = () => { (P || modeRow).scrollIntoView({ block: 'center' }); };   // a jump, never a smooth scroll
        const ap = bar.querySelector('[data-apply]');
        if (ap) { const go = () => { const g = bar.querySelector('[data-proxy="gap"]'), s2 = bar.querySelector('[data-proxy="seg"]'); if (g && gap) gap.value = g.value; if (s2 && seg) seg.value = s2.value; if (typeof updateGapFill === 'function') updateGapFill(); };
            ap.onclick = go; bar.querySelectorAll('[data-proxy]').forEach((x) => x.onkeydown = (e) => { if (e.key === 'Enter') go(); }); }
    }
    render();
    // the Scripts card's chips and boxes changed (a click there, a page reset): show the same state here
    const watch = [...(palRow ? palRow.querySelectorAll('.preset-chip') : []), ...(modeRow ? modeRow.querySelectorAll('.preset-chip') : []), ...document.querySelectorAll('#presets-pale .preset-chip, #presets-contacts .preset-chip')];
    let pend = 0; const later = () => { if (!pend) pend = requestAnimationFrame(() => { pend = 0; if (!bar.contains(document.activeElement) || document.activeElement.type !== 'number') render(); }); };
    for (const o of watch) new MutationObserver(later).observe(o, { attributes: true, attributeFilter: ['class'] });
    const mc = document.getElementById('multichain-presets'); if (mc) new MutationObserver(later).observe(mc, { attributes: true, childList: true, subtree: true, attributeFilter: ['class', 'style'] });
    for (const x of [gap, seg]) if (x) x.addEventListener('change', later);
}

// base: the page's color components; chains: [{ chain, label, uniprot: [...], ted: [...] }], domains as the contact maps use
// them ({ name, start, end, segments? } in the chain's residue numbers). Chains hidden by the chain toggles stay hidden.
function applyViewer3dDomains(base, structText, fmt, chains) {
    const key = document.getElementById('viewer3d-domkey');
    if (!viewer3dDomainMode || !structText) { if (key) key.hidden = true; return base; }
    const res = parseBfactorsPerResidue(structText, fmt === 'mmcif' || fmt === 'cif' ? 'cif' : 'pdb'), per = new Map();
    for (const k of res.keys()) { const i = k.lastIndexOf(':'), ch = k.slice(0, i), rn = +k.slice(i + 1); if (!per.has(ch)) per.set(ch, new Set()); per.get(ch).add(rn); }
    const drawn = base.filter(c => !c.isIon && !c.stick);
    let shown = [...new Set(drawn.map(c => c.chain))];
    if (!shown.length) shown = (chains || []).map(c => c.chain);   // a page that drew nothing yet (a monomer without regions): every chain it names
    // the residues the page draws on each chain (LIR only unless Complete structure is on): the domains recolor those, never more
    const inView = new Map();
    for (const c of drawn) { if (!inView.has(c.chain)) inView.set(c.chain, new Set());
        const s = inView.get(c.chain); if (s.all) continue;
        if (!c.ranges) { s.all = true; continue; }
        for (const g of c.ranges) for (let r = +g.start; r <= +g.end; r++) s.add(r); }
    const out = base.filter(c => c.isIon || c.stick), items = [];
    const runs = (list) => { const r = []; for (const n of list) { const last = r[r.length - 1]; if (last && n === last.end + 1) last.end = n; else r.push({ start: n, end: n }); } return r; };
    let n = 0;
    for (const ch of shown) {
        const spec = (chains || []).find(c => c.chain === ch) || { chain: ch, label: ch };
        const doms = [...((viewer3dDomainMode === 'ted' ? spec.ted : spec.uniprot) || [])].filter(d => d && Number.isFinite(+d.start)).sort((a, b) => a.start - b.start);
        const owner = new Map();
        doms.forEach((d, di) => { const col = DOM3D_COLORS[(n + di) % DOM3D_COLORS.length];
            items.push({ col, chain: spec.label || ch, name: d.cath && !String(d.name || '').includes(d.cath) ? `${d.name || 'TED'} (${d.cath})` : d.name, start: d.start, end: d.end });
            for (const sg of (d.segments && d.segments.length ? d.segments : [{ start: d.start, end: d.end }])) for (let r = +sg.start; r <= +sg.end; r++) if (!owner.has(r)) owner.set(r, col); });
        n += doms.length;
        const vis = inView.get(ch), resid = [...(per.get(ch) || [])].filter(r => !vis || vis.all || vis.has(r)).sort((a, b) => a - b), byCol = new Map();
        for (const r of resid) { const c = owner.get(r) || '#FFFFFF'; if (!byCol.has(c)) byCol.set(c, []); byCol.get(c).push(r); }
        for (const [c, list] of byCol) out.push({ chain: ch, ranges: runs(list), color: c });
    }
    if (key) {
        key.hidden = false;
        key.innerHTML = `<b>${DOM3D_NAMES[viewer3dDomainMode]}</b>` + (items.length ? items.map(d => `<span><i style="background:${d.col}"></i>${d.chain} · ${d.name} (${d.start}–${d.end})</span>`).join('')
            : `<span class="muted">none for the chains shown</span>`) + `<span><i style="background:#fff"></i>not in a domain</span><span class="muted">3D view only; the scripts keep their colors</span>`;
    }
    return out;
}
function _ensureViewer3dFindCss() {
    if (document.getElementById('v3d-find-css')) return;
    const st = document.createElement('style'); st.id = 'v3d-find-css';
    st.textContent = `.v3d-find { display:flex; flex-direction:column; gap:0.35rem; margin:0 0 0.55rem; font-size:0.82rem; color:#555; }
.v3d-find-row { display:flex; flex-wrap:wrap; align-items:center; gap:0.35rem 0.45rem; }
.v3d-find-lab { font-size:0.7rem; font-weight:600; letter-spacing:0.06em; text-transform:uppercase; color:#7A8899; min-width:5.6rem; }
.v3d-find-in { flex:1 1 240px; max-width:360px; height:28px; box-sizing:border-box; padding:0 0.55rem; border:1px solid #DCE3EB; border-radius:6px; font-family:inherit; font-size:0.8rem; color:#17263A; background:#fff; }
.v3d-find-in:focus { outline:2px solid #2471A3; outline-offset:1px; border-color:#2471A3; }
.v3d-find-btn { height:28px; padding:0 0.65rem; border-radius:6px; font-family:inherit; font-size:0.78rem; font-weight:600; cursor:pointer; border:1px solid #DCE3EB; background:#fff; color:#17263A; }
.v3d-find-btn:hover { border-color:#2471A3; color:#2471A3; }
.v3d-find-btn.pri { background:#1A5276; border-color:#1A5276; color:#fff; }
.v3d-find-btn.pri:hover { background:#2471A3; border-color:#2471A3; color:#fff; }
.v3d-find-res { display:flex; flex-wrap:wrap; align-items:center; gap:0.3rem 0.4rem; padding-left:6.05rem; font-size:0.78rem; }
.v3d-find-res:empty { display:none; }
.v3d-hit { display:inline-flex; align-items:center; gap:0.3rem; padding:0.18rem 0.5rem; border-radius:999px; border:1px solid #F3B3CB; background:#FDF2F6; color:#9C1143; font-family:inherit; font-size:0.76rem; font-weight:600; cursor:pointer; }
.v3d-hit:hover, .v3d-hit.on { background:#E91E63; border-color:#E91E63; color:#fff; }
.v3d-hit code { font-family:'IBM Plex Mono',ui-monospace,SFMono-Regular,Menlo,monospace; font-size:0.74rem; font-weight:500; opacity:0.88; }
.v3d-find-msg { color:#6B7785; font-size:0.78rem; }
.v3d-find-msg.miss { color:#B03A2E; }
@media (max-width:600px) { .v3d-find-lab { min-width:0; width:100%; } .v3d-find-res { padding-left:0; } }`;
    document.head.appendChild(st);
}
function mountViewer3dFindBar(frameId) {
    const id = frameId || 'viewer3d-frame';
    if (id !== 'viewer3d-frame') return;
    const frame = document.getElementById(id);
    if (!frame || document.getElementById('viewer3d-findbar')) return;
    _ensureViewer3dFindCss();
    const wrap = document.createElement('div');
    wrap.id = 'viewer3d-findbar';
    wrap.className = 'v3d-find';
    wrap.innerHTML = `<div class="v3d-find-row"><span class="v3d-find-lab">Find in 3D</span>`
        + `<input type="search" id="v3d-find-in" class="v3d-find-in" spellcheck="false" autocomplete="off" placeholder="Sequence (e.g. PPXY, TNCESL) or residue (A:120-135)" aria-label="Search sequence motif or residue range in 3D structure">`
        + `<button type="button" id="v3d-find-go" class="v3d-find-btn pri">Highlight</button>`
        + `<button type="button" id="v3d-find-clr" class="v3d-find-btn" hidden>Clear</button>`
        + `<button type="button" id="v3d-find-rst" class="v3d-find-btn" title="Reset 3D camera to the full structure">Reset view</button></div>`
        + `<div id="v3d-find-res" class="v3d-find-res" aria-live="polite"></div>`;
    const colorBar = document.getElementById('viewer3d-colorbar');
    const host = frame.parentElement && frame.parentElement.classList.contains('viewer3d') ? frame.parentElement : frame;
    if (colorBar && colorBar.parentElement === host.parentElement) colorBar.insertAdjacentElement('afterend', wrap);
    else if (host.parentElement) host.parentElement.insertBefore(wrap, host);
    const inp = wrap.querySelector('#v3d-find-in');
    const runSearch = () => {
        const q = (inp.value || '').trim();
        if (!q) { clearMolstarHighlight(id); return; }
        const st = _getV3dState(id);
        st.query = q;
        const seqs = parseChainSequences(st.structData, st.fmt);
        const res = searchStructure3d(seqs, q);
        st.hits = res.hits;
        st.activeHitIdx = res.hits.length === 1 ? 0 : -1;
        if (!res.hits.length) {
            st.highlightSpans = [];
            st.autoCompleteNote = '';
            if (st.baseComponents && st.baseComponents.length) applyColorsToMolstarFrame(id, st.baseComponents, st.fmt, null);
            _syncViewer3dFindUi(id, res);
            return;
        }
        const spans = res.hits.map(h => ({ chain: h.chain, start: h.start, end: h.end }));
        highlightSpansInMolstar(id, spans, spans);
        _syncViewer3dFindUi(id, res);
    };
    wrap.querySelector('#v3d-find-go').onclick = runSearch;
    wrap.querySelector('#v3d-find-clr').onclick = () => clearMolstarHighlight(id);
    wrap.querySelector('#v3d-find-rst').onclick = () => resetMolstarCamera(id);
    inp.addEventListener('keydown', (e) => {
        if (e.key === 'Enter') { e.preventDefault(); runSearch(); }
        else if (e.key === 'Escape') { e.preventDefault(); clearMolstarHighlight(id); }
    });
    inp.addEventListener('search', () => { if (!inp.value.trim()) clearMolstarHighlight(id); });
    _syncViewer3dFindUi(id);
}
function _syncViewer3dFindUi(frameId, lastRes) {
    const id = frameId || 'viewer3d-frame';
    if (id !== 'viewer3d-frame') return;
    if (!document.getElementById('viewer3d-findbar')) mountViewer3dFindBar(id);
    const resEl = document.getElementById('v3d-find-res');
    const clrBtn = document.getElementById('v3d-find-clr');
    if (!resEl) return;
    const st = _getV3dState(id);
    const hasHl = !!(st.highlightSpans && st.highlightSpans.length);
    if (clrBtn) clrBtn.hidden = !hasHl && !st.query;
    const esc = (t) => String(t).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
    if (!st.hits.length) {
        if (st.query) {
            const seqs = parseChainSequences(st.structData, st.fmt);
            const chTxt = seqs.chains.length ? ` in chain${seqs.chains.length === 1 ? '' : 's'} ${seqs.chains.join(', ')}` : ' (no structure loaded yet)';
            const msg = lastRes && lastRes.kind === 'invalid'
                ? 'Enter amino acid letters (X or . for any residue) or residue numbers (e.g. A:120-135).'
                : `No match for "${st.query}"${chTxt}.`;
            resEl.innerHTML = `<span class="v3d-find-msg miss">${esc(msg)}</span>`;
        } else {
            resEl.innerHTML = '';
        }
        return;
    }
    const pills = st.hits.slice(0, 12).map((h, i) => {
        const sub = h.seq ? (h.seq.length > 10 ? h.seq.slice(0, 8) + '\u2026' : h.seq) : '';
        return `<button type="button" class="v3d-hit${st.activeHitIdx === i ? ' on' : ''}" data-hi="${i}" title="Focus ${esc(h.label)}${h.seq ? ' (' + esc(h.seq) + ')' : ''} in 3D">${esc(h.label)}${sub ? ` <code>${esc(sub)}</code>` : ''}</button>`;
    }).join('');
    const allPill = st.hits.length > 1 ? `<button type="button" class="v3d-hit${st.activeHitIdx === -1 ? ' on' : ''}" data-hi="-1">All (${st.hits.length})</button>` : '';
    const extra = st.hits.length > 12 ? `<span class="v3d-find-msg">+${st.hits.length - 12} more</span>` : '';
    const note = st.autoCompleteNote ? `<span class="v3d-find-msg">\u00b7 ${esc(st.autoCompleteNote)}</span>` : '';
    const prefix = st.query ? `<span class="v3d-find-msg">${st.hits.length} match${st.hits.length === 1 ? '' : 'es'}:</span>` : `<span class="v3d-find-msg">Highlighted:</span>`;
    resEl.innerHTML = `${prefix}${allPill}${pills}${extra}${note}`;
    resEl.querySelectorAll('[data-hi]').forEach((b) => {
        b.onclick = () => {
            const idx = +b.dataset.hi;
            st.activeHitIdx = idx;
            const spans = st.hits.map(h => ({ chain: h.chain, start: h.start, end: h.end }));
            const focus = idx >= 0 && st.hits[idx] ? [{ chain: st.hits[idx].chain, start: st.hits[idx].start, end: st.hits[idx].end }] : spans;
            highlightSpansInMolstar(id, spans, focus, { keepFindUi: true });
            resEl.querySelectorAll('[data-hi]').forEach(x => x.classList.toggle('on', +x.dataset.hi === idx));
        };
    });
}
function mountViewer3dPngButton() {
    const frame = document.getElementById('viewer3d-frame');
    if (!frame || document.getElementById('viewer3d-png-btn')) return;
    const fsBtn = frame.parentElement && frame.parentElement.querySelector('button[onclick*="requestFullscreen"]');
    let bar = fsBtn ? fsBtn.parentElement : null;
    if (!bar || bar === frame.parentElement) {
        bar = document.createElement('div');
        bar.style.cssText = 'text-align:right; margin-top:0.3rem;';
        const after = document.getElementById('viewer3d-domkey') || frame;
        after.insertAdjacentElement('afterend', bar);
    }
    const btn = document.createElement('button');
    btn.type = 'button';
    btn.id = 'viewer3d-png-btn';
    btn.textContent = '\u2193 3D PNG';
    btn.title = 'Download a transparent-background PNG snapshot of the current 3D view';
    btn.style.cssText = 'padding:0.25rem 0.6rem; border:1px solid #2471A3; border-radius:4px; background:#fff; cursor:pointer; font-size:0.78rem; color:#2471A3; font-weight:600; margin-right:0.35rem;';
    btn.onclick = () => {
        const stem = (document.getElementById('result-title')?.textContent || document.getElementById('title')?.textContent || 'livia_3d')
            .replace(/Interaction Analysis Summary\s*[—:-]?\s*/i, '').trim() || 'livia_3d';
        downloadMolstarPng('viewer3d-frame', `${stem}_3d.png`);
    };
    bar.insertBefore(btn, bar.firstChild);
}
if (typeof document !== 'undefined') {
    const _initV3dUi = () => { mountViewer3dColorBar(); mountViewer3dFindBar(); mountViewer3dPngButton(); };
    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', _initV3dUi); else _initV3dUi();
}
