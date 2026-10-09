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
function parseBfactorsPerResidue(text, format) {
    const m = new Map();
    // Auto-detect format if not specified
    if (!format) format = text.includes('_atom_site.') ? 'cif' : 'pdb';
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
        let inA = false; const cols = [];
        for (const line of lines) {
            if (line.startsWith('_atom_site.')) { inA = true; cols.push(line.trim().split('.')[1]); continue; }
            if (inA && !line.startsWith('_atom_site.') && !line.startsWith('#') && line.trim()) {
                if (line.startsWith('loop_') || line.startsWith('_')) { inA = false; continue; }
                const p = line.trim().split(/\s+/);
                if (p.length < cols.length) continue;
                const g = (n) => { const i = cols.indexOf(n); return i >= 0 ? p[i] : ''; };
                const group = g('group_PDB');
                const atom = g('label_atom_id');
                const ch = g('label_asym_id');
                const bf = parseFloat(g('B_iso_or_equiv'));
                if (isNaN(bf)) continue;
                if (group === 'ATOM' && atom === 'CA') {
                    const rn = parseInt(g('label_seq_id'));
                    if (!isNaN(rn)) m.set(`${ch}:${rn}`, bf);
                } else if (group === 'HETATM') {
                    const rn = parseInt(g('label_seq_id'));
                    m.set(`${ch}:${isNaN(rn) ? 1 : rn}`, bf);
                }
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
function applyColorsToMolstarFrame(frameId, colorComponents, fmt) {
    const frame = document.getElementById(frameId);
    if (!frame || !frame.contentWindow) return false;
    const mvsJson = buildMvsJson(colorComponents, fmt);
    frame.contentWindow.postMessage({
        type: 'updateColors',
        mvsStr: JSON.stringify(mvsJson),
    }, '*');
    return true;
}

// ── Swap the structure inside an already-booted Mol* iframe (warm-up pattern) ──
// The iframe must have been built with buildMolstarPage. Use this to replace the
// placeholder structure loaded at page-entry warmup with the real prediction structure
// once analysis finishes, without rebuilding the iframe (avoids re-downloading Mol* JS).
function swapMolstarStructure(frameId, structData, colorComponents, fmt, label) {   // label: a readable name for Mol*'s structure panel, in place of the blob URL
    const frame = document.getElementById(frameId);
    if (!frame || !frame.contentWindow) return false;
    const mvsJson = buildMvsJson(colorComponents, fmt);
    frame.contentWindow.postMessage({
        type: 'loadStructure',
        structData: structData,
        mvsStr: JSON.stringify(mvsJson),
        fmt: fmt,
        label: label || '',
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
const _molstarWarmState = new Map(); // frameId → { warmedUp, blobUrl }
function warmupMolstarFrame(frameId, parentBaseUrl) {
    const frame = document.getElementById(frameId);
    if (!frame) return false;
    const prev = _molstarWarmState.get(frameId);
    if (prev && prev.warmedUp) return true;
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
function isMolstarFrameWarmedUp(frameId) {
    const s = _molstarWarmState.get(frameId);
    return !!(s && s.warmedUp);
}

// ── Build complete Mol* viewer HTML page with MVS coloring ──
// structData: raw structure text (PDB or mmCIF)
// fmt: 'pdb' or 'mmcif'
// colorComponents: array of { chain, ranges: [{start, end}], color: '#hex' }
// parentBaseUrl: absolute URL of the parent LIVIA page (used for self-hosted Mol* fallback)
const MOLSTAR_VERSION = '5.9.0';
function buildMolstarPage(structData, fmt, colorComponents, parentBaseUrl) {
    const mvsJson = buildMvsJson(colorComponents, fmt);
    const selfHostBase = parentBaseUrl || '';

    return `<!DOCTYPE html>
<html><head>
<style>
#viewer1 { position:absolute; top:0; left:0; right:0; bottom:0; }
/* Mol*'s state-snapshot picker ("[1/1] <load time>" + play): every structure loads as one MVS snapshot, so the picker
   only ever offers that one entry, stamped with the load time. Trajectory controls in the same corner stay. */
.msp-state-snapshot-viewport-controls { display: none !important; }
@media (max-width: 768px) {
  .msp-layout-right { display: none !important; }
  .msp-viewport-controls { display: none !important; }
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
var _structLabel = '', _relabelT = 0, _lastWH = '';
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

// Listen for parent postMessage: structure swap (warm-up pattern) + color update
window.addEventListener('message', function(ev) {
    if (!ev.data) return;
    if (ev.data.type === 'loadStructure' && ev.data.structData && ev.data.mvsStr) {
        if (!_ready || !_viewer) { _pendingStructure = ev.data; return; }
        _structLabel = ev.data.label || '';
        try {
            if (_structUrl) { try { URL.revokeObjectURL(_structUrl); } catch(_) {} }
            var newBlob = new Blob([ev.data.structData], { type: 'text/plain' });
            _structUrl = URL.createObjectURL(newBlob);
            var newMvs = ev.data.mvsStr.replace('__STRUCT_BLOB_URL__', _structUrl);
            // Structure changed → must explicitly refit camera, otherwise Mol* keeps the
            // canvas3d state from the previous load (e.g. the tight zoom around the
            // warm-up placeholder atom at the origin) and the new structure ends up
            // mostly off-screen.
            var p = _viewer.loadMvsData(newMvs, 'mvsj');
            var refit = function() {
                _applyIllustrativeStyle(_viewer);
                var resetCamera = function() {
                    try {
                        var c3d = _viewer && _viewer.plugin && _viewer.plugin.canvas3d;
                        // Force Mol* to recompute viewport size: when the iframe just
                        // transitioned from display:none → visible (warm-up → real
                        // results), ResizeObserver may not have fired yet and the
                        // camera would fit using stale 0x0 dimensions (looks broken
                        // on mobile in particular).
                        var wh = window.innerWidth + 'x' + window.innerHeight;   // the retries below resize only when the frame's size moved (each resize can reallocate the render buffers)
                        if (c3d && typeof c3d.handleResize === 'function' && wh !== _lastWH) { c3d.handleResize(); _lastWH = wh; }
                        if (c3d && typeof c3d.requestCameraReset === 'function') c3d.requestCameraReset();
                        else if (_viewer && _viewer.plugin && _viewer.plugin.managers && _viewer.plugin.managers.camera && typeof _viewer.plugin.managers.camera.reset === 'function') _viewer.plugin.managers.camera.reset();
                    } catch(_e) {}
                };
                // Retry across a few frames — Mol*'s structure ingestion may finish
                // a tick after loadMvsData's promise resolves.
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
        _pendingColorMvs = ev.data.mvsStr;
        return;
    }
    try {
        var newMvs = ev.data.mvsStr.replace('__STRUCT_BLOB_URL__', _structUrl);
        newMvs = _injectCurrentCameraIntoMvs(newMvs);
        // Snapshot the FULL camera state (incl. radius/radiusMax for ortho zoom) before reload.
        // MVS only carries target/position/up, so we need this to keep zoom stable.
        var savedCamera = null;
        try {
            var c3d = _viewer.plugin && _viewer.plugin.canvas3d;
            if (c3d && c3d.camera && c3d.camera.state) {
                // Clone so Mol* can't mutate it after we save.
                savedCamera = JSON.parse(JSON.stringify(c3d.camera.state));
            }
        } catch(_e) {}
        var restore = function() {
            try {
                var c3dr = _viewer && _viewer.plugin && _viewer.plugin.canvas3d;
                if (savedCamera && c3dr && c3dr.camera) c3dr.camera.setState(savedCamera, 0);
            } catch(_e2) {}
        };
        var p = _viewer.loadMvsData(newMvs, 'mvsj');
        var reapply = function() {
            _applyIllustrativeStyle(_viewer);
            // Restore at multiple ticks to win against Mol*'s post-load auto-fit which
            // can run asynchronously after loadMvsData's promise resolves.
            restore();
            requestAnimationFrame(function() { restore(); requestAnimationFrame(restore); });
            setTimeout(restore, 100);
            setTimeout(restore, 300);
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
        // Process queued structure swap (warm-up → real prediction transition)
        try {
            if (_structUrl) { try { URL.revokeObjectURL(_structUrl); } catch(_) {} }
            var psBlob = new Blob([_pendingStructure.structData], { type: 'text/plain' });
            _structUrl = URL.createObjectURL(psBlob);
            var psMvs = _pendingStructure.mvsStr.replace('__STRUCT_BLOB_URL__', _structUrl);
            var psP = _viewer.loadMvsData(psMvs, 'mvsj');
            var psRefit = function() {
                _applyIllustrativeStyle(_viewer);
                var resetCamera = function() {
                    try {
                        var c3d = _viewer && _viewer.plugin && _viewer.plugin.canvas3d;
                        // Force Mol* to recompute viewport size: when the iframe just
                        // transitioned from display:none → visible (warm-up → real
                        // results), ResizeObserver may not have fired yet and the
                        // camera would fit using stale 0x0 dimensions (looks broken
                        // on mobile in particular).
                        var wh = window.innerWidth + 'x' + window.innerHeight;   // the retries below resize only when the frame's size moved (each resize can reallocate the render buffers)
                        if (c3d && typeof c3d.handleResize === 'function' && wh !== _lastWH) { c3d.handleResize(); _lastWH = wh; }
                        if (c3d && typeof c3d.requestCameraReset === 'function') c3d.requestCameraReset();
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
        _pendingColorMvs = null; // colors are baked into MVS that came with the structure
        return;
    }
    if (_pendingColorMvs && _viewer && _structUrl) {
        try {
            var newMvs = _pendingColorMvs.replace('__STRUCT_BLOB_URL__', _structUrl);
            var p = _viewer.loadMvsData(newMvs, 'mvsj');
            var reapply = function() { _applyIllustrativeStyle(_viewer); };
            if (p && typeof p.then === 'function') p.then(reapply, reapply);
            else setTimeout(reapply, 50);
        } catch (e) { console.warn('Pending color update failed:', e); }
        _pendingColorMvs = null;
    }
}

async function init() {
    // Mobile detection: the iframe boots inside a display:none parent during warm-up,
    // which collapses window.innerWidth to 0. Fall back to the parent window's width
    // (same-origin blob iframe → accessible) so the original desktop-vs-mobile choice
    // is restored. CSS @media @max-width:768px also hides the panel visually as a
    // belt-and-suspenders measure.
    // When the page is opened via file:// the browser treats the parent as a null
    // origin and blocks parent.innerWidth access; if iframe's own width is also 0
    // (warm-up phase), default to desktop so the Structure Tools panel still shows.
    var pw = 0;
    try { pw = (parent && parent !== window && parent.innerWidth) || 0; } catch(_e) {}
    var iw = window.innerWidth;
    var isMobile = pw > 0 ? pw < 768 : (iw > 0 ? iw < 768 : false);
    var viewer = await molstar.Viewer.create('viewer1', {
        layoutIsExpanded: false,
        layoutShowControls: !isMobile,
        layoutShowRemoteState: false,
        layoutShowSequence: false,
        layoutShowLog: false,
        layoutShowLeftPanel: false,
        viewportShowExpand: false,
        viewportShowSelectionMode: false,
        viewportShowAnimation: false,
    });
    _viewer = viewer;
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
@media (max-width:600px) { .v3d-lab { min-width:0; width:100%; } }`;
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
    const isSolid = () => { const mc = mcRow(); if (mc) { const a = mc.querySelector('.mode-chip.active'); return !!a && /solid/i.test(a.textContent); }
        return !!val('color-clir-a') && val('color-lir-a') === val('color-clir-a') && val('color-lir-b') === val('color-clir-b'); };
    const sameDark = (o) => { const c = cols(o); return c.length === 4 && c[1] === val('color-clir-a') && c[2] === val('color-clir-b'); };
    const setTwo = (la, ca, lb, cb) => {
        const twin = [...document.querySelectorAll('.presets-row .preset-chip')].find((x) => { const c = cols(x); return c.length === 4 && c[0] === la && c[1] === ca && c[2] === cb && c[3] === lb; });
        if (twin) { twin.click(); return; }
        document.querySelectorAll('.preset-chip.active').forEach((x) => x.classList.remove('active'));
        if (typeof applyPreset === 'function') applyPreset(la, ca, lb, cb, null);
    };
    const lighten = (h) => '#' + [1, 3, 5].map((i) => Math.round(parseInt(h.slice(i, i + 2), 16) * 0.45 + 255 * 0.55).toString(16).padStart(2, '0')).join('');
    const pickPal = (o) => { viewer3dDomainMode = '';
        if (mcRow() || !isSolid()) { o.click(); return; }
        const c = cols(o); if (c.length === 4) setTwo(c[1], c[1], c[2], c[2]); else o.click(); };
    const setShade = (solid) => { const mc = mcRow();
        if (mc) { const o = [...mc.querySelectorAll('.mode-chip')].find((x) => /solid/i.test(x.textContent) === solid); if (o) o.click(); return; }
        const ca = val('color-clir-a'), cb = val('color-clir-b'); if (!ca || !cb) return;
        if (solid) { setTwo(ca, ca, cb, cb); return; }
        const grad = [...document.querySelectorAll('.presets-row .preset-chip')].find((x) => { const c = cols(x); return c.length === 4 && c[0] !== c[1] && c[1] === ca && c[2] === cb; });
        if (grad) grad.click(); else setTwo(lighten(ca), ca, lighten(cb), cb); };
    function render() {
        const cur = current(), P = visPal();
        let h = `<div class="v3d-row"><span class="v3d-lab">Color by</span><span class="v3d-seg" role="group" aria-label="Color by">${MODES.map(([m, l]) => `<button type="button" data-mode="${m}" class="${cur === m ? 'on' : ''}" aria-pressed="${cur === m}">${l}</button>`).join('')}</span></div>`;
        if (cur === 'interface' && P) { const solid = isSolid(), multi = !!mcRow();
            h += `<div class="v3d-row"><span class="v3d-lab">Palette</span><span class="v3d-seg" role="group" aria-label="Shade">${[[false, 'Gradient', 'LIR light, cLIR dark'], [true, 'Solid', 'one color per chain']].map(([v, l, t]) => `<button type="button" data-shade="${v ? 1 : 0}" class="${solid === v ? 'on' : ''}" aria-pressed="${solid === v}" title="${t}">${l}</button>`).join('')}</span><span class="v3d-sep" aria-hidden="true"></span>`
              + pals().map((o, i) => { const s = o.querySelector('.preset-strip'), t = (o.textContent || '').trim(), c = cols(o);
                const on = multi || !solid ? o.classList.contains('active') : sameDark(o);
                const strip = !multi && solid && c.length === 4 ? `<div class="preset-strip"><span style="background:${c[1]}"></span><span style="background:${c[2]}"></span></div>` : s ? s.outerHTML : '';
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
        bar.querySelectorAll('[data-shade]').forEach((b) => b.onclick = () => { setShade(b.dataset.shade === '1'); render(); });
        const more = bar.querySelector('[data-more]'); if (more) more.onclick = () => { (P || modeRow).scrollIntoView({ block: 'center' }); };   // a jump, never a smooth scroll
        const ap = bar.querySelector('[data-apply]');
        if (ap) { const go = () => { const g = bar.querySelector('[data-proxy="gap"]'), s2 = bar.querySelector('[data-proxy="seg"]'); if (g && gap) gap.value = g.value; if (s2 && seg) seg.value = s2.value; if (typeof updateGapFill === 'function') updateGapFill(); };
            ap.onclick = go; bar.querySelectorAll('[data-proxy]').forEach((x) => x.onkeydown = (e) => { if (e.key === 'Enter') go(); }); }
    }
    render();
    // the Scripts card's chips and boxes changed (a click there, a page reset): show the same state here
    const watch = [...(palRow ? palRow.querySelectorAll('.preset-chip') : []), ...(modeRow ? modeRow.querySelectorAll('.preset-chip') : [])];
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
    let shown = [...new Set(base.filter(c => !c.isIon && !c.stick).map(c => c.chain))];
    if (!shown.length) shown = (chains || []).map(c => c.chain);   // a page that drew nothing yet (a monomer without regions): every chain it names
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
        const resid = [...(per.get(ch) || [])].sort((a, b) => a - b), byCol = new Map();
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
if (typeof document !== 'undefined') {
    if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', mountViewer3dColorBar); else mountViewer3dColorBar();
}
