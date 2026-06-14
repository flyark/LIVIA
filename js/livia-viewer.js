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

// ── Build the MVS structureChildren (per-component representations + colors) ──
function _buildMvsStructureChildren(colorComponents) {
    const structureChildren = [];
    for (const comp of colorComponents) {
        if (comp.isIon) {
            structureChildren.push({
                kind: 'component',
                params: { selector: 'ion' },
                children: [
                    { kind: 'representation', params: { type: 'ball_and_stick' },
                      children: [{ kind: 'color', params: { color: comp.color } }] }
                ]
            });
            structureChildren.push({
                kind: 'component',
                params: { selector: 'ligand' },
                children: [
                    { kind: 'representation', params: { type: 'ball_and_stick' },
                      children: [{ kind: 'color', params: { color: comp.color } }] }
                ]
            });
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

// ── Send a color-only update to a Mol* iframe (no structure reload) ──
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
function swapMolstarStructure(frameId, structData, colorComponents, fmt) {
    const frame = document.getElementById(frameId);
    if (!frame || !frame.contentWindow) return false;
    const mvsJson = buildMvsJson(colorComponents, fmt);
    frame.contentWindow.postMessage({
        type: 'loadStructure',
        structData: structData,
        mvsStr: JSON.stringify(mvsJson),
        fmt: fmt,
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
#loading-overlay .loading-source { font-size: 11px; color: #aaa; margin-top: 8px; }
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
                        if (c3d && typeof c3d.handleResize === 'function') c3d.handleResize();
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
                        if (c3d && typeof c3d.handleResize === 'function') c3d.handleResize();
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
