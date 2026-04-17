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

// ── Build complete Mol* viewer HTML page with MVS coloring ──
// structData: raw structure text (PDB or mmCIF)
// fmt: 'pdb' or 'mmcif'
// colorComponents: array of { chain, ranges: [{start, end}], color: '#hex' }
function buildMolstarPage(structData, fmt, colorComponents) {
    const mvsJson = buildMvsJson(colorComponents, fmt);

    return `<!DOCTYPE html>
<html><head>
<script src="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.js"><\/script>
<link rel="stylesheet" type="text/css" href="https://cdn.jsdelivr.net/npm/molstar@latest/build/viewer/molstar.css" />
<style>
#viewer1 { position:absolute; top:0; left:0; right:0; bottom:0; }
@media (max-width: 768px) {
  .msp-layout-right { display: none !important; }
  .msp-viewport-controls { display: none !important; }
}
</style>
</head><body>
<div id="viewer1"></div>
<script>
var structData = ${JSON.stringify(structData)};
var fmt = "${fmt}";
var mvsTemplate = ${JSON.stringify(JSON.stringify(mvsJson))};
var _viewer = null;
var _structUrl = null;
var _ready = false;
var _pendingColorMvs = null;

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
        var c3d = plugin.canvas3d;
        if (c3d) {
            var pp = c3d.props.postprocessing || {};
            c3d.setProps({
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

// Listen for parent postMessage color updates
window.addEventListener('message', function(ev) {
    if (!ev.data || ev.data.type !== 'updateColors' || !ev.data.mvsStr) return;
    if (!_ready || !_viewer || !_structUrl) {
        _pendingColorMvs = ev.data.mvsStr;
        return;
    }
    try {
        var newMvs = ev.data.mvsStr.replace('__STRUCT_BLOB_URL__', _structUrl);
        var p = _viewer.loadMvsData(newMvs, 'mvsj');
        // Re-apply illustrative style after MVS reload (loadMvsData resets postprocessing)
        var reapply = function() { _applyIllustrativeStyle(_viewer); };
        if (p && typeof p.then === 'function') p.then(reapply, reapply);
        else setTimeout(reapply, 50);
    } catch (e) { console.warn('Color update failed:', e); }
});

function _notifyReady() {
    _ready = true;
    try { parent.postMessage({ type: 'molstarReady' }, '*'); } catch(e) {}
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
    var isMobile = window.innerWidth < 768;
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

    try {
        viewer.loadMvsData(mvsStr, 'mvsj');
    } catch(e) {
        console.warn('MVS failed, falling back:', e);
        viewer.loadStructureFromData(structData, fmt);
    }

    // Apply illustrative style AFTER loading
    setTimeout(function() {
        _applyIllustrativeStyle(viewer);
        // Default screenshot to transparent background
        try {
            var sh = viewer.plugin.helpers.viewportScreenshot;
            if (sh) {
                var cur = sh.behaviors.values.value;
                sh.behaviors.values.next(Object.assign({}, cur, { transparent: true }));
            }
        } catch(e3) {}
        _notifyReady();
    }, 2000);
}
init().catch(function(e) { console.error('Mol* error:', e); });
<\/script>
</body></html>`;
}
