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

// ── Build complete Mol* viewer HTML page with MVS coloring ──
// structData: raw structure text (PDB or mmCIF)
// fmt: 'pdb' or 'mmcif'
// colorComponents: array of { chain, ranges: [{start, end}], color: '#hex' }
function buildMolstarPage(structData, fmt, colorComponents) {
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

    const mvsJson = {
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
                            children: structureChildren
                        }]
                    }]
                }
            ]
        },
        metadata: { version: '1.6' }
    };

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

    var structBlob = new Blob([structData], { type: 'text/plain' });
    var structUrl = URL.createObjectURL(structBlob);
    var mvsStr = mvsTemplate.replace('__STRUCT_BLOB_URL__', structUrl);

    try {
        viewer.loadMvsData(mvsStr, 'mvsj');
    } catch(e) {
        console.warn('MVS failed, falling back:', e);
        viewer.loadStructureFromData(structData, fmt);
    }

    // Apply illustrative style AFTER loading
    setTimeout(function() {
        try {
            if (viewer.plugin.managers && viewer.plugin.managers.canvas3dContext) {
                viewer.plugin.managers.canvas3dContext.setProps({ style: { name: 'illustrative' } });
            }
        } catch(e1) {}
        try {
            var c3d = viewer.plugin.canvas3d;
            if (c3d) {
                var rp = Object.assign({}, c3d.props.renderer);
                rp.style = { name: 'flat-shaded', params: {} };
                c3d.setProps({
                    renderer: rp,
                    postprocessing: {
                        occlusion: { name: 'on', params: { multiScale: { name: 'off', params: {} }, radius: 5, bias: 0.8, blurKernelSize: 15, resolutionScale: 1, color: 0x000000 } },
                        outline: { name: 'on', params: { scale: 1, threshold: 0.33, color: 0x000000, includeTransparent: true } },
                    },
                });

            }
        } catch(e2) { console.warn('Illustrative style error:', e2); }
        // Default screenshot to transparent background
        try {
            var sh = viewer.plugin.helpers.viewportScreenshot;
            if (sh) {
                var cur = sh.behaviors.values.value;
                sh.behaviors.values.next(Object.assign({}, cur, { transparent: true }));
            }
        } catch(e3) {}
    }, 2000);
}
init().catch(function(e) { console.error('Mol* error:', e); });
<\/script>
</body></html>`;
}
