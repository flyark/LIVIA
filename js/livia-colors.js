/**
 * livia-colors.js — Color picker and preset management for LIVIA tool pages
 *
 * Provides:
 *   cxcColorMode           — current ChimeraX color command preset ('' | 'plddt' | 'bychain' | 'bypolymer')
 *   showComplete           — display flag: true = show all residues, false = show LIR only
 *   onColorChange          — callback: page sets this to its reload function
 *   initColorPickers()     — bind color input sync (call after DOM ready)
 *   updateColorStrip()     — sync color-strip swatches with picker values
 *   applyPreset()          — set 4 LIR/cLIR colors and trigger reload
 *   applyCxcPreset()       — set CXC color mode (plddt, bychain, etc.)
 *   toggleShowComplete()   — toggle complete structure display
 *   swapColors()           — swap A/B color pairs
 *
 * Dependencies: none (self-contained)
 *
 * Pages should set onColorChange to their reload function:
 *   onColorChange = () => { if (parsed && currentRank) loadRank(currentRank); };
 */

// ── Color mode state ──
let cxcColorMode = '';

// ── Display mode: show all residues vs LIR only ──
let showComplete = false;

// ── Gray non-LIR: color non-LIR residues gray when showComplete is on ──
let grayNonLir = false;

// ── Callback for page-specific reload after color change ──
let onColorChange = null;

// ── Validate a hex color string (returns normalized "#rrggbb" or null) ──
function _normalizeHex(s) {
    if (!s) return null;
    let h = String(s).trim();
    if (h[0] !== '#') h = '#' + h;
    if (/^#[0-9a-fA-F]{6}$/.test(h)) return h.toLowerCase();
    if (/^#[0-9a-fA-F]{3}$/.test(h)) {
        // expand #rgb → #rrggbb
        return ('#' + h[1] + h[1] + h[2] + h[2] + h[3] + h[3]).toLowerCase();
    }
    return null;
}

// ── A picker's hex field: a text input once initColorPickers() has upgraded it (an input shows .value, not its
// text), the original <div> before that. Presets and Swap A↔B write through here so the field follows the swatch. ──
function setColorHex(id, hex) {
    const el = document.getElementById(id);
    if (!el) return;
    if (el.tagName === 'INPUT') el.value = String(hex).toLowerCase(); else el.textContent = hex;
}

// ── Bind color picker input events (call once after DOM ready) ──
// Converts read-only .color-hex divs into editable text inputs and wires
// bidirectional sync with the paired color picker. Triggers onColorChange()
// when the user commits a new color (picker close or hex input blur/Enter).
function initColorPickers() {
    document.querySelectorAll('input[type="color"]').forEach(input => {
        const hexId = input.id.replace('color-', 'hex-');
        let hexEl = document.getElementById(hexId);

        // Upgrade <div class="color-hex"> → <input type="text" class="color-hex">
        if (hexEl && hexEl.tagName !== 'INPUT') {
            const newInput = document.createElement('input');
            newInput.type = 'text';
            newInput.id = hexId;
            newInput.className = hexEl.className;
            newInput.value = input.value;
            newInput.spellcheck = false;
            newInput.maxLength = 7;
            newInput.setAttribute('aria-label', 'Hex color code');
            hexEl.parentNode.replaceChild(newInput, hexEl);
            hexEl = newInput;
        } else if (hexEl) {
            hexEl.value = input.value;
        }

        // Picker live drag → update hex display + strip (no onColorChange yet)
        input.addEventListener('input', () => {
            if (hexEl) hexEl.value = input.value;
            updateColorStrip();
        });
        // Picker commit (close) → trigger onColorChange so 3D viewer + scripts update
        input.addEventListener('change', () => {
            if (hexEl) hexEl.value = input.value;
            updateColorStrip();
            if (onColorChange) onColorChange();
        });

        if (hexEl) {
            // User typing a hex → live update strip; commit on blur/Enter
            hexEl.addEventListener('input', () => {
                const norm = _normalizeHex(hexEl.value);
                if (norm) {
                    input.value = norm;
                    updateColorStrip();
                }
            });
            const commit = () => {
                const norm = _normalizeHex(hexEl.value);
                if (norm) {
                    hexEl.value = norm;
                    input.value = norm;
                    updateColorStrip();
                    if (onColorChange) onColorChange();
                } else {
                    // invalid input → revert to picker's current value
                    hexEl.value = input.value;
                }
            };
            hexEl.addEventListener('change', commit);
            hexEl.addEventListener('blur', commit);
            hexEl.addEventListener('keydown', (e) => {
                if (e.key === 'Enter') { e.preventDefault(); hexEl.blur(); }
            });
        }
    });
}

// ── Sync the color-strip swatch bar with current picker values ──
function updateColorStrip() {
    const strip = document.getElementById('color-strip');
    if (!strip) return;
    const s = strip.children;
    s[0].style.background = document.getElementById('color-lir-a').value;
    s[1].style.background = document.getElementById('color-clir-a').value;
    s[2].style.background = document.getElementById('color-clir-b').value;
    s[3].style.background = document.getElementById('color-lir-b').value;
}

// ── Apply a 4-color preset (LIR A, cLIR A, LIR B, cLIR B) ──
function applyPreset(lirA, clirA, lirB, clirB, el) {
    cxcColorMode = '';
    document.getElementById('color-lir-a').value = lirA;
    document.getElementById('color-clir-a').value = clirA;
    document.getElementById('color-lir-b').value = lirB;
    document.getElementById('color-clir-b').value = clirB;
    setColorHex('hex-lir-a', lirA);
    setColorHex('hex-clir-a', clirA);
    setColorHex('hex-lir-b', lirB);
    setColorHex('hex-clir-b', clirB);
    updateColorStrip();
    // Update active state
    if (el) {
        document.querySelectorAll('.preset-chip, .preset').forEach(p => p.classList.remove('active'));
        el.classList.add('active');
    }
    // Regenerate with new colors via page-specific callback
    if (onColorChange) onColorChange();
}

// ── Apply a ChimeraX color command preset (plddt, bychain, bypolymer) ──
function applyCxcPreset(mode, el) {
    cxcColorMode = mode;
    document.querySelectorAll('.preset-chip').forEach(p => p.classList.remove('active'));
    if (el) el.classList.add('active');
    if (onColorChange) onColorChange();
}

// ── Toggle complete structure display ──
function toggleShowComplete(cb) {
    showComplete = cb.checked;
    // Sync all checkboxes with same class on the page
    document.querySelectorAll('.show-complete-cb').forEach(el => { el.checked = showComplete; });
    // If turning off showComplete, also turn off grayNonLir
    if (!showComplete && grayNonLir) {
        grayNonLir = false;
        document.querySelectorAll('.gray-nonlir-cb').forEach(el => { el.checked = false; });
    }
    if (onColorChange) onColorChange();
}

// ── Toggle gray non-LIR display ──
function toggleGrayNonLir(cb) {
    grayNonLir = cb.checked;
    document.querySelectorAll('.gray-nonlir-cb').forEach(el => { el.checked = grayNonLir; });
    // Gray non-LIR implies showComplete
    if (grayNonLir && !showComplete) {
        showComplete = true;
        document.querySelectorAll('.show-complete-cb').forEach(el => { el.checked = true; });
    }
    if (onColorChange) onColorChange();
}

// ── Swap A and B color pairs ──
function swapColors() {
    const la = document.getElementById('color-lir-a').value;
    const ca = document.getElementById('color-clir-a').value;
    const lb = document.getElementById('color-lir-b').value;
    const cb = document.getElementById('color-clir-b').value;
    applyPreset(lb, cb, la, ca, null);
}

// Lighten a color toward white by amt (0..1) → the LIR (light) shade paired with a cLIR (dark) color.
function _lightenHex(hex, amt) {
    const h = _normalizeHex(hex); if (!h) return hex;
    const m = (c) => Math.round(c + (255 - c) * amt).toString(16).padStart(2, '0');
    return '#' + m(parseInt(h.slice(1, 3), 16)) + m(parseInt(h.slice(3, 5), 16)) + m(parseInt(h.slice(5, 7), 16));
}
// Resolve any CSS color (#hex, #rgb, or a name like "teal") to #rrggbb, or null if invalid.
function _colorToHex(c) {
    const norm = _normalizeHex(c); if (norm) return norm;
    try {
        const ctx = (_colorToHex._ctx || (_colorToHex._ctx = document.createElement('canvas').getContext('2d')));
        ctx.fillStyle = '#000'; ctx.fillStyle = String(c); const r = ctx.fillStyle;
        return /^#[0-9a-f]{6}$/i.test(r) ? r.toLowerCase() : null;
    } catch (e) { return null; }
}

// CSV/paste color upload for a 2-chain (A/B) figure. Accepts "chain,color" (also 3-col "chain,name,color";
// the name is ignored here). Sets each chain's cLIR (main) color and derives the light LIR shade, then
// applies through applyPreset() (which re-renders on every page). Match by chain letter A/B (or 1/2), or an
// optional gene name via opts.geneOf → { A:'sym1', B:'sym2' }.
function attachChainColorUpload(container, opts) {
    opts = opts || {};
    if (!container || !window.LiviaMaps || container.dataset.ccu) return;
    container.dataset.ccu = '1';
    const geneOf = opts.geneOf || (() => ({}));
    const val = (id) => (document.getElementById(id) || {}).value;
    window.LiviaMaps.attachColorUpload(container, {
        label: 'chain colors', keyHeader: 'chain',
        placeholder: 'chain,color\nA,#00897B\nB,#E64A19',
        currentRows: () => {
            const g = geneOf();
            return [
                { key: g.A || 'A', color: val('color-clir-a') || '#00897B' },
                { key: g.B || 'B', color: val('color-clir-b') || '#E64A19' },
            ];
        },
        apply: (rows) => {
            const g = geneOf();
            let lirA = val('color-lir-a'), clirA = val('color-clir-a'), lirB = val('color-lir-b'), clirB = val('color-clir-b');
            let applied = 0;
            for (const { key, color } of rows) {
                const hex = _colorToHex(color); if (!hex) continue;
                const k = String(key).toLowerCase();
                if (k === 'a' || k === '1' || k === String(g.A || '').toLowerCase()) { clirA = hex; lirA = _lightenHex(hex, 0.55); applied++; }
                else if (k === 'b' || k === '2' || k === String(g.B || '').toLowerCase()) { clirB = hex; lirB = _lightenHex(hex, 0.55); applied++; }
            }
            if (applied) applyPreset(lirA, clirA, lirB, clirB, null);
            return applied;
        },
    });
}

// ── Map color scales: PAE, LIS, cLIS ──
// Every PAE / LIS / cLIS map LIVIA draws from numeric PAE (universal after a full upload, dimer, monomer), with
// its color bar and exports, takes its colors from mapColorFn(kind). FlyPredictome data and lightweight bundles
// carry these maps only as images, which cannot be recolored, so those views never show mapScaleControl().
// A scale runs from the low end of the map's range to the high end — PAE low (confident) → high (uncertain);
// LIS / cLIS 0 → 1 (confident) — and each map's choice is remembered per browser.
function _mapHex6(h) {   // '#abc' / 'ABC' / '#aabbcc' → '#aabbcc'
    h = String(h || '').replace('#', '').toLowerCase();
    return '#' + (h.length === 3 ? h.split('').map(c => c + c).join('') : h.padEnd(6, '0').slice(0, 6));
}
function _mapRamp(hexStops) {   // piecewise-linear across 2+ color stops, low → high
    const stops = hexStops.map(h => { const x = _mapHex6(h); return [parseInt(x.slice(1, 3), 16), parseInt(x.slice(3, 5), 16), parseInt(x.slice(5, 7), 16)]; });
    const seg = stops.length - 1;
    return (value, vmin, vmax) => {
        const t = Math.max(0, Math.min(1, (value - vmin) / (vmax - vmin)));
        const f = t * seg, i = Math.min(seg - 1, Math.floor(f)), u = f - i, a = stops[i], b = stops[i + 1];
        return [Math.round(a[0] + (b[0] - a[0]) * u), Math.round(a[1] + (b[1] - a[1]) * u), Math.round(a[2] + (b[2] - a[2]) * u)];
    };
}
// The defaults are the pages' original bwrColor / bluesColor / greensColor, computed exactly as before so
// default maps stay pixel-identical.
function _mapBwr(value, vmin, vmax) {
    const t = Math.max(0, Math.min(1, (value - vmin) / (vmax - vmin)));
    let r, g, b;
    if (t < 0.5) { const s = t / 0.5; r = Math.round(s * 255); g = Math.round(s * 255); b = 255; }
    else { const s = (t - 0.5) / 0.5; r = 255; g = Math.round((1 - s) * 255); b = Math.round((1 - s) * 255); }
    return [r, g, b];
}
function _mapBlues(value, vmin, vmax) {   // matplotlib Blues: (247,251,255) → (107,174,214) → (8,48,107)
    const t = Math.max(0, Math.min(1, (value - vmin) / (vmax - vmin)));
    if (t < 0.5) { const s = t / 0.5; return [Math.round(247 - s * 140), Math.round(251 - s * 77), Math.round(255 - s * 41)]; }
    const s = (t - 0.5) / 0.5; return [Math.round(107 - s * 99), Math.round(174 - s * 126), Math.round(214 - s * 107)];
}
function _mapGreens(value, vmin, vmax) {   // matplotlib Greens: (247,252,245) → (116,196,118) → (0,68,27)
    const t = Math.max(0, Math.min(1, (value - vmin) / (vmax - vmin)));
    if (t < 0.5) { const s = t / 0.5; return [Math.round(247 - s * 131), Math.round(252 - s * 56), Math.round(245 - s * 127)]; }
    const s = (t - 0.5) / 0.5; return [Math.round(116 - s * 116), Math.round(196 - s * 128), Math.round(118 - s * 91)];
}
const _MAP_SEQ = {   // LIS / cLIS ramps, 0 (white) → 1 (dark); ColorBrewer anchors
    purples: ['#fcfbfd', '#9e9ac8', '#3f007d'],
    oranges: ['#fff5eb', '#fd8d3c', '#7f2704'],
    reds: ['#fff5f0', '#fb6a4a', '#67000d'],
    greys: ['#ffffff', '#969696', '#000000'],
    viridis: ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725'],
};
const MAP_SCALES = {
    pae: { label: 'PAE', low: 'low PAE (confident)', high: 'high PAE (uncertain)', key: 'livia.paeScale', def: 'bwr', custom: ['#00897b', '#ffffff', '#ff0000'],
        oldDefaultCustom: ['#00897b', '#ffffff'],   // the 2-color default before 3 became the default — a saved copy of it was never customized
        options: [['bwr', 'Blue–white–red (default)'], ['alphafold', 'AlphaFold (green)'], ['blues', 'Blues'], ['viridis', 'Viridis'], ['greys', 'Grays'], ['custom', 'Custom']],
        fns: { bwr: _mapBwr },
        stops: {
            alphafold: ['#00441b', '#74c476', '#f7fcf5'],      // dark green = confident, as in AlphaFold DB / AlphaFold 3 PAE plots
            blues: ['#08306b', '#6baed6', '#f7fbff'],
            viridis: ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725'],
            greys: ['#000000', '#ffffff'],
        } },
    lis: { label: 'LIS', low: 'LIS 0', high: 'LIS 1 (confident)', key: 'livia.lisScale', def: 'blues', custom: ['#ffffff', '#2471a3'],
        options: [['blues', 'Blues (default)'], ['greens', 'Greens'], ['purples', 'Purples'], ['oranges', 'Oranges'], ['reds', 'Reds'], ['greys', 'Grays'], ['viridis', 'Viridis'], ['custom', 'Custom']],
        fns: { blues: _mapBlues, greens: _mapGreens }, stops: _MAP_SEQ },
    clis: { label: 'cLIS', low: 'cLIS 0', high: 'cLIS 1 (confident)', key: 'livia.clisScale', def: 'greens', custom: ['#ffffff', '#00897b'],
        options: [['greens', 'Greens (default)'], ['blues', 'Blues'], ['purples', 'Purples'], ['oranges', 'Oranges'], ['reds', 'Reds'], ['greys', 'Grays'], ['viridis', 'Viridis'], ['custom', 'Custom']],
        fns: { blues: _mapBlues, greens: _mapGreens }, stops: _MAP_SEQ },
};
const mapScale = {};
for (const k of Object.keys(MAP_SCALES)) {
    const c = MAP_SCALES[k];
    mapScale[k] = { name: c.def, custom: c.custom.slice() };
    try {
        const saved = JSON.parse(localStorage.getItem(c.key) || 'null');
        if (saved && (saved.name === 'custom' || c.fns[saved.name] || c.stops[saved.name])) {
            const own = Array.isArray(saved.custom) && saved.custom.length >= 2
                && !(c.oldDefaultCustom && saved.custom.join() === c.oldDefaultCustom.join());
            mapScale[k] = { name: saved.name, custom: own ? saved.custom : c.custom.slice() };
        }
    } catch (e) { /* storage blocked: keep the default */ }
}
let onMapScaleChange = null;   // page sets this: (kind) => redraw that kind's maps and color bars

function mapColorFn(kind) {
    const c = MAP_SCALES[kind], s = mapScale[kind];
    if (s.name === 'custom') return _mapRamp(s.custom);
    return c.fns[s.name] || (c.stops[s.name] ? _mapRamp(c.stops[s.name]) : c.fns[c.def]);
}
function setMapScale(kind, name, custom) {
    const c = MAP_SCALES[kind];
    mapScale[kind] = { name: name || c.def, custom: custom && custom.length >= 2 ? custom.map(_mapHex6) : mapScale[kind].custom };
    try { localStorage.setItem(c.key, JSON.stringify(mapScale[kind])); } catch (e) {}
    document.querySelectorAll('.map-scale-ctl[data-kind="' + kind + '"]').forEach(_syncMapControl);
    if (typeof onMapScaleChange === 'function') onMapScaleChange(kind);
}
// "<map> colors" control: the map's preset list plus Custom (2 or 3 picked colors, or pasted hex codes).
function mapScaleControl(kind) {
    const c = MAP_SCALES[kind];
    const el = document.createElement('div');
    el.className = 'map-scale-ctl';
    el.dataset.kind = kind;
    el.style.cssText = 'display:inline-flex; align-items:center; gap:6px; flex-wrap:wrap; font-size:0.78rem; color:#555;';
    const hexHint = 'or paste hex: ' + c.custom.join(', ');   // the box is sized to show it whole
    const pick = (i, title) => '<input type="color" class="map-stop" data-i="' + i + '" title="' + title + '" style="width:24px; height:18px; padding:0; border:1px solid #ccd; border-radius:3px; cursor:pointer;">';
    el.innerHTML = '<span>' + c.label + ' colors</span>'
        + '<select class="map-scale-sel" title="Color scale for the ' + c.label + ' maps, ' + c.low + ' to ' + c.high + '" style="font-size:0.78rem; padding:1px 3px; border:1px solid #ccd; border-radius:4px;">'
        + c.options.map(([v, l]) => '<option value="' + v + '">' + l + '</option>').join('') + '</select>'
        + '<span class="map-custom" style="display:none; align-items:center; gap:4px;">'
        + pick(0, c.low) + '<span style="color:#aaa;">&rarr;</span>'
        + '<span class="map-mid" style="display:none; align-items:center; gap:4px;">' + pick(1, 'middle') + '<span style="color:#aaa;">&rarr;</span></span>'
        + pick(2, c.high)
        + '<label style="display:inline-flex; align-items:center; gap:3px; cursor:pointer; margin:0;"><input type="checkbox" class="map-3" style="margin:0;"> 3 colors</label>'
        + '<input type="text" class="map-hex" placeholder="' + hexHint + '" title="2 or 3 hex colors, ' + c.low + ' to ' + c.high + '" style="width:' + Math.round(hexHint.length * 6.2 + 10) + 'px; font-size:0.72rem; padding:1px 4px; border:1px solid #ccd; border-radius:4px;">'
        + '</span>';
    const stop = (i) => el.querySelector('.map-stop[data-i="' + i + '"]').value;
    const custom = () => el.querySelector('.map-3').checked ? [stop(0), stop(1), stop(2)] : [stop(0), stop(2)];
    el.querySelector('.map-scale-sel').onchange = (e) => setMapScale(kind, e.target.value, mapScale[kind].custom);
    el.querySelectorAll('.map-stop').forEach(x => { x.onchange = () => setMapScale(kind, 'custom', custom()); });   // 'change' = on release, not every drag step
    el.querySelector('.map-3').onchange = () => setMapScale(kind, 'custom', custom());
    el.querySelector('.map-hex').oninput = (e) => {
        const toks = (String(e.target.value).match(/#?[0-9a-fA-F]+/g) || []).map(h => h.replace('#', '')).filter(h => h.length === 3 || h.length === 6);
        if (toks.length >= 2) setMapScale(kind, 'custom', toks.slice(0, 3));
    };
    _syncMapControl(el);
    return el;
}
// One row of controls for several maps (dimer / monomer, whose PAE, LIS and cLIS panels share a row);
// it sits above the panels, so opening Custom never widens a panel. opts.align: 'start' | 'center'.
function mapScaleRow(kinds, opts) {
    const row = document.createElement('div');
    row.className = 'map-scale-row';
    row.style.cssText = 'display:flex; align-items:center; justify-content:' + (opts && opts.align === 'center' ? 'center' : 'flex-start') + '; gap:6px 18px; flex-wrap:wrap; margin:0.35rem 0;';
    for (const k of kinds) row.appendChild(mapScaleControl(k));
    return row;
}
function _syncMapControl(el) {
    const s = mapScale[el.dataset.kind], c = s.custom, three = c.length >= 3;
    el.querySelector('.map-scale-sel').value = s.name;
    el.querySelector('.map-custom').style.display = s.name === 'custom' ? 'inline-flex' : 'none';
    el.querySelector('.map-mid').style.display = three ? 'inline-flex' : 'none';
    el.querySelector('.map-3').checked = three;
    el.querySelector('.map-stop[data-i="0"]').value = _mapHex6(c[0]);
    el.querySelector('.map-stop[data-i="1"]').value = _mapHex6(three ? c[1] : '#ffffff');
    el.querySelector('.map-stop[data-i="2"]').value = _mapHex6(c[c.length - 1]);
}
