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
    document.getElementById('hex-lir-a').textContent = lirA;
    document.getElementById('hex-clir-a').textContent = clirA;
    document.getElementById('hex-lir-b').textContent = lirB;
    document.getElementById('hex-clir-b').textContent = clirB;
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

// ── PAE color scale ──
// Every PAE map LIVIA draws from numeric PAE (universal after a full upload, dimer, monomer), with its
// color bar and exports, takes its colors from paeColorFn(). FlyPredictome data and lightweight bundles
// carry PAE only as images, which cannot be recolored, so those views never show paeScaleControl().
// Scales run from low PAE (confident) to high PAE (uncertain); the choice is remembered per browser.
const PAE_SCALES = [
    ['bwr', 'Blue–white–red (default)'],
    ['alphafold', 'AlphaFold (green)'],
    ['blues', 'Blues'],
    ['viridis', 'Viridis'],
    ['greys', 'Greys'],
    ['custom', 'Custom'],
];
const _PAE_STOPS = {
    alphafold: ['#00441b', '#74c476', '#f7fcf5'],      // dark green = confident, as in AlphaFold DB / AlphaFold 3 PAE plots
    blues: ['#08306b', '#6baed6', '#f7fbff'],
    viridis: ['#440154', '#3b528b', '#21918c', '#5ec962', '#fde725'],
    greys: ['#000000', '#ffffff'],
};
let paeScale = { name: 'bwr', custom: ['#00897b', '#ffffff'] };
try {
    const saved = JSON.parse(localStorage.getItem('livia.paeScale') || 'null');
    if (saved && typeof saved.name === 'string') paeScale = { name: saved.name, custom: Array.isArray(saved.custom) && saved.custom.length >= 2 ? saved.custom : paeScale.custom };
} catch (e) { /* storage blocked: keep the default */ }
let onPaeScaleChange = null;   // page sets this to redraw its PAE maps and color bars

function _paeHex6(h) {   // '#abc' / 'ABC' / '#aabbcc' → '#aabbcc'
    h = String(h || '').replace('#', '').toLowerCase();
    return '#' + (h.length === 3 ? h.split('').map(c => c + c).join('') : h.padEnd(6, '0').slice(0, 6));
}
// The default is the original blue-white-red, computed exactly as before so default maps stay pixel-identical.
function _paeBwr(value, vmin, vmax) {
    const t = Math.max(0, Math.min(1, (value - vmin) / (vmax - vmin)));
    let r, g, b;
    if (t < 0.5) { const s = t / 0.5; r = Math.round(s * 255); g = Math.round(s * 255); b = 255; }
    else { const s = (t - 0.5) / 0.5; r = 255; g = Math.round((1 - s) * 255); b = Math.round((1 - s) * 255); }
    return [r, g, b];
}
function _paeRamp(hexStops) {   // piecewise-linear across 2+ color stops, low PAE → high PAE
    const stops = hexStops.map(h => { const x = _paeHex6(h); return [parseInt(x.slice(1, 3), 16), parseInt(x.slice(3, 5), 16), parseInt(x.slice(5, 7), 16)]; });
    const seg = stops.length - 1;
    return (value, vmin, vmax) => {
        const t = Math.max(0, Math.min(1, (value - vmin) / (vmax - vmin)));
        const f = t * seg, i = Math.min(seg - 1, Math.floor(f)), u = f - i, a = stops[i], b = stops[i + 1];
        return [Math.round(a[0] + (b[0] - a[0]) * u), Math.round(a[1] + (b[1] - a[1]) * u), Math.round(a[2] + (b[2] - a[2]) * u)];
    };
}
function paeColorFn() {
    const stops = paeScale.name === 'custom' ? paeScale.custom : _PAE_STOPS[paeScale.name];
    return stops && stops.length >= 2 ? _paeRamp(stops) : _paeBwr;
}
function setPaeScale(name, custom) {
    paeScale = { name: name || 'bwr', custom: custom && custom.length >= 2 ? custom.map(_paeHex6) : paeScale.custom };
    try { localStorage.setItem('livia.paeScale', JSON.stringify(paeScale)); } catch (e) {}
    document.querySelectorAll('.pae-scale-ctl').forEach(_syncPaeControl);
    if (typeof onPaeScaleChange === 'function') onPaeScaleChange();
}
// "PAE colors" control: a preset list plus Custom (2 or 3 picked colors, or pasted hex codes).
// opts.align: 'start' (default) or 'center', to line up with the maps it sits above.
function paeScaleControl(opts) {
    const el = document.createElement('div');
    el.className = 'pae-scale-ctl';
    el.style.cssText = 'display:flex; align-items:center; justify-content:' + (opts && opts.align === 'center' ? 'center' : 'flex-start') + '; gap:6px; flex-wrap:wrap; font-size:0.78rem; color:#555; margin:0.35rem 0;';
    const pick = (i, title) => '<input type="color" class="pae-stop" data-i="' + i + '" title="' + title + '" style="width:24px; height:18px; padding:0; border:1px solid #ccd; border-radius:3px; cursor:pointer;">';
    el.innerHTML = '<span>PAE colors</span>'
        + '<select class="pae-scale-sel" title="Color scale for the PAE maps, low (confident) to high (uncertain)" style="font-size:0.78rem; padding:1px 3px; border:1px solid #ccd; border-radius:4px;">'
        + PAE_SCALES.map(([v, l]) => '<option value="' + v + '">' + l + '</option>').join('') + '</select>'
        + '<span class="pae-custom" style="display:none; align-items:center; gap:4px;">'
        + pick(0, 'low PAE (confident)') + '<span style="color:#aaa;">&rarr;</span>'
        + '<span class="pae-mid" style="display:none; align-items:center; gap:4px;">' + pick(1, 'middle') + '<span style="color:#aaa;">&rarr;</span></span>'
        + pick(2, 'high PAE (uncertain)')
        + '<label style="display:inline-flex; align-items:center; gap:3px; cursor:pointer; margin:0;"><input type="checkbox" class="pae-3" style="margin:0;"> 3 colors</label>'
        + '<input type="text" class="pae-hex" placeholder="or paste hex: #00897b, #fff, #b2182b" title="2 or 3 hex colors, low to high PAE" style="width:200px; font-size:0.72rem; padding:1px 4px; border:1px solid #ccd; border-radius:4px;">'
        + '</span>';
    const stop = (i) => el.querySelector('.pae-stop[data-i="' + i + '"]').value;
    const custom = () => el.querySelector('.pae-3').checked ? [stop(0), stop(1), stop(2)] : [stop(0), stop(2)];
    el.querySelector('.pae-scale-sel').onchange = (e) => setPaeScale(e.target.value, paeScale.custom);
    el.querySelectorAll('.pae-stop').forEach(x => { x.onchange = () => setPaeScale('custom', custom()); });   // 'change' = on release, not every drag step
    el.querySelector('.pae-3').onchange = () => setPaeScale('custom', custom());
    el.querySelector('.pae-hex').oninput = (e) => {
        const toks = (String(e.target.value).match(/#?[0-9a-fA-F]+/g) || []).map(h => h.replace('#', '')).filter(h => h.length === 3 || h.length === 6);
        if (toks.length >= 2) setPaeScale('custom', toks.slice(0, 3));
    };
    _syncPaeControl(el);
    return el;
}
function _syncPaeControl(el) {
    const c = paeScale.custom, three = c.length >= 3;
    el.querySelector('.pae-scale-sel').value = paeScale.name;
    el.querySelector('.pae-custom').style.display = paeScale.name === 'custom' ? 'inline-flex' : 'none';
    el.querySelector('.pae-mid').style.display = three ? 'inline-flex' : 'none';
    el.querySelector('.pae-3').checked = three;
    el.querySelector('.pae-stop[data-i="0"]').value = _paeHex6(c[0]);
    el.querySelector('.pae-stop[data-i="1"]').value = _paeHex6(three ? c[1] : '#ffffff');
    el.querySelector('.pae-stop[data-i="2"]').value = _paeHex6(c[c.length - 1]);
}
