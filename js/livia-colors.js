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

// Lighten a colour toward white by amt (0..1) → the LIR (light) shade paired with a cLIR (dark) colour.
function _lightenHex(hex, amt) {
    const h = _normalizeHex(hex); if (!h) return hex;
    const m = (c) => Math.round(c + (255 - c) * amt).toString(16).padStart(2, '0');
    return '#' + m(parseInt(h.slice(1, 3), 16)) + m(parseInt(h.slice(3, 5), 16)) + m(parseInt(h.slice(5, 7), 16));
}
// Resolve any CSS colour (#hex, #rgb, or a name like "teal") to #rrggbb, or null if invalid.
function _colorToHex(c) {
    const norm = _normalizeHex(c); if (norm) return norm;
    try {
        const ctx = (_colorToHex._ctx || (_colorToHex._ctx = document.createElement('canvas').getContext('2d')));
        ctx.fillStyle = '#000'; ctx.fillStyle = String(c); const r = ctx.fillStyle;
        return /^#[0-9a-f]{6}$/i.test(r) ? r.toLowerCase() : null;
    } catch (e) { return null; }
}

// CSV/paste colour upload for a 2-chain (A/B) figure. Accepts "chain,color" (also 3-col "chain,name,color";
// the name is ignored here). Sets each chain's cLIR (main) colour and derives the light LIR shade, then
// applies through applyPreset() (which re-renders on every page). Match by chain letter A/B (or 1/2), or an
// optional gene name via opts.geneOf → { A:'sym1', B:'sym2' }.
function attachChainColorUpload(container, opts) {
    opts = opts || {};
    if (!container || !window.LiviaMaps || container.dataset.ccu) return;
    container.dataset.ccu = '1';
    const geneOf = opts.geneOf || (() => ({}));
    const val = (id) => (document.getElementById(id) || {}).value;
    window.LiviaMaps.attachColorUpload(container, {
        label: 'chain colours', keyHeader: 'chain',
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
