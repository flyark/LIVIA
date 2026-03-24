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

// ── Callback for page-specific reload after color change ──
let onColorChange = null;

// ── Bind color picker input events (call once after DOM ready) ──
function initColorPickers() {
    document.querySelectorAll('input[type="color"]').forEach(input => {
        input.addEventListener('input', () => {
            const hexEl = document.getElementById(input.id.replace('color-', 'hex-'));
            if (hexEl) hexEl.textContent = input.value;
            updateColorStrip();
        });
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
