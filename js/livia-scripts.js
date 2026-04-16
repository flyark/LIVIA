/**
 * livia-scripts.js — Script generation helpers for LIVIA tool pages
 *
 * Provides:
 *   hexToRgbFloat()    — hex color to PyMOL [r, g, b] float string
 *   pmlResSpec()       — PyMOL residue spec builder
 *   highlightCxc()     — syntax highlighting for CXC/PML script preview
 *   scriptMode         — current mode: 'chimerax' or 'pymol'
 *   switchScriptTab()  — toggle between ChimeraX and PyMOL tabs
 *
 * Dependencies: livia-core.js (uses esc())
 *
 * Note: generateCxc() and generatePml() are page-specific (different data
 * structures per tool) and remain in each page's inline script.
 *
 * Pages should set onScriptTabSwitch to their updateScriptPreview function:
 *   onScriptTabSwitch = updateScriptPreview;
 */

// ── Callback for page-specific script preview update ──
let onScriptTabSwitch = null;

// ── Convert hex color to PyMOL RGB float array string ──
function hexToRgbFloat(hex) {
    hex = hex.replace('#', '');
    const r = parseInt(hex.substring(0, 2), 16) / 255;
    const g = parseInt(hex.substring(2, 4), 16) / 255;
    const b = parseInt(hex.substring(4, 6), 16) / 255;
    return `[${r.toFixed(3)}, ${g.toFixed(3)}, ${b.toFixed(3)}]`;
}

// ── Build PyMOL residue spec from positions Set ──
function pmlResSpec(positions, chain) {
    if (positions.size === 0) return '';
    const sorted = [...positions].sort((a, b) => a - b);
    const ranges = [];
    let start = sorted[0], end = sorted[0];
    for (let i = 1; i < sorted.length; i++) {
        if (sorted[i] === end + 1) { end = sorted[i]; }
        else { ranges.push(start === end ? `${start}` : `${start}-${end}`); start = end = sorted[i]; }
    }
    ranges.push(start === end ? `${start}` : `${start}-${end}`);
    return `chain ${chain} and resi ${ranges.join('+')}`;
}

// ── Syntax highlighting for ChimeraX/PyMOL script preview ──
function highlightCxc(text) {
    return text.split('\n').map(line => {
        if (line.startsWith('#')) return `<span class="comment">${esc(line)}</span>`;
        line = esc(line).replace(/(#[0-9A-Fa-f]{6})/g, '<span class="color-val">$1</span>');
        const m = line.match(/^(\w+)/);
        if (m) line = `<span class="command">${m[1]}</span>` + line.slice(m[1].length);
        return line;
    }).join('\n');
}

// ── Script mode state ──
let scriptMode = 'chimerax';

// ── Toggle between ChimeraX and PyMOL script tabs ──
function switchScriptTab(mode) {
    scriptMode = mode;
    document.querySelectorAll('.input-tab[data-tab="chimerax"], .input-tab[data-tab="pymol"]').forEach(t => t.classList.remove('active'));
    document.querySelector('.input-tab[data-tab="' + mode + '"]').classList.add('active');
    // Update download button labels
    var ext = mode === 'pymol' ? '.pml' : '.cxc';
    var dlExt = document.getElementById('dl-script-ext');
    if (dlExt) dlExt.textContent = ext;
    var dlAll = document.getElementById('dl-all-chains');
    if (dlAll) dlAll.textContent = '\u2193 Download All Chains ' + ext;
    if (onScriptTabSwitch) onScriptTabSwitch();
}
