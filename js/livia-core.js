/**
 * livia-core.js — Shared core utilities for LIVIA tool pages
 *
 * Provides:
 *   esc()                  — HTML escaping
 *   ilisColor()            — iLIS FPR coloring with thresholds
 *   parseResidueRanges()   — parse residue range strings to Set of ints
 *   resSpec()              — ChimeraX residue spec builder
 *   fillGaps()             — fill gaps in residue position sets
 *   filterSmallSegments()  — remove LIR segments below minimum size
 *   getGapFillValue()      — read gap-fill-input element value
 *   getMinSegmentValue()   — read min-segment-input element value
 *   rangesStr()            — compact range string for display
 *   hexToRgba()            — hex color to rgba string
 *   lightenColor()         — lighten a hex color by amount
 *   CORS_PROXIES           — array of CORS proxy URL builders
 *   fetchTextViaProxy()    — generic proxy-based fetch with direct-first fallback
 *   AA3TO1                 — amino acid 3-letter to 1-letter mapping
 */

// ── HTML escaping ──
function esc(s) { return s.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;'); }

// ── iLIS FPR coloring (purple/green/amber/gray thresholds) ──
function ilisColor(val, isAvg) {
    const color = isAvg
        ? (val >= 0.303 ? '#6B21A8' : val >= 0.120 ? '#0e8a6e' : val >= 0.073 ? '#bf8700' : '#8b949e')
        : (val >= 0.551 ? '#6B21A8' : val >= 0.339 ? '#0e8a6e' : val >= 0.223 ? '#bf8700' : '#8b949e');
    return '<span style="color:' + color + '; font-weight:700;">' + val.toFixed(3) + '</span>';
}

// ── Parse residue ranges: '["107-114","117","120-121"]' → Set of ints ──
function parseResidueRanges(raw) {
    if (!raw || raw === '[]') return new Set();
    let items;
    try { items = JSON.parse(raw); } catch {
        items = raw.match(/[\d]+-[\d]+|\d+/g) || [];
    }
    const positions = new Set();
    for (const item of items) {
        const s = String(item).trim();
        if (s.includes('-')) {
            const [a, b] = s.split('-').map(Number);
            for (let i = a; i <= b; i++) positions.add(i);
        } else {
            const n = parseInt(s);
            if (!isNaN(n)) positions.add(n);
        }
    }
    return positions;
}

// ── Build ChimeraX residue spec ──
function resSpec(positions, chain) {
    if (positions.size === 0) return '';
    const sorted = [...positions].sort((a, b) => a - b);
    const ranges = [];
    let start = sorted[0], end = sorted[0];
    for (let i = 1; i < sorted.length; i++) {
        if (sorted[i] === end + 1) { end = sorted[i]; }
        else { ranges.push(start === end ? `${start}` : `${start}-${end}`); start = end = sorted[i]; }
    }
    ranges.push(start === end ? `${start}` : `${start}-${end}`);
    return `/${chain}:${ranges.join(',')}`;
}

// ── Read gap-fill input value (default 10) ──
function getGapFillValue() {
    const el = document.getElementById('gap-fill-input');
    return el ? parseInt(el.value) || 10 : 10;
}

// ── Read min-segment input value (default 3) ──
function getMinSegmentValue() {
    const el = document.getElementById('min-segment-input');
    return el ? parseInt(el.value) || 3 : 3;
}

// ── Remove LIR segments shorter than minSize residues ──
function filterSmallSegments(positions, minSize) {
    if (minSize === undefined) minSize = getMinSegmentValue();
    if (minSize <= 1 || positions.size === 0) return positions;
    const sorted = [...positions].sort((a, b) => a - b);
    const segments = [];
    let seg = [sorted[0]];
    for (let i = 1; i < sorted.length; i++) {
        if (sorted[i] === seg[seg.length - 1] + 1) {
            seg.push(sorted[i]);
        } else {
            segments.push(seg);
            seg = [sorted[i]];
        }
    }
    segments.push(seg);
    const result = new Set();
    for (const s of segments) {
        if (s.length >= minSize) {
            for (const r of s) result.add(r);
        }
    }
    return result;
}

// ── Fill gaps for continuous cartoon ──
function fillGaps(positions, maxGap) {
    if (maxGap === undefined) maxGap = getGapFillValue();
    if (positions.size === 0) return new Set();
    const sorted = [...positions].sort((a, b) => a - b);
    const filled = new Set(sorted);
    for (let i = 0; i < sorted.length - 1; i++) {
        const gap = sorted[i + 1] - sorted[i] - 1;
        if (gap > 0 && gap <= maxGap) {
            for (let j = sorted[i] + 1; j < sorted[i + 1]; j++) filled.add(j);
        }
    }
    return filled;
}

// ── Compact range string for display (merges within 15 residues) ──
function rangesStr(positions) {
    if (positions.size === 0) return 'none';
    const sorted = [...positions].sort((a, b) => a - b);
    const ranges = [];
    let start = sorted[0], end = sorted[0];
    for (let i = 1; i < sorted.length; i++) {
        if (sorted[i] <= end + 15) { end = sorted[i]; }
        else { ranges.push(`${start}-${end}`); start = end = sorted[i]; }
    }
    ranges.push(`${start}-${end}`);
    return ranges.join(', ');
}

// ── Hex color to rgba string ──
function hexToRgba(hex, alpha) {
    hex = hex.replace('#', '');
    if (hex.length === 3) hex = hex[0]+hex[0]+hex[1]+hex[1]+hex[2]+hex[2];
    const r = parseInt(hex.substring(0, 2), 16);
    const g = parseInt(hex.substring(2, 4), 16);
    const b = parseInt(hex.substring(4, 6), 16);
    return `rgba(${r},${g},${b},${alpha})`;
}

// ── Lighten a hex color by amount (0–1) ──
function lightenColor(hex, amount) {
    hex = hex.replace('#', '');
    if (hex.length === 3) hex = hex[0]+hex[0]+hex[1]+hex[1]+hex[2]+hex[2];
    const r = parseInt(hex.substring(0, 2), 16);
    const g = parseInt(hex.substring(2, 4), 16);
    const b = parseInt(hex.substring(4, 6), 16);
    return `rgb(${Math.round(r+(255-r)*amount)},${Math.round(g+(255-g)*amount)},${Math.round(b+(255-b)*amount)})`;
}

// ── CORS proxy URL builders ──
const CORS_PROXIES = [
    url => `https://api.codetabs.com/v1/proxy?quest=${encodeURIComponent(url)}`,
    url => `https://api.allorigins.win/raw?url=${encodeURIComponent(url)}`,
    url => `https://corsproxy.io/?${encodeURIComponent(url)}`,
];

// ── Generic proxy fetch (tries direct first, then each proxy) ──
async function fetchTextViaProxy(url) {
    // Try direct first
    try {
        const resp = await fetch(url);
        if (resp.ok) return await resp.text();
    } catch {}
    // Try proxies
    for (const makeProxy of CORS_PROXIES) {
        try {
            const resp = await fetch(makeProxy(url), { signal: AbortSignal.timeout(10000) });
            if (resp.ok) {
                const text = await resp.text();
                if (text.length > 100) return text;
            }
        } catch {}
    }
    return null;
}

// ── Amino acid 3-letter to 1-letter code mapping ──
const AA3TO1 = {ALA:'A',ARG:'R',ASN:'N',ASP:'D',CYS:'C',GLN:'Q',GLU:'E',GLY:'G',HIS:'H',ILE:'I',LEU:'L',LYS:'K',MET:'M',PHE:'F',PRO:'P',SER:'S',THR:'T',TRP:'W',TYR:'Y',VAL:'V'};
