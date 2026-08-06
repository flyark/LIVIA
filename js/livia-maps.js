/*
 * livia-maps.js — shared figure export for the LIVIA map pages.
 *
 * Vector SVG export works by re-running a figure's own canvas draw function against a
 * canvas2svg (C2S) context instead of the 2D bitmap context, then serializing the result.
 * That keeps one exporter for every page rather than a copy per file.
 *
 * Requires js/canvas2svg.js to be loaded first (it provides window.C2S).
 *
 * Usage:
 *   LiviaMaps.attachExportBar('contact-canvas', {
 *       name: 'linear_contact_map',
 *       svg: () => buildDimerContactMap(),   // a no-arg redraw of THIS canvas
 *       png: true,
 *   });
 */
(function (global) {
    'use strict';

    let EXPORT_FONT = 1;   // font-size multiplier applied to the serialized SVG
    let EXPORT_W = 0;      // target SVG width in px (0 = keep the drawn width / aspect from height)
    let EXPORT_H = 0;      // target SVG height in px (0 = keep the drawn height / aspect from width)

    function setExportOpts(o) {
        if (!o) return;
        if (o.font != null) EXPORT_FONT = +o.font;
        if (o.width != null) EXPORT_W = +o.width;
        if (o.height != null) EXPORT_H = +o.height;
    }

    // Add a viewBox (canvas2svg omits it) and optionally rescale width/height/fonts for publication.
    // Width and height each default to "auto": set one → the other follows aspect; set both → exact.
    function applyExportOpts(svg) {
        if (EXPORT_FONT && EXPORT_FONT !== 1) {
            svg = svg.replace(/font-size="([\d.]+)px"/g, (m, s) => 'font-size="' + (parseFloat(s) * EXPORT_FONT).toFixed(2) + 'px"');
        }
        const m = svg.match(/<svg\b[^>]*\bwidth="([\d.]+)"[^>]*\bheight="([\d.]+)"/);
        if (m) {
            const ow = parseFloat(m[1]), oh = parseFloat(m[2]);
            if (!/viewBox=/.test(svg)) svg = svg.replace(/<svg\b/, '<svg viewBox="0 0 ' + ow + ' ' + oh + '"');
            if (EXPORT_W || EXPORT_H) {
                let tw = EXPORT_W, th = EXPORT_H;
                if (tw && !th) th = Math.round(oh * tw / ow);        // width set → height by aspect
                else if (!tw && th) tw = Math.round(ow * th / oh);   // height set → width by aspect
                svg = svg.replace(/(<svg\b[^>]*?)\bwidth="[\d.]+"/, '$1width="' + tw + '"')
                         .replace(/(<svg\b[^>]*?)\bheight="[\d.]+"/, '$1height="' + th + '"');
            }
        }
        return svg;
    }

    // For draw functions that accept the canvas (or a canvas-like object) as an argument.
    function svgFromDraw(drawClosure) {
        const C = global.C2S;
        if (!C) throw new Error('SVG library not loaded (needs internet)');
        let ctx = null;
        const fake = {
            style: {}, _w: 300, _h: 150,
            set width(w) { this._w = w; }, get width() { return this._w; },
            set height(h) { this._h = h; }, get height() { return this._h; },
            getContext() { ctx = new C(this._w, this._h); return ctx; },
        };
        drawClosure(fake);
        return ctx ? ctx.getSerializedSvg() : '';
    }

    // For figures that draw into a fixed <canvas> element: swap its getContext for a C2S
    // context, re-run the draw, then restore the element and repaint the on-screen bitmap.
    function svgFromFixedCanvas(canvasId, redraw) {
        const C = global.C2S;
        if (!C) throw new Error('SVG library not loaded (needs internet)');
        const cv = document.getElementById(canvasId);
        if (!cv) return '';
        const real = cv.getContext;
        let c2s = null;
        cv.getContext = function () { c2s = new C(this.width || 300, this.height || 150); return c2s; };
        try { redraw(); } finally { cv.getContext = real; }
        const svg = c2s ? c2s.getSerializedSvg() : '';
        redraw();   // restore the on-screen bitmap render
        return svg;
    }

    function safeName(name) { return String(name || 'figure').replace(/[^\w.-]+/g, '_'); }

    function downloadSVGFile(name, svg) {
        const a = document.createElement('a');
        a.href = URL.createObjectURL(new Blob([applyExportOpts(svg)], { type: 'image/svg+xml' }));
        a.download = safeName(String(name).replace(/\.svg$/i, '')) + '.svg';
        a.click();
        setTimeout(() => URL.revokeObjectURL(a.href), 1000);
    }

    function downloadSVGFromCanvas(canvasId, name, redraw) {
        try { downloadSVGFile(name, svgFromFixedCanvas(canvasId, redraw)); }
        catch (e) { alert('SVG export failed: ' + e.message); }
    }

    function downloadCanvasPNG(canvasId, name) {
        const cv = document.getElementById(canvasId);
        if (!cv) return;
        const a = document.createElement('a');
        a.href = cv.toDataURL('image/png');
        a.download = safeName(String(name).replace(/\.png$/i, '')) + '.png';
        a.click();
    }

    // PNG at a custom width/height (blank/0 = native; one set → the other follows aspect).
    function downloadCanvasPNGScaled(canvasId, name, w, h) {
        const cv = document.getElementById(canvasId);
        if (!cv) return;
        let out = cv; const ow = cv.width, oh = cv.height; let tw = w || 0, th = h || 0;
        if (tw || th) {
            if (tw && !th) th = Math.round(oh * tw / ow); else if (!tw && th) tw = Math.round(ow * th / oh);
            const tc = document.createElement('canvas'); tc.width = tw; tc.height = th;
            const cx = tc.getContext('2d'); cx.imageSmoothingEnabled = true; cx.imageSmoothingQuality = 'high'; cx.drawImage(cv, 0, 0, tw, th); out = tc;
        }
        const a = document.createElement('a'); a.href = out.toDataURL('image/png'); a.download = safeName(String(name).replace(/\.png$/i, '')) + '.png'; a.click();
    }

    // Shared datalists backing the export size widgets: "auto" + suggested values (custom typing allowed).
    function ensureExportDatalists() {
        if (typeof document === 'undefined' || document.getElementById('lm-exp-dim')) return;
        const mk = (id, vals) => { const dl = document.createElement('datalist'); dl.id = id; dl.innerHTML = vals.map((v) => '<option value="' + v + '"></option>').join(''); document.body.appendChild(dl); };
        mk('lm-exp-dim', ['auto', '600', '900', '1200', '1600']);
        mk('lm-exp-font', ['1', '1.25', '1.5', '2', '0.8']);
    }

    // Idempotent "↓ SVG · ↓ PNG" bar (plus custom W / H / font× entry) directly under a canvas.
    // Pass opts.dims === false to omit the size controls. W/H are px (blank = auto, aspect-preserving);
    // font× is a free multiplier applied to the SVG text.
    function attachExportBar(canvasId, opts) {
        opts = opts || {};
        const cv = document.getElementById(canvasId);
        if (!cv || !cv.parentElement) return;
        const barId = canvasId + '-export-bar';
        const existing = document.getElementById(barId);
        if (existing) existing.remove();          // rebuild so the redraw closure stays current
        const showDims = opts.dims !== false;
        const bar = document.createElement('div');
        bar.id = barId;
        bar.style.cssText = showDims
            ? 'display:flex; flex-wrap:wrap; gap:5px 8px; align-items:center; justify-content:center; margin-top:3px; font-size:0.7rem; color:#888;'
            : 'text-align:center; margin-top:3px;';
        let wIn, hIn, fIn;
        const readOpts = () => ({ font: (fIn && +fIn.value) || 1, width: (wIn && +wIn.value) || 0, height: (hIn && +hIn.value) || 0 });
        const mk = (label, fn) => {
            const b = document.createElement('button');
            b.type = 'button'; b.textContent = label;
            b.style.cssText = 'font-size:0.72rem; color:#2471A3; background:none; border:none; cursor:pointer; font-weight:600; padding:0 6px;';
            b.onclick = fn;
            bar.appendChild(b);
        };
        if (opts.svg) mk('↓ SVG', () => { if (showDims) setExportOpts(readOpts()); downloadSVGFromCanvas(canvasId, opts.name, opts.svg); });
        if (opts.png) mk('↓ PNG', () => { const o = showDims ? readOpts() : {}; downloadCanvasPNGScaled(canvasId, opts.name, o.width || 0, o.height || 0); });
        if (showDims && bar.children.length) {
            ensureExportDatalists();
            const mkIn = (val, ph, list, w) => { const i = document.createElement('input'); i.type = 'text'; i.setAttribute('list', list); if (val) i.value = val; if (ph) i.placeholder = ph; i.style.cssText = 'width:' + w + 'px; font-size:0.7rem; padding:1px 3px; color:#555;'; return i; };
            const lbl = (t) => { const s = document.createElement('span'); s.textContent = t; return s; };
            wIn = mkIn('', 'auto', 'lm-exp-dim', 58); hIn = mkIn('', 'auto', 'lm-exp-dim', 58); fIn = mkIn('1', '', 'lm-exp-font', 52);
            bar.append(lbl('W'), wIn, lbl('H'), hIn, lbl('font ×'), fIn);
        }
        bar.__opts = opts;                        // exposed so callers/tests can reuse the redraw closure
        if (bar.children.length) cv.parentElement.appendChild(bar);
    }

    // Merge a Set of 1-based residue numbers into [start, end] runs, so contiguous residues
    // paint as one block instead of hairline-gapped per-residue slivers.
    function setToRanges(set) {
        if (!set || !set.size) return [];
        const xs = [...set].map(Number).filter(n => !isNaN(n)).sort((a, b) => a - b);
        const out = [];
        let s = xs[0], p = xs[0];
        for (let i = 1; i < xs.length; i++) {
            if (xs[i] === p + 1) { p = xs[i]; continue; }
            out.push([s, p]); s = p = xs[i];
        }
        out.push([s, p]);
        return out;
    }

    // A chord arc drawn as three bands over its residue span: non-LIR gray base,
    // LIR in the light shade, cLIR in the dark shade. `len` is the arc's residue count and
    // lirSet/clirSet hold 1-based indices within that span.
    function paintChordArcBand(ctx, cx, cy, outerR, innerR, sA, eA, len, lirSet, clirSet, lightCol, darkCol) {
        const donut = (a0, a1) => {
            ctx.beginPath();
            ctx.arc(cx, cy, outerR, a0, a1);
            ctx.arc(cx, cy, innerR, a1, a0, true);
            ctx.closePath();
        };
        donut(sA, eA); ctx.fillStyle = '#e8e8e8'; ctx.fill();
        const paint = (set, color) => {
            if (!set || !set.size) return;
            ctx.fillStyle = color;
            for (const [s, e] of setToRanges(set)) {
                const a0 = sA + (eA - sA) * ((s - 1) / Math.max(1, len));
                const a1 = sA + (eA - sA) * (e / Math.max(1, len));
                donut(a0, a1); ctx.fill();
            }
        };
        paint(lirSet, lightCol);
        paint(clirSet, darkCol);
        donut(sA, eA); ctx.strokeStyle = '#fff'; ctx.lineWidth = 1; ctx.stroke();
    }

    // Domain segments as a thin ring band over their residue spans on a chord arc. rIn/rOut = ring radii,
    // sA/eA = the arc's angular span, len = arc residue count, domains = [{name,start,end}] (1-based).
    // Each domain gets a fill from `fills` (cycled); `muted` greys a hidden chain. Mirrors the linear
    // map's UniProt/Pfam + TED tracks as concentric rings just inside the arc.
    function paintChordDomainRing(ctx, cx, cy, rIn, rOut, sA, eA, len, domains, fills, muted) {
        const hits = [];   // returned so the caller can build hover hit-regions
        if (!domains || !domains.length) return hits;
        const ang = (res) => sA + (eA - sA) * ((res - 1) / Math.max(1, len - 1));
        for (let i = 0; i < domains.length; i++) {
            const d = domains[i], a0 = ang(d.start), a1 = ang(d.end);
            ctx.beginPath();
            ctx.arc(cx, cy, rOut, a0, a1);
            ctx.arc(cx, cy, rIn, a1, a0, true);
            ctx.closePath();
            ctx.fillStyle = muted ? '#dddddd' : fills[i % fills.length];
            ctx.fill();
            ctx.strokeStyle = 'rgba(255,255,255,0.9)'; ctx.lineWidth = 0.75; ctx.stroke();
            hits.push({ rIn, rOut, a0, a1, name: d.name, start: d.start, end: d.end });
        }
        return hits;
    }

    // Hover tooltips over the chord domain rings. Each drawChord* stashes per-canvas hit-regions on
    // `canvas._domHits` (each {rIn,rOut,a0,a1,name,start,end,chain,src}) and geometry on `canvas._domGeom`
    // ({cx,cy,size}). Attaching is idempotent per canvas; call it whenever the chord is (re)built.
    let _cdTip = null;
    const _cdEsc = (s) => String(s == null ? '' : s).replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
    function _cdAt(cv, ev) {
        const hits = cv._domHits, geom = cv._domGeom;
        if (!hits || !hits.length || !geom) return null;
        const rect = cv.getBoundingClientRect(); if (!rect.width) return null;
        const scale = geom.size / rect.width;
        const x = (ev.clientX - rect.left) * scale - geom.cx, y = (ev.clientY - rect.top) * scale - geom.cy;
        const r = Math.hypot(x, y), a = Math.atan2(y, x);
        for (const h of hits) { if (r < h.rIn - 1 || r > h.rOut + 1) continue; let aa = a; while (aa < h.a0) aa += 2 * Math.PI; if (aa <= h.a1) return h; }
        return null;
    }
    function attachChordDomainHover(ids) {
        for (const id of ids) {
            const cv = document.getElementById(id);
            if (!cv || cv._domHoverAttached) continue;
            cv._domHoverAttached = true;
            cv.addEventListener('mousemove', (ev) => {
                const h = _cdAt(cv, ev);
                if (h) {
                    if (!_cdTip) { _cdTip = document.createElement('div'); _cdTip.style.cssText = 'position:fixed;z-index:9999;pointer-events:none;background:rgba(30,30,30,0.92);color:#fff;padding:5px 9px;border-radius:5px;font-size:0.78rem;line-height:1.45;max-width:280px;box-shadow:0 2px 10px rgba(0,0,0,0.35)'; document.body.appendChild(_cdTip); }
                    _cdTip.innerHTML = '<b>' + _cdEsc(h.name) + '</b><br>' + _cdEsc(h.chain) + ' &middot; ' + h.start + '–' + h.end + ' <span style="color:#9ecbff;">(' + _cdEsc(h.src) + ')</span>';
                    _cdTip.style.left = Math.round(ev.clientX + 13) + 'px'; _cdTip.style.top = Math.round(ev.clientY + 13) + 'px'; _cdTip.style.display = ''; cv.style.cursor = 'help';
                } else { if (_cdTip) _cdTip.style.display = 'none'; cv.style.cursor = ''; }
            });
            cv.addEventListener('mouseleave', () => { if (_cdTip) _cdTip.style.display = 'none'; cv.style.cursor = ''; });
        }
    }

    // Nice-number major-tick spacing (1/2/5 × 10^k) aiming for ~5 major ticks across the span.
    function chordRulerStep(span) {
        const raw = span / 5, mag = Math.pow(10, Math.floor(Math.log10(raw))), n = raw / mag;
        return Math.max(1, (n < 1.5 ? 1 : n < 3 ? 2 : n < 7 ? 5 : 10) * mag);
    }
    // Circos-style residue ruler around one arc spanning [sA,eA] over `len` residues: minor + major
    // ticks radiating outward from outerR, radial numbers at the majors (at nice round residue
    // coordinates). Same ctx as the display canvas AND the canvas2svg export, so it exports too.
    // opts: {muted, type, minResidues, offset, color, textColor, font}. `offset` shifts the printed
    // number to an absolute coordinate (residue = within-arc index + offset) — 0 for a full chain,
    // start-1 for a subdomain. Short (< minResidues, default 20) or non-polymer arcs get no ruler.
    function drawChordRuler(ctx, cx, cy, outerR, sA, eA, len, opts) {
        opts = opts || {};
        if (opts.muted || len < (opts.minResidues || 20)) return;
        const t = opts.type; if (t === 'ion' || t === 'glycan' || t === 'ligand') return;
        const off = opts.offset || 0, lo = off + 1, hi = off + len;     // absolute residue range on this arc
        const angOf = abs => sA + ((abs - off - 1) / Math.max(1, len - 1)) * (eA - sA);
        const major = chordRulerStep(len), minor = Math.max(1, Math.round(major / 5));
        ctx.save();
        ctx.strokeStyle = opts.color || '#9aa3ac'; ctx.lineWidth = 0.8;
        for (let a = Math.ceil(lo / minor) * minor; a <= hi; a += minor) {   // minor ticks at round coords (skip majors)
            if (major && a % major === 0) continue;
            const ang = angOf(a), c = Math.cos(ang), s = Math.sin(ang);
            ctx.beginPath(); ctx.moveTo(cx + c * outerR, cy + s * outerR); ctx.lineTo(cx + c * (outerR + 3), cy + s * (outerR + 3)); ctx.stroke();
        }
        ctx.fillStyle = opts.textColor || '#666'; ctx.font = opts.font || '9px sans-serif'; ctx.textBaseline = 'middle';
        for (let a = Math.ceil(lo / major) * major; a <= hi; a += major) {   // major ticks + radial numbers
            const ang = angOf(a), c = Math.cos(ang), s = Math.sin(ang);
            ctx.beginPath(); ctx.moveTo(cx + c * outerR, cy + s * outerR); ctx.lineTo(cx + c * (outerR + 6), cy + s * (outerR + 6)); ctx.stroke();
            ctx.save();
            ctx.translate(cx + c * (outerR + 8), cy + s * (outerR + 8));
            if (c < 0) { ctx.rotate(ang + Math.PI); ctx.textAlign = 'right'; } else { ctx.rotate(ang); ctx.textAlign = 'left'; }
            ctx.fillText(String(a), 0, 0);
            ctx.restore();
        }
        ctx.restore();
    }

    // Chord arc label (subunit name + residue count) in either orientation.
    //   radial=false : tangential — text follows the arc curve (the original look).
    //   radial=true  : text points outward along the radius, so labels fan out and don't collide
    //                  when there are many subunits. Name then length stack outward; flipped on the
    //                  left half so it always reads outward.
    // opts: { cx, cy, outerR, angle, ticks, radial, name, showName, length, showLength, hidden }
    function drawChordArcLabel(ctx, opts) {
        const cx = opts.cx, cy = opts.cy, R0 = opts.outerR, a = opts.angle, ticks = opts.ticks;
        const nameFont = 'bold 13px sans-serif', lenFont = '11px sans-serif';
        const nameColor = opts.hidden ? '#bbb' : '#333', lenColor = opts.hidden ? '#ccc' : '#888';
        const at = (r) => [cx + Math.cos(a) * r, cy + Math.sin(a) * r];
        if (!opts.radial) {
            let ta = a + Math.PI / 2;
            const n = ((ta % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
            if (n > Math.PI / 2 && n < 3 * Math.PI / 2) ta += Math.PI;
            const draw = (txt, r, font, color) => {
                const p = at(r); ctx.save(); ctx.translate(p[0], p[1]); ctx.rotate(ta);
                ctx.fillStyle = color; ctx.font = font; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
                ctx.fillText(txt, 0, 0); ctx.restore();
            };
            if (opts.showName) draw(opts.name, R0 + (ticks ? 34 : 18), nameFont, nameColor);
            if (opts.showLength) draw(opts.length, R0 + (ticks ? 50 : 35), lenFont, lenColor);
            return;
        }
        const norm = ((a % (2 * Math.PI)) + 2 * Math.PI) % (2 * Math.PI);
        const flip = norm > Math.PI / 2 && norm < 3 * Math.PI / 2;
        const ra = flip ? a + Math.PI : a;
        const rStart = R0 + (ticks ? 26 : 8);               // clear the ruler numbers
        const nm = opts.showName ? (opts.name || '') : '', ln = opts.showLength ? (opts.length || '') : '';
        ctx.font = nameFont; const nameW = nm ? ctx.measureText(nm).width : 0;
        ctx.font = lenFont; const lenW = ln ? ctx.measureText(ln).width : 0;
        const maxW = Math.max(nameW, lenW);
        if (maxW <= 0) return;
        const both = nameW > 0 && lenW > 0;                 // length goes on a second line BELOW the name (not further out — avoids clipping)
        const p = at(rStart + maxW / 2); ctx.save(); ctx.translate(p[0], p[1]); ctx.rotate(ra);
        ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
        if (nm) { ctx.font = nameFont; ctx.fillStyle = nameColor; ctx.fillText(nm, 0, both ? -8 : 0); }
        if (ln) { ctx.font = lenFont; ctx.fillStyle = lenColor; ctx.fillText(ln, 0, both ? 7 : 0); }
        ctx.restore();
    }

    // Max radial extent (px beyond outerR) the labels need, so a page can shrink outerR to keep them
    // on-canvas. labels: [{name, length}]. Mirrors drawChordArcLabel's geometry.
    function chordLabelExtent(ctx, labels, opts) {
        opts = opts || {};
        const ticks = opts.ticks;
        if (!opts.radial) return (ticks ? 50 : 35) + 14;    // tangential: length ring + a line of text
        let maxW = 0;
        for (const L of (labels || [])) {
            ctx.font = 'bold 13px sans-serif'; const nw = L.name ? ctx.measureText(L.name).width : 0;
            ctx.font = '11px sans-serif'; const lw = L.length ? ctx.measureText(L.length).width : 0;
            maxW = Math.max(maxW, nw, lw);
        }
        return (ticks ? 26 : 8) + maxW + 4;                 // ruler gap + longest label + a little padding
    }

    // Highlight the clicked button in a 2-option segmented control (used for the chord label-direction
    // picker); un-highlights its siblings. The caller sets the state + redraws.
    function segPick(btn) {
        for (const b of btn.parentElement.children) {
            const on = b === btn;
            b.style.background = on ? '#2471A3' : 'transparent';
            b.style.color = on ? '#fff' : '#555';
        }
    }

    // PNG bar for a canvas ELEMENT (charts are built inside a container and carry no id).
    function attachPngBtnFor(cv, name) {
        if (!cv || !cv.parentElement) return;
        if (cv.nextElementSibling && cv.nextElementSibling.classList.contains('livia-png-bar')) return;
        const bar = document.createElement('div');
        bar.className = 'livia-png-bar';
        bar.style.cssText = 'text-align:center; margin-top:2px;';
        const b = document.createElement('button');
        b.type = 'button'; b.textContent = '↓ PNG';
        b.style.cssText = 'font-size:0.72rem; color:#2471A3; background:none; border:none; cursor:pointer; font-weight:600; padding:0 6px;';
        b.onclick = () => {
            const a = document.createElement('a');
            a.href = cv.toDataURL('image/png');
            a.download = safeName(String(name).replace(/\.png$/i, '')) + '.png';
            a.click();
        };
        bar.appendChild(b);
        cv.parentElement.insertBefore(bar, cv.nextSibling);
    }

    // Give every canvas inside a chart container its own ↓PNG button.
    function attachChartPngButtons(containerId, prefix) {
        const box = document.getElementById(containerId);
        if (!box) return;
        box.querySelectorAll('canvas').forEach((cv, i) => attachPngBtnFor(cv, (prefix || 'chart') + '_' + (i + 1)));
    }

    // Parse a colour CSV. Two forms, auto-detected per row:
    //   "<key>,<color>"          — key (chain/node/cluster) + colour
    //   "<key>,<name>,<color>"   — key + display-label override + colour (chain relabel)
    // Tolerant of an optional header, quotes, and comma/tab/semicolon separators. Colour = #RGB / #RRGGBB / CSS name.
    function parseColorCSV(text) {
        const out = [];
        const isColor = (c) => /^#([0-9a-fA-F]{3}|[0-9a-fA-F]{6})$/.test(c) || /^[a-zA-Z]{3,20}$/.test(c);
        for (const raw of String(text || '').split(/\r?\n/)) {
            const line = raw.trim(); if (!line) continue;
            const parts = line.split(/[,\t;]+/).map((s) => s.trim().replace(/^["']|["']$/g, ''));
            if (parts.length < 2 || !parts[0]) continue;
            if (/^(chain|node|protein|cluster|community|id|name|key|label)$/i.test(parts[0]) && parts.some((p) => /^colou?r$/i.test(p))) continue;   // header row
            if (parts.length >= 3 && parts[2] && isColor(parts[2])) { out.push({ key: parts[0], name: parts[1] || '', color: parts[2] }); continue; }   // key,name,color
            if (parts[1] && isColor(parts[1])) out.push({ key: parts[0], color: parts[1] });                                                            // key,color
        }
        return out;
    }

    // Compact "custom colours" widget: upload a CSV file, paste rows, apply, and download the current
    // set (round-trip). opts = { label, keyHeader, currentRows: () => [{key,color}], apply: rows => appliedCount }.
    function attachColorUpload(container, opts) {
        opts = opts || {};
        const kh = opts.keyHeader || 'chain';
        const wrap = document.createElement('div');
        wrap.style.cssText = 'display:flex; flex-wrap:wrap; align-items:center; gap:6px; font-size:0.75rem; color:#666; margin-top:6px;';
        const mkBtn = (t) => { const b = document.createElement('button'); b.type = 'button'; b.textContent = t; b.style.cssText = 'font-size:0.72rem; color:#2471A3; background:none; border:1px solid #ccd6e0; border-radius:5px; padding:2px 7px; cursor:pointer; font-weight:600;'; return b; };
        const status = document.createElement('span'); status.style.cssText = 'font-size:0.72rem; color:#888;';
        const setStatus = (m, err) => { status.textContent = m; status.style.color = err ? '#c0392b' : '#888'; };
        const doApply = (text) => {
            const rows = parseColorCSV(text);
            if (!rows.length) { setStatus('no valid rows (need "' + kh + ',#hex")', true); return; }
            let applied = 0; try { applied = opts.apply(rows) || 0; } catch (e) { setStatus('error: ' + e.message, true); return; }
            setStatus('applied ' + applied + ' / ' + rows.length + (applied < rows.length ? ' (some keys not matched)' : ''), applied === 0);
        };
        const file = document.createElement('input'); file.type = 'file'; file.accept = '.csv,.txt'; file.style.display = 'none';
        file.onchange = () => { const f = file.files && file.files[0]; if (!f) return; const r = new FileReader(); r.onload = () => doApply(r.result); r.readAsText(f); file.value = ''; };
        const upBtn = mkBtn('↑ CSV file'); upBtn.onclick = () => file.click();
        const ta = document.createElement('textarea'); ta.placeholder = opts.placeholder || (kh + ',color\nA,#1f77b4\nB,#ff7f0e'); ta.rows = 2;
        ta.style.cssText = 'display:none; font:inherit; font-size:0.72rem; width:180px; padding:2px 4px; border:1px solid #ccd6e0; border-radius:5px;';
        const pasteBtn = mkBtn('paste'); pasteBtn.onclick = () => { ta.style.display = ta.style.display === 'none' ? 'inline-block' : 'none'; if (ta.style.display !== 'none') ta.focus(); };
        const applyBtn = mkBtn('Apply'); applyBtn.onclick = () => { if (ta.value.trim()) doApply(ta.value); };
        const dlBtn = mkBtn('↓ CSV'); dlBtn.onclick = () => {
            const rows = (opts.currentRows && opts.currentRows()) || [];
            const has3 = rows.some((r) => 'name' in r);
            const csv = (has3 ? kh + ',name,color\n' : kh + ',color\n')
                + rows.map((r) => has3 ? (r.key + ',' + (r.name || '') + ',' + r.color) : (r.key + ',' + r.color)).join('\n') + '\n';
            const a = document.createElement('a'); a.href = URL.createObjectURL(new Blob([csv], { type: 'text/csv' })); a.download = kh + '_colors.csv'; a.click(); setTimeout(() => URL.revokeObjectURL(a.href), 1000);
        };
        const lbl = document.createElement('span'); lbl.textContent = '🎨 custom ' + (opts.label || 'colours') + ':'; lbl.style.fontWeight = '600';
        wrap.append(lbl, upBtn, pasteBtn, ta, applyBtn, dlBtn, status);
        container.appendChild(wrap);
        return wrap;
    }

    // ── Two-colour contact curves ─────────────────────────────────────────────
    //
    // A contact line fades smoothly from strokeA at one end to strokeB at the other
    // via a canvas linear gradient — the look users expect. Caller passes ready
    // stroke styles (usually hexToRgba(col, alpha)) so the colour logic stays at the
    // call site, where it differs per figure. In the SVG export canvas2svg emits one
    // <linearGradient> per line; the coordinate rounding in canvas2svg.js keeps each
    // compact. A follow-up
    // pass can consolidate the per-line defs into a few shared ones if size matters.
    function strokeFadeQuad(ctx, ax, ay, cx, cy, bx, by, strokeA, strokeB, lw) {
        const g = ctx.createLinearGradient(ax, ay, bx, by);
        g.addColorStop(0, strokeA); g.addColorStop(1, strokeB);
        ctx.strokeStyle = g; ctx.lineWidth = lw;
        ctx.beginPath(); ctx.moveTo(ax, ay); ctx.quadraticCurveTo(cx, cy, bx, by); ctx.stroke();
    }
    function strokeFadeCubic(ctx, ax, ay, c1x, c1y, c2x, c2y, bx, by, strokeA, strokeB, lw) {
        const g = ctx.createLinearGradient(ax, ay, bx, by);
        g.addColorStop(0, strokeA); g.addColorStop(1, strokeB);
        ctx.strokeStyle = g; ctx.lineWidth = lw;
        ctx.beginPath(); ctx.moveTo(ax, ay); ctx.bezierCurveTo(c1x, c1y, c2x, c2y, bx, by); ctx.stroke();
    }

    // ── Sequence viewer ──────────────────────────────────────────────────────
    //
    // Renders residues as ONE continuous, reflowing string so that find-in-page
    // (Cmd-F) matches a motif even where it spans a line break, and so the
    // sequence re-wraps to whatever width the card happens to have.
    //
    // The rule that makes this work: nothing but residues may enter the text
    // stream. Element boundaries are fine — browsers match straight through
    // inline elements — but a single stray character breaks a match. So each of
    // the three things a sequence viewer normally interleaves with the residues
    // is moved out of the text:
    //
    //   gap every 10   -> CSS margin on .seq-g   (not a space)
    //   line break     -> <wbr>                  (a zero-width break opportunity)
    //   position label -> .seq-g::before { content: attr(data-n) }  (generated)
    //
    // The predecessor chunked residues into fixed 80-column rows with the start
    // and end numbers as real text nodes, which put digits mid-sequence and made
    // any motif crossing a row boundary unfindable.
    //
    //   residues  [{ aa, resnum }, ...]  (aa may be a multi-letter ion label)
    //   decorate  (res, i) -> { bg, clir, tip } | null
    //             bg   background colour, or null/undefined for none
    //             clir true for the bold/white cLIR treatment
    //             tip  hover text; omit for no tooltip on that residue
    //   opts      { group: 10 }
    //
    // Returns a <div class="seq-flow">. Tooltips are delegated from that div, so
    // a 9-chain complex costs two listeners rather than thousands.
    // Size a sequence host to the largest whole number of residue blocks that fit,
    // then centre it. Measured from the live box rather than computed from font
    // metrics: letter-spacing, the ch unit and the inter-block margin all feed the
    // real block pitch, and guessing any of them puts the grid one block out.
    function fitSeqHost(host) {
        if (!host) return;
        host.style.width = '';
        host.style.marginInline = '';
        const gs = host.querySelectorAll('.seq-g');
        if (gs.length < 2) return;                  // one block: nothing to align to

        // Let it lay out naturally, then read how many blocks the browser actually put
        // on the first line, and the real pitch between two of them. Deriving this from
        // font metrics does not work: blockWidth + marginRight underestimates what a
        // line costs, because Chrome charges the trailing margin-right of the LAST box
        // on a line too. That arithmetic sized the host to fit 10 blocks when only 9
        // fit, leaving exactly enough slack for a short trailing group to squeeze in -
        // the ragged line this whole thing exists to prevent. Measure, do not derive.
        const top0 = gs[0].getBoundingClientRect().top;
        const first = [];
        for (const g of gs) {
            if (Math.abs(g.getBoundingClientRect().top - top0) > 1) break;
            first.push(g);
        }
        const n = first.length;
        if (n < 2) return;
        const pitch = first[1].getBoundingClientRect().left - first[0].getBoundingClientRect().left;
        if (!(pitch > 0)) return;

        const w = Math.ceil(n * pitch);             // n blocks INCLUDING the nth trailing gap
        const avail = host.clientWidth;
        if (w > 0 && w <= avail) { host.style.width = w + 'px'; host.style.marginInline = 'auto'; }
    }

    const _seqObserved = new WeakSet();
    function observeSeqHost(host) {
        if (!host || !host.parentElement || _seqObserved.has(host)) return;
        _seqObserved.add(host);
        // Watch the PARENT: the host's own width is what we set, so observing it
        // would retrigger on our own write and loop.
        const ro = new ResizeObserver(() => requestAnimationFrame(() => fitSeqHost(host)));
        ro.observe(host.parentElement);
    }

    function seqFlow(residues, decorate, opts) {
        opts = opts || {};
        const group = opts.group || 10;
        const div = document.createElement('div');
        div.className = 'seq-flow';
        if (!residues || !residues.length) return div;

        let block = null;
        residues.forEach((res, i) => {
            if (i % group === 0) {
                if (block) div.appendChild(document.createElement('wbr'));
                block = document.createElement('span');
                block.className = 'seq-g';
                // Ruler reads the LAST residue number in the block, so it stays
                // correct for constructs that do not start at 1. A trailing group of
                // one or two residues sits so close to the previous marker that the
                // two numbers collide, and the chain length is already in the label.
                const last = Math.min(i + group, residues.length) - 1;
                if (i === 0 || last - i + 1 >= 3) block.dataset.n = residues[last].resnum;
                div.appendChild(block);
            }
            const d = (decorate && decorate(res, i)) || {};
            const el = document.createElement('span');
            el.className = 'seq-res' + (d.clir ? ' clir' : '');
            if (d.bg) el.style.background = d.bg;
            if (d.tip) el.dataset.t = d.tip;
            el.textContent = res.aa;
            block.appendChild(el);
        });

        // Pad a short trailing group out to a full group's width. Without this it is
        // narrower than the rest and squeezes onto a line that is otherwise full, so
        // the sequence runs past the right margin instead of wrapping - a 205-residue
        // chain put residues 201-205 beyond the edge of an otherwise 10-block grid.
        // The pad is a sibling (so the group's ruler stays over its real last residue)
        // and carries no text, so find-in-page is unaffected. There is no <wbr> before
        // it, so it cannot be separated from the group it belongs to.
        const rem = residues.length % group;
        if (rem) {
            const pad = document.createElement('span');
            pad.className = 'seq-pad';
            pad.style.width = (group - rem) + 'ch';
            div.appendChild(pad);
        }

        // Centre the grid, keep the residues left-aligned inside it. Sizing the host
        // to a whole number of blocks is what makes both true at once: the leftover
        // slack moves to the margins instead of piling up on the right, and every
        // chain keeps the same left edge because they share the host. Sizing each
        // chain to its own content instead would centre a short chain differently
        // from a long one and they would stop lining up.
        requestAnimationFrame(() => fitSeqHost(div.parentElement));
        if (!opts.noObserve && div.parentElement && typeof ResizeObserver !== 'undefined') observeSeqHost(div.parentElement);

        let tip = document.querySelector('.seq-tooltip');
        if (!tip) { tip = document.createElement('div'); tip.className = 'seq-tooltip'; document.body.appendChild(tip); }
        div.addEventListener('mouseover', (e) => {
            const t = e.target.closest && e.target.closest('.seq-res');
            if (!t || !t.dataset.t) { tip.style.display = 'none'; return; }
            tip.textContent = t.dataset.t;
            tip.style.display = '';
            tip.style.left = (e.clientX + 8) + 'px';
            tip.style.top = (e.clientY - 20) + 'px';
        });
        div.addEventListener('mouseleave', () => { tip.style.display = 'none'; });
        return div;
    }

    global.LiviaMaps = {
        setExportOpts, applyExportOpts,
        svgFromDraw, svgFromFixedCanvas,
        downloadSVGFile, downloadSVGFromCanvas, downloadCanvasPNG,
        attachExportBar, attachPngBtnFor, attachChartPngButtons,
        setToRanges, paintChordArcBand, paintChordDomainRing, attachChordDomainHover, drawChordRuler, drawChordArcLabel, chordLabelExtent, segPick,
        parseColorCSV, attachColorUpload,
        seqFlow, fitSeqHost,
        strokeFadeQuad, strokeFadeCubic,
    };
})(window);
