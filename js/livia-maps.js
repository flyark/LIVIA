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

    // Idempotent "↓ SVG · ↓ PNG" bar directly under a canvas.
    function attachExportBar(canvasId, opts) {
        opts = opts || {};
        const cv = document.getElementById(canvasId);
        if (!cv || !cv.parentElement) return;
        const barId = canvasId + '-export-bar';
        const existing = document.getElementById(barId);
        if (existing) existing.remove();          // rebuild so the redraw closure stays current
        const bar = document.createElement('div');
        bar.id = barId;
        bar.style.cssText = 'text-align:center; margin-top:3px;';
        const mk = (label, fn) => {
            const b = document.createElement('button');
            b.type = 'button'; b.textContent = label;
            b.style.cssText = 'font-size:0.72rem; color:#2471A3; background:none; border:none; cursor:pointer; font-weight:600; padding:0 6px;';
            b.onclick = fn;
            bar.appendChild(b);
        };
        if (opts.svg) mk('↓ SVG', () => downloadSVGFromCanvas(canvasId, opts.name, opts.svg));
        if (opts.png) mk('↓ PNG', () => downloadCanvasPNG(canvasId, opts.name));
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

    global.LiviaMaps = {
        setExportOpts, applyExportOpts,
        svgFromDraw, svgFromFixedCanvas,
        downloadSVGFile, downloadSVGFromCanvas, downloadCanvasPNG,
        attachExportBar, attachPngBtnFor, attachChartPngButtons,
        setToRanges, paintChordArcBand,
    };
})(window);
