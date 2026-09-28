/**
 * livia-subnav.js — a sticky "jump to section" bar for LIVIA's long report pages.
 *
 * One button per visible report card: the card's <h2>, labeled by its data-nav attribute (data-nav="" leaves the
 * card out) or else by the heading's own leading text. The bar rebuilds itself as cards appear or hide, stays
 * pinned to the top of the window while the report scrolls, and marks the section currently in view. It shows
 * only once the report has at least three visible sections.
 *
 * Usage: <script src="js/livia-subnav.js"></script> and a data-subnav attribute on the report container
 * (<div id="results" data-subnav>). It mounts itself; LiviaSubnav.mount(el) mounts one by hand.
 *
 * News: LiviaSubnav.notify(el, badge, title) marks the section holding `el` as having new information (a chain
 * identified, reported interactions found). Its button pulses three times and keeps a small badge, which stays while the
 * reader scrolls or opens the section: it tells what the section holds. Reduced-motion readers get the badge without the pulse.
 * Dependencies: none (injects its own CSS).
 */
(function (root) {
    'use strict';
    const CSS = [
        '.livia-subnav { position: sticky; top: 0; z-index: 900; display: flex; flex-wrap: wrap; justify-content: center; gap: 2px 2px;',   // tight enough for 13 sections on one line in the 1052px page
        '  margin: 0 0 12px; padding: 6px 8px; background: rgba(248, 250, 252, 0.94); -webkit-backdrop-filter: blur(8px); backdrop-filter: blur(8px);',
        '  border: 1px solid #dde3ea; border-radius: 10px; box-shadow: 0 2px 10px rgba(26, 82, 118, 0.07); }',
        '.livia-subnav[hidden] { display: none; }',
        '.livia-subnav button { all: unset; cursor: pointer; font-size: 12.5px; font-weight: 500; line-height: 1.2; color: #4a5a6a; padding: 5px 6px; border-radius: 7px; white-space: nowrap; }',
        '.livia-subnav button:hover { background: #e8f1f8; color: #2471a3; }',
        '.livia-subnav button.on { background: #1a5276; color: #fff; }',
        '.livia-subnav button:focus-visible { outline: 2px solid #2471a3; outline-offset: 1px; }',
        '[data-subnav] .card, [data-subnav] h2 { scroll-margin-top: 64px; }   /* scrollIntoView leaves room for the pinned bar */',
        '.livia-subnav button.news { position: relative; color: #7a3208; background: #fdf0e6; }',
        '.livia-subnav button.news.on { color: #fff; background: #1a5276; }',
        // a corner badge, like an app's: it takes no room in the row, so the bar stays on one line
        '.livia-subnav button.news::after { content: attr(data-badge); position: absolute; top: -5px; right: -4px; padding: 0 4px; min-width: 7px;',
        '  border-radius: 7px; background: #d35400; color: #fff; font-size: 9.5px; font-weight: 700; line-height: 13px; text-align: center; pointer-events: none; }',
        '.livia-subnav button.pulse { animation: livia-news 0.9s ease-out 3; }',
        '@keyframes livia-news { 0% { box-shadow: 0 0 0 0 rgba(211, 84, 0, 0.6); } 100% { box-shadow: 0 0 0 10px rgba(211, 84, 0, 0); } }',
        '@media (prefers-reduced-motion: reduce) { .livia-subnav button.pulse { animation: none; } }',
        '@media (max-width: 720px) { .livia-subnav { flex-wrap: nowrap; justify-content: flex-start; overflow-x: auto; } }',
        '@media print { .livia-subnav { display: none !important; } }',
    ].join('\n');

    const visible = (el) => el.getClientRects().length > 0;
    function labelOf(h) {
        if (h.hasAttribute('data-nav')) return h.getAttribute('data-nav').trim();
        let t = '';
        for (const n of h.childNodes) { if (n.nodeType === 3) t += n.textContent; else if (t.trim()) break; }
        return t.replace(/\s+/g, ' ').trim();
    }
    const esc = (s) => s.replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));

    function mount(box, opts) {
        box = typeof box === 'string' ? document.getElementById(box) : box;
        if (!box || box.__liviaSubnav) return box && box.__liviaSubnav;
        opts = opts || {};
        if (!document.getElementById('livia-subnav-css')) {
            const st = document.createElement('style'); st.id = 'livia-subnav-css'; st.textContent = CSS; document.head.appendChild(st);
        }
        const bar = document.createElement('nav');
        bar.className = 'livia-subnav'; bar.setAttribute('aria-label', 'Report sections'); bar.hidden = true;
        box.insertBefore(bar, box.firstChild);
        let items = [], key = null, cur = -2, raf = 0, spyRaf = 0;
        const target = (it) => it.h.closest('.card') || it.h;
        const news = new Map(), toPulse = new Set();                     // section (card) → { badge, title }; sections still to pulse
        function paintNews() {
            items.forEach((it, i) => {
                const b = bar.children[i], card = target(it); if (!b) return;
                const has = news.has(card);
                b.classList.toggle('news', has);
                if (has) { b.dataset.badge = news.get(card).badge; b.title = news.get(card).title; } else { delete b.dataset.badge; b.removeAttribute('title'); }
                if (has && toPulse.has(card) && !bar.hidden) {
                    toPulse.delete(card); b.classList.remove('pulse'); void b.offsetWidth; b.classList.add('pulse');
                    b.addEventListener('animationend', () => b.classList.remove('pulse'), { once: true });
                }
            });
        }

        function spy() {
            spyRaf = 0;
            if (bar.hidden) return;
            const lim = bar.getBoundingClientRect().bottom + 24;
            let c = -1;
            items.forEach((it, i) => { if (target(it).getBoundingClientRect().top <= lim) c = i; });
            if (c === cur) return;
            cur = c;
            Array.prototype.forEach.call(bar.children, (b, i) => b.classList.toggle('on', i === c));
        }
        function rebuild() {
            raf = 0;
            const hs = Array.prototype.filter.call(box.querySelectorAll('h2'), (h) => !bar.contains(h) && visible(h) && labelOf(h));
            const k = hs.map(labelOf).join('\u0001');
            const show = hs.length >= (opts.min || 3);
            if (k !== key) {
                key = k; cur = -2;
                items = hs.map((h) => ({ h, label: labelOf(h) }));
                bar.innerHTML = items.map((it, i) => '<button type="button" data-i="' + i + '">' + esc(it.label) + '</button>').join('');
            }
            if (bar.hidden === show) bar.hidden = !show;
            paintNews();
            spy();
        }
        bar.addEventListener('click', (e) => {
            const b = e.target.closest('button'); if (!b) return;
            const it = items[+b.dataset.i]; if (!it) return;
            const y = target(it).getBoundingClientRect().top + window.scrollY - bar.offsetHeight - 12;
            window.scrollTo({ top: Math.max(0, y), behavior: 'auto' });
        });
        // Cards appear, hide and retitle as a report fills in; rebuild on the next frame, ignoring the bar's own edits.
        new MutationObserver((muts) => {
            if (muts.every((m) => bar.contains(m.target))) return;
            if (!raf) raf = requestAnimationFrame(rebuild);
        }).observe(box, { subtree: true, childList: true, characterData: true, attributes: true, attributeFilter: ['style', 'class', 'hidden'] });
        window.addEventListener('scroll', () => { if (!spyRaf) spyRaf = requestAnimationFrame(spy); }, { passive: true });
        window.addEventListener('resize', () => { if (!raf) raf = requestAnimationFrame(rebuild); });
        function notify(el, badge, title) {
            const card = (el && el.closest && el.closest('.card')) || el; if (!card) return;
            const b = badge == null || badge === '' ? 'new' : String(badge);
            news.set(card, { badge: b, title: 'New: ' + (title || b) }); toPulse.add(card);
            if (!raf) raf = requestAnimationFrame(rebuild);                 // the card may only now be showing
        }
        rebuild();
        return (box.__liviaSubnav = { bar, rebuild, notify });
    }

    function auto() { Array.prototype.forEach.call(document.querySelectorAll('[data-subnav]'), (el) => mount(el)); }
    if (typeof document !== 'undefined') {
        if (document.readyState === 'loading') document.addEventListener('DOMContentLoaded', auto); else auto();
    }
    function notify(el, badge, title) {                                   // the bar of the report holding `el`
        const host = el && el.closest && el.closest('[data-subnav]');
        const nav = host && (host.__liviaSubnav || mount(host));
        if (nav) nav.notify(el, badge, title);
    }
    root.LiviaSubnav = { mount, notify };
})(typeof window !== 'undefined' ? window : this);
