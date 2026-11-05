// GoatCounter counts one visit per page per session (8 hours), so a reader opening fifty models from the Atlas or a
// collaborator's data site is one visit to universal.html. Each arrival from another site also counts as an event named by
// the page and the site it came from, every time: 'from/livia-atlas/universal'. Never the file, pair or query.
(function () {
    function siteOf(url) {
        const u = new URL(url);
        if (/^(localhost|127\.|\[::1\])/.test(u.hostname)) return 'offline';
        return u.hostname === 'flyark.github.io' ? (u.pathname.split('/')[1] || u.hostname) : u.hostname;
    }
    function arrivalSource() {
        const nav = performance.getEntriesByType('navigation')[0];
        if (nav && nav.type !== 'navigate') return '';   // a reload or the back button is not another arrival
        const here = siteOf(location.href), data = new URLSearchParams(location.search).get('data');
        try {
            if (document.referrer) {   // the Atlas and the data sites share this origin, so their referrer keeps the path that names them
                const s = siteOf(document.referrer), own = new URL(document.referrer).hostname === location.hostname;
                if (s === here) return '';
                return own || location.search ? s : '';   // a plain page reached from a search engine is an ordinary visit
            }
            if (!location.search) return '';
            return data ? siteOf(new URL(data, location.href).href) : 'direct';   // a link with no referrer: the site its data comes from
        } catch (e) { return location.search ? 'direct' : ''; }
    }
    window.addEventListener('load', () => setTimeout(() => {
        const source = arrivalSource();
        if (!source || !(window.goatcounter && window.goatcounter.count)) return;
        const page = (location.pathname.split('/').pop() || 'index.html').replace(/\.html$/, '');
        window.goatcounter.count({ path: `from/${source}/${page}`, title: `${page}.html opened from ${source}`, event: true, no_session: true });
    }, 0), { once: true });
})();
