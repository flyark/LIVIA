/**
 * livia-reported.js — the "Reported Interactions" card: what MIST and BioGRID already report for a
 * complex's subunit pairs, merged into ONE set of per-type heatmaps (a paper reported by both counts once).
 *
 * Shared by universal (any number of subunits), dimer, FlyPredictome and Ortholog Predictome (one pair).
 * Species-independent: MIST is keyed by NCBI Gene ID (one get_interactions call per subunit, through the
 * page's fetchTextViaProxy); BioGRID is asked through LIVIA's own proxy Worker, which holds the access key
 * server-side (tools/cloudflare-worker/worker.js) — the browser only ever sends Gene IDs and taxon ids.
 *
 * A subunit is { key, sym, acc?, geneId?, taxId?, copies? }:
 *   key     identifies it in pair keys ("keyA|keyB", sorted) — the UniProt accession where there is one
 *   acc     UniProt accession; when geneId/taxId are missing they are looked up from it (one UniProt call)
 *   copies  ≥ 2 marks a homomeric subunit: its self-pair is a real contact, drawn on the heatmap diagonal
 *
 * Needs livia-core.js (fetchTextViaProxy). Exposes window.LiviaReported.
 */
(function () {
'use strict';

const TYPES = [['ppi', 'PPI'], ['interolog', 'Interolog PPI'], ['genetic', 'Genetic (GI)'], ['interolog-genetic', 'Interolog GI']];
const BIOGRID_PROXY = 'https://livia-proxy.flyark.workers.dev/?biogrid=1';
const esc = (s) => String(s == null ? '' : s).replace(/[&<>"]/g, c => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' }[c]));
const pairKey = (a, b) => [a.key, b.key].sort().join('|');

// MIST answers JSON, sometimes behind an HTML preamble ending in <br>.
function fpJson(text){ if (!text) return null; const m = text.match(/<br>\s*([\[{])/); const i = m ? text.indexOf(m[1], m.index) : text.search(/[\[{]/); if (i < 0) return null; try { return JSON.parse(text.slice(i)); } catch (e) { return null; } }
async function mistInteractions(geneId){   // all MIST interactions for one gene id (or [])
    const resp = await fetchTextViaProxy('https://www.flyrnai.org/tools/mist/api/get_interactions/' + geneId + '/ppi,interolog,genetic,interolog-genetic');
    const data = fpJson(resp); return (data && data.mist_data && data.mist_data.results) || [];
}
async function accToXrefs(acc){   // NCBI Gene ID (for MIST) + BioGRID id/interactor-count + taxon id (for the pairwise check) — one combined UniProt call
    try { const d = await (await fetch('https://rest.uniprot.org/uniprotkb/' + encodeURIComponent(acc) + '?fields=xref_geneid,xref_biogrid,organism_id&format=json')).json();
          const xrefs = d.uniProtKBCrossReferences || [];
          const gid = xrefs.find(r => r.database === 'GeneID');
          const bg = xrefs.find(r => r.database === 'BioGRID');
          const bgN = bg && (bg.properties || []).find(p => p.key === 'Interactions');
          return { geneId: gid ? gid.id : null, bioGridId: bg ? bg.id : null, bioGridN: bgN ? bgN.value : null, taxId: (d.organism && d.organism.taxonId) || null }; }
    catch (e){ return { geneId: null, bioGridId: null, bioGridN: null, taxId: null }; }
}
// Which of these subunit pairs does BioGRID report as interacting? Gene IDs rather than UniProt's BioGRID
// links: only reviewed UniProt entries link to BioGRID, so most subunits of e.g. a fly complex would drop
// out. The Worker needs ≥2 distinct genes and leaves self-interactions out. See worker.js handleBiogrid().
async function biogridPairCheck(subs){
    const ids = [...new Set(subs.map(s => s.geneId).filter(Boolean).map(String))];
    const taxa = [...new Set(subs.map(s => s.taxId).filter(Boolean).map(String))];
    if (ids.length < 2 || !taxa.length) return null;
    try {
        const resp = await fetch(BIOGRID_PROXY + '&ids=' + encodeURIComponent(ids.join(',')) + '&taxId=' + encodeURIComponent(taxa.join('|')));
        const data = await resp.json();
        if (!resp.ok || data.error) return { error: data.error || ('proxy returned ' + resp.status) };
        return { pairs: new Map(Object.entries(data.pairs || {})), biogridIds: data.biogridIds || {} };
    } catch (e){ return { error: String(e.message || e) }; }
}

let uid = 0;
// opts: { types, source, intro, introExtra, cellNote } — defaults render MIST alone; query() passes the
// merged MIST + BioGRID map with its source name, an intro sentence, and a per-source cell tooltip.
function render(withId, subs, pairMap, opts){
    const types = (opts && opts.types) || TYPES, source = (opts && opts.source) || 'MIST', intro = !(opts && opts.intro === false);
    const introExtra = (opts && opts.introExtra) || '', cellNote = (opts && opts.cellNote) || null;
    const pmidLink = (p) => '<a href="https://pubmed.ncbi.nlm.nih.gov/' + p + '/" target="_blank" style="color:#2471A3;">' + p + '</a>';
    const pubmedSet = (ids) => 'https://pubmed.ncbi.nlm.nih.gov/?term=' + ids.join(',');   // comma-joined PMIDs → PubMed shows exactly that set
    const pmidRefs = (list) => {                                        // inline list; ">5" collapses behind a "+N more" toggle so every PMID is viewable
        if (list.length <= 5) return list.map(pmidLink).join(', ');
        const id = 'mpm' + (uid++);
        const n = list.length - 5;
        return list.slice(0, 5).map(pmidLink).join(', ')
            + '<span id="' + id + '" style="display:none;">, ' + list.slice(5).map(pmidLink).join(', ') + '</span>'
            + ' <a href="javascript:void(0)" onclick="var e=document.getElementById(\'' + id + '\'),o=e.style.display===\'none\';e.style.display=o?\'\':\'none\';this.textContent=o?\' less\':\' +' + n + ' more\'" style="color:#888; font-size:0.72rem;">+' + n + ' more</a>';
    };
    const warn = (msg) => '<div style="background:#fff3cd; border-left:4px solid #ffc107; padding:0.6rem 0.8rem; border-radius:0 6px 6px 0; font-size:0.85rem; color:#856404;">' + msg + '</div>';
    const selfOk = (s) => (s.copies || 1) >= 2;                         // a homomeric subunit's self-pair is a real contact
    if (withId.length === 1 && !selfOk(withId[0])){                    // single resolved protein, one copy → its reported-partner count
        const s = withId[0], sg = String(s.geneId);
        const partners = new Set((s.results || []).map(r => { const ga = String(r.GeneA), gb = String(r.GeneB); return ga === sg ? gb : (gb === sg ? ga : null); }).filter(g => g && g !== sg));
        const note = subs.length > 1 ? ' <span style="color:#999;">(other chains not resolved to a Gene ID)</span>' : '';
        return '<div style="font-size:0.85rem; color:#444;"><strong>' + esc(s.sym) + '</strong> has <strong>' + partners.size + '</strong> reported interaction partner' + (partners.size === 1 ? '' : 's') + ' in MIST.' + note + ' Single resolved protein — browse them via the link below.</div>';
    }
    const pairs = [];
    for (let i = 0; i < withId.length; i++) for (let j = i; j < withId.length; j++){
        const a = withId[i], b = withId[j];
        if (i === j && !selfOk(a)) continue;
        pairs.push({ a, b, byType: pairMap.get(pairKey(a, b)) || null });
    }
    const pairName = (p) => p.a === p.b ? esc(p.a.sym) + ' (self)' : esc(p.a.sym) + '–' + esc(p.b.sym);
    const reported = pairs.filter(p => p.byType && Object.keys(p.byType).length);
    const unreported = pairs.filter(p => !(p.byType && Object.keys(p.byType).length));
    if (!reported.length) return warn(pairs.length === 1
        ? (pairs[0].a === pairs[0].b
            ? '<strong>No reported self-interaction</strong> for ' + esc(pairs[0].a.sym) + ' in ' + source + '.'
            : '<strong>No reported interaction</strong> between ' + esc(pairs[0].a.sym) + ' and ' + esc(pairs[0].b.sym) + ' in ' + source + '.')
        : 'None of the ' + pairs.length + ' subunit pairs is reported in ' + source + '.');
    // Four heatmaps — one per interaction type. N×N subunit grid; each cell = PMID count for that
    // pair, colored by magnitude (number shown), clickable → PubMed. Scales to large complexes
    // far better than a per-pair row list (a 14-mer is 91 rows but a 14×14 grid).
    const n = withId.length;
    const cnt = (a, b, t) => { const bt = pairMap.get(pairKey(a, b)); const s = bt && bt[t]; return s ? s.size : 0; };
    const idsOf = (a, b, t) => { const bt = pairMap.get(pairKey(a, b)); const s = bt && bt[t]; return s ? [...s] : []; };
    const maxOf = {};                                                  // per-type max, for color normalisation
    for (const [t] of types){ let m = 0; for (const p of pairs) m = Math.max(m, cnt(p.a, p.b, t)); maxOf[t] = m; }
    const cs = n <= 6 ? 30 : (n <= 10 ? 24 : 19);                      // cell size (px), adaptive to subunit count
    const fs = n <= 10 ? '0.72rem' : '0.64rem';
    const maxLen = withId.reduce((m, s) => Math.max(m, String(s.sym).length), 0);
    const hdrH = Math.min(96, 18 + Math.round(maxLen * 4.6));          // header height fits the rotated labels
    const heat = (c, max) => { if (!c) return 'background:#fbfcfd; color:#d0d0d0;'; const t = max ? c / max : 1; return 'background:rgba(36,113,163,' + (0.14 + t * 0.78).toFixed(3) + '); color:' + (t > 0.5 ? '#fff' : '#14324a') + ';'; };
    function heatmap(t, label){
        const max = maxOf[t];
        if (!max) return '';                                           // no evidence of this type → omit its heatmap
        let h = '<div data-mist-hm data-n="' + n + '" style="margin:0 ' + Math.round(hdrH * 0.85) + 'px 0.7rem 0;"><div style="font-size:0.8rem; font-weight:700; color:#2471A3; margin-bottom:0.3rem;">' + esc(label) + ' <span style="color:#aaa; font-weight:400;">(max ' + max + ')</span></div>';
        h += '<table style="border-collapse:collapse; font-size:' + fs + ';"><thead><tr><th style="padding:0;"></th>';
        for (let j = 0; j < n; j++) h += '<th style="padding:0; vertical-align:bottom; height:' + hdrH + 'px;"><div style="width:var(--cs); height:' + hdrH + 'px; position:relative;"><span style="position:absolute; left:50%; bottom:2px; transform-origin:left bottom; transform:rotate(-45deg); white-space:nowrap; font-weight:600; color:#444;">' + esc(withId[j].sym) + '</span></div></th>';
        h += '</tr></thead><tbody>';
        for (let i = 0; i < n; i++){
            h += '<tr><th style="padding:1px 6px 1px 0; text-align:right; white-space:nowrap; font-weight:600; color:#444;">' + esc(withId[i].sym) + '</th>';
            for (let j = 0; j < n; j++){
                const a = withId[i], b = withId[j];
                if (i === j && !selfOk(a)){ h += '<td style="width:var(--cs); height:var(--cs); background:#eceff2; border:1px solid #fff;"></td>'; continue; }
                const c = cnt(a, b, t);
                const who = i === j ? esc(a.sym) + ' (self)' : esc(a.sym) + ' – ' + esc(b.sym);
                if (c) h += '<td style="padding:0; border:1px solid #fff;"><a href="' + pubmedSet(idsOf(a, b, t)) + '" target="_blank" title="' + who + ': ' + c + ' paper' + (c > 1 ? 's' : '') + (cellNote ? cellNote(a, b, t) : '') + ' — open in PubMed" style="display:flex; align-items:center; justify-content:center; width:var(--cs); height:var(--cs); font-size:min(0.72rem, calc(var(--cs) * 0.5)); text-decoration:none; font-weight:700; ' + heat(c, max) + '">' + c + '</a></td>';
                else h += '<td style="width:var(--cs); height:var(--cs); border:1px solid #fff; text-align:center; ' + heat(0, max) + '">·</td>';
            }
            h += '</tr>';
        }
        return h + '</tbody></table></div>';
    }
    // overflow:hidden (not auto) — we never want a horizontal scrollbar; layout() shrinks --cs so even a
    // single very wide heatmap fits the card. flex-wrap picks 4-up / 2×2 by width.
    // One statement per line: what a cell is, where the numbers come from, then what to do with it.
    let html = intro ? '<div style="font-size:0.78rem; color:#777; margin:0 0 0.5rem; line-height:1.55;">'
        + '<div>Each cell is the number of reported papers for that subunit pair (darker = more)' + (pairs.some(p => p.a === p.b) ? '; the diagonal is a subunit with itself.' : '.') + '</div>'
        + (introExtra ? '<div>' + introExtra.trim() + '</div>' : '')
        + '<div><b style="color:#2471A3;">Click a cell to open those papers in PubMed</b>, or expand “Show PMIDs” below to copy the IDs.</div></div>' : '';
    html += '<div data-mist-wrap data-natcs="' + cs + '" style="display:flex; flex-wrap:wrap; align-items:flex-start; overflow:hidden; --cs:' + cs + 'px;">' + types.map(([t, label]) => heatmap(t, label)).join('') + '</div>';
    // PMIDs — collapsed by default so the card stays compact; click a count above to open all papers in
    // PubMed, or expand here to read/copy the raw PMIDs (each links to its paper; Copy grabs them all).
    html += '<details style="margin-top:0.5rem;"><summary style="cursor:pointer; color:#2471A3; font-size:0.8rem;">Show PMIDs</summary><div style="margin-top:0.3rem; font-size:0.78rem; color:#555;">';
    const copyLink = (ids) => '<a href="javascript:void(0)" onclick="navigator.clipboard&&navigator.clipboard.writeText(\'' + ids.join(', ') + '\');this.textContent=\'copied\';setTimeout(()=>this.textContent=\'copy\',1200)" style="color:#888; font-size:0.72rem;">copy</a>';
    for (const p of reported){
        html += '<div style="margin:0.4rem 0;"><b>' + (p.a === p.b ? esc(p.a.sym) + ' (self)' : esc(p.a.sym) + ' — ' + esc(p.b.sym)) + '</b>';
        for (const [t, label] of types){ const set = p.byType[t]; if (!set || !set.size) continue; const list = [...set];   // one line per evidence category: "copy — Label: pmids", each separately copyable
            html += '<div style="margin:0.12rem 0 0.12rem 1rem;">' + copyLink(list) + ' — <a href="' + pubmedSet(list) + '" target="_blank" style="color:#2471A3;">' + label + '</a>: ' + pmidRefs(list) + '</div>';
        }
        html += '</div>';
    }
    html += '</div></details>';
    if (unreported.length) html += '<div style="margin-top:0.4rem; font-size:0.78rem; color:#999;">Not reported in ' + source + ': ' + unreported.map(pairName).join(', ') + '</div>';
    return html;
}

// Fit the heatmaps to the card width without ever producing a horizontal scrollbar. flex-wrap already
// lays them out 4-up (wide) or 2×2 (narrow); this only kicks in for the rare case where a SINGLE heatmap
// is still wider than the card (very large complexes) — then it shrinks the shared --cs cell size just
// enough to fit, instead of scrolling or clipping.
function fitWrap(wrap){
    const boxes = [...wrap.querySelectorAll('[data-mist-hm]')]; if (!boxes.length) return;
    const avail = wrap.clientWidth; if (avail < 60) return;                 // card hidden / not laid out yet
    const natCs = parseFloat(wrap.dataset.natcs) || 19;
    const nSub = parseInt(boxes[0].dataset.n, 10) || 1;
    wrap.style.setProperty('--cs', natCs + 'px');                           // measure the natural size first
    const mr = parseFloat(getComputedStyle(boxes[0]).marginRight) || 0;
    let widest = 0; for (const b of boxes) widest = Math.max(widest, b.getBoundingClientRect().width);
    if (widest + mr <= avail) return;                                       // natural fits → flex-wrap handles columns
    const fixed = widest - nSub * natCs;                                    // row-label column + rotated-header overhang + borders
    let cs = Math.floor((avail - fixed - mr - 4) / nSub);
    cs = Math.max(9, Math.min(natCs, cs));                                  // never below 9px, never above natural
    wrap.style.setProperty('--cs', cs + 'px');
}
function layout(){ document.querySelectorAll('#mist-card [data-mist-wrap]').forEach(fitWrap); }
let resizeBound = false;
function bindResize(){
    if (resizeBound) return; resizeBound = true;
    let t; window.addEventListener('resize', () => { clearTimeout(t); t = setTimeout(layout, 120); });
}

// Compact one-liner for a pair shown elsewhere on the page (universal's Linear Contact Map): per-type
// counts, each linking to PubMed, and a jump up to the card, which holds the full copyable PMID lists.
function pairHtml(byType, sources){
    if (!byType || !Object.keys(byType).length) return '';
    const pubmedSet = (ids) => 'https://pubmed.ncbi.nlm.nih.gov/?term=' + ids.join(',');
    const parts = [];
    for (const [t, label] of TYPES){ const set = byType[t]; if (set && set.size){ const list = [...set];
        parts.push(esc(label) + ' <a href="' + pubmedSet(list) + '" target="_blank" title="Open all ' + list.length + ' papers in PubMed" style="color:#2471A3; font-weight:700; text-decoration:none;">' + list.length + '</a>'); } }
    return '<div style="margin-top:6px; font-size:0.78rem; color:#555; background:#f6f9fc; border-left:3px solid #2471A3; padding:6px 10px; border-radius:0 5px 5px 0;">'
        + '<b style="color:#1a5276;">Reported in ' + (sources || 'MIST').split(' + ').map(s => s === 'BioGRID'
            ? '<a href="https://thebiogrid.org/" target="_blank" style="color:#1a5276;">BioGRID</a>'
            : '<a href="https://fgrtools.hms.harvard.edu/MIST/" target="_blank" style="color:#1a5276;">MIST</a>').join(' + ') + '</b>: '
        + parts.join(' &middot; ')
        + ' &nbsp;<a href="javascript:void(0)" onclick="var c=document.getElementById(\'mist-card\'); if(c){ c.style.display=\'\'; c.scrollIntoView({behavior:\'smooth\', block:\'center\'}); }" style="color:#2471A3; font-size:0.74rem; white-space:nowrap;">full PMIDs &uarr;</a></div>';
}

// Fill the card for these subunits: MIST (one call per subunit) and BioGRID (one call for all pairs)
// answer independently; the heatmaps are drawn when MIST settles and redrawn if BioGRID lands later.
// o: { subs, content, linkEl?, isStale?() → true once the page moved on, onDrawn?(merged, sourcesLabel) }
async function query(o){
    const content = o.content, linkEl = o.linkEl || null, subs = o.subs || [];
    const stale = o.isStale || (() => false);
    if (!content || !subs.length) return;
    if (linkEl) linkEl.innerHTML = '';
    content.innerHTML = '<span style="color:#888; font-size:0.85rem;">Querying MIST for ' + subs.map(s => esc(s.sym)).join(', ') + '…</span>';
    await Promise.all(subs.map(async s => {                            // accession → NCBI Gene ID (MIST) + BioGRID xref + taxon
        if (!s.acc) return;
        const x = await accToXrefs(s.acc);
        s.geneId = s.geneId || x.geneId; s.taxId = s.taxId || x.taxId; s.bioGridId = x.bioGridId; s.bioGridN = x.bioGridN;
    }));
    if (stale()) return;
    const withId = subs.filter(s => s.geneId);
    // BioGRID's physical evidence joins MIST's PPI, its genetic evidence MIST's Genetic (GI).
    const rep = { mist: null, mistError: null, mistDone: false, bg: null };
    const drawReported = () => {
        if (stale() || !rep.mistDone) return;
        if (!rep.mist && !rep.bg){ content.innerHTML = '<span style="color:#999; font-size:0.85rem;">Could not query MIST: ' + esc(rep.mistError || 'no answer') + '</span>'; return; }
        const merged = new Map(), src = new Map();                     // src: pair → type → { MIST: n, BioGRID: n } for the cell tooltips
        const add = (map, name) => { if (!map) return; for (const [key, byType] of map){
            const m = merged.get(key) || {}, s = src.get(key) || {}; merged.set(key, m); src.set(key, s);
            for (const [t, set] of Object.entries(byType)){ const u = m[t] || (m[t] = new Set()); set.forEach(p => u.add(p)); (s[t] || (s[t] = {}))[name] = set.size; } } };
        add(rep.mist, 'MIST'); add(rep.bg, 'BioGRID');
        const names = rep.bg ? (rep.mist ? 'MIST + BioGRID' : 'BioGRID') : 'MIST';
        content.innerHTML = render(withId, subs, merged, {
            source: rep.bg && rep.mist ? 'MIST or BioGRID' : names,
            introExtra: rep.bg ? (rep.mist ? ' MIST and BioGRID are combined; a paper reported by both counts once.' : ' From BioGRID (MIST did not answer).') : '',
            cellNote: rep.bg ? (a, b, t) => { const s = (src.get(pairKey(a, b)) || {})[t] || {}; return ' (' + Object.entries(s).map(([k, n]) => k + ' ' + n).join(' · ') + ')'; } : null,
        });
        bindResize(); requestAnimationFrame(layout);                    // fit heatmaps to card width, no h-scroll
        if (o.onDrawn) o.onDrawn(merged, names);
    };
    if (linkEl){                                                        // BioGRID cross-reference — shown independent of MIST outcome, since MIST's coverage skews toward a few reference organisms
        // A subunit is linked when UniProt cross-references its BioGRID entry (with UniProt's interactor
        // count) or when the Gene-ID pair check below finds it in BioGRID.
        const bgLine = (found) => {
            const shown = subs.filter(s => s.bioGridId || (s.geneId && found[s.geneId]));
            return shown.length ? 'Also on <a href="https://thebiogrid.org/" target="_blank" style="color:#2471A3;">BioGRID</a>: '
                + shown.map(s => '<a href="https://thebiogrid.org/' + encodeURIComponent(s.bioGridId || found[s.geneId]) + '" target="_blank" style="color:#2471A3;">' + esc(s.sym) + '</a>'
                    + (s.bioGridN != null ? ' (' + s.bioGridN + ' interactor' + (s.bioGridN == 1 ? '' : 's') + ')' : '')).join(', ') : '';
        };
        linkEl.innerHTML = '<div id="bg-also" style="color:#888;">' + bgLine({}) + '</div>';
        if (new Set(withId.map(s => String(s.geneId))).size >= 2){
            linkEl.innerHTML += '<div id="bg-pairs" style="color:#888; margin-top:2px;">Checking BioGRID for reported pairs among these subunits…</div>';
            biogridPairCheck(withId).then(res => {
                if (stale()) return;
                const el = document.getElementById('bg-pairs'); if (el) el.remove();
                if (!res || res.error) return;                          // proxy unavailable/unconfigured — fail quiet, MIST and the UniProt link-outs still stand
                const also = document.getElementById('bg-also'); if (also) also.innerHTML = bgLine(res.biogridIds);
                const bgMap = new Map();                                // same keys and type names as MIST's map, so the two merge cell for cell
                for (let i = 0; i < withId.length; i++) for (let j = i + 1; j < withId.length; j++){
                    const a = withId[i], b = withId[j];
                    const p = res.pairs.get([String(a.geneId), String(b.geneId)].sort((x, y) => x - y).join('|'));
                    if (!p) continue;
                    const byType = { ppi: new Set(p.pmidsPhysical || []), genetic: new Set(p.pmidsGenetic || []) };
                    for (const t of Object.keys(byType)) if (!byType[t].size) delete byType[t];
                    if (Object.keys(byType).length) bgMap.set(pairKey(a, b), byType);
                }
                rep.bg = bgMap; drawReported();
            });
        }
    }
    if (!withId.length){ content.innerHTML = '<span style="color:#999; font-size:0.85rem;">No NCBI Gene ID for these chains — cannot query MIST.</span>'; return; }
    try {
        await Promise.all(withId.map(async s => { s.results = await mistInteractions(s.geneId); }));   // one query per subunit
        if (stale()) return;
        const gidToSub = new Map(); withId.forEach(s => gidToSub.set(String(s.geneId), s));
        const pairMap = new Map();                                     // pairKey → { type → Set(pmids) }
        for (const s of withId){
            const sg = String(s.geneId);
            for (const r of (s.results || [])){
                const ga = String(r.GeneA), gb = String(r.GeneB);
                if (ga !== sg && gb !== sg) continue;                  // row must actually involve this subunit
                const partner = gidToSub.get(ga === sg ? gb : ga);
                if (!partner) continue;                                // only contacts between subunits of THIS complex
                if (partner === s && !((s.copies || 1) >= 2)) continue;   // a self-interaction counts only for a homomeric subunit
                const key = pairKey(s, partner);
                if (!pairMap.has(key)) pairMap.set(key, {});
                const bt = pairMap.get(key); const t = r.Interaction_type; if (!bt[t]) bt[t] = new Set();
                if (r.Reference) String(r.Reference).split(';').map(p => p.trim()).filter(Boolean).forEach(p => bt[t].add(p));
            }
        }
        rep.mist = pairMap;
        // (no per-protein deep link: MIST is a JS app with no stable gene-page URL. The header links to
        //  the MIST tool, and each count links to PubMed for that pair's papers.)
    } catch (e){ rep.mistError = e.message; }
    if (stale()) return;
    rep.mistDone = true; drawReported();
}

window.LiviaReported = { TYPES, esc, query, render, layout, pairHtml, pairKey, accToXrefs, mistInteractions, biogridPairCheck, fpJson };
})();
