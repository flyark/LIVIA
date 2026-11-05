// Helpers shared by the FlyPredictome and Ortholog Interactome pages, the same code on both before they moved here.

// pDockQ: faithful original (Bryant et al. 2022) -- pure Cb-contact geometric interface, no PAE at
// all, freshly computed from the structure (O(lenA*lenB) pairwise distance, trivial at these sizes).
function parsePdbCbCoordsAB(pdbText) {
    const out = { A: [], B: [] };
    const seen = { A: new Set(), B: new Set() };
    for (const line of pdbText.split('\n')) {
        if (!line.startsWith('ATOM') || line.length < 54) continue;
        const chain = line.substring(21, 22).trim();
        if (chain !== 'A' && chain !== 'B') continue;
        const atomName = line.substring(12, 16).trim();
        const resName = line.substring(17, 20).trim();
        const resnum = parseInt(line.substring(22, 26));
        if (isNaN(resnum)) continue;
        const x = parseFloat(line.substring(30, 38)), y = parseFloat(line.substring(38, 46)), z = parseFloat(line.substring(46, 54));
        if (atomName === 'CB') {
            const idx = out[chain].findIndex(c => c.resnum === resnum);
            if (idx >= 0) out[chain][idx] = { resnum, x, y, z }; else out[chain].push({ resnum, x, y, z });
            seen[chain].add(resnum);
        } else if (atomName === 'CA' && resName === 'GLY' && !seen[chain].has(resnum)) {
            out[chain].push({ resnum, x, y, z });
            seen[chain].add(resnum);
        }
    }
    return out;
}

function combinePlddt(pA, pB, wA, wB) {
    if (!isNaN(pA) && !isNaN(pB) && (wA + wB) > 0) return (pA * wA + pB * wB) / (wA + wB);
    if (!isNaN(pA)) return pA;
    if (!isNaN(pB)) return pB;
    return NaN;
}

// Per-rank pDockQ/LIpDockQ: unlike the partner-summary table above (aggregated scalars only), a
// selected pair's actual structure IS fetched (flyPdbPerRank) and its per-rank cLIR residue sets
// ARE already parsed (r._cLIR_A/_cLIR_B, from parseResidueRanges on the scraped detail table) --
// so both metrics are fully computable here, not just LIpDockQ. pDockQ2/LIpDockQ2 still aren't: no
// numeric PAE is ever exposed, only a PAE image.
function weightedPlddt(bfs, set, chain) {
    let s = 0, n = 0;
    for (const res of set) { const v = bfs.get(`${chain}:${res}`); if (v != null) { s += v; n++; } }
    return n > 0 ? s / n : NaN;
}

function pDockqForRank(pdbText, cbCutoff) {
    if (!pdbText) return null;
    cbCutoff = cbCutoff || 8;
    const cut2 = cbCutoff * cbCutoff;
    const coords = parsePdbCbCoordsAB(pdbText);
    const geomA = new Set(), geomB = new Set();
    let nContacts = 0;
    for (const a of coords.A) {
        for (const b of coords.B) {
            const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
            if (dx * dx + dy * dy + dz * dz <= cut2) { geomA.add(a.resnum); geomB.add(b.resnum); nContacts++; }
        }
    }
    const bfs = parseBfactorsPerResidue(pdbText, 'pdb');
    const plddt = combinePlddt(weightedPlddt(bfs, geomA, 'A'), weightedPlddt(bfs, geomB, 'B'), geomA.size, geomB.size);
    // N is the contact PAIR count, matching contacts.shape[0] in the reference implementation --
    // NOT geomA.size+geomB.size (a residue count).
    return calcPdockq(plddt, nContacts);
}

// LIpDockQ: the page's own cLIR interface (already-known residue sets, r._cLIR_A/_cLIR_B) + a
// FRESH Cb contact count restricted to those residues, computed directly from the fetched
// structure -- not FlyPredictome's own cLIA (which sums both PAE directions, double-counting a
// cell confident in both, and isn't on the same scale as pDockQ's own N). Restricting a real
// contact matrix to (row in cLIR_A, col in cLIR_B) recovers a cell count very close to the true
// either_contact.sum() (validated against live-PAE ground truth on real multi-chain predictions: N
// off by 0-6 cells out of 40-250, LIpDockQ off by <0.04, usually <0.003) -- still an approximation (a
// residue's *specific* qualifying partner isn't recoverable without PAE), but much closer than cLIA.
function LIpDockqForRank(r, pdbText) {
    if (!pdbText || r['cLIR Indice A'] == null || r['cLIR Indice B'] == null) return null;
    // A rank with no cLIR has no interface, as in lis.py: 0, a value the Avg row counts.
    if (r._cLIR_A.size === 0 || r._cLIR_B.size === 0) return 0;
    const bfs = parseBfactorsPerResidue(pdbText, 'pdb');
    const plddt = combinePlddt(weightedPlddt(bfs, r._cLIR_A, 'A'), weightedPlddt(bfs, r._cLIR_B, 'B'), r._cLIR_A.size, r._cLIR_B.size);
    const coords = parsePdbCbCoordsAB(pdbText);
    const cut2 = 8 * 8;
    let n = 0;
    for (const a of coords.A) {
        if (!r._cLIR_A.has(a.resnum)) continue;
        for (const b of coords.B) {
            if (!r._cLIR_B.has(b.resnum)) continue;
            const dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
            if (dx * dx + dy * dy + dz * dz <= cut2) n++;
        }
    }
    return calcPdockq(plddt, n);
}

function pdqFmt(v) { return v == null || isNaN(v) ? '—' : v.toFixed(3); }

function ranksHint() {   // the summary table is wider than a phone: say so, as the Atlas does
    const w = document.getElementById('ranks-wrap'), h = document.getElementById('ranks-hint'); if (!w || !h) return;
    if (!w._ro && window.ResizeObserver) { w._ro = new ResizeObserver(ranksHint); w._ro.observe(w); const t = document.getElementById('ranks-table'); if (t) w._ro.observe(t); }   // also when the card first shows
    h.hidden = !(w.scrollWidth > w.clientWidth + 2);
}

function showError(msg) {
    const el = document.getElementById('error-msg');
    el.textContent = msg;
    el.classList.add('active');
}

// PAE plot is fetched as a same-origin blob, so a plain <a download> works (no re-fetch).
function addPaeDownloadBtn(img) {
    const id = 'pae-dl-bar';
    const old = document.getElementById(id); if (old) old.remove();
    const bar = document.createElement('div'); bar.id = id;
    bar.style.cssText = 'text-align:center; margin-top:4px;';
    const a = document.createElement('a');
    a.textContent = '\u2193 PNG (all ranks)';
    a.href = img.src; a.download = 'pae_plot.png';
    a.style.cssText = 'font-size:0.72rem; color:#2471A3; text-decoration:none; font-weight:600;';
    bar.appendChild(a);
    if (img.parentElement) img.parentElement.appendChild(bar);
}

// One lis.py-style CSV line from a parsed rank row. name uses `_vs_` so cLIP's parseProteinNames splits it correctly.
function flyRowToCsvLine(r) {
    const sym1 = r['Symbol 1'] || '';
    const sym2 = r['Symbol 2'] || '';
    return [
        `${sym1}_vs_${sym2}`, r['Rank'] || '', r['Rank'] || '', sym1, sym2,
        r['iLIS'] != null ? parseFloat(r['iLIS']).toFixed(6) : '',
        (() => { const la=parseFloat(r['LIA']||0)||(parseFloat(r['LIA (AB)']||0)+parseFloat(r['LIA (BA)']||0)),cl=parseFloat(r['cLIA']||0)||(parseFloat(r['cLIA (AB)']||0)+parseFloat(r['cLIA (BA)']||0)); return la>0&&cl>0?Math.sqrt(la*cl).toFixed(1):''; })(),
        (() => { const la=parseFloat(r['LIA']||0)||(parseFloat(r['LIA (AB)']||0)+parseFloat(r['LIA (BA)']||0)),cl=parseFloat(r['cLIA']||0)||(parseFloat(r['cLIA (AB)']||0)+parseFloat(r['cLIA (BA)']||0)),il=parseFloat(r['iLIS']||0); const ilia=Math.sqrt(la*cl); return il>0&&ilia>0?(il*ilia).toFixed(1):''; })(),
        r['LIS'] != null ? parseFloat(r['LIS']).toFixed(6) : '',
        r['cLIS'] != null ? parseFloat(r['cLIS']).toFixed(6) : '',
        (parseFloat(r['LIA']||0)||(parseFloat(r['LIA (AB)']||0)+parseFloat(r['LIA (BA)']||0))).toFixed(1),
        (parseFloat(r['cLIA']||0)||(parseFloat(r['cLIA (AB)']||0)+parseFloat(r['cLIA (BA)']||0))).toFixed(1),
        r['ipTM'] != null ? parseFloat(r['ipTM']).toFixed(6) : '',
        r['pLDDT'] != null ? parseFloat(r['pLDDT']).toFixed(1) : '',
        r['pTM'] != null ? parseFloat(r['pTM']).toFixed(6) : '',
        r._LIR_A.size, r._LIR_B.size, r._cLIR_A.size, r._cLIR_B.size,
        r['Len A'] || '', r['Len B'] || '',
        `"[${[...r._LIR_A].sort((a,b)=>a-b).join(', ')}]"`,
        `"[${[...r._LIR_B].sort((a,b)=>a-b).join(', ')}]"`,
        `"[${[...r._cLIR_A].sort((a,b)=>a-b).join(', ')}]"`,
        `"[${[...r._cLIR_B].sort((a,b)=>a-b).join(', ')}]"`,
    ].join(',');
}

function buildGatherOverlay(total) {
    const ov = document.createElement('div');
    ov.style.cssText = 'position:fixed; inset:0; background:rgba(0,0,0,0.55); display:flex; align-items:center; justify-content:center; z-index:99999;';
    const box = document.createElement('div');
    box.style.cssText = 'background:#fff; border-radius:12px; padding:1.6rem 2rem; min-width:330px; text-align:center; box-shadow:0 8px 30px rgba(0,0,0,0.3); font-family:inherit;';
    box.innerHTML =
        '<div style="font-weight:600; font-size:1rem; margin-bottom:0.3rem;">Gathering interactome for cLIP…</div>' +
        '<div style="font-size:0.8rem; color:#6B6B6B; margin-bottom:0.9rem;">Fetching contact residues (cLIR) from each partner prediction page</div>' +
        '<div style="height:8px; background:#eee; border-radius:4px; overflow:hidden; margin-bottom:0.5rem;"><div id="gather-bar" style="height:100%; width:0%; background:#2471A3; transition:width 0.25s;"></div></div>' +
        '<div id="gather-stat" style="font-size:0.85rem; color:#555; margin-bottom:1rem;">0 / ' + total + '</div>' +
        '<button id="gather-stop" style="padding:0.4rem 1rem; border:1px solid #ccc; border-radius:8px; background:#fff; cursor:pointer; font-size:0.85rem;">Stop &amp; send gathered so far</button>';
    ov.appendChild(box);
    document.body.appendChild(ov);
    const api = {
        stopped: false, total,
        update(d, okc) { const b=document.getElementById('gather-bar'); if(b) b.style.width = Math.round(d/total*100)+'%'; const s=document.getElementById('gather-stat'); if(s) s.textContent = d+' / '+total+'  ('+okc+' with data)'; },
        close() { ov.remove(); }
    };
    box.querySelector('#gather-stop').onclick = () => { api.stopped = true; box.querySelector('#gather-stop').textContent = 'Stopping…'; };
    return api;
}

function handleGeneAcKey(e) {
    const list = document.getElementById('gene-ac-list');
    const items = list.querySelectorAll('div');
    if (list.style.display !== 'none' && items.length > 0) {
        if (e.key === 'ArrowDown') {
            e.preventDefault();
            geneAcIdx = Math.min(geneAcIdx + 1, items.length - 1);
            items.forEach((it, i) => { it.style.background = i === geneAcIdx ? '#eaf2f8' : ''; });
            items[geneAcIdx].scrollIntoView({ block: 'nearest' });
            return;
        }
        if (e.key === 'ArrowUp') {
            e.preventDefault();
            geneAcIdx = Math.max(geneAcIdx - 1, 0);
            items.forEach((it, i) => { it.style.background = i === geneAcIdx ? '#eaf2f8' : ''; });
            items[geneAcIdx].scrollIntoView({ block: 'nearest' });
            return;
        }
        if (e.key === 'Enter') {
            e.preventDefault();
            if (geneAcIdx >= 0 && geneAcIdx < items.length) {
                items[geneAcIdx].click();
            } else {
                hideGeneAc();
                searchGene();
            }
            return;
        }
        if (e.key === 'Escape') {
            hideGeneAc();
            return;
        }
    } else if (e.key === 'Enter') {
        hideGeneAc();
        searchGene();
    }
}

function hideGeneAc() { document.getElementById('gene-ac-list').style.display = 'none'; geneAcIdx = -1; }

// pDockQ (Bryant, Pozzati & Elofsson 2022, author-corrected constants), fed this file's own cLIR
// interface to give LIpDockQ -- matches lis.py's calc_pdockq() and universal.html's calcPdockq()
// exactly. Verified against the authors' own reference implementation (ElofssonLab/FoldDock,
// src/pdockq.py): nIfContacts must be the COUNT OF CONTACT PAIRS (cells), NOT a residue count and
// not LIA/cLIA-style double-counted sums -- x = avgIfPlddt * log10(n), no +1. pDockQ2/LIpDockQ2 are
// NOT computable here: FlyPredictome's summary table carries no numeric PAE, only a PAE image.
function calcPdockq(avgIfPlddt, nIfContacts) {
    if (nIfContacts == null || nIfContacts < 0) return null;
    // Matches the reference implementation's own explicit special case (contacts.shape[0]<1 ->
    // pdockq=0) -- a literal 0, not the sigmoid evaluated at x=0.
    if (nIfContacts === 0) return 0;
    if (avgIfPlddt == null || isNaN(avgIfPlddt)) return null;
    const x = avgIfPlddt * Math.log10(nIfContacts);
    return 0.724 / (1 + Math.exp(-0.052 * (x - 152.611))) + 0.018;
}

// ── Parse the FlyPredictome details page's rank table ──
// The redesigned page renders ranks in a DataTables table (#rank_table_test). Its <thead> prepends a
// compact 15-column display header before the real ~60-column data header, so each row's cells align
// with the LAST N headers, not the first N. Columns are read by name, so the rest is layout-agnostic.
function parseTable(html) {
    const parser = new DOMParser();
    const doc = parser.parseFromString(html, 'text/html');
    const table = doc.getElementById('rank_table_test') || doc.getElementById('rank_table');   // new id, old fallback
    if (!table) return null;

    const headers = [];
    table.querySelectorAll('thead th').forEach(th => headers.push(th.textContent.trim()));

    const rows = [];
    table.querySelectorAll('tbody tr').forEach(tr => {
        const cells = [];
        tr.querySelectorAll('td').forEach(td => cells.push(td.textContent.trim()));
        if (cells.length >= 10) {                                             // skip DataTables placeholder/empty rows
            const hdr = headers.length > cells.length ? headers.slice(headers.length - cells.length) : headers;
            const row = {};
            cells.forEach((val, i) => { if (i < hdr.length) row[hdr[i]] = val; });
            // The rank column is now "Rank/Model" ("1 / 5"); keep row['Rank'] = the rank number for back-compat.
            if (!row['Rank'] && row['Rank/Model']) row['Rank'] = String(row['Rank/Model']).split('/')[0].trim();
            if (!/^\d+$/.test(String(row['Rank'] || '').trim())) return;   // skip FlyPredictome's own "Average" summary row — LIVIA computes its own Avg
            row._cLIR_A = parseResidueRanges(row['cLIR Indice A'] || '');
            row._cLIR_B = parseResidueRanges(row['cLIR Indice B'] || '');
            row._LIR_A = parseResidueRanges(row['LIR Indice A'] || '');
            row._LIR_B = parseResidueRanges(row['LIR Indice B'] || '');
            rows.push(row);
        }
    });

    // Extract PDB URL — prefer the full /tools/fly_predictome/web/ path (the page also has a bare
    // /colabfold-output/ reference that 404s), matching how the PAE image is resolved below.
    let pdbUrl = '';
    const pdbMatch = html.match(/["'](\/tools\/fly_predictome\/web\/colabfold-output\/[^"']*\.pdb)["']/)
                  || html.match(/["']([^"']*colabfold-output[^"']*\.pdb)["']/);
    if (pdbMatch) {
        pdbUrl = pdbMatch[1];
        if (!pdbUrl.startsWith('http')) pdbUrl = 'https://www.flyrnai.org' + pdbUrl;
    }

    // NCBI gene IDs are now table columns (ncbi_gene_id1/2); fall back to the old class/regex probe.
    let ncbiGeneId1 = (rows[0] && (rows[0]['ncbi_gene_id1'] || '').trim()) || '';
    let ncbiGeneId2 = (rows[0] && (rows[0]['ncbi_gene_id2'] || '').trim()) || '';
    if (!ncbiGeneId1) {
        const ncbi1El = doc.querySelector('.ncbi_gene_id1'); ncbiGeneId1 = ncbi1El ? ncbi1El.textContent.trim() : '';
        if (!ncbiGeneId1) { const m1 = html.match(/ncbi_gene_id1[^>]*>\s*(\d+)/); if (m1) ncbiGeneId1 = m1[1]; }
    }
    if (!ncbiGeneId2) {
        const ncbi2El = doc.querySelector('.ncbi_gene_id2'); ncbiGeneId2 = ncbi2El ? ncbi2El.textContent.trim() : '';
        if (!ncbiGeneId2) { const m2 = html.match(/ncbi_gene_id2[^>]*>\s*(\d+)/); if (m2) ncbiGeneId2 = m2[1]; }
    }

    // Extract PAE image URL (prefer colabfold-output path which is publicly accessible)
    let paeImgUrl = '';
    const paeImgMatch = html.match(/["'](\/tools\/fly_predictome\/web\/colabfold-output\/[^"']*_pae\.png)["']/);
    if (paeImgMatch) {
        paeImgUrl = 'https://www.flyrnai.org' + paeImgMatch[1];
    } else {
        // Fallback: try data/extracted path
        const paeImgMatch2 = html.match(/["'](\/tools\/fly_predictome\/web\/data\/extracted\/png\/[^"']*_pae\.png)["']/);
        if (paeImgMatch2) paeImgUrl = 'https://www.flyrnai.org' + paeImgMatch2[1];
    }

    return { headers, rows, pdbUrl, ncbiGeneId1, ncbiGeneId2, paeImgUrl };
}

// ── Chain visibility toggle (3D viewer only — fly page has no chord) ──
function renderChainToggleBar() {
    const bar = document.getElementById('chain-toggle-bar');
    const btnHolder = document.getElementById('chain-toggle-buttons');
    if (!bar || !btnHolder || !flyParsed) { if (bar) bar.style.display = 'none'; return; }
    bar.style.display = 'flex';
    btnHolder.innerHTML = '';
    for (const name of ['A', 'B']) {
        const btn = document.createElement('button');
        btn.dataset.chain = name;
        btn.textContent = name;
        btn.onclick = () => toggleChainVisibility(name);
        applyChainToggleStyle(btn, !flyHiddenChains.has(name));
        btnHolder.appendChild(btn);
    }
}

function getChainToggleColor(chainName) {
    if (chainName === 'A') return document.getElementById('color-clir-a')?.value || '#00897B';
    if (chainName === 'B') return document.getElementById('color-clir-b')?.value || '#E64A19';
    return '#888';
}

function toggleChainVisibility(name) {
    if (flyHiddenChains.has(name)) flyHiddenChains.delete(name);
    else flyHiddenChains.add(name);
    const btn = document.querySelector(`#chain-toggle-buttons button[data-chain="${name}"]`);
    if (btn) applyChainToggleStyle(btn, !flyHiddenChains.has(name));
    if (flyCurrentRank) loadRank(flyCurrentRank);
}

function resetChainVisibility() {
    if (flyHiddenChains.size === 0) return;
    flyHiddenChains = new Set();
    document.querySelectorAll('#chain-toggle-buttons button').forEach(b => applyChainToggleStyle(b, true));
    if (flyCurrentRank) loadRank(flyCurrentRank);
}

async function fetchViaProxy(url) {
    // Race all proxies in parallel — use whichever responds first with valid data
    return new Promise(resolve => {
        let resolved = false;
        let pending = CORS_PROXIES.length;

        CORS_PROXIES.forEach(makeProxy => {
            const controller = new AbortController();
            const timeout = setTimeout(() => controller.abort(), 12000);

            fetch(makeProxy(url), { signal: controller.signal })
                .then(resp => resp.ok ? resp.text() : '')
                .then(text => {
                    clearTimeout(timeout);
                    if (!resolved && text.includes('rank_table')) {
                        resolved = true;
                        resolve(text);
                    }
                })
                .catch(() => {})
                .finally(() => {
                    pending--;
                    if (pending === 0 && !resolved) resolve(null);
                });
        });
    });
}

async function showPaeMap() {
    const card = document.getElementById('pae-card');
    const img = document.getElementById('pae-img');
    if (!card || !img || !flyParsed?.paeImgUrl) { if (card) card.style.display = 'none'; return; }

    card.style.display = '';
    // Fetch image as blob via proxy (same approach as PDB download)
    try {
        let blob = null;
        // Try direct first
        try {
            const resp = await fetch(flyParsed.paeImgUrl);
            if (resp.ok) blob = await resp.blob();
        } catch {}
        // Try proxies
        if (!blob) {
            for (const makeProxy of CORS_PROXIES) {
                try {
                    const resp = await fetch(makeProxy(flyParsed.paeImgUrl), { signal: AbortSignal.timeout(10000) });
                    if (resp.ok) {
                        blob = await resp.blob();
                        if (blob.size > 100) break;
                        blob = null;
                    }
                } catch {}
            }
        }
        if (blob && blob.size > 100) {
            img.src = URL.createObjectURL(blob);
            addPaeDownloadBtn(img);
        } else {
            card.style.display = 'none';
        }
    } catch {
        card.style.display = 'none';
    }
}

function readAtlasScores() {
    flyAtlas = null; flyAtlasOpen = null;
    const m = /(?:^#|&)atlas=([^&]*)/.exec(location.hash); if (!m || !flyParsed || !flyParsed.rows) return;
    const cur = ((document.getElementById('url-input') || {}).value || '').trim(); if (!ATLAS_URL0 || (cur && cur !== ATLAS_URL0)) return;   // another pair opened since
    let d; try { d = JSON.parse(decodeURIComponent(m[1])); } catch (e) { return; }
    if (!d || !Array.isArray(d.k) || !d.r || typeof d.r !== 'object') return;
    const by = new Map();
    for (const [rk, vals] of Object.entries(d.r)) { if (!Array.isArray(vals)) return; const o = {}; d.k.forEach((k, i) => { if (vals[i] != null && Number.isFinite(+vals[i])) o[k] = +vals[i]; }); by.set(+rk, o); }
    for (const r of flyParsed.rows) { const a = by.get(parseInt(r['Rank'])), fi = parseFloat(r['iLIS']);
        if (!a || !Number.isFinite(a.iLIS) || !Number.isFinite(fi) || Math.abs(a.iLIS - fi) > 0.002) return; }
    flyAtlas = by;
    const o = parseInt(d.o); if (by.has(o) && flyParsed.rows.some((r) => parseInt(r['Rank']) === o)) flyAtlasOpen = o;   // the Atlas link names the row's rank
}

function atlasHead() {
    const tr = document.getElementById('ranks-head'); if (!tr) return;
    if (flyRanksHead0 == null) flyRanksHead0 = tr.innerHTML;
    const r0 = (flyParsed && flyParsed.rows && flyParsed.rows[0]) || {}, e = (x) => String(x || '').replace(/[&<>"]/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;' })[c]);
    const ab = r0['Symbol 1'] && r0['Symbol 2'] ? `(${e(r0['Symbol 1'])} / ${e(r0['Symbol 2'])})` : '(A/B)';   // the two proteins, chain A first, in place of "A/B"
    tr.innerHTML = (flyAtlas ? '<th>Rank</th><th>Symbol 1</th><th>Symbol 2</th><th>iLIS</th><th>iLIA</th><th>iLISA</th><th>ipTM</th><th>ipSAE</th><th title="lis.py\'s approximation of actifpTM (Varga et al. 2025)">actifpTM*</th><th>pDockQ</th><th>LIpDockQ</th><th>pDockQ2</th><th>LIpDockQ2</th><th>LIR (A/B)</th><th>cLIR (A/B)</th><th>Len (A/B)</th>' : flyRanksHead0).replace(/\(A\/B\)/g, ab);
    const nt = document.getElementById('atlas-scores-note'); if (nt) nt.hidden = !flyAtlas;
    const pn = document.getElementById('pdq2-note');   // with the Atlas scores, pDockQ2 and LIpDockQ2 are shown: their rows appear and the note keeps only its caveat
    if (pn) { pn.dataset.base = pn.dataset.base || pn.innerHTML; pn.innerHTML = flyAtlas ? 'pDockQ2 and LIpDockQ2 come with the LIVIA Atlas scores (lis.py on the numeric PAE). Because their constants were fitted to DockQ on the original interfaces, LIpDockQ and LIpDockQ2 rank predictions; they are not DockQ estimates.' : pn.dataset.base; }
    document.querySelectorAll('.pdq2-row').forEach(r => { r.hidden = !flyAtlas; });
}

function atlasCells(rk) {   // the Atlas's six scores for one rank, or their mean over the ranks for the Avg row
    let a = {};
    if (rk === 'avg') { for (const k of ATLAS_COLS) { const xs = [...flyAtlas.values()].map((o) => o[k]).filter(Number.isFinite); if (xs.length) a[k] = xs.reduce((x, y) => x + y, 0) / xs.length; } }
    else a = flyAtlas.get(parseInt(rk)) || {};
    return ATLAS_COLS.map((k) => `<td>${Number.isFinite(a[k]) ? a[k].toFixed(3) : '—'}</td>`).join('');
}

async function fetchAllRankPdbs() {
    if (!flyParsed || !flyParsed.rows || flyParsed.rows.length <= 1) return;
    for (const r of flyParsed.rows) {
        const rankNum = parseInt(r['Rank']);
        if (flyPdbPerRank[rankNum]) continue;
        const fetched = await fetchRankPdb(rankNum);
        if (fetched) {
            // Update combined chord diagram with newly available data
            updateCombinedChord();
            // Fill in this rank's pDockQ/LIpDockQ cells now that its structure has arrived
            const pdbText = flyPdbPerRank[rankNum];
            const pdqCell = document.getElementById('pdq-rank-' + rankNum);
            const cpdqCell = document.getElementById('cpdq-rank-' + rankNum);
            if (pdqCell) pdqCell.textContent = pdqFmt(pDockqForRank(pdbText));
            if (cpdqCell) cpdqCell.textContent = pdqFmt(LIpDockqForRank(r, pdbText));
        }
    }
    // Rank 1's structure was already available at first render, but the Avg row's LIpDockQ needs
    // every rank's cLIA (mean across ranks) — recompute it now that all ranks have been fetched.
    if (flyAverageRow && !flyAtlas) {   // the Avg row: the mean of each rank's own pDockQ and LIpDockQ (a rank with no value is left out)
        const pd = [], lp = [];
        flyParsed.rows.forEach((rr, k) => { const t = flyPdbPerRank[k + 1]; if (!t) return;
            const a = pDockqForRank(t), b = LIpDockqForRank(rr, t); if (Number.isFinite(a)) pd.push(a); if (Number.isFinite(b)) lp.push(b); });
        const mean = (xs) => (xs.length ? xs.reduce((x, y) => x + y, 0) / xs.length : null);
        const pc = document.getElementById('pdq-rank-avg'), lc = document.getElementById('cpdq-rank-avg');
        if (pc) { pc.textContent = pd.length ? pdqFmt(mean(pd)) : '—'; pc.title = `mean over ${pd.length} rank${pd.length === 1 ? '' : 's'}`; }
        if (lc) { lc.textContent = lp.length ? pdqFmt(mean(lp)) : '—'; lc.title = `mean over the ${lp.length} rank${lp.length === 1 ? '' : 's'} with a value`; }
    }
    // Rebuild pLDDT charts now that all rank PDBs are available
    if (flyCurrentRank) {
        const colors = {
            lir_a: document.getElementById('color-lir-a').value,
            lir_b: document.getElementById('color-lir-b').value,
            clir_a: document.getElementById('color-clir-a').value,
            clir_b: document.getElementById('color-clir-b').value,
        };
        buildFlyPlddtCharts(flyCurrentRank, colors);
    }
}

function updateCombinedChord() {
    if (!flyParsed || !flyAverageRow) return;
    const canvasAll = document.getElementById('chord-canvas-all');
    const chordAllWrap = document.getElementById('chord-all-wrap');
    if (!canvasAll || !chordAllWrap) return;
    LiviaMaps.attachExportBar('chord-canvas-all', { name: 'chord_diagram_all_ranks', svg: () => updateCombinedChord(), png: true });

    const row0 = flyParsed.rows[0];
    const sym1 = row0['Symbol 1'] || '?';
    const sym2 = row0['Symbol 2'] || '?';
    const lenA = parseInt(row0['Len A'] || row0['protein_1_size']) || 0;
    const lenB = parseInt(row0['Len B'] || row0['protein_2_size']) || 0;
    if (lenA === 0 || lenB === 0) return;

    const colors = {
        lir_a: document.getElementById('color-lir-a').value,
        lir_b: document.getElementById('color-lir-b').value,
        clir_a: document.getElementById('color-clir-a').value,
        clir_b: document.getElementById('color-clir-b').value,
    };

    // Collect contacts from all available rank PDBs
    const contactDist = 8;
    const allContacts = new Map();
    for (const r of flyParsed.rows) {
        const rankNum = parseInt(r['Rank']);
        const rankPdb = flyPdbPerRank[rankNum];
        if (!rankPdb) continue;
        const rankCoords = parsePdbCbCoords(rankPdb);
        const rCoordA = rankCoords.get('A');
        const rCoordB = rankCoords.get('B');
        if (!rCoordA || !rCoordB) continue;
        for (const resA of r._cLIR_A) {
            const cA = rCoordA.get(resA);
            if (!cA) continue;
            for (const resB of r._cLIR_B) {
                const key = `${resA}-${resB}`;
                if (allContacts.has(key)) continue;
                const cB = rCoordB.get(resB);
                if (!cB) continue;
                const dx = cA.x - cB.x, dy = cA.y - cB.y, dz = cA.z - cB.z;
                const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
                if (dist <= contactDist) allContacts.set(key, { resA, resB, dist });
            }
        }
    }

    const unionContacts = [...allContacts.values()];
    if (unionContacts.length > 0) {
        chordAllWrap.style.display = '';
        drawFlyChordFromContacts(canvasAll, unionContacts, lenA, lenB, sym1, sym2, colors, flyAverageRow._LIR_A, flyAverageRow._cLIR_A, flyAverageRow._LIR_B, flyAverageRow._cLIR_B);
    }
}

// Populate the per-rank selector chip strips (one per viewer card; all identical because
// rank selection is global). Each chip shows "Rank N" + that rank's iLIS; the chip matching
// flyCurrentRank is highlighted. Clicking routes through loadRank, which also highlights the
// matching All-Ranks table row — so chips, table, and all viewers stay in sync. No "Average"
// chip (the All-Ranks table keeps its Avg row); the Visualization Script card has no selector.
function renderModelSelectors() {
    const selectors = document.querySelectorAll('.model-selector');
    if (!selectors.length) return;
    if (!flyParsed || !flyParsed.rows || !flyParsed.rows.length) {
        selectors.forEach(s => { s.innerHTML = ''; s.style.display = 'none'; });
        return;
    }
    const esc = (s) => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
    const html = flyParsed.rows.map(r => {
        const rk = parseInt(r['Rank']);
        if (isNaN(rk)) return '';
        const sel = (flyCurrentRank !== 'avg' && rk === parseInt(flyCurrentRank)) ? ' selected' : '';
        const ilis = parseFloat(r['iLIS']);
        const ilisStr = isNaN(ilis) ? '—' : ilis.toFixed(3);
        const iptm = parseFloat(r['ipTM']);
        const iptmStr = isNaN(iptm) ? '—' : iptm.toFixed(2);
        const label = 'Rank ' + rk;
        return `<button type="button" class="model-chip${sel}" onclick="loadRank(${rk})" title="Select ${esc(label)}: FlyPredictome's ipTM-based rank (iLIS ${ilisStr}, ipTM ${iptmStr})"><span class="model-chip-name">${esc(label)}</span><span class="model-chip-ilis"><span class="model-chip-metric">iLIS</span>${ilisStr}</span><span class="model-chip-ilis"><span class="model-chip-metric">ipTM</span>${iptmStr}</span></button>`;
    }).join('');
    selectors.forEach(s => { s.innerHTML = html; s.style.display = 'flex'; });
}

// ============================================================================
// PER-RESIDUE pLDDT CHARTS (across ranks)
// ============================================================================
function parsePdbPlddt(pdbText) {
    // Returns Map<chain, Map<resnum, bfactor>>
    const chains = new Map();
    for (const line of pdbText.split('\n')) {
        if (!line.startsWith('ATOM')) continue;
        if (line.substring(12, 16).trim() !== 'CA') continue;
        const chain = line.substring(21, 22).trim();
        const resnum = parseInt(line.substring(22, 26));
        const bf = parseFloat(line.substring(60, 66));
        if (isNaN(resnum) || isNaN(bf)) continue;
        if (!chains.has(chain)) chains.set(chain, new Map());
        chains.get(chain).set(resnum, bf);
    }
    return chains;
}

function parsePdbSequence(pdbText) {
    // Extract per-chain sequences from PDB CA atoms
    const chains = new Map(); // chain -> [{resnum, aa}]
    for (const line of pdbText.split('\n')) {
        if (!line.startsWith('ATOM')) continue;
        const atomName = line.substring(12, 16).trim();
        if (atomName !== 'CA') continue;
        const chain = line.substring(21, 22).trim();
        const resName = line.substring(17, 20).trim();
        const resSeq = parseInt(line.substring(22, 26).trim());
        if (isNaN(resSeq)) continue;
        if (!chains.has(chain)) chains.set(chain, []);
        const arr = chains.get(chain);
        if (!arr.find(r => r.resnum === resSeq)) {
            arr.push({ resnum: resSeq, aa: AA3TO1[resName] || 'X' });
        }
    }
    return chains;
}

async function fetchDomainsFromPfam(uniprotId) {
    if (domainCache.has(uniprotId)) return domainCache.get(uniprotId);
    try {
        const resp = await fetch(`https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/${uniprotId}?format=json`);
        if (!resp.ok) { domainCache.set(uniprotId, []); return []; }
        const data = await resp.json();
        const domains = [];
        for (const entry of (data.results || [])) {
            const name = entry.metadata?.name || 'Pfam domain';
            for (const prot of (entry.proteins || [])) {
                for (const loc of (prot.entry_protein_locations || [])) {
                    for (const frag of (loc.fragments || [])) {
                        if (frag.start != null && frag.end != null) {
                            domains.push({ name, start: frag.start, end: frag.end });
                        }
                    }
                }
            }
        }
        domains.sort((a, b) => a.start - b.start);
        domainCache.set(uniprotId, domains);
        return domains;
    } catch (e) {
        domainCache.set(uniprotId, []);
        return [];
    }
}

function parsePdbCbCoords(pdbText) {
    // Extract Cβ (or CA for GLY) coordinates per chain
    const chains = new Map(); // chain -> Map(resnum -> {x,y,z})
    for (const line of pdbText.split('\n')) {
        if (!line.startsWith('ATOM')) continue;
        const atomName = line.substring(12, 16).trim();
        const resName = line.substring(17, 20).trim();
        const chain = line.substring(21, 22).trim();
        const resSeq = parseInt(line.substring(22, 26).trim());
        const x = parseFloat(line.substring(30, 38));
        const y = parseFloat(line.substring(38, 46));
        const z = parseFloat(line.substring(46, 54));
        if (isNaN(resSeq) || isNaN(x)) continue;
        if (!chains.has(chain)) chains.set(chain, new Map());
        const m = chains.get(chain);
        if (atomName === 'CB') {
            m.set(resSeq, {x, y, z});
        } else if (atomName === 'CA' && resName === 'GLY' && !m.has(resSeq)) {
            m.set(resSeq, {x, y, z});
        }
    }
    return chains;
}

async function handoffToClip(cap) {
    const query = window._searchQuery;
    let partners = (window._searchPartners || []).filter(partnerPassesFilters);
    const uEl = document.getElementById('clip-unique');
    if (!uEl || uEl.checked) {   // dedup to unique partner: keep highest iLIS(best); tie -> highest iLIS(avg)
        const best = {};
        for (const p of partners) { const k = p.partnerSymbol, c = best[k]; if (!c || p.ilisBest > c.ilisBest || (p.ilisBest === c.ilisBest && p.ilisAvg > c.ilisAvg)) best[k] = p; }
        partners = Object.values(best);
    }
    partners.sort((a,b) => (b.ilisBest - a.ilisBest) || (b.ilisAvg - a.ilisAvg));   // best iLIS first
    if (partners.length < 2) { alert('Need at least 2 unique partners passing the current filters to cluster in cLIP.'); return; }
    partners = partners.slice(0, cap || 100);

    const ov = buildGatherOverlay(partners.length);
    const CONC = 10;                              // concurrent Worker fetches
    const cache = window._detailCache = window._detailCache || {};   // avoid re-fetching within a session
    const lines = [FLY_CSV_HEADER];
    let done = 0, ok = 0, idx = 0;
    const worker = async () => {
        while (idx < partners.length && !ov.stopped) {
            const p = partners[idx++];
            if (p.detailUrl) {
                try {
                    let html = cache[p.detailUrl];
                    if (!html) { html = await fetchViaProxy(p.detailUrl); if (html) cache[p.detailUrl] = html; }
                    if (html) { const parsed = parseTable(html); if (parsed && parsed.rows && parsed.rows.length) { parsed.rows.forEach(r => lines.push(flyRowToCsvLine(r))); ok++; } }
                } catch (e) { /* skip this partner */ }
            }
            ov.update(++done, ok);
        }
    };
    await Promise.all(Array.from({ length: CONC }, worker));
    ov.close();
    if (ok < 2) { alert('cLIP needs at least 2 partners with contact data — only ' + ok + ' gathered. Try another gene or a larger cap.'); return; }
    // Derive the query's symbol AS IT APPEARS in the gathered CSV — flypredictome may label e.g. Akt1 as "Akt",
    // so the searched name won't match cLIP's #gene option. The query is the symbol common to (almost) every pair.
    const symCount = {};
    for (let i = 1; i < lines.length; i++) { const c = lines[i].split(','); symCount[c[2]] = (symCount[c[2]]||0)+1; symCount[c[3]] = (symCount[c[3]]||0)+1; }
    const querySym = Object.keys(symCount).sort((a,b)=>symCount[b]-symCount[a])[0] || query;
    localStorage.setItem('livia_net2clip', JSON.stringify({ csv: lines.join('\n'), query, querySym }));
    window.open('clip.html', '_blank');
}

function buildFlyColorComponents(row, colors, pdbText) {
    const components = [];
    const LIR_A_filtered = filterSmallSegments(row._LIR_A);
    const LIR_B_filtered = filterSmallSegments(row._LIR_B);
    const cLIR_A = new Set([...row._cLIR_A].filter(r => LIR_A_filtered.has(r)));
    const cLIR_B = new Set([...row._cLIR_B].filter(r => LIR_B_filtered.has(r)));
    let lirAFilled = fillGaps(LIR_A_filtered);
    let lirBFilled = fillGaps(LIR_B_filtered);

    // When showComplete, expand to all residues so color mode applies to full chains
    if (showComplete) {
        const lenA = parseInt(row['Len A'] || row['protein_1_size']) || 0;
        const lenB = parseInt(row['Len B'] || row['protein_2_size']) || 0;
        lirAFilled = new Set(); for (let i = 1; i <= lenA; i++) lirAFilled.add(i);
        lirBFilled = new Set(); for (let i = 1; i <= lenB; i++) lirBFilled.add(i);
    }

    // Helper: convert Set to array of {start, end} ranges
    function toRanges(positions) {
        const sorted = [...positions].sort((a, b) => a - b);
        if (sorted.length === 0) return [];
        const ranges = [];
        let s = sorted[0], e = sorted[0];
        for (let i = 1; i <= sorted.length; i++) {
            if (i < sorted.length && sorted[i] === e + 1) { e = sorted[i]; continue; }
            ranges.push({ start: s, end: e });
            if (i < sorted.length) { s = e = sorted[i]; }
        }
        return ranges;
    }

    if (cxcColorMode === 'plddt') {
        const bfactors = parseBfactorsPerResidue(pdbText);
        const addPlddt = (chain, positions) => {
            const byColor = {};
            for (const resnum of positions) {
                const bf = bfactors.get(`${chain}:${resnum}`);
                const color = bf !== undefined ? plddtColor(bf) : '#65CBF3';
                if (!byColor[color]) byColor[color] = new Set();
                byColor[color].add(resnum);
            }
            for (const [color, posSet] of Object.entries(byColor)) {
                components.push({ chain, ranges: toRanges(posSet), color });
            }
        };
        addPlddt('A', lirAFilled);
        addPlddt('B', lirBFilled);
    } else if (cxcColorMode === 'bychain' || cxcColorMode === 'bypolymer') {
        components.push({ chain: 'A', ranges: toRanges(lirAFilled), color: colors.clir_a });
        components.push({ chain: 'B', ranges: toRanges(lirBFilled), color: colors.clir_b });
    } else {
        // Gradient: LIR-only in light, cLIR in dark (separate components)
        if (grayNonLir && showComplete) {
            // Gray non-LIR: separate LIR from non-LIR
            const LIR_A_orig = fillGaps(filterSmallSegments(row._LIR_A));
            const LIR_B_orig = fillGaps(filterSmallSegments(row._LIR_B));
            const nonLirA = new Set([...lirAFilled].filter(r => !LIR_A_orig.has(r)));
            const nonLirB = new Set([...lirBFilled].filter(r => !LIR_B_orig.has(r)));
            components.push({ chain: 'A', ranges: toRanges(nonLirA), color: '#B0B0B0' });
            components.push({ chain: 'B', ranges: toRanges(nonLirB), color: '#B0B0B0' });
            const lirOnlyA = new Set([...LIR_A_orig].filter(r => !cLIR_A.has(r)));
            const lirOnlyB = new Set([...LIR_B_orig].filter(r => !cLIR_B.has(r)));
            components.push({ chain: 'A', ranges: toRanges(lirOnlyA), color: colors.lir_a });
            components.push({ chain: 'B', ranges: toRanges(lirOnlyB), color: colors.lir_b });
        } else {
            // Default: all shown residues (lirAFilled) minus cLIR get LIR color
            const lirOnlyA = new Set([...lirAFilled].filter(r => !cLIR_A.has(r)));
            const lirOnlyB = new Set([...lirBFilled].filter(r => !cLIR_B.has(r)));
            components.push({ chain: 'A', ranges: toRanges(lirOnlyA), color: colors.lir_a });
            components.push({ chain: 'B', ranges: toRanges(lirOnlyB), color: colors.lir_b });
        }
        components.push({ chain: 'A', ranges: toRanges(cLIR_A), color: colors.clir_a });
        components.push({ chain: 'B', ranges: toRanges(cLIR_B), color: colors.clir_b });
    }
    return components.filter(c => c.ranges.length > 0 && !flyHiddenChains.has(c.chain));
}

function updateViewer3d(row, colors, pdbText) {
    const card = document.getElementById('viewer3d-card');
    const frame = document.getElementById('viewer3d-frame');
    if (!pdbText || !frame) { if (card) card.style.display = 'none'; return; }

    card.style.display = '';

    // Build separate components per colored group (non-LIR = hidden)
    const colorComponents = applyViewer3dDomains(buildFlyColorComponents(row, colors, pdbText), pdbText, 'pdb', ['A', 'B'].map((c, k) => ({ chain: c, label: (row && row[`Symbol ${k + 1}`]) || c, uniprot: flyDomains[c], ted: flyTedDomains[c] })));   // domain coloring (the 3D color bar), else the page's colors

    // Fast path: same structure → only colors changed → postMessage (no iframe reload)
    if (frame.src && viewer3dStructText === pdbText) {
        applyColorsToMolstarFrame('viewer3d-frame', colorComponents, 'pdb');
        return;
    }

    // Warm-up path: iframe pre-booted at page entry with placeholder → swap in place
    if (isMolstarFrameWarmedUp('viewer3d-frame')) {
        viewer3dStructText = pdbText;
        swapMolstarStructure('viewer3d-frame', pdbText, colorComponents, 'pdb', `${row['Symbol 1'] || 'A'} × ${row['Symbol 2'] || 'B'}, ${row['Rank'] === 'Avg' ? 'average row' : 'rank ' + (row['Rank'] || '?')}`);
        return;
    }

    // Slow path: structure changed AND no warm-up → full iframe rebuild
    viewer3dStructText = pdbText;
    const page = buildMolstarPage(pdbText, 'pdb', colorComponents);
    const blob = new Blob([page], { type: 'text/html' });
    if (viewer3dBlobUrl) URL.revokeObjectURL(viewer3dBlobUrl);
    viewer3dBlobUrl = URL.createObjectURL(blob);
    const parent = frame.parentNode;
    const newFrame = frame.cloneNode(false);
    newFrame.src = viewer3dBlobUrl;
    parent.replaceChild(newFrame, frame);
}

function buildFlyPlddtCharts(selectedRank, colors) {
    const card = document.getElementById('plddt-card');
    const container = document.getElementById('plddt-charts');
    if (!card || !container || !flyParsed || !flyParsed.rows) return;

    const availableRanks = Object.keys(flyPdbPerRank).map(Number).sort((a, b) => a - b);
    if (availableRanks.length === 0) { card.style.display = 'none'; return; }

    // Parse pLDDT from all available rank PDBs
    const allPlddt = []; // [{rank, chains: Map<chain, Map<resnum, bf>>}]
    for (const rk of availableRanks) {
        allPlddt.push({ rank: rk, chains: parsePdbPlddt(flyPdbPerRank[rk]) });
    }

    // Determine chain names and sizes from rank 1 (or first available)
    const firstPlddt = allPlddt[0].chains;
    const chainNames = [...firstPlddt.keys()];
    if (chainNames.length === 0) { card.style.display = 'none'; return; }

    const row0 = flyParsed.rows[0];
    const sym1 = row0?.['Symbol 1'] || 'A';
    const sym2 = row0?.['Symbol 2'] || 'B';
    const chainLabels = { A: sym1, B: sym2 };

    const isAvg = selectedRank === 'avg';
    const selRankNum = isAvg ? null : parseInt(selectedRank);
    const singleRank = availableRanks.length === 1;

    container.innerHTML = '';
    card.style.display = '';
    container.offsetWidth;

    // (Rank label now shown by the per-rank selector chips via renderModelSelectors.)

    // Get LIR/cLIR for current rank
    const target = isAvg ? flyAverageRow : (flyParsed.rows.find(r => parseInt(r['Rank']) === selRankNum) || flyParsed.rows[0]);

    for (const chainName of chainNames) {
        const chainSize = firstPlddt.get(chainName)?.size || 0;
        if (chainSize === 0) continue;
        const maxRes = Math.max(...firstPlddt.get(chainName).keys());

        // Build per-rank pLDDT arrays
        const rankArrays = [];
        for (const rp of allPlddt) {
            const cm = rp.chains.get(chainName);
            const arr = [];
            for (let r = 1; r <= maxRes; r++) {
                arr.push(cm?.get(r) ?? null);
            }
            rankArrays.push({ rank: rp.rank, plddt: arr });
        }

        // Compute average and std
        const n = maxRes;
        const avgArr = [], stdArr = [];
        for (let i = 0; i < n; i++) {
            let sum = 0, cnt = 0;
            for (const ra of rankArrays) { if (ra.plddt[i] !== null) { sum += ra.plddt[i]; cnt++; } }
            const avg = cnt > 0 ? sum / cnt : null;
            avgArr.push(avg);
            let sqSum = 0;
            if (cnt > 1) {
                for (const ra of rankArrays) { if (ra.plddt[i] !== null) sqSum += (ra.plddt[i] - avg) ** 2; }
                stdArr.push(Math.sqrt(sqSum / cnt));
            } else { stdArr.push(0); }
        }

        // Get LIR/cLIR
        const lirKey = chainName === 'A' ? '_LIR_A' : '_LIR_B';
        const clirKey = chainName === 'A' ? '_cLIR_A' : '_cLIR_B';
        const lir = target?.[lirKey] || new Set();
        const clir = target?.[clirKey] || new Set();
        const lirColor = chainName === 'A' ? colors.lir_a : colors.lir_b;
        const clirColor = chainName === 'A' ? colors.clir_a : colors.clir_b;

        // Selected rank data
        const selData = selRankNum ? rankArrays.find(ra => ra.rank === selRankNum)?.plddt : avgArr;

        // Avg pLDDT for legend
        let selSum = 0, selN = 0, avgS = 0, avgN2 = 0;
        for (let i = 0; i < n; i++) {
            if (selData && selData[i] !== null) { selSum += selData[i]; selN++; }
            if (avgArr[i] !== null) { avgS += avgArr[i]; avgN2++; }
        }
        const selAvg = selN > 0 ? (selSum / selN).toFixed(1) : '—';
        const totalAvg = avgN2 > 0 ? (avgS / avgN2).toFixed(1) : '—';

        // Chain header
        const header = document.createElement('div');
        header.style.cssText = 'font-size:0.9rem; font-weight:700; color:#1a5276; margin-top:1rem; margin-bottom:0.2rem; padding-bottom:0.3rem; border-bottom:2px solid #e0e0e0;';
        { const _g = chainLabels[chainName]; header.textContent = 'Chain ' + chainName + (_g && _g !== chainName ? ' (' + _g + ')' : '') + ' · ' + n + ' residues'; }
        container.appendChild(header);

        // Legend
        const legendEl = document.createElement('div');
        legendEl.style.cssText = 'display:flex; align-items:center; gap:1.5rem; flex-wrap:wrap; margin-top:0.5rem; margin-bottom:0.3rem;';
        let legendHtml = '<span style="font-size:0.85rem; font-weight:600; color:#1a5276;">pLDDT</span>';
        if (!singleRank) {
            legendHtml += '<span style="display:flex; align-items:center; gap:0.3rem; font-size:0.75rem; color:#333;"><span style="display:inline-block; width:20px; height:2px; background:#C62828;"></span> ' + (isAvg ? 'Average' : 'Rank ' + selRankNum) + ' (' + selAvg + ')</span>';
            legendHtml += '<span style="display:flex; align-items:center; gap:0.3rem; font-size:0.75rem; color:#333;"><span style="display:inline-block; width:20px; height:0; border-top:2px dashed #1565C0;"></span> Avg \u00b1 SD (' + totalAvg + ')</span>';
        } else {
            legendHtml += '<span style="font-size:0.75rem; color:#555;">(' + selAvg + ')</span>';
        }
        legendEl.innerHTML = legendHtml;
        container.appendChild(legendEl);

        // Canvas
        const W = Math.max(Math.min(container.clientWidth || 900, 1200), 400);
        const H = 145;
        const PAD = { left: 50, right: 20, top: 14, bottom: 35 };
        const plotW = W - PAD.left - PAD.right;
        const plotH = H - PAD.top - PAD.bottom;

        const canvas = document.createElement('canvas');
        canvas.width = W * 2; canvas.height = H * 2;
        canvas.style.cssText = `width:${W}px; height:${H}px;`;
        container.appendChild(canvas);

        const ctx = canvas.getContext('2d');
        ctx.scale(2, 2);
        const xScale = (r) => PAD.left + (r / n) * plotW;
        const yPlddt = (v) => PAD.top + plotH - (Math.max(0, Math.min(100, v)) / 100) * plotH;

        // Background
        ctx.fillStyle = '#fafafa';
        ctx.fillRect(PAD.left, PAD.top, plotW, plotH);

        // LIR/cLIR bar
        const barY = H - PAD.bottom + 2, barH = 6;
        for (let r = 1; r <= n; r++) {
            const x0 = xScale(r - 1), x1 = xScale(r);
            if (clir.has(r)) { ctx.fillStyle = clirColor; ctx.fillRect(x0, barY, x1 - x0 + 0.5, barH); }
            else if (lir.has(r)) { ctx.fillStyle = lirColor; ctx.fillRect(x0, barY, x1 - x0 + 0.5, barH); }
        }

        // Y-axis gridlines
        ctx.strokeStyle = '#e0e0e0'; ctx.lineWidth = 0.5;
        ctx.fillStyle = '#999'; ctx.font = '10px sans-serif'; ctx.textAlign = 'right';
        for (const tick of [0, 25, 50, 75, 100]) {
            const y = yPlddt(tick);
            ctx.beginPath(); ctx.moveTo(PAD.left, y); ctx.lineTo(PAD.left + plotW, y); ctx.stroke();
            ctx.fillText(tick, PAD.left - 5, y + 3);
        }

        // X-axis labels
        ctx.textAlign = 'center'; ctx.fillStyle = '#999';
        const step = n > 800 ? 200 : n > 400 ? 100 : n > 200 ? 50 : 25;
        for (let r = step; r <= n; r += step) ctx.fillText(r, xScale(r), H - PAD.bottom + 22);

        // Avg ± SD shading
        if (!singleRank) {
            ctx.fillStyle = 'rgba(21,101,192,0.12)';
            ctx.beginPath();
            for (let i = 0; i < n; i++) {
                if (avgArr[i] === null) continue;
                ctx.lineTo(xScale(i + 0.5), yPlddt(avgArr[i] + stdArr[i]));
            }
            for (let i = n - 1; i >= 0; i--) {
                if (avgArr[i] === null) continue;
                ctx.lineTo(xScale(i + 0.5), yPlddt(avgArr[i] - stdArr[i]));
            }
            ctx.closePath();
            ctx.fill();

            // Average line (dashed)
            ctx.strokeStyle = '#1565C0'; ctx.lineWidth = 1;
            ctx.setLineDash([4, 3]);
            ctx.beginPath();
            let started = false;
            for (let i = 0; i < n; i++) {
                if (avgArr[i] === null) { started = false; continue; }
                const x = xScale(i + 0.5), y = yPlddt(avgArr[i]);
                if (!started) { ctx.moveTo(x, y); started = true; } else ctx.lineTo(x, y);
            }
            ctx.stroke();
            ctx.setLineDash([]);
        }

        // Selected rank (solid)
        const drawArr = selData || avgArr;
        ctx.strokeStyle = '#C62828'; ctx.lineWidth = 1;
        ctx.beginPath();
        let started2 = false;
        for (let i = 0; i < n; i++) {
            if (drawArr[i] === null) { started2 = false; continue; }
            const x = xScale(i + 0.5), y = yPlddt(drawArr[i]);
            if (!started2) { ctx.moveTo(x, y); started2 = true; } else ctx.lineTo(x, y);
        }
        ctx.stroke();

        // Y-axis label
        ctx.save();
        ctx.translate(12, PAD.top + plotH / 2);
        ctx.rotate(-Math.PI / 2);
        ctx.fillStyle = '#666'; ctx.font = '11px sans-serif'; ctx.textAlign = 'center';
        ctx.fillText('pLDDT', 0, 0);
        ctx.restore();
    }
    LiviaMaps.attachChartPngButtons('plddt-charts', 'plddt_chart');
}

function buildConnectionDiagram(row, sym1, sym2, colors, pdbText) {
    const card = document.getElementById('connection-diagram-card');
    const canvas = document.getElementById('connection-canvas');
    if (!card || !canvas || !pdbText) { if (card) card.style.display = 'none'; return; }
    // Panel options (parity with universal.html): swap chains, domain show/hide, adjustable x-ticks.
    const swap = !!(document.getElementById('contact-swap-ab') && document.getElementById('contact-swap-ab').checked);
    const showDom = document.getElementById('contact-show-domains') ? document.getElementById('contact-show-domains').checked : true;   // UniProt/Pfam
    const showTed = document.getElementById('contact-show-ted') ? document.getElementById('contact-show-ted').checked : true;           // TED
    const showPlddt = document.getElementById('contact-show-plddt') ? document.getElementById('contact-show-plddt').checked : true;   // per-residue pLDDT strip
    // The export closure re-runs this function; hand it the ORIGINAL args so the swap isn't applied twice.
    const _xp = { row, sym1, sym2, colors, pdbText };
    LiviaMaps.attachExportBar('connection-canvas', { name: 'linear_contact_map', svg: () => buildConnectionDiagram(_xp.row, _xp.sym1, _xp.sym2, _xp.colors, _xp.pdbText), png: true });
    window.__redrawConn = () => buildConnectionDiagram(_xp.row, _xp.sym1, _xp.sym2, _xp.colors, _xp.pdbText);   // panel-control re-render
    if (swap) { const _s = sym1; sym1 = sym2; sym2 = _s; colors = { lir_a: colors.lir_b, clir_a: colors.clir_b, lir_b: colors.lir_a, clir_b: colors.clir_a }; }

    const cLIR_A = swap ? row._cLIR_B : row._cLIR_A, cLIR_B = swap ? row._cLIR_A : row._cLIR_B;
    const LIR_A = swap ? row._LIR_B : row._LIR_A, LIR_B = swap ? row._LIR_A : row._LIR_B;
    const lenA = parseInt((swap ? (row['Len B'] || row['protein_2_size']) : (row['Len A'] || row['protein_1_size']))) || 0;
    const lenB = parseInt((swap ? (row['Len A'] || row['protein_1_size']) : (row['Len B'] || row['protein_2_size']))) || 0;
    // Per-side domain data honoring swap + show/hide (empty when hidden → no domain rows / legend)
    const _dmA = showDom ? (swap ? flyDomains.B : flyDomains.A) : [], _dmB = showDom ? (swap ? flyDomains.A : flyDomains.B) : [];
    const _tdA = showTed ? (swap ? flyTedDomains.B : flyTedDomains.A) : [], _tdB = showTed ? (swap ? flyTedDomains.A : flyTedDomains.B) : [];
    const _diA = (showDom || showTed) ? (swap ? flyDomainInfo.B : flyDomainInfo.A) : null, _diB = (showDom || showTed) ? (swap ? flyDomainInfo.A : flyDomainInfo.B) : null;
    if (lenA === 0 || lenB === 0 || cLIR_A.size === 0 || cLIR_B.size === 0) {
        card.style.display = 'none';
        return;
    }

    // Parse 3D coordinates
    const coordChains = parsePdbCbCoords(pdbText);
    const coordA = coordChains.get('A');
    const coordB = coordChains.get('B');
    if (!coordA || !coordB) { card.style.display = 'none'; return; }

    // Per-residue pLDDT (B-factors), honouring swap: top bar = real B when swapped
    const _plBfs = parseBfactorsPerResidue(pdbText);
    const _plArr = (ch, len) => { const a = new Array(len).fill(null); for (let r = 1; r <= len; r++) { const val = _plBfs.get(`${ch}:${r}`); if (val != null) a[r - 1] = val; } return a; };
    const plddtA = _plArr(swap ? 'B' : 'A', lenA), plddtB = _plArr(swap ? 'A' : 'B', lenB);
    const plStripH = 5, plBarGap = 2, plGap = 2;
    const plOnA = showPlddt && plddtA.some(x => x != null), plOnB = showPlddt && plddtB.some(x => x != null);
    const plBlockA = plOnA ? (plBarGap + plStripH + plGap) : 0, plBlockB = plOnB ? (plBarGap + plStripH + plGap) : 0;

    // Find cLIR-cLIR contact pairs (Cβ distance ≤ 8 Å)
    const contactDist = 8;
    const contacts = [];
    for (const resA of row._cLIR_A) {            // detection uses the REAL chain-A/B residues + coords (display swap is applied later)
        const cA = coordA.get(resA);
        if (!cA) continue;
        for (const resB of row._cLIR_B) {
            const cB = coordB.get(resB);
            if (!cB) continue;
            const dx = cA.x - cB.x, dy = cA.y - cB.y, dz = cA.z - cB.z;
            const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
            if (dist <= contactDist) {
                contacts.push({ resA, resB, dist });
            }
        }
    }

    if (contacts.length === 0) { card.style.display = 'none'; return; }
    card.style.display = '';

    // Drawing
    const S = 2; // HiDPI
    const hasDomains = (_dmA.length > 0 || _dmB.length > 0);
    const hasTed = (_tdA.length > 0 || _tdB.length > 0);
    const domRowH = hasDomains ? 20 : 0;
    const tedRowH = hasTed ? 20 : 0;
    const topAnno = (_dmA.length ? 20 : 0) + (_tdA.length ? 20 : 0);   // top-chain rows (above the ticks)
    const botAnno = (_dmB.length ? 20 : 0) + (_tdB.length ? 20 : 0);   // bottom-chain rows (below the ticks)
    const tickH = 14;                                    // room for residue-number ticks (kept off the ribbons)
    const width = Math.min(900, (card.clientWidth || 900) - 32);
    const barH = 12;                                     // compact bars (was 14)
    const barY_A = 10 + topAnno + plBlockA + tickH;                   // [margin][domains][ticks][bar]
    const barY_B = barY_A + barH + 60;                   // compact ribbon gap — ribbons only (was ~186)
    const height = barY_B + barH + tickH + plBlockB + botAnno + 8;   // count text now lives in the HTML below, not the canvas
    const barRight = width - 20;

    canvas.width = width * S;
    canvas.height = height * S;
    // Responsive display: scale to fit the container (capped at the drawn width) so the
    // map never overflows on narrow windows and rescales on resize.
    canvas.style.width = '100%';
    canvas.style.height = 'auto';
    canvas.style.maxWidth = width + 'px';
    const ctx = canvas.getContext('2d');
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.save();
    ctx.scale(S, S);

    // Size the left margin to fit the (gene-symbol) chain labels so they aren't clipped
    // off the left edge; the font shrinks as a fallback for very long names when drawn.
    ctx.font = 'bold 13px sans-serif';
    const labelA = `${sym1}${window.showConnIDs !== false ? ` (${swap ? 'B' : 'A'})` : ''}`, labelB = `${sym2}${window.showConnIDs !== false ? ` (${swap ? 'A' : 'B'})` : ''}`;
    const barLeft = Math.max(60, Math.min(170, Math.ceil(Math.max(ctx.measureText(labelA).width, ctx.measureText(labelB).width)) + 12));
    const barW = barRight - barLeft;

    // Helper: residue to x position
    const xA = (res) => barLeft + ((res - 1) / Math.max(1, lenA - 1)) * barW;
    const xB = (res) => barLeft + ((res - 1) / Math.max(1, lenB - 1)) * barW;

    // Paint contiguous LIR/cLIR runs as solid blocks. One rect per residue leaves hairline
    // gray gaps between adjacent residues (the "lines between residues" barcode artifact);
    // merging runs via setToRanges gives clean blocks like the chord map.
    function paintBarRanges(yBar, set, xFn, len, color) {
        const resW = barW / Math.max(1, len - 1);
        ctx.fillStyle = color;
        for (const [s, e] of setToRanges(set)) {
            const x0 = Math.max(barLeft, xFn(s) - resW / 2);
            const x1 = Math.min(barRight, xFn(e) + resW / 2);
            ctx.fillRect(x0, yBar, Math.max(1.5, x1 - x0), barH);
        }
    }

    // Contiguous LIR/cLIR runs merged into solid blocks (via setToRanges in paintBarRanges),
    // matching the circular map. cLIR is NOT gap-filled, so it stays exact.
    // Draw Chain A bar
    ctx.fillStyle = '#f0f0f0';
    ctx.fillRect(barLeft, barY_A, barW, barH);
    paintBarRanges(barY_A, LIR_A, xA, lenA, colors.lir_a);
    paintBarRanges(barY_A, cLIR_A, xA, lenA, colors.clir_a);
    ctx.strokeStyle = '#999';
    ctx.lineWidth = 0.5;
    ctx.strokeRect(barLeft, barY_A, barW, barH);

    // Draw Chain B bar
    ctx.fillStyle = '#f0f0f0';
    ctx.fillRect(barLeft, barY_B, barW, barH);
    paintBarRanges(barY_B, LIR_B, xB, lenB, colors.lir_b);
    paintBarRanges(barY_B, cLIR_B, xB, lenB, colors.clir_b);
    ctx.strokeStyle = '#999';
    ctx.lineWidth = 0.5;
    ctx.strokeRect(barLeft, barY_B, barW, barH);

    // pLDDT strips (per-residue AlphaFold confidence bands), hugging each bar with a small gap; runs RLE-merged
    const _plBand = (val) => (val == null || isNaN(val)) ? -1 : (val > 90 ? 0 : val > 70 ? 1 : val > 50 ? 2 : 3);
    const _plCOL = ['#0053D6', '#65CBF3', '#FFDB13', '#FF7D45'];
    const _drawPlStrip = (arr, len, xFn, yStrip) => {
        const n = Math.min(len, arr.length), resW = barW / Math.max(1, len - 1);
        const seg = (a, b, bd) => { if (bd < 0) return; const x0 = Math.max(barLeft, xFn(a) - resW / 2), x1 = Math.min(barRight, xFn(b) + resW / 2); ctx.fillStyle = _plCOL[bd]; ctx.fillRect(x0, yStrip, Math.max(1, x1 - x0), plStripH); };
        if (!n) return;
        let r0 = 1, cur = _plBand(arr[0]);
        for (let r = 2; r <= n; r++) { const b = _plBand(arr[r - 1]); if (b !== cur) { seg(r0, r - 1, cur); r0 = r; cur = b; } }
        seg(r0, n, cur);
        ctx.strokeStyle = '#ccc'; ctx.lineWidth = 0.5; ctx.strokeRect(barLeft, yStrip, barW, plStripH);
    };
    if (plOnA) _drawPlStrip(plddtA, lenA, xA, barY_A - plBarGap - plStripH);
    if (plOnB) _drawPlStrip(plddtB, lenB, xB, barY_B + barH + plBarGap);

    // Labels — right-aligned in the left margin; shrink font to fit if the name is long.
    const drawChainLabel = (text, y, color) => {
        let fs = 13; ctx.font = `bold ${fs}px sans-serif`;
        while (ctx.measureText(text).width > barLeft - 4 && fs > 8) { fs--; ctx.font = `bold ${fs}px sans-serif`; }
        ctx.fillStyle = color; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
        ctx.fillText(text, barLeft - 5, y);
    };
    drawChainLabel(labelA, barY_A + barH / 2, colors.clir_a);
    drawChainLabel(labelB, barY_B + barH / 2, colors.clir_b);

    // Residue number ticks
    ctx.fillStyle = '#888';
    ctx.font = '10px sans-serif';
    ctx.textAlign = 'center';
    ctx.textBaseline = 'bottom';
    const nTicks = Math.max(2, Math.min(30, parseInt(document.getElementById('contact-xticks')?.value) || 5));
    const tickPos = (len) => { const n = Math.max(2, Math.min(nTicks, len)); const s = new Set(); for (let i = 0; i < n; i++) s.add(Math.max(1, Math.min(len, Math.round(1 + (len - 1) * i / (n - 1))))); return [...s]; };
    for (const t of tickPos(lenA)) ctx.fillText(t, xA(t), barY_A - plBlockA - 3);            // top-bar ticks ABOVE
    ctx.textBaseline = 'top';
    for (const t of tickPos(lenB)) ctx.fillText(t, xB(t), barY_B + barH + plBlockB + 3);     // bottom-bar ticks BELOW

    // Draw connection lines
    const minDist = Math.min(...contacts.map(c => c.dist));
    const maxDist = Math.max(...contacts.map(c => c.dist));

    // Parse colors to rgba
    const colA = colors.clir_a || '#00897B';
    const colB = colors.clir_b || '#E64A19';

    for (const c of contacts) {
        const ax = xA(swap ? c.resB : c.resA);
        const bx = xB(swap ? c.resA : c.resB);
        const t = maxDist > minDist ? 1 - (c.dist - minDist) / (maxDist - minDist) : 1;
        const alpha = 0.08 + t * 0.2;

        const midY = (barY_A + barH + barY_B) / 2;
        LiviaMaps.strokeFadeCubic(ctx, ax, barY_A + barH, ax, midY, bx, midY, bx, barY_B, hexToRgba(colA, alpha), hexToRgba(colB, alpha), 1);
    }

    // Contact-count line moved to the HTML below (kept out of the exported figure)
    { const _cc = document.getElementById('contact-count-text'); if (_cc) _cc.textContent = `${contacts.length} cLIR contact pairs (C\u03B2 \u2264 ${contactDist} \u00C5)`; }

    // Domain annotations
    const domPalette = ['rgba(143,169,196,0.45)', 'rgba(117,147,179,0.45)'];
    const domBorderPalette = ['rgba(75,85,99,0.8)', 'rgba(75,85,99,0.8)'];
    if (hasDomains) {
        const domBlockH = 14, domGap = 2;
        // Chain A domains above bar
        for (let di = 0; di < _dmA.length; di++) {
            const d = _dmA[di];
            const x1 = xA(Math.max(1, d.start)), x2 = xA(Math.min(lenA, d.end));
            const dw = x2 - x1, dy = barY_A - plBlockA - tickH - domGap - domBlockH;
            ctx.fillStyle = domPalette[di % domPalette.length];
            ctx.fillRect(x1, dy, dw, domBlockH);
            ctx.strokeStyle = domBorderPalette[di % domBorderPalette.length];
            ctx.lineWidth = 0.8;
            ctx.strokeRect(x1, dy, dw, domBlockH);
            ctx.fillStyle = '#444'; ctx.font = '9px sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
            const label = d.name, labelW = ctx.measureText(label).width;
            if (labelW < dw - 2) ctx.fillText(label, x1 + dw/2, dy + domBlockH/2);
            else { const abbr = label.length > 6 ? label.substring(0,5)+'..' : label; if (ctx.measureText(abbr).width < dw-2) ctx.fillText(abbr, x1+dw/2, dy+domBlockH/2); }
        }
        // Chain B domains below bar
        for (let di = 0; di < _dmB.length; di++) {
            const d = _dmB[di];
            const x1 = xB(Math.max(1, d.start)), x2 = xB(Math.min(lenB, d.end));
            const dw = x2 - x1, dy = barY_B + barH + plBlockB + tickH + domGap;
            ctx.fillStyle = domPalette[di % domPalette.length];
            ctx.fillRect(x1, dy, dw, domBlockH);
            ctx.strokeStyle = domBorderPalette[di % domBorderPalette.length];
            ctx.lineWidth = 0.8;
            ctx.strokeRect(x1, dy, dw, domBlockH);
            ctx.fillStyle = '#444'; ctx.font = '9px sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
            const label = d.name, labelW = ctx.measureText(label).width;
            if (labelW < dw - 2) ctx.fillText(label, x1 + dw/2, dy + domBlockH/2);
            else { const abbr = label.length > 6 ? label.substring(0,5)+'..' : label; if (ctx.measureText(abbr).width < dw-2) ctx.fillText(abbr, x1+dw/2, dy+domBlockH/2); }
        }
    }

    // Draw TED domain annotations
    if (hasTed) {
        const TED_PALETTE = ['#A898C8', '#8E7DB8'];
        const tedBlockH = 14, tedGap = 2;

        // Chain A TED domains — above UniProt domains (further above chain A bar)
        for (let di = 0; di < _tdA.length; di++) {
            const d = _tdA[di];
            for (const seg of d.segments) {
                const x1 = xA(Math.max(1, seg.start)), x2 = xA(Math.min(lenA, seg.end));
                const dw = x2 - x1;
                const dy = barY_A - plBlockA - tickH - tedGap - tedBlockH - (_dmA.length ? tedBlockH + tedGap : 0);
                ctx.fillStyle = TED_PALETTE[di % TED_PALETTE.length];
                ctx.globalAlpha = 0.6;
                ctx.fillRect(x1, dy, dw, tedBlockH);
                ctx.globalAlpha = 1;
                ctx.strokeStyle = TED_PALETTE[di % TED_PALETTE.length];
                ctx.lineWidth = 0.8;
                ctx.strokeRect(x1, dy, dw, tedBlockH);
            }
            const mainSeg = d.segments.reduce((a, b) => (b.end - b.start > a.end - a.start ? b : a), d.segments[0]);
            const x1 = xA(Math.max(1, mainSeg.start)), x2 = xA(Math.min(lenA, mainSeg.end));
            const dw = x2 - x1;
            const dy = barY_A - plBlockA - tickH - tedGap - tedBlockH - (_dmA.length ? tedBlockH + tedGap : 0);
            ctx.fillStyle = '#333'; ctx.font = '9px sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
            const label = d.cath ? d.cath : `D${d.tedNo}`;
            const labelW = ctx.measureText(label).width;
            if (labelW < dw - 2) ctx.fillText(label, x1 + dw / 2, dy + tedBlockH / 2);
            else { const abbr = label.length > 6 ? label.substring(0, 5) + '..' : label; if (ctx.measureText(abbr).width < dw - 2) ctx.fillText(abbr, x1 + dw / 2, dy + tedBlockH / 2); }
        }

        // Chain B TED domains — below UniProt domains (further below chain B bar)
        for (let di = 0; di < _tdB.length; di++) {
            const d = _tdB[di];
            for (const seg of d.segments) {
                const x1 = xB(Math.max(1, seg.start)), x2 = xB(Math.min(lenB, seg.end));
                const dw = x2 - x1;
                const dy = barY_B + barH + plBlockB + tickH + tedGap + (_dmB.length ? tedBlockH + tedGap : 0);
                ctx.fillStyle = TED_PALETTE[di % TED_PALETTE.length];
                ctx.globalAlpha = 0.6;
                ctx.fillRect(x1, dy, dw, tedBlockH);
                ctx.globalAlpha = 1;
                ctx.strokeStyle = TED_PALETTE[di % TED_PALETTE.length];
                ctx.lineWidth = 0.8;
                ctx.strokeRect(x1, dy, dw, tedBlockH);
            }
            const mainSeg = d.segments.reduce((a, b) => (b.end - b.start > a.end - a.start ? b : a), d.segments[0]);
            const x1 = xB(Math.max(1, mainSeg.start)), x2 = xB(Math.min(lenB, mainSeg.end));
            const dw = x2 - x1;
            const dy = barY_B + barH + plBlockB + tickH + tedGap + (_dmB.length ? tedBlockH + tedGap : 0);
            ctx.fillStyle = '#333'; ctx.font = '9px sans-serif'; ctx.textAlign = 'center'; ctx.textBaseline = 'middle';
            const label = d.cath ? d.cath : `D${d.tedNo}`;
            const labelW = ctx.measureText(label).width;
            if (labelW < dw - 2) ctx.fillText(label, x1 + dw / 2, dy + tedBlockH / 2);
            else { const abbr = label.length > 6 ? label.substring(0, 5) + '..' : label; if (ctx.measureText(abbr).width < dw - 2) ctx.fillText(abbr, x1 + dw / 2, dy + tedBlockH / 2); }
        }
    }

    ctx.restore();

    // Domain legend — each chain on its own row, with partial protein note
    const legendEl = document.getElementById('fly-domain-legend');
    if (legendEl) {
        legendEl.innerHTML = '';
        legendEl.style.flexDirection = 'column';
        legendEl.style.gap = '0.3rem';
        const TED_PALETTE_LEG = ['#A898C8', '#8E7DB8'];
        const chainLabels = [
            { key: 'A', label: `${sym1} (${swap ? 'B' : 'A'})`, domains: _dmA, tedDomains: _tdA, info: _diA },
            { key: 'B', label: `${sym2} (${swap ? 'A' : 'B'})`, domains: _dmB, tedDomains: _tdB, info: _diB },
        ];
        for (const cl of chainLabels) {
            if (cl.domains.length === 0 && cl.tedDomains.length === 0 && !cl.info) continue;
            const row = document.createElement('div');
            row.style.cssText = 'display:flex; gap:0.5rem; flex-wrap:wrap; align-items:center;';
            const chainLabel = document.createElement('span');
            chainLabel.style.cssText = 'font-weight:600; color:#444;';
            chainLabel.textContent = cl.label + ':';
            row.appendChild(chainLabel);
            // Partial protein note
            if (cl.info && cl.info.isPartial) {
                const note = document.createElement('span');
                note.style.cssText = 'font-size:0.68rem; color:#c0392b; font-style:italic;';
                note.textContent = `partial (aa ${cl.info.startInFull}\u2013${cl.info.endInFull} of ${cl.info.fullLength})`;
                row.appendChild(note);
            }
            // Domain items with full-length positions
            for (let di = 0; di < cl.domains.length; di++) {
                const d = cl.domains[di];
                const item = document.createElement('span');
                // Show predicted positions, and full-length positions if partial
                let posText = `${d.start}\u2013${d.end}`;
                if (cl.info && cl.info.offset > 0) {
                    const fullStart = d.start + cl.info.offset;
                    const fullEnd = d.end + cl.info.offset;
                    posText += ` [${fullStart}\u2013${fullEnd} in full]`;
                }
                item.innerHTML = `<span style="display:inline-block; width:10px; height:10px; background:${domPalette[di % domPalette.length]}; border:1px solid ${domBorderPalette[di % domBorderPalette.length]}; vertical-align:middle; margin-right:3px; border-radius:2px;"></span>${d.name} (${posText})`;
                row.appendChild(item);
            }
            // TED domain legend items
            for (let di = 0; di < cl.tedDomains.length; di++) {
                const d = cl.tedDomains[di];
                const item = document.createElement('span');
                const color = TED_PALETTE_LEG[di % TED_PALETTE_LEG.length];
                let posText = `${d.start}\u2013${d.end}`;
                if (cl.info && cl.info.offset > 0) {
                    posText += ` [${d.start + cl.info.offset}\u2013${d.end + cl.info.offset} in full]`;
                }
                item.innerHTML = `<span style="display:inline-block; width:10px; height:10px; background:${color}99; border:1px solid ${color}; vertical-align:middle; margin-right:3px; border-radius:2px;"></span>${d.name} (${posText})`;
                row.appendChild(item);
            }
            // UniProt link
            if (cl.info && cl.info.accession) {
                const link = document.createElement('a');
                link.href = `https://www.uniprot.org/uniprot/${cl.info.accession}`;
                link.target = '_blank';
                link.style.cssText = 'font-size:0.68rem; color:#2471A3;';
                link.textContent = `UniProt: ${cl.info.accession}`;
                row.appendChild(link);
            }
            legendEl.appendChild(row);
        }
    }
    if (showPlddt && (plOnA || plOnB) && legendEl) legendEl.innerHTML = '<div style="width:100%; display:flex; flex-wrap:wrap; align-items:center; gap:0.5rem;">' + LiviaMaps.plddtLegendHtml() + '</div>' + legendEl.innerHTML;
}

function drawFlyChord(canvas, cLIR_A, cLIR_B, lenA, lenB, sym1, sym2, colors, coordA, coordB, skipDistFilter, LIR_A, LIR_B) {
    const S = 2;
    const size = 350;
    canvas.width = size * S; canvas.height = size * S;
    canvas.style.width = size + 'px'; canvas.style.height = size + 'px';
    const ctx = canvas.getContext('2d');
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.save();
    ctx.scale(S, S);

    const cx = size / 2, cy = size / 2;
    const _plddtOn = flyShowChordPlddt && ((flyPlddt.A && flyPlddt.A.some(x => x != null)) || (flyPlddt.B && flyPlddt.B.some(x => x != null)));   // pLDDT ring only when data present
    const _nDomRings = (_plddtOn ? 1 : 0) + (flyShowChordDom ? 1 : 0) + (flyShowChordTed ? 1 : 0), _domSpace = _nDomRings ? (_nDomRings * 7 + 3) : 0, _domHits = [];   // pLDDT + domain rings sit OUTSIDE the arc; reserve adaptive outer space
    const outerR = Math.max(40, Math.min(size / 2 - (flyShowChordTicks ? 62 : 48) - _domSpace, size / 2 - 4 - _domSpace - LiviaMaps.chordLabelExtent(ctx, [
        { name: flyChordLabel(sym1, 'A'), length: flyShowChordLength ? `${lenA}` : '' },
        { name: flyChordLabel(sym2, 'B'), length: flyShowChordLength ? `${lenB}` : '' },
    ], { ticks: flyShowChordTicks, radial: flyShowChordRadial }))), innerR = outerR - 12, chordR = innerR - 4;   // shrink to keep long radial labels on-canvas
    ctx.canvas._domHits = _domHits; ctx.canvas._domGeom = { cx, cy, size };   // hover
    const gap = 0.04;
    const totalLen = lenA + lenB;
    const totalAngle = 2 * Math.PI - 2 * gap;
    const angleA = (lenA / totalLen) * totalAngle;
    const angleB = (lenB / totalLen) * totalAngle;

    const startA = -Math.PI / 2;
    const endA = startA + angleA;
    const startB = endA + gap;
    const endB = startB + angleB;

    const colA = colors.clir_a || '#00897B';
    const colB = colors.clir_b || '#E64A19';

    const resAngleA = (res) => startA + ((res - 1) / Math.max(1, lenA - 1)) * angleA;
    const resAngleB = (res) => startB + ((res - 1) / Math.max(1, lenB - 1)) * angleB;

    // Draw contact lines
    if (skipDistFilter) {
        // All cLIR pairs without distance filtering
        for (const resA of cLIR_A) {
            for (const resB of cLIR_B) {
                const aA = resAngleA(resA), aB = resAngleB(resB);
                const xA = cx + Math.cos(aA) * chordR, yA = cy + Math.sin(aA) * chordR;
                const xB = cx + Math.cos(aB) * chordR, yB = cy + Math.sin(aB) * chordR;
                LiviaMaps.strokeFadeQuad(ctx, xA, yA, cx, cy, xB, yB, hexToRgba(colA, 0.12), hexToRgba(colB, 0.12), 1);
            }
        }
    } else {
        const contactDist = 8;
        const contacts = [];
        for (const resA of cLIR_A) {
            const cA = coordA.get(resA);
            if (!cA) continue;
            for (const resB of cLIR_B) {
                const cB = coordB.get(resB);
                if (!cB) continue;
                const dx = cA.x - cB.x, dy = cA.y - cB.y, dz = cA.z - cB.z;
                const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
                if (dist <= contactDist) contacts.push({ resA, resB, dist });
            }
        }
        const minDist = contacts.length > 0 ? Math.min(...contacts.map(c => c.dist)) : 0;
        const maxDist = contacts.length > 0 ? Math.max(...contacts.map(c => c.dist)) : 1;
        for (const c of contacts) {
            const aA = resAngleA(c.resA), aB = resAngleB(c.resB);
            const pA = { x: cx + Math.cos(aA) * chordR, y: cy + Math.sin(aA) * chordR };
            const pB = { x: cx + Math.cos(aB) * chordR, y: cy + Math.sin(aB) * chordR };
            const t = maxDist > minDist ? 1 - (c.dist - minDist) / (maxDist - minDist) : 1;
            const alpha = 0.08 + t * 0.2;
            LiviaMaps.strokeFadeQuad(ctx, pA.x, pA.y, cx, cy, pB.x, pB.y, hexToRgba(colA, alpha), hexToRgba(colB, alpha), 0.7);
        }
    }

    // Draw arcs (gray base + LIR light + cLIR dark, matching seq viewer / linear map)
    const lightA = colors.lir_a || colA, lightB = colors.lir_b || colB;
    paintChordArcBand(ctx, cx, cy, outerR, innerR, startA, endA, lenA, LIR_A, cLIR_A, lightA, colA);
    paintChordArcBand(ctx, cx, cy, outerR, innerR, startB, endB, lenB, LIR_B, cLIR_B, lightB, colB);
    { let _rA = outerR + 3, _rB = outerR + 3;   // pLDDT ring first, then domain rings (UniProt/Pfam, TED) stacked outward, with hover hit-regions
      if (_plddtOn) { for (const _h of LiviaMaps.paintChordPlddtRing(ctx, cx, cy, _rA, _rA + 6, startA, endA, lenA, flyPlddt.A, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym1, 'A'), src: 'pLDDT' }, _h)); _rA += 7;
                     for (const _h of LiviaMaps.paintChordPlddtRing(ctx, cx, cy, _rB, _rB + 6, startB, endB, lenB, flyPlddt.B, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym2, 'B'), src: 'pLDDT' }, _h)); _rB += 7; }
      if (flyShowChordDom) { for (const _h of LiviaMaps.paintChordDomainRing(ctx, cx, cy, _rA, _rA + 6, startA, endA, lenA, flyDomains.A, CHORD_DOM_FILL, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym1, 'A'), src: 'UniProt/Pfam' }, _h)); _rA += 7;
                             for (const _h of LiviaMaps.paintChordDomainRing(ctx, cx, cy, _rB, _rB + 6, startB, endB, lenB, flyDomains.B, CHORD_DOM_FILL, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym2, 'B'), src: 'UniProt/Pfam' }, _h)); _rB += 7; }
      if (flyShowChordTed) { for (const _h of LiviaMaps.paintChordDomainRing(ctx, cx, cy, _rA, _rA + 6, startA, endA, lenA, flyTedDomains.A, CHORD_TED_FILL, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym1, 'A'), src: 'TED' }, _h));
                             for (const _h of LiviaMaps.paintChordDomainRing(ctx, cx, cy, _rB, _rB + 6, startB, endB, lenB, flyTedDomains.B, CHORD_TED_FILL, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym2, 'B'), src: 'TED' }, _h)); } }
    if (flyShowChordTicks) { LiviaMaps.drawChordRuler(ctx, cx, cy, outerR + _domSpace, startA, endA, lenA, {}); LiviaMaps.drawChordRuler(ctx, cx, cy, outerR + _domSpace, startB, endB, lenB, {}); }

    // Labels
    for (const arc of [
        { mid: (startA + endA) / 2, name: flyChordLabel(sym1, 'A'), size: lenA },
        { mid: (startB + endB) / 2, name: flyChordLabel(sym2, 'B'), size: lenB }
    ]) {
        LiviaMaps.drawChordArcLabel(ctx, {
            cx, cy, outerR: outerR + _domSpace, angle: arc.mid, ticks: flyShowChordTicks, radial: flyShowChordRadial, hidden: false,
            name: arc.name, showName: true,
            length: `${arc.size}`, showLength: flyShowChordLength,
        });
    }
    ctx.restore();
}

function drawFlyChordFromContacts(canvas, contacts, lenA, lenB, sym1, sym2, colors, LIR_A, cLIR_A, LIR_B, cLIR_B) {
    const S = 2, size = 350;
    canvas.width = size * S; canvas.height = size * S;
    canvas.style.width = size + 'px'; canvas.style.height = size + 'px';
    const ctx = canvas.getContext('2d');
    ctx.clearRect(0, 0, canvas.width, canvas.height);
    ctx.save(); ctx.scale(S, S);

    const cx = size/2, cy = size/2;
    const _plddtOn = flyShowChordPlddt && ((flyPlddt.A && flyPlddt.A.some(x => x != null)) || (flyPlddt.B && flyPlddt.B.some(x => x != null)));   // pLDDT ring only when data present
    const _nDomRings = (_plddtOn ? 1 : 0) + (flyShowChordDom ? 1 : 0) + (flyShowChordTed ? 1 : 0), _domSpace = _nDomRings ? (_nDomRings * 7 + 3) : 0, _domHits = [];   // pLDDT + domain rings sit OUTSIDE the arc; reserve adaptive outer space
    const outerR = Math.max(40, Math.min(size/2 - (flyShowChordTicks ? 62 : 48) - _domSpace, size/2 - 4 - _domSpace - LiviaMaps.chordLabelExtent(ctx, [
        { name: flyChordLabel(sym1, 'A'), length: flyShowChordLength ? `${lenA}` : '' },
        { name: flyChordLabel(sym2, 'B'), length: flyShowChordLength ? `${lenB}` : '' },
    ], { ticks: flyShowChordTicks, radial: flyShowChordRadial }))), innerR = outerR - 12, chordR = innerR - 4;   // shrink to keep long radial labels on-canvas
    ctx.canvas._domHits = _domHits; ctx.canvas._domGeom = { cx, cy, size };   // hover
    const gap = 0.04, totalLen = lenA + lenB;
    const totalAngle = 2*Math.PI - 2*gap;
    const angleA = (lenA/totalLen)*totalAngle, angleB = (lenB/totalLen)*totalAngle;
    const startA = -Math.PI/2, endA = startA + angleA;
    const startB = endA + gap, endB = startB + angleB;
    const colA = colors.clir_a || '#00897B', colB = colors.clir_b || '#E64A19';
    const resAngleA = r => startA + ((r-1)/Math.max(1,lenA-1)) * angleA;
    const resAngleB = r => startB + ((r-1)/Math.max(1,lenB-1)) * angleB;

    // Draw contact lines
    const minDist = Math.min(...contacts.map(c => c.dist));
    const maxDist = Math.max(...contacts.map(c => c.dist));
    for (const c of contacts) {
        const aA = resAngleA(c.resA), aB = resAngleB(c.resB);
        const pA = {x: cx+Math.cos(aA)*chordR, y: cy+Math.sin(aA)*chordR};
        const pB = {x: cx+Math.cos(aB)*chordR, y: cy+Math.sin(aB)*chordR};
        const t = maxDist > minDist ? 1-(c.dist-minDist)/(maxDist-minDist) : 1;
        const alpha = 0.08 + t*0.2;
        LiviaMaps.strokeFadeQuad(ctx, pA.x, pA.y, cx, cy, pB.x, pB.y, hexToRgba(colA, alpha), hexToRgba(colB, alpha), 0.7);
    }

    // Arcs (gray base + LIR light + cLIR dark)
    const lightA = colors.lir_a || colA, lightB = colors.lir_b || colB;
    paintChordArcBand(ctx, cx, cy, outerR, innerR, startA, endA, lenA, LIR_A, cLIR_A, lightA, colA);
    paintChordArcBand(ctx, cx, cy, outerR, innerR, startB, endB, lenB, LIR_B, cLIR_B, lightB, colB);
    { let _rA = outerR + 3, _rB = outerR + 3;   // pLDDT ring first, then domain rings (UniProt/Pfam, TED) stacked outward, with hover hit-regions
      if (_plddtOn) { for (const _h of LiviaMaps.paintChordPlddtRing(ctx, cx, cy, _rA, _rA + 6, startA, endA, lenA, flyPlddt.A, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym1, 'A'), src: 'pLDDT' }, _h)); _rA += 7;
                     for (const _h of LiviaMaps.paintChordPlddtRing(ctx, cx, cy, _rB, _rB + 6, startB, endB, lenB, flyPlddt.B, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym2, 'B'), src: 'pLDDT' }, _h)); _rB += 7; }
      if (flyShowChordDom) { for (const _h of LiviaMaps.paintChordDomainRing(ctx, cx, cy, _rA, _rA + 6, startA, endA, lenA, flyDomains.A, CHORD_DOM_FILL, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym1, 'A'), src: 'UniProt/Pfam' }, _h)); _rA += 7;
                             for (const _h of LiviaMaps.paintChordDomainRing(ctx, cx, cy, _rB, _rB + 6, startB, endB, lenB, flyDomains.B, CHORD_DOM_FILL, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym2, 'B'), src: 'UniProt/Pfam' }, _h)); _rB += 7; }
      if (flyShowChordTed) { for (const _h of LiviaMaps.paintChordDomainRing(ctx, cx, cy, _rA, _rA + 6, startA, endA, lenA, flyTedDomains.A, CHORD_TED_FILL, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym1, 'A'), src: 'TED' }, _h));
                             for (const _h of LiviaMaps.paintChordDomainRing(ctx, cx, cy, _rB, _rB + 6, startB, endB, lenB, flyTedDomains.B, CHORD_TED_FILL, false)) _domHits.push(Object.assign({ chain: flyChordLabel(sym2, 'B'), src: 'TED' }, _h)); } }
    if (flyShowChordTicks) { LiviaMaps.drawChordRuler(ctx, cx, cy, outerR + _domSpace, startA, endA, lenA, {}); LiviaMaps.drawChordRuler(ctx, cx, cy, outerR + _domSpace, startB, endB, lenB, {}); }

    // Labels
    for (const arc of [{mid:(startA+endA)/2, name:flyChordLabel(sym1,'A'), size:lenA}, {mid:(startB+endB)/2, name:flyChordLabel(sym2,'B'), size:lenB}]) {
        LiviaMaps.drawChordArcLabel(ctx, {
            cx, cy, outerR: outerR + _domSpace, angle: arc.mid, ticks: flyShowChordTicks, radial: flyShowChordRadial, hidden: false,
            name: arc.name, showName: true, length: `${arc.size}`, showLength: flyShowChordLength,
        });
    }
    ctx.restore();
}

function buildChordDiagram(row, sym1, sym2, colors, pdbText) {
    const card = document.getElementById('chord-card');
    const canvas = document.getElementById('chord-canvas');
    const canvasAll = document.getElementById('chord-canvas-all');
    const chordAllWrap = document.getElementById('chord-all-wrap');
    if (!card || !canvas || !pdbText) { if (card) card.style.display = 'none'; return; }
    window.__redrawFlyChord = () => buildChordDiagram(row, sym1, sym2, colors, pdbText);   // chain-ID toggle re-render
    LiviaMaps.attachExportBar('chord-canvas', { name: 'chord_diagram', svg: () => buildChordDiagram(row, sym1, sym2, colors, pdbText), png: true });
    LiviaMaps.attachChordDomainHover(['chord-canvas', 'chord-canvas-all']);
    if (canvasAll) LiviaMaps.attachExportBar('chord-canvas-all', { name: 'chord_diagram_all_ranks', svg: () => buildChordDiagram(row, sym1, sym2, colors, pdbText), png: true });

    const cLIR_A = row._cLIR_A, cLIR_B = row._cLIR_B;
    const lenA = parseInt(row['Len A'] || row['protein_1_size']) || 0;
    const lenB = parseInt(row['Len B'] || row['protein_2_size']) || 0;
    if (lenA === 0 || lenB === 0 || cLIR_A.size === 0 || cLIR_B.size === 0) {
        card.style.display = 'none'; return;
    }

    const coordChains = parsePdbCbCoords(pdbText);
    const coordA = coordChains.get('A');
    const coordB = coordChains.get('B');
    if (!coordA || !coordB) { card.style.display = 'none'; return; }

    // Per-residue pLDDT per chain (B-factors) → chord pLDDT ring. Indexed plddt[r-1] = residue r.
    flyPlddt = { A: [], B: [] };
    try {
        const _bfs = parseBfactorsPerResidue(pdbText);
        for (const ch of ['A', 'B']) { const L = ch === 'A' ? lenA : lenB, a = new Array(L).fill(null); for (let r = 1; r <= L; r++) { const val = _bfs.get(`${ch}:${r}`); if (val != null) a[r - 1] = val; } flyPlddt[ch] = a; }
    } catch (e) { flyPlddt = { A: [], B: [] }; }
    { const _lg = document.getElementById('chord-plddt-legend'); if (_lg) { const _has = flyShowChordPlddt && ((flyPlddt.A && flyPlddt.A.some(v => v != null)) || (flyPlddt.B && flyPlddt.B.some(v => v != null))); _lg.innerHTML = LiviaMaps.plddtLegendHtml(); _lg.style.display = _has ? 'flex' : 'none'; } }

    card.style.display = '';

    // Draw selected rank chord
    drawFlyChord(canvas, cLIR_A, cLIR_B, lenA, lenB, sym1, sym2, colors, coordA, coordB, false, row._LIR_A, row._LIR_B);

    // Draw combined all-ranks chord — union of contacts from each rank's own PDB
    const hasMultipleRanks = flyParsed && flyParsed.rows && flyParsed.rows.length > 1;
    if (chordAllWrap) chordAllWrap.style.display = hasMultipleRanks ? '' : 'none';

    if (hasMultipleRanks && canvasAll && flyAverageRow) {
        // Collect contacts from each rank using its own PDB structure
        const contactDist = 8;
        const allContacts = new Map(); // "resA-resB" → {resA, resB, dist}
        for (const r of flyParsed.rows) {
            const rankNum = parseInt(r['Rank']);
            const rankPdb = flyPdbPerRank[rankNum];
            if (!rankPdb) continue;
            const rankCoords = parsePdbCbCoords(rankPdb);
            const rCoordA = rankCoords.get('A');
            const rCoordB = rankCoords.get('B');
            if (!rCoordA || !rCoordB) continue;
            for (const resA of r._cLIR_A) {
                const cA = rCoordA.get(resA);
                if (!cA) continue;
                for (const resB of r._cLIR_B) {
                    const key = `${resA}-${resB}`;
                    if (allContacts.has(key)) continue;
                    const cB = rCoordB.get(resB);
                    if (!cB) continue;
                    const dx = cA.x - cB.x, dy = cA.y - cB.y, dz = cA.z - cB.z;
                    const dist = Math.sqrt(dx*dx + dy*dy + dz*dz);
                    if (dist <= contactDist) allContacts.set(key, { resA, resB, dist });
                }
            }
        }
        // Draw using rank 1 coordinates for positioning, but with union of all contacts
        const unionContacts = [...allContacts.values()];
        if (unionContacts.length > 0) {
            // Use a modified draw that takes pre-computed contacts
            drawFlyChordFromContacts(canvasAll, unionContacts, lenA, lenB, sym1, sym2, colors, flyAverageRow._LIR_A, flyAverageRow._cLIR_A, flyAverageRow._LIR_B, flyAverageRow._cLIR_B);
        }
    }
}

// ── Generate PyMOL script ──
function generatePml(row, pdbFilename, colors) {
    const sym1 = row['Symbol 1'] || 'Protein1';
    const sym2 = row['Symbol 2'] || 'Protein2';
    const cLIR_A = row._cLIR_A, cLIR_B = row._cLIR_B;
    const LIR_A = row._LIR_A, LIR_B = row._LIR_B;
    const iLIS = row['iLIS'] || 'N/A';
    const ipTM = row['ipTM'] || 'N/A';

    const LIR_A_filtered = filterSmallSegments(LIR_A);
    const LIR_B_filtered = filterSmallSegments(LIR_B);
    const cLIR_A_filtered = new Set([...cLIR_A].filter(r => LIR_A_filtered.has(r)));
    const cLIR_B_filtered = new Set([...cLIR_B].filter(r => LIR_B_filtered.has(r)));
    const LIR_A_filled = fillGaps(LIR_A_filtered);
    const LIR_B_filled = fillGaps(LIR_B_filtered);

    const lirASpec = pmlResSpec(LIR_A_filled, 'A');
    const lirBSpec = pmlResSpec(LIR_B_filled, 'B');
    const clirASpec = pmlResSpec(cLIR_A_filtered, 'A');
    const clirBSpec = pmlResSpec(cLIR_B_filtered, 'B');

    let lines = [
        `# PyMOL: ${sym1} \u2014 ${sym2} (iLIS=${iLIS}, ipTM=${ipTM})`,
        `# Generated by LIVIA (Local Interaction Visualization and Analysis), build ${window.LIVIA_BUILD || "?"}, ${new Date().toLocaleDateString("en-CA")}`,
        ``,
        `reinitialize`,
        `load ${pdbFilename}`,
        ``,
        `# Define colors`,
        `set_color lir_a, ${hexToRgbFloat(colors.lir_a)}`,
        `set_color lir_b, ${hexToRgbFloat(colors.lir_b)}`,
        `set_color clir_a, ${hexToRgbFloat(colors.clir_a)}`,
        `set_color clir_b, ${hexToRgbFloat(colors.clir_b)}`,
        ``,
        ``,
        `bg_color white`,
        `set cartoon_fancy_helices, 1`,
        `set cartoon_fancy_sheets, 1`,
        `set antialias, 2`,
        `set ray_shadows, 1`,
        `set ray_trace_mode, 1`,
        `hide everything, all`,
    ];

    // ── Show cartoons (showComplete controls scope) ──
    if (showComplete) {
        lines.push('show cartoon, all');
    } else {
        if (lirASpec) lines.push(`show cartoon, ${lirASpec}`);
        if (lirBSpec) lines.push(`show cartoon, ${lirBSpec}`);
    }

    // ── Color mode (independent of showComplete) ──
    if (cxcColorMode === 'plddt') {
        lines.push('set_color af_vhigh, [0.051, 0.341, 0.827]');
        lines.push('set_color af_high, [0.416, 0.796, 0.945]');
        lines.push('set_color af_low, [0.996, 0.851, 0.212]');
        lines.push('set_color af_vlow, [0.992, 0.490, 0.302]');
        lines.push('color af_vlow, all and not b > 50');
        lines.push('color af_low, all and b > 50 and not b > 70');
        lines.push('color af_high, all and b > 70 and not b > 90');
        lines.push('color af_vhigh, all and b > 90');
    } else if (cxcColorMode === 'bychain' || cxcColorMode === 'bypolymer') {
        lines.push('util.cbc');
    } else {
        if (grayNonLir && showComplete) {
            lines.push(`set_color nonlir_gray, [0.69, 0.69, 0.69]`);
            lines.push(`color nonlir_gray, all`);
        }
        if (lirASpec) lines.push(`color lir_a, ${lirASpec}`);
        if (lirBSpec) lines.push(`color lir_b, ${lirBSpec}`);
        if (clirASpec) lines.push(`color clir_a, ${clirASpec}`);
        if (clirBSpec) lines.push(`color clir_b, ${clirBSpec}`);
    }

    lines.push('');
    lines.push('orient');
    lines.push('');
    lines.push('# Save commands (uncomment to use):');
    lines.push(`# bg_color white`);
    lines.push(`# set ray_opaque_background, 0`);
    lines.push(`# ray 2400, 2400`);
    lines.push(`# png ${sym1}_${sym2}.png, dpi=300`);
    lines.push(`# save ${sym1}_${sym2}_session.pse`);

    return lines.join('\n');
}

// ── Gene Pair search: the parts both predictome pages share ──
let pairAcDebounce = null;
let pairAcIdx = -1;

function hidePairAc() {
    document.getElementById('pair-ac-a').style.display = 'none';
    document.getElementById('pair-ac-b').style.display = 'none';
    pairAcIdx = -1;
}

function handlePairAcKey(e, inputId, listId) {
    const list = document.getElementById(listId);
    const items = list.querySelectorAll('div');
    if (list.style.display !== 'none' && items.length > 0) {
        if (e.key === 'ArrowDown') {
            e.preventDefault();
            pairAcIdx = Math.min(pairAcIdx + 1, items.length - 1);
            items.forEach((it, i) => { it.style.background = i === pairAcIdx ? '#eaf2f8' : ''; });
            items[pairAcIdx].scrollIntoView({ block: 'nearest' });
            return;
        }
        if (e.key === 'ArrowUp') {
            e.preventDefault();
            pairAcIdx = Math.max(pairAcIdx - 1, 0);
            items.forEach((it, i) => { it.style.background = i === pairAcIdx ? '#eaf2f8' : ''; });
            items[pairAcIdx].scrollIntoView({ block: 'nearest' });
            return;
        }
        if (e.key === 'Enter') {
            e.preventDefault();
            if (pairAcIdx >= 0 && pairAcIdx < items.length) {
                items[pairAcIdx].click();
            } else {
                list.style.display = 'none';
            }
            return;
        }
        if (e.key === 'Escape') { list.style.display = 'none'; pairAcIdx = -1; return; }
    } else if (e.key === 'Enter') {
        searchGenePair();
    }
}

function setPairAndSearch(a, b) {
    const inputA = document.getElementById('pair-gene-a');
    const inputB = document.getElementById('pair-gene-b');
    inputA.value = a; inputA.dataset.fbgn = '';
    inputB.value = b; inputB.dataset.fbgn = '';
    searchGenePair();
}

async function offerPairMatches(matches, labelA, labelB, statusEl) {
    // FlyPredictome's own pair index carries stale entries — a listed SET can return "Not found"
    // on their server even though the summary page still references it. Check each candidate's
    // detail page actually exists before offering it, so a click never dead-ends on a broken link.
    statusEl.textContent = labelA + ' vs ' + labelB + ': checking ' + matches.length +
        (matches.length === 1 ? ' prediction set' : ' prediction sets') + ' against FlyPredictome...';
    const liveFlags = await Promise.all(matches.map(p => p.detailUrl ? fetchViaProxy(p.detailUrl).then(Boolean) : false));
    const liveMatches = matches.filter((p, i) => liveFlags[i]);
    const deadCount = matches.length - liveMatches.length;
    const deadNote = deadCount ? (' (' + deadCount + ' other listed set' + (deadCount === 1 ? '' : 's') + ' skipped — broken on FlyPredictome)') : '';

    if (liveMatches.length === 0) {
        statusEl.textContent = labelA + ' vs ' + labelB + ': FlyPredictome lists ' + matches.length +
            (matches.length === 1 ? ' prediction set' : ' prediction sets') +
            ' for this pair, but every linked page currently returns nothing on their end (stale index entries) — nothing to analyze right now.';
        statusEl.style.color = '#c0392b';
        return;
    }

    // Single match → load directly
    if (liveMatches.length === 1 && liveMatches[0].detailUrl) {
        statusEl.textContent = labelA + ' vs ' + labelB + ': loading analysis...' + deadNote;
        statusEl.style.color = '#27ae60';
        document.getElementById('url-input').value = liveMatches[0].detailUrl;
        switchTab('url');
        generate();
        return;
    }

    // Multiple matches → show selection table
    statusEl.textContent = labelA + ' vs ' + labelB + ': ' + liveMatches.length + ' prediction sets found' + deadNote + '. Click a row to analyze.';
    statusEl.style.color = '#27ae60';
    displayPairResults(liveMatches);
}

function displayPairResults(matches) {
    const container = document.getElementById('pair-results');
    container.innerHTML = '';
    container.style.display = '';

    const tableWrap = document.createElement('div');
    tableWrap.style.cssText = 'overflow-x:auto; max-height:300px; overflow-y:auto; border:1px solid #e0e0e0; border-radius:8px;';

    const table = document.createElement('table');
    table.className = 'data-table';
    table.style.fontSize = '0.8rem';
    table.innerHTML =
        '<thead style="position:sticky; top:0; z-index:1;">' +
            '<tr>' +
                '<th>Protein A</th>' +
                '<th>Protein B</th>' +
                '<th>iLIS (best)</th>' +
                '<th>ipTM (best)</th>' +
                '<th>Len A</th>' +
                '<th>Len B</th>' +
                '<th>MIST</th>' +
            '</tr>' +
        '</thead>' +
        '<tbody></tbody>';

    const tbody = table.querySelector('tbody');
    for (const p of matches) {
        const tr = document.createElement('tr');
        tr.className = 'clickable';
        tr.style.cursor = 'pointer';
        const mistStr = [p.mistPPI ? 'PPI' : '', p.mistGI ? 'GI' : ''].filter(Boolean).join(', ') || '\u2014';
        const partial1 = p.isPartial1 ? ' <span style="color:#c0392b; font-size:0.7rem; font-weight:600;">(Partial)</span>' : '';
        const partial2 = p.isPartial2 ? ' <span style="color:#c0392b; font-size:0.7rem; font-weight:600;">(Partial)</span>' : '';
        tr.innerHTML =
            '<td><strong>' + esc(p.sym1) + '</strong> <span style="color:#707070; font-size:0.72rem;">' + esc(p.fbgn1 || p.prot1) + '</span>' + partial1 + '</td>' +
            '<td><strong>' + esc(p.sym2) + '</strong> <span style="color:#707070; font-size:0.72rem;">' + esc(p.fbgn2 || p.prot2) + '</span>' + partial2 + '</td>' +
            '<td>' + ilisColor(p.ilisBest) + '</td>' +
            '<td>' + p.iptmBest.toFixed(2) + '</td>' +
            '<td>' + p.lenA + '</td>' +
            '<td>' + p.lenB + '</td>' +
            '<td>' + mistStr + '</td>';
        if (p.detailUrl) {
            tr.addEventListener('click', ((url) => (e) => {
                if (e.ctrlKey || e.metaKey) {
                    window.open(window.location.pathname + '?url=' + encodeURIComponent(url), '_blank');
                } else {
                    document.getElementById('url-input').value = url;
                    switchTab('url');
                    generate();
                }
            })(p.detailUrl));
        }
        tbody.appendChild(tr);
    }

    tableWrap.appendChild(table);
    container.appendChild(tableWrap);
}
