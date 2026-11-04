// Helpers shared by the FlyPredictome and Ortholog Interactome pages, byte-identical on both before they moved here.

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
    if (!pdbText || !r._cLIR_A || !r._cLIR_B || r._cLIR_A.size === 0 || r._cLIR_B.size === 0) return null;
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
