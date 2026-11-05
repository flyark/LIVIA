// Small helpers shared by the prediction pages (universal, dimer, monomer, flypredictome, ortholog_predictome). They take what they need as arguments, or call a page function by name at call time.

function copyScript() {
    const text = document.getElementById('script-preview').textContent;
    navigator.clipboard.writeText(text).then(() => {
        const btn = document.getElementById('btn-copy');
        btn.textContent = 'Copied!';
        setTimeout(() => { btn.textContent = 'Copy to Clipboard'; }, 1500);
    });
}

// A download link's file, made in memory: the link's previous file is freed first, so a session of clicks does not keep a copy
// of the structure per click.
function setBlobHref(el, text) {
    if (el._blobUrl) { try { URL.revokeObjectURL(el._blobUrl); } catch (e) {} }
    el._blobUrl = URL.createObjectURL(new Blob([text], { type: 'text/plain' }));
    el.href = el._blobUrl;
}

// Weighted mean of two per-chain numeric values (e.g. cLIpLDDT_i/j weighted by interface-residue
// count), matching lis.py's _wmean_plddt. Returns NaN if either value or the total weight is missing.
function wmeanNum(vi, vj, wi, wj) {
    if (vi == null || vj == null || isNaN(vi) || isNaN(vj)) return NaN;
    const w = (wi || 0) + (wj || 0);
    if (w <= 0) return NaN;
    return (wi * vi + wj * vj) / w;
}

// pDockQ2/LIpDockQ2 are inherently per-chain-direction (Zhu et al. 2023's own definition), unlike
// every other metric on these pages, which is naturally symmetric. Mean-of-two-directions collapses
// them to one displayed number -- verified against real DockQ-labeled benchmark data to perform
// equivalently to max/min (ColabFold 1.6.3's own console log uses max, for reference), and chosen
// over max/min because it matches how every other metric here already symmetrizes (ipTM, ipSAE,
// actifpTM are all (fwd+rev)/2). Returns null only when both inputs are null (e.g. lightweight
// bundle with no numeric PAE at all); a single missing side falls back to the other.
function mean2Nullable(a, b) {
    if (a == null && b == null) return null;
    if (a == null) return b;
    if (b == null) return a;
    return (a + b) / 2;
}

function roundRectPath(ctx, x, y, w, h, r) {
    r = Math.min(r, w / 2, h / 2);
    ctx.beginPath();
    ctx.moveTo(x + r, y); ctx.arcTo(x + w, y, x + w, y + h, r); ctx.arcTo(x + w, y + h, x, y + h, r);
    ctx.arcTo(x, y + h, x, y, r); ctx.arcTo(x, y, x + w, y, r); ctx.closePath();
}

// ── Confident-region structure trimming for ChimeraX / PyMOL post-processing ──
// keepMap[chain] = null (whole chain) | [minResnum, maxResnum]; chains absent are dropped.
// Original residue numbering is preserved so the generated script's coloring still works.
function trimCif(cifText, keepMap){
  const lines=cifText.split('\n'), out=[]; let i=0;
  while(i<lines.length){
    const line=lines[i];
    if(line.trim()==='loop_'){
      let j=i+1; const headers=[];
      while(j<lines.length && lines[j].trim().startsWith('_')){ headers.push(lines[j].trim()); j++; }
      if(headers.length && headers[0].startsWith('_atom_site.')){
        const col=n=>headers.indexOf('_atom_site.'+n);
        let cc=col('auth_asym_id'); if(cc<0)cc=col('label_asym_id');
        let rc=col('auth_seq_id'); if(rc<0)rc=col('label_seq_id');
        out.push(line); headers.forEach(h=>out.push(h));
        let k=j;
        while(k<lines.length){
          const t=lines[k].trim();
          if(t===''||t==='#'||t==='loop_'||t.startsWith('_')||t.startsWith('data_'))break;
          const toks=t.split(/\s+/); const ch=toks[cc], rn=parseInt(toks[rc],10); const kr=keepMap[ch];
          if(keepMap.hasOwnProperty(ch)&&(kr===null||(kr instanceof Set?kr.has(rn):(rn>=kr[0]&&rn<=kr[1])))) out.push(lines[k]);
          k++;
        }
        i=k; continue;
      }
    }
    out.push(line); i++;
  }
  return out.join('\n');
}

function applyChainToggleStyle(btn, visible) {
    const baseColor = getChainToggleColor(btn.dataset.chain);
    btn.style.cssText = visible
        ? `padding:0.2rem 0.6rem; border:1px solid ${LiviaMaps.chipFill(baseColor).fill}; background:${LiviaMaps.chipFill(baseColor).fill}; color:${LiviaMaps.chipFill(baseColor).ink}; border-radius:4px; cursor:pointer; font-size:0.78rem; font-weight:600; min-width:1.6rem;`
        : `padding:0.2rem 0.6rem; border:1px solid #ccc; background:#fff; color:#6B6B6B; border-radius:4px; cursor:pointer; font-size:0.78rem; min-width:1.6rem; text-decoration:line-through;`;
}

// Fill and text color for a filled name chip: white text on a fill dark enough for it (at least 4.5:1),
// the fill darkened just enough when it is not (js/livia-maps.js chipFill).
function contrastText(hex) { return LiviaMaps.chipFill(hex).ink; }

// ============================================================================
// TED DOMAIN ANNOTATIONS (from AlphaFold DB)
// ============================================================================
// The page keeps the results in its own tedCache.
async function fetchTedDomains(uniprotId) {
    if (!uniprotId) return [];
    if (tedCache.has(uniprotId)) return tedCache.get(uniprotId);
    let domains = [];
    try {
        const resp = await fetch(`https://alphafold.ebi.ac.uk/api/domains/${uniprotId}`);
        if (resp.ok) {
            const data = await resp.json();
            for (const ann of (data.annotations || [])) {
                const segments = (ann.segments || []).map(s => ({ start: s.af_start, end: s.af_end }));
                if (segments.length === 0) continue;
                const start = Math.min(...segments.map(s => s.start));
                const end = Math.max(...segments.map(s => s.end));
                const cath = ann.cath_label || '';
                const cathLevel = ann.cath_assignment_level || '';
                const name = cath ? `TED ${ann.ted_domain_no} (${cath})` : `TED ${ann.ted_domain_no}`;
                domains.push({
                    name, start, end, segments,
                    cath, cathLevel,
                    plddt: ann.plddt || 0,
                    qscore: ann.qscore || 0,
                    nres: ann.nres_domain || 0,
                    tedNo: ann.ted_domain_no
                });
            }
            domains.sort((a, b) => a.start - b.start);
        }
    } catch (e) {
        console.warn(`Failed to fetch TED domains for ${uniprotId}:`, e);
    }
    tedCache.set(uniprotId, domains);
    return domains;
}

// ── Confident-region structure trimming (PDB) for ChimeraX / PyMOL post-processing ──
// keepMap[chain] = null (whole chain) | [minResnum,maxResnum]; chains absent are dropped. Residue numbers preserved.
function trimPdb(pdbText, keepMap){
  return pdbText.split('\n').filter(l=>{
    const rec=l.substr(0,6);
    if(!(rec.startsWith('ATOM')||rec.startsWith('HETATM')||rec.startsWith('ANISOU')||rec.startsWith('TER'))) return true;
    const ch=l[21], rn=parseInt(l.substr(22,4),10), kr=keepMap[ch];
    if(rec.startsWith('TER')) return keepMap.hasOwnProperty(ch);
    return keepMap.hasOwnProperty(ch) && (kr===null || (kr instanceof Set ? kr.has(rn) : (rn>=kr[0]&&rn<=kr[1])));
  }).join('\n');
}

function makeSortable(table) {
    if (!table) return;
    const headers = table.querySelectorAll('thead th');
    headers.forEach((th, colIdx) => {
        th.classList.add('sortable');
        th.addEventListener('click', () => {
            const tbody = table.querySelector('tbody');
            if (!tbody) return;
            const rows = Array.from(tbody.querySelectorAll('tr'));
            // Separate pinned rows (Average, selected) from sortable rows
            const sortableRows = rows.filter(r => !r.dataset.pinned && r.textContent.trim() !== '');
            const pinnedRows = rows.filter(r => r.dataset.pinned);

            const isAsc = th.classList.contains('sort-asc');
            headers.forEach(h => { h.classList.remove('sort-asc', 'sort-desc'); });
            th.classList.add(isAsc ? 'sort-desc' : 'sort-asc');

            sortableRows.sort((a, b) => {
                const aVal = a.cells[colIdx]?.textContent.trim() || '';
                const bVal = b.cells[colIdx]?.textContent.trim() || '';
                const aNum = parseFloat(aVal);
                const bNum = parseFloat(bVal);
                if (!isNaN(aNum) && !isNaN(bNum)) {
                    return isAsc ? bNum - aNum : aNum - bNum;
                }
                return isAsc ? bVal.localeCompare(aVal) : aVal.localeCompare(bVal);
            });

            // Re-append: sortable rows first, then pinned rows at bottom
            sortableRows.forEach(r => tbody.appendChild(r));
            pinnedRows.forEach(r => tbody.appendChild(r));
        });
    });
}
