/*
 * livia-xlsx.js — read an .xlsx workbook in the browser with no extra dependency.
 *
 * An .xlsx file is a zip of XML parts, and every LIVIA page that accepts a .zip already
 * loads fflate, so the whole reader is fflate + DOMParser. Cells are either numbers or
 * indices into the shared string table (t="s"); inline strings are handled too.
 *
 * LIVIA-specific helpers:
 *   pickLisSheet(wb)   → the sheet that carries lis.py columns, preferring one that also
 *                        knows the protein lengths (they let cLIP place cLIR on the structure).
 *   sequenceSheet(wb)  → a sheet with Seq_1/Seq_2, which we can hand to cLIP as FASTA.
 *
 * Requires: fflate (window.fflate).
 */
(function (global) {
    'use strict';

    const REL_NS = 'http://schemas.openxmlformats.org/officeDocument/2006/relationships';

    // Collaborator sheets carry a Score column (one point per metric passing its 10% FPR
    // cutoff). It is an upstream QC number, not something LIVIA reasons about — drop it.
    const DROP_COLS = ['Score'];

    function colIndex(ref) {                         // "BA12" → 52
        let n = 0;
        for (const ch of String(ref).replace(/\d+/g, '')) n = n * 26 + (ch.charCodeAt(0) - 64);
        return n - 1;
    }

    function readSheet(doc, shared) {
        const grid = [];
        for (const row of doc.getElementsByTagName('row')) {
            const cells = [];
            for (const c of row.getElementsByTagName('c')) {
                const i = colIndex(c.getAttribute('r') || 'A1');
                const t = c.getAttribute('t');
                let v = '';
                if (t === 'inlineStr') {
                    const is = c.getElementsByTagName('t')[0];
                    v = is ? is.textContent : '';
                } else {
                    const ve = c.getElementsByTagName('v')[0];
                    v = ve ? ve.textContent : '';
                    if (t === 's') v = shared[+v] || '';
                }
                cells[i] = v;
            }
            grid.push(cells);
        }
        if (!grid.length) return { header: [], rows: [] };

        const header = (grid[0] || []).map((x) => String(x == null ? '' : x).trim());
        const rows = [];
        for (let r = 1; r < grid.length; r++) {
            const cur = grid[r];
            if (!cur) continue;
            const o = {};
            let any = false;
            for (let i = 0; i < header.length; i++) {
                const h = header[i];
                if (!h || DROP_COLS.indexOf(h) >= 0) continue;
                const v = cur[i];
                o[h] = v == null ? '' : v;
                if (v != null && v !== '') any = true;
            }
            if (any) rows.push(o);
        }
        return { header: header.filter((h) => h && DROP_COLS.indexOf(h) < 0), rows };
    }

    // → { sheetName: {header, rows} } in workbook order
    function parseXlsx(bytes) {
        if (typeof global.fflate === 'undefined') throw new Error('xlsx support needs fflate (no internet connection?)');
        const files = global.fflate.unzipSync(bytes instanceof Uint8Array ? bytes : new Uint8Array(bytes));
        const dec = new TextDecoder();
        const dom = (s) => new DOMParser().parseFromString(s, 'application/xml');

        let shared = [];
        if (files['xl/sharedStrings.xml']) {
            const d = dom(dec.decode(files['xl/sharedStrings.xml']));
            shared = [...d.getElementsByTagName('si')].map((si) =>
                [...si.getElementsByTagName('t')].map((t) => t.textContent).join(''));
        }

        const wb = dom(dec.decode(files['xl/workbook.xml']));
        const relMap = {};
        if (files['xl/_rels/workbook.xml.rels']) {
            const rd = dom(dec.decode(files['xl/_rels/workbook.xml.rels']));
            for (const r of rd.getElementsByTagName('Relationship')) relMap[r.getAttribute('Id')] = r.getAttribute('Target');
        }

        const out = {};
        for (const s of wb.getElementsByTagName('sheet')) {
            const name = s.getAttribute('name');
            const rid = s.getAttributeNS(REL_NS, 'id') || s.getAttribute('r:id');
            let target = relMap[rid];
            if (!target) continue;
            // Targets come in three flavours: relative ("worksheets/sheet1.xml"), package-absolute
            // ("/xl/worksheets/sheet1.xml") and already-prefixed. Normalize before the lookup —
            // prefixing an absolute path yields "xl/xl/..." and silently loses every sheet.
            target = target.replace(/^\//, '').replace(/^\.\//, '');
            if (!target.startsWith('xl/')) target = 'xl/' + target;
            const part = files[target];
            if (!part) continue;
            out[name] = readSheet(dom(dec.decode(part)), shared);
        }
        return out;
    }

    const has = (h, c) => h.indexOf(c) >= 0;

    // Detect by columns, not by sheet name — collaborators rename sheets.
    function sheetKind(header) {
        const clip = has(header, 'Symbol_1') && has(header, 'Symbol_2') &&
                     has(header, 'cLIR_indice_A') && has(header, 'cLIR_indice_B') && has(header, 'iLIS');
        if (clip) return has(header, 'Protein_Len_A') ? 'lis-full' : 'lis';
        if (has(header, 'Seq_1') && has(header, 'Seq_2') && has(header, 'Symbol_1')) return 'sequence';
        if (has(header, 'Symbol_1') && has(header, 'Symbol_2') && has(header, 'iLIS')) return 'pairs';   // network-only
        return null;
    }

    // Every sheet LIVIA could use, richest first.
    function candidateSheets(wb) {
        const rank = { 'lis-full': 0, 'lis': 1, 'pairs': 2, 'sequence': 3 };
        return Object.keys(wb)
            .map((name) => ({ name, kind: sheetKind(wb[name].header), rows: wb[name].rows.length }))
            .filter((s) => s.kind)
            .sort((a, b) => rank[a.kind] - rank[b.kind] || b.rows - a.rows);
    }

    // The sheet cLIP should cluster: needs cLIR indices, prefers known protein lengths.
    function pickLisSheet(wb) {
        const c = candidateSheets(wb).find((s) => s.kind === 'lis-full' || s.kind === 'lis');
        return c ? Object.assign({}, c, { rows: wb[c.name].rows, header: wb[c.name].header }) : null;
    }

    // The sheet the PPI network builds from. The network applies its own iLIS cutoff, so the
    // widest sheet wins — a positives-only sheet would silently hide every edge below it.
    function pickPairSheet(wb) {
        const c = candidateSheets(wb).filter((s) => s.kind !== 'sequence').sort((a, b) => b.rows - a.rows)[0];
        return c ? Object.assign({}, c, { rows: wb[c.name].rows, header: wb[c.name].header }) : null;
    }

    // Sheet → lis.py-style CSV. A workbook has no prediction name or chain columns, and the
    // network keys its pair records on them, so synthesize both. `name` leads the row because
    // the cLIP hand-off matches predictions on the first field.
    //
    // Residue-index cells arrive as pretty-printed numpy arrays that wrap across lines. A
    // quoted newline is legal CSV, but the hand-off splits the text line-by-line, so a wrapped
    // cell would tear one prediction into several broken rows. Indices are whitespace-separated
    // either way — flatten the cell to a single line.
    function sheetToCsv(sheet, lengths) {
        const esc = (v) => {
            const s = String(v == null ? '' : v).replace(/\s*\r?\n\s*/g, ' ');
            return /[",]/.test(s) ? '"' + s.replace(/"/g, '""') + '"' : s;
        };
        const cols = sheet.header.slice();
        const addLen = lengths && lengths.size && !has(cols, 'Protein_Len_A');
        const head = 'name,chain_i,chain_j,' + cols.map(esc).join(',') + (addLen ? ',Protein_Len_A,Protein_Len_B' : '');
        const out = [head];
        for (const r of sheet.rows) {
            const row = [esc(r.Symbol_1 + '_vs_' + r.Symbol_2), 'A', 'B'].concat(cols.map((c) => esc(r[c])));
            if (addLen) row.push(lengths.get(String(r.Symbol_1 || '').trim()) || '', lengths.get(String(r.Symbol_2 || '').trim()) || '');
            out.push(row.join(','));
        }
        return out.join('\n');
    }

    // symbol → protein length, gathered from wherever the workbook happens to state it. The
    // widest sheet often omits Protein_Len, and resolveStructure scores candidates on length —
    // guessing it from the largest contact index can land on the wrong protein.
    function lengthMap(wb) {
        const len = new Map();
        for (const name of Object.keys(wb)) {
            const { header, rows } = wb[name];
            const pairs = [];
            if (has(header, 'Symbol_1') && has(header, 'Protein_Len_A')) pairs.push(['Symbol_1', 'Protein_Len_A']);
            if (has(header, 'Symbol_2') && has(header, 'Protein_Len_B')) pairs.push(['Symbol_2', 'Protein_Len_B']);
            for (const r of rows) for (const [s, l] of pairs) {
                const sym = String(r[s] || '').trim(), n = parseInt(r[l], 10);
                if (sym && n > 0 && !len.has(sym)) len.set(sym, n);
            }
            if (sheetKind(header) === 'sequence') {
                for (const r of rows) for (const [s, q] of [['Symbol_1', 'Seq_1'], ['Symbol_2', 'Seq_2']]) {
                    const sym = String(r[s] || '').trim(), seq = String(r[q] || '').trim();
                    if (sym && seq && !len.has(sym)) len.set(sym, seq.length);
                }
            }
        }
        return len;
    }

    function sequenceSheet(wb) {
        const c = candidateSheets(wb).find((s) => s.kind === 'sequence');
        return c ? Object.assign({}, c, { rows: wb[c.name].rows }) : null;
    }

    // Sequence sheet → FASTA text, one record per distinct symbol.
    function sequenceSheetToFasta(sheet) {
        const seen = new Map();
        for (const r of (sheet.rows || [])) {
            for (const [sym, seq] of [[r.Symbol_1, r.Seq_1], [r.Symbol_2, r.Seq_2]]) {
                const s = String(sym || '').trim(), q = String(seq || '').trim();
                if (s && q && !seen.has(s)) seen.set(s, q);
            }
        }
        let out = '';
        for (const [sym, seq] of seen) out += '>' + sym + '\n' + seq + '\n';
        return { fasta: out, count: seen.size };
    }

    // UniProt entry-name suffixes and accession prefixes both name an organism, and the
    // Protein_*/Uniprot_* columns carry one even when Symbol_* has been reduced to bare
    // gene names. Without this the resolver falls back on symbol+length and can land on
    // the wrong species' ortholog.
    const ORG_SUFFIX = { DROME: '7227', CAEEL: '6239', HUMAN: '9606', MOUSE: '10090', RAT: '10116',
                         YEAST: '559292', DANRE: '7955', ARATH: '3702', XENLA: '8355', CHICK: '9031',
                         BOVIN: '9913', PIG: '9823', SCHPO: '284812', ECOLI: '83333' };
    const ORG_PREFIX = [[/^FB(gn|pp|tr)\d/i, '7227'], [/^WBGene\d/i, '6239'], [/^ENSMUS[GPT]\d/i, '10090'],
                        [/^ENSDAR[GPT]\d/i, '7955'], [/^ENSRNO[GPT]\d/i, '10116'], [/^ENS[GPT]\d/i, '9606'],
                        [/^AT\dG\d{5}/i, '3702'], [/^Y[A-P][LR]\d{3}[WC]/, '559292']];

    function speciesHint(rows) {
        const votes = {};
        const vote = (tax) => { if (tax) votes[tax] = (votes[tax] || 0) + 1; };
        for (const r of (rows || []).slice(0, 400)) {
            for (const col of ['Protein_1', 'Protein_2', 'Uniprot_1', 'Uniprot_2']) {
                const v = String(r[col] || '').trim();
                if (!v) continue;
                const suf = v.match(/_([A-Z][A-Z0-9]{2,5})$/);
                if (suf && ORG_SUFFIX[suf[1]]) { vote(ORG_SUFFIX[suf[1]]); continue; }
                for (const [re, tax] of ORG_PREFIX) if (re.test(v)) { vote(tax); break; }
            }
        }
        const best = Object.keys(votes).sort((a, b) => votes[b] - votes[a])[0];
        return best || '';
    }

    global.LiviaXlsx = { parseXlsx, candidateSheets, pickLisSheet, pickPairSheet, sequenceSheet, sequenceSheetToFasta, speciesHint, sheetToCsv, lengthMap };
})(window);
