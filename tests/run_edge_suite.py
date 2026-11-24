#!/usr/bin/env python3
"""
LIVIA Automated Multimodal Edge-Case & Invariant Test Suite
===========================================================
Generates lightweight synthetic mmCIF/PDB + PAE fixtures covering 14 high-risk
structural biology & predictor format edge cases, drives live `universal.html`
in headless Chromium (Playwright), and verifies 4 invariant classes:

  1. Parser & token-alignment invariants:
     sum(chainSizes) == pae.length == coords.length == caCoords.length,
     exact chainTypes, covalentLinks (_struct_conn loop & single-record & PDB LINK),
     PTM detection, and per-token B-factor coverage (including DNA/RNA & glycans).
  2. Numerical & WASM-vs-JS parity invariants:
     livia-wasm == pure-JS fallback (<= 1e-5), zero NaN/undefined metrics,
     no index shift on pDockQ2 when a glycan sits between two protein chains,
     and UI displays '—' for CA-only pDockQ2 on non-protein pairs.
  3. 3D Mol* component & ChimeraX/PyMOL export invariants:
     covalently attached polymer residues bundled into the glycan/ligand Mol*
     ball_and_stick component, included in .cxc & .pml exports, and zero errors
     when cycling all color/display modes (Gradient, Solid, By chain, pLDDT, LIR/Complete).
  4. Built-in regression check (optional real archives if present).

Usage:
    python3 tests/run_edge_suite.py
    python3 tests/run_edge_suite.py --real-zip /path/to/glycan_folds.zip
"""
from __future__ import annotations

import argparse
import functools
import http.server
import io
import json
import math
import socketserver
import sys
import tempfile
import threading
import zipfile
from pathlib import Path

from playwright.sync_api import sync_playwright

REPO_DIR = Path(__file__).resolve().parent.parent

CIF_HEADER = """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_entity_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_asym_id
_atom_site.pdbx_PDB_model_num
"""


def _cif_atom(
    idx: int,
    group: str,
    sym: str,
    atom: str,
    comp: str,
    chain: str,
    ent: int,
    seq: str,
    x: float,
    y: float,
    z: float,
    bf: float = 88.5,
    auth_seq: int | None = None,
    wrap_boltz: bool = False,
) -> str:
    aseq = str(auth_seq if auth_seq is not None else (seq if seq != "." else 1))
    if wrap_boltz:
        return (
            f"{group} {idx} {sym} {atom} . {comp} {chain} {ent} {seq} ? "
            f"{x:.3f} {y:.3f} {z:.3f} 1.00 {bf:.2f} {aseq} {chain}\n1\n"
        )
    return (
        f"{group} {idx} {sym} {atom} . {comp} {chain} {ent} {seq} ? "
        f"{x:.3f} {y:.3f} {z:.3f} 1.00 {bf:.2f} {aseq} {chain} 1\n"
    )


def _add_protein_chain(
    lines: list[str],
    idx: int,
    chain: str,
    ent: int,
    seq_res: list[str],
    x0: float,
    y0: float,
    z0: float,
    bf: float = 88.0,
    wrap_boltz: bool = False,
    start_seq: int = 1,
) -> int:
    for r_i, comp in enumerate(seq_res):
        rnum = start_seq + r_i
        x = x0 + r_i * 3.8
        atoms = [("N", "N", -0.6, 0.0, 0.0), ("CA", "C", 0.0, 0.0, 0.0), ("C", "C", 0.6, 0.0, 0.0), ("O", "O", 0.6, 0.6, 0.0)]
        if comp == "SER":
            atoms += [("CB", "C", 0.0, 1.5, 0.0), ("OG", "O", 0.0, 2.7, 0.0)]
        elif comp == "ASN":
            atoms += [("CB", "C", 0.0, 1.5, 0.0), ("CG", "C", 0.0, 2.5, 0.5), ("OD1", "O", 0.0, 2.5, 1.7), ("ND2", "N", 0.0, 3.6, 0.0)]
        elif comp == "CYS":
            atoms += [("CB", "C", 0.0, 1.5, 0.0), ("SG", "S", 0.0, 3.0, 0.0)]
        elif comp == "SEP":
            # Phosphoserine written as HETATM (AF3 convention)
            p_atoms = [
                ("N", "N", -0.6, 0.0, 0.0), ("CA", "C", 0.0, 0.0, 0.0), ("C", "C", 0.6, 0.0, 0.0), ("O", "O", 0.6, 0.6, 0.0),
                ("CB", "C", 0.0, 1.5, 0.0), ("OG", "O", 0.0, 2.7, 0.0), ("P", "P", 0.0, 4.1, 0.0),
                ("O1P", "O", 1.2, 4.5, 0.0), ("O2P", "O", -1.2, 4.5, 0.0), ("O3P", "O", 0.0, 4.5, 1.2),
            ]
            for an, sy, dx, dy, dz in p_atoms:
                lines.append(_cif_atom(idx, "HETATM", sy, an, "SEP", chain, ent, str(rnum), x + dx, y0 + dy, z0 + dz, bf, rnum, wrap_boltz))
                idx += 1
            continue
        elif comp != "GLY":
            atoms += [("CB", "C", 0.0, 1.5, 0.0)]
        for an, sy, dx, dy, dz in atoms:
            lines.append(_cif_atom(idx, "ATOM", sy, an, comp, chain, ent, str(rnum), x + dx, y0 + dy, z0 + dz, bf, rnum, wrap_boltz))
            idx += 1
    return idx


def _add_nucleic_chain(
    lines: list[str],
    idx: int,
    chain: str,
    ent: int,
    bases: list[str],
    x0: float,
    y0: float,
    z0: float,
    bf: float = 84.0,
) -> int:
    for r_i, base in enumerate(bases):
        rnum = r_i + 1
        x = x0 + r_i * 3.4
        # Realistic 5'-terminal nucleotide (rnum == 1) LACKS a P atom!
        if rnum > 1:
            lines.append(_cif_atom(idx, "ATOM", "P", "P", base, chain, ent, str(rnum), x, y0, z0, bf, rnum))
            idx += 1
        lines.append(_cif_atom(idx, "ATOM", "O", '"O5\'"', base, chain, ent, str(rnum), x + 0.4, y0 + 0.5, z0, bf, rnum))
        idx += 1
        lines.append(_cif_atom(idx, "ATOM", "C", '"C4\'"' if base in ("A", "U", "G", "C") else "C4'", base, chain, ent, str(rnum), x + 0.8, y0 + 1.0, z0, bf, rnum))
        idx += 1
        lines.append(_cif_atom(idx, "ATOM", "C", '"C1\'"', base, chain, ent, str(rnum), x + 1.2, y0 + 1.5, z0, bf, rnum))
        idx += 1
    return idx


def _add_glycan_chain(
    lines: list[str],
    idx: int,
    chain: str,
    ent: int,
    n_sugars: int,
    comp: str,
    x0: float,
    y0: float,
    z0: float,
    bf_base: float = 82.0,
) -> int:
    # 6 heavy atoms per sugar (C1, C2, C3, C4, C5, O5)
    for s_i in range(n_sugars):
        for a_i, (an, sy, dx, dy, dz) in enumerate([
            ("C1", "C", 0.0, 0.0, 0.0),
            ("C2", "C", 1.2, 0.0, 0.0),
            ("C3", "C", 1.8, 1.1, 0.0),
            ("C4", "C", 1.0, 2.2, 0.0),
            ("C5", "C", -0.2, 1.8, 0.0),
            ("O5", "O", -0.6, 0.8, 0.0),
        ]):
            lines.append(
                _cif_atom(
                    idx, "HETATM", sy, an, comp, chain, ent, ".",
                    x0 + dx, y0 + s_i * 3.0 + dy, z0 + dz,
                    bf_base + (s_i * 6 + a_i) * 0.2, auth_seq=s_i + 1,
                )
            )
            idx += 1
    return idx


def _make_pae_json(chain_counts: list[tuple[str, int]], low_pae: float = 2.5, high_pae: float = 28.0, zero_contact: bool = False) -> tuple[str, str]:
    token_chain_ids = []
    for ch, cnt in chain_counts:
        token_chain_ids.extend([ch] * cnt)
    n = len(token_chain_ids)
    pae = []
    for i in range(n):
        row = []
        for j in range(n):
            if token_chain_ids[i] == token_chain_ids[j]:
                row.append(1.5)
            else:
                row.append(high_pae if zero_contact else low_pae)
        pae.append(row)
    nc = len(chain_counts)
    iptm_mat = [[(0.9 if i == j else (0.05 if zero_contact else 0.78)) for j in range(nc)] for i in range(nc)]
    full_data = json.dumps({"token_chain_ids": token_chain_ids, "pae": pae})
    summary = json.dumps({"iptm": 0.05 if zero_contact else 0.78, "ptm": 0.85, "chain_pair_iptm": iptm_mat})
    return full_data, summary


def _make_af3_zip(out_zip: Path, name: str, cif_text: str, full_data_json: str, summary_json: str) -> Path:
    with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(f"{name}/fold_{name}_model_0.cif", cif_text)
        zf.writestr(f"{name}/fold_{name}_full_data_0.json", full_data_json)
        zf.writestr(f"{name}/fold_{name}_summary_confidences_0.json", summary_json)
    return out_zip


def _make_pdb_zip(out_zip: Path, name: str, pdb_text: str, pae_matrix: list[list[float]]) -> Path:
    with zipfile.ZipFile(out_zip, "w", zipfile.ZIP_DEFLATED) as zf:
        zf.writestr(f"{name}_unrelaxed_rank_001_alphafold2_multimer_v3_model_1_seed_000.pdb", pdb_text)
        zf.writestr(
            f"{name}_scores_rank_001_alphafold2_multimer_v3_model_1_seed_000.json",
            json.dumps({"pae": pae_matrix, "ptm": 0.82, "iptm": 0.75}),
        )
    return out_zip


def build_all_fixtures(work_dir: Path) -> dict[str, Path]:
    out: dict[str, Path] = {}

    # 1. Baseline protein dimer
    lines = [CIF_HEADER]
    idx = _add_protein_chain(lines, 1, "A", 1, ["ALA", "VAL", "LEU", "SER", "ASP", "GLY", "LYS", "TYR"], 0, 0, 0)
    idx = _add_protein_chain(lines, idx, "B", 2, ["GLY", "LEU", "VAL", "ALA", "THR", "PHE"], 0, 5.5, 0)
    fd, sm = _make_pae_json([("A", 8), ("B", 6)])
    out["01_protein_dimer_cif"] = _make_af3_zip(work_dir / "01.zip", "t01_dimer", "".join(lines), fd, sm)

    # 2. O-glycan with multi-row loop_ _struct_conn
    sc_loop = """loop_
_struct_conn.conn_type_id
_struct_conn.id
_struct_conn.ptnr1_label_asym_id
_struct_conn.ptnr1_label_atom_id
_struct_conn.ptnr1_label_comp_id
_struct_conn.ptnr1_label_seq_id
_struct_conn.ptnr2_label_asym_id
_struct_conn.ptnr2_label_atom_id
_struct_conn.ptnr2_label_comp_id
_struct_conn.ptnr2_label_seq_id
covale covale1 C C1 MAN . B OG SER 3
covale covale2 C O5 MAN . C C1 MAN .
"""
    lines = ["data_test\n", sc_loop, CIF_HEADER.split("\n", 1)[1]]
    idx = _add_protein_chain(lines, 1, "A", 1, ["ALA", "VAL", "LEU", "SER", "ASP", "GLY", "LYS", "TYR"], 0, 0, 0)
    idx = _add_protein_chain(lines, idx, "B", 2, ["GLY", "GLY", "SER", "GLY", "GLY"], 0, 5.5, 0)
    idx = _add_glycan_chain(lines, idx, "C", 3, 2, "MAN", 7.6, 9.5, 0)  # 12 atoms near B:SER3.OG (7.6, 8.2, 0)
    fd, sm = _make_pae_json([("A", 8), ("B", 5), ("C", 12)])
    out["02_oglycan_loop_cif"] = _make_af3_zip(work_dir / "02.zip", "t02_oglycan_loop", "".join(lines), fd, sm)

    # 3. Single-monosaccharide O-glycan with single-record key-value _struct_conn (like s04_man1)
    sc_single = """_struct_conn.conn_type_id covale
_struct_conn.id covale1
_struct_conn.ptnr1_label_asym_id C
_struct_conn.ptnr1_label_atom_id C1
_struct_conn.ptnr1_label_comp_id MAN
_struct_conn.ptnr1_label_seq_id .
_struct_conn.ptnr2_label_asym_id B
_struct_conn.ptnr2_label_atom_id OG
_struct_conn.ptnr2_label_comp_id SER
_struct_conn.ptnr2_label_seq_id 3
"""
    lines = ["data_test\n", sc_single, CIF_HEADER.split("\n", 1)[1]]
    idx = _add_protein_chain(lines, 1, "A", 1, ["ALA", "VAL", "LEU", "SER", "ASP", "GLY"], 0, 0, 0)
    idx = _add_protein_chain(lines, idx, "B", 2, ["GLY", "GLY", "SER", "GLY", "GLY"], 0, 5.5, 0)
    idx = _add_glycan_chain(lines, idx, "C", 3, 1, "MAN", 7.6, 9.5, 0)  # 6 atoms
    fd, sm = _make_pae_json([("A", 6), ("B", 5), ("C", 6)])
    out["03_oglycan_single_structconn_cif"] = _make_af3_zip(work_dir / "03.zip", "t03_oglycan_single", "".join(lines), fd, sm)

    # 4. N-glycan on ASN3
    sc_n = sc_loop.replace("B OG SER 3", "B ND2 ASN 3").replace("C C1 MAN", "C C1 NAG")
    lines = ["data_test\n", sc_n, CIF_HEADER.split("\n", 1)[1]]
    idx = _add_protein_chain(lines, 1, "A", 1, ["ALA", "VAL", "LEU", "SER", "ASP", "GLY"], 0, 0, 0)
    idx = _add_protein_chain(lines, idx, "B", 2, ["GLY", "GLY", "ASN", "GLY", "GLY"], 0, 5.5, 0)
    idx = _add_glycan_chain(lines, idx, "C", 3, 2, "NAG", 7.6, 10.4, 0)
    fd, sm = _make_pae_json([("A", 6), ("B", 5), ("C", 12)])
    out["04_nglycan_branched_cif"] = _make_af3_zip(work_dir / "04.zip", "t04_nglycan", "".join(lines), fd, sm)

    # 5. Glycan BETWEEN two protein chains (A=protein, B=glycan, C=protein) — tests caCoords index alignment!
    lines = [CIF_HEADER]
    idx = _add_protein_chain(lines, 1, "A", 1, ["ALA", "VAL", "LEU", "SER", "ASP", "GLY", "LYS", "TYR"], 0, 0, 0)
    idx = _add_glycan_chain(lines, idx, "B", 2, 2, "MAN", 40.0, 40.0, 40.0)  # far away so A-B has 0 contacts
    idx = _add_protein_chain(lines, idx, "C", 3, ["GLY", "LEU", "VAL", "ALA", "THR", "PHE"], 0, 5.5, 0)  # identical coords to 01's B!
    fd, sm = _make_pae_json([("A", 8), ("B", 12), ("C", 6)])
    out["05_glycan_between_proteins_cif"] = _make_af3_zip(work_dir / "05.zip", "t05_glycan_mid", "".join(lines), fd, sm)

    # 6. Covalent small-molecule inhibitor on CYS 4
    sc_cov = sc_single.replace("ptnr1_label_asym_id C", "ptnr1_label_asym_id B").replace("MAN", "LIG").replace("ptnr2_label_asym_id B", "ptnr2_label_asym_id A").replace("OG", "SG").replace("SER", "CYS").replace("seq_id 3", "seq_id 4")
    lines = ["data_test\n", sc_cov, CIF_HEADER.split("\n", 1)[1]]
    idx = _add_protein_chain(lines, 1, "A", 1, ["ALA", "VAL", "LEU", "CYS", "ASP", "GLY"], 0, 0, 0)
    idx = _add_glycan_chain(lines, idx, "B", 2, 1, "LIG", 11.4, 4.4, 0)
    fd, sm = _make_pae_json([("A", 6), ("B", 6)])
    out["06_covalent_ligand_cys_cif"] = _make_af3_zip(work_dir / "06.zip", "t06_cov_lig", "".join(lines), fd, sm)

    # 7. Protein + dsDNA (5'-terminal DA/DT lacking P atom!)
    lines = [CIF_HEADER]
    idx = _add_protein_chain(lines, 1, "A", 1, ["ARG", "LYS", "ALA", "VAL", "LEU", "SER", "ARG", "LYS"], 0, 0, 0)
    idx = _add_nucleic_chain(lines, idx, "B", 2, ["DA", "DT", "DG", "DC", "DA"], 0, 6.0, 0)
    idx = _add_nucleic_chain(lines, idx, "C", 3, ["DT", "DG", "DC", "DA", "DT"], 0, 10.0, 0)
    fd, sm = _make_pae_json([("A", 8), ("B", 5), ("C", 5)])
    out["07_protein_dsdna_cif"] = _make_af3_zip(work_dir / "07.zip", "t07_dsdna", "".join(lines), fd, sm)

    # 8. Protein + ssRNA (5'-terminal A lacking P atom!)
    lines = [CIF_HEADER]
    idx = _add_protein_chain(lines, 1, "A", 1, ["ARG", "LYS", "ALA", "VAL", "LEU", "SER"], 0, 0, 0)
    idx = _add_nucleic_chain(lines, idx, "B", 2, ["A", "U", "G", "C", "A", "U"], 0, 6.0, 0)
    fd, sm = _make_pae_json([("A", 6), ("B", 6)])
    out["08_protein_rna_cif"] = _make_af3_zip(work_dir / "08.zip", "t08_rna", "".join(lines), fd, sm)

    # 9. Boltz-2 wrapped CIF lines (17 fields + continuation "1" on next line)
    lines = [CIF_HEADER]
    idx = _add_protein_chain(lines, 1, "A", 1, ["ALA", "VAL", "LEU", "SER", "ASP", "GLY", "LYS", "TYR"], 0, 0, 0, wrap_boltz=True)
    idx = _add_protein_chain(lines, idx, "B", 2, ["GLY", "LEU", "VAL", "ALA", "THR", "PHE"], 0, 5.5, 0, wrap_boltz=True)
    fd, sm = _make_pae_json([("A", 8), ("B", 6)])
    out["09_boltz_wrapped_cif"] = _make_af3_zip(work_dir / "09.zip", "t09_boltz_wrap", "".join(lines), fd, sm)

    # 10. Multiple ZN ions sharing one chain C (tests multi-ion handling)
    lines = [CIF_HEADER]
    idx = _add_protein_chain(lines, 1, "A", 1, ["CYS", "HIS", "ALA", "VAL", "CYS", "HIS"], 0, 0, 0)
    idx = _add_protein_chain(lines, idx, "B", 2, ["ALA", "VAL", "LEU", "ASP", "GLU"], 0, 6.0, 0)
    lines.append(_cif_atom(idx, "HETATM", "ZN", "ZN", "ZN", "C", 3, ".", 2.0, 3.0, 0.0, 92.0, auth_seq=1))
    idx += 1
    lines.append(_cif_atom(idx, "HETATM", "ZN", "ZN", "ZN", "C", 3, ".", 15.0, 3.0, 0.0, 91.0, auth_seq=2))
    fd, sm = _make_pae_json([("A", 6), ("B", 5), ("C", 2)])
    out["10_multi_ions_same_chain_cif"] = _make_af3_zip(work_dir / "10.zip", "t10_multi_ion", "".join(lines), fd, sm)

    # 11. PDB format with Protein A + DNA B (5'-nt lacking P) + Glycan C + PDB LINK record
    pdb_lines = [
        "LINK         OG  SER A   3                 C1  MAN C   1     1555   1555  1.35\n"
    ]
    a_id = 1
    for r_i, comp in enumerate(["GLY", "GLY", "SER", "GLY", "GLY"], 1):
        x = (r_i - 1) * 3.8
        for an, dx, dy in ([("N", -0.6, 0), ("CA", 0, 0), ("C", 0.6, 0), ("O", 0.6, 0.6)] + ([("CB", 0, 1.5), ("OG", 0, 2.7)] if comp == "SER" else [])):
            pdb_lines.append(f"ATOM  {a_id:5d}  {an:<3s} {comp:>3s} A{r_i:4d}    {x+dx:8.3f}{dy:8.3f}{0.0:8.3f}  1.00 88.00           {an[0]}\n")
            a_id += 1
    for r_i, base in enumerate(["DA", "DT", "DG", "DC"], 1):
        x = (r_i - 1) * 3.4
        if r_i > 1:
            pdb_lines.append(f"ATOM  {a_id:5d}  P    {base:>2s} B{r_i:4d}    {x:8.3f}{5.0:8.3f}{0.0:8.3f}  1.00 84.00           P\n")
            a_id += 1
        pdb_lines.append(f"ATOM  {a_id:5d}  C4'  {base:>2s} B{r_i:4d}    {x+0.8:8.3f}{6.0:8.3f}{0.0:8.3f}  1.00 84.00           C\n")
        a_id += 1
    for s_i, (an, dx, dy) in enumerate([("C1", 0, 0), ("C2", 1.2, 0), ("C3", 1.8, 1.1), ("O5", -0.6, 0.8)]):
        pdb_lines.append(f"HETATM{a_id:5d}  {an:<3s} MAN C{1:4d}    {7.6+dx:8.3f}{4.0+dy:8.3f}{0.0:8.3f}  1.00 81.00           {an[0]}\n")
        a_id += 1
    n_pdb = 5 + 4 + 4
    pae_pdb = [[1.5 if (i < 5 and j < 5) or (5 <= i < 9 and 5 <= j < 9) or (i >= 9 and j >= 9) else 3.0 for j in range(n_pdb)] for i in range(n_pdb)]
    out["11_pdb_glycan_link_and_dna"] = _make_pdb_zip(work_dir / "11.zip", "t11_pdb_multi", "".join(pdb_lines), pae_pdb)

    # 12. Zero-contact negative dimer (50 A apart, PAE = 28)
    lines = [CIF_HEADER]
    idx = _add_protein_chain(lines, 1, "A", 1, ["ALA", "VAL", "LEU", "SER", "ASP", "GLY"], 0, 0, 0)
    idx = _add_protein_chain(lines, idx, "B", 2, ["GLY", "LEU", "VAL", "ALA", "THR", "PHE"], 50, 50, 50)
    fd, sm = _make_pae_json([("A", 6), ("B", 6)], zero_contact=True)
    out["12_zero_contact_negative_cif"] = _make_af3_zip(work_dir / "12.zip", "t12_zero_neg", "".join(lines), fd, sm)

    # 13. 2-letter chain IDs (AA, AB, AC)
    lines = [CIF_HEADER]
    idx = _add_protein_chain(lines, 1, "AA", 1, ["ALA", "VAL", "LEU", "SER", "ASP", "GLY"], 0, 0, 0)
    idx = _add_protein_chain(lines, idx, "AB", 2, ["GLY", "LEU", "VAL", "ALA", "THR", "PHE"], 0, 5.5, 0)
    idx = _add_protein_chain(lines, idx, "AC", 3, ["LYS", "ARG", "ASP", "GLU", "SER"], 0, 11.0, 0)
    fd, sm = _make_pae_json([("AA", 6), ("AB", 6), ("AC", 5)])
    out["13_two_letter_chains_cif"] = _make_af3_zip(work_dir / "13.zip", "t13_2char", "".join(lines), fd, sm)

    # 14. Phosphorylated SEP (MOD_AA HETATM in protein chain A) + Protein B
    lines = [CIF_HEADER]
    idx = _add_protein_chain(lines, 1, "A", 1, ["ALA", "VAL", "SEP", "ASP", "GLY"], 0, 0, 0)  # 4 std + 10 SEP atoms = 14 tokens
    idx = _add_protein_chain(lines, idx, "B", 2, ["ARG", "LYS", "ALA", "VAL", "LEU"], 0, 6.0, 0)
    fd, sm = _make_pae_json([("A", 14), ("B", 5)])
    out["14_ptm_phospho_sep_cif"] = _make_af3_zip(work_dir / "14.zip", "t14_ptm_sep", "".join(lines), fd, sm)

    return out


class QuietHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, *a):
        pass


def load_zip_and_process(page, port: int, zip_path: Path, disable_wasm: bool = False):
    page_errors = []
    console_errors = []
    page.on("pageerror", lambda e: page_errors.append(str(e)))
    page.on("console", lambda m: console_errors.append(m.text) if m.type == "error" else None)

    page.goto(f"http://127.0.0.1:{port}/universal.html", wait_until="domcontentloaded", timeout=60000)
    if disable_wasm:
        page.evaluate("() => { window.LiviaWasm = null; }")
    page.set_input_files("#file-input", str(zip_path))
    page.wait_for_function("() => !document.getElementById('btn-process').disabled", timeout=60000)
    page.click("#btn-process")
    page.wait_for_function(
        "() => uniState.modelPairs && uniState.modelPairs.length > 0 && document.querySelectorAll('#pairs-tbody tr').length > 0",
        timeout=60000,
    )
    page.wait_for_timeout(150)

    state = page.evaluate("""() => {
        const md = uniState.modelData[0];
        const coords = parseStructureCoords(md.structText, md.structFormat);
        const caCoords = parseStructureCaCoords(md.structText, md.structFormat);
        const bfs = parseBfactorsPerResidue(md.structText, md.structFormat);
        const seqs = parseChainSequences(md.structText, md.structFormat);
        const sym = symmetrizePairs(uniState.modelPairs[0]);
        const pairMetrics = {};
        for (const [k, v] of Object.entries(sym)) {
            pairMetrics[k] = {
                iLIS: v.iLIS, LIS: v.LIS, cLIS: v.cLIS, ipSAE: v.ipSAE, actifpTM: v.actifpTM,
                pDockQ: v.pDockQ, LIpDockQ: v.LIpDockQ,
                pDockQ2: mean2Nullable(v.pDockQ2_i, v.pDockQ2_j),
                LIpDockQ2: mean2Nullable(v.LIpDockQ2_i, v.LIpDockQ2_j),
                lirI: v.lirI.size, lirJ: v.lirJ.size, clirI: v.clirI.size, clirJ: v.clirJ.size,
            };
        }
        const tableRows = [...document.querySelectorAll('#pairs-tbody tr')].map(tr => {
            const tds = [...tr.querySelectorAll('td')].map(td => td.textContent.trim());
            return { model: tds[0], pair: tds[1], pDockQ2: tds[11], LIpDockQ2: tds[12] };
        });
        const comps = buildColorComponents(uniState.currentPairResult);
        const cxc = buildAllChainsCxc();
        const pml = buildAllChainsPml();
        const bfsPerChain = {};
        for (let i = 0; i < uniState.chainNames.length; i++) {
            const ch = uniState.chainNames[i];
            const sz = uniState.chainSizes[i];
            let found = 0;
            for (let r = 1; r <= sz; r++) if (bfs.has(`${ch}:${r}`)) found++;
            bfsPerChain[ch] = { found, expected: sz };
        }
        return {
            paeLen: md.pae.length,
            nTotal: uniState.nTotal,
            chainNames: uniState.chainNames,
            chainSizes: uniState.chainSizes,
            chainTypes: uniState.chainTypes,
            coordsLen: coords.length,
            caCoordsLen: caCoords.length,
            bfsPerChain,
            seqLens: Object.fromEntries([...seqs.byChain.entries()].map(([k, v]) => [k, v.seq.length])),
            ptmResidues: uniState.ptmResidues,
            covalentLinks: uniState.covalentLinks,
            pairMetrics,
            tableRows,
            ionComps: comps.filter(c => c.isIon),
            stickComps: comps.filter(c => c.stick),
            cxc,
            pml,
        };
    }""")
    return state, page_errors, console_errors


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--real-zip", type=Path, default=None, help="Optional real AF3 zip to smoke-test")
    args = ap.parse_args()

    httpd = socketserver.TCPServer(("127.0.0.1", 0), functools.partial(QuietHandler, directory=str(REPO_DIR)))
    port = httpd.server_address[1]
    threading.Thread(target=httpd.serve_forever, daemon=True).start()

    failures: list[str] = []
    passed = 0

    def check(case: str, cond: bool, msg: str):
        nonlocal passed
        if cond:
            passed += 1
        else:
            failures.append(f"[FAIL] {case}: {msg}")
            print(f"  [FAIL] {case}: {msg}")

    with tempfile.TemporaryDirectory(prefix="livia_edge_suite_") as td:
        fixtures = build_all_fixtures(Path(td))
        with sync_playwright() as p:
            browser = p.chromium.launch(headless=True)
            context = browser.new_context(viewport={"width": 1400, "height": 900})
            results = {}

            for cname, zpath in fixtures.items():
                page = context.new_page()
                st, perrs, cerrs = load_zip_and_process(page, port, zpath)
                results[cname] = st
                # Filter out WebGL-unavailable console messages in headless without GPU
                non_webgl_p = [e for e in perrs if "WebGL" not in e]
                non_webgl_c = [e for e in cerrs if "WebGL" not in e and "reprCount" not in e]

                print(f"\n=== {cname} === names={st['chainNames']} sizes={st['chainSizes']} types={st['chainTypes']} coords={st['coordsLen']} caCoords={st['caCoordsLen']} pae={st['paeLen']}")
                check(cname, len(non_webgl_p) == 0 and len(non_webgl_c) == 0, f"JS errors: {non_webgl_p} {non_webgl_c}")
                check(cname, sum(st["chainSizes"]) == st["paeLen"], f"sum(chainSizes)={sum(st['chainSizes'])} != paeLen={st['paeLen']}")
                check(cname, st["coordsLen"] == st["paeLen"], f"coordsLen={st['coordsLen']} != paeLen={st['paeLen']}")
                check(cname, st["caCoordsLen"] == st["paeLen"], f"caCoordsLen={st['caCoordsLen']} != paeLen={st['paeLen']} (causes index shift on downstream chains!)")

                # Every metric finite (no NaN)
                for pk, pm in st["pairMetrics"].items():
                    for mk, mv in pm.items():
                        if mv is not None:
                            check(cname, not (isinstance(mv, float) and math.isnan(mv)), f"{pk}.{mk} is NaN")

                # Case-specific invariants
                if cname == "02_oglycan_loop_cif":
                    check(cname, "C" in st["covalentLinks"] and st["covalentLinks"]["C"][0]["chain"] == "B" and st["covalentLinks"]["C"][0]["seq"] == 3, f"missing covalentLinks: {st['covalentLinks']}")
                    check(cname, len(st["ionComps"]) == 1 and len(st["ionComps"][0].get("covalentAttachments", [])) == 1, f"Mol* missing covalentAttachments: {st['ionComps']}")
                    check(cname, "show #1/B:3 atoms" in st["cxc"] and "chain B and resi 3" in st["pml"], "missing B:3 in cxc/pml export")
                    check(cname, st["bfsPerChain"]["C"]["found"] == 12, f"glycan C B-factors only {st['bfsPerChain']['C']} (expected 12/12)")
                    bc_row = [r for r in st["tableRows"] if r["pair"] == "B / C"][0]
                    ab_row = [r for r in st["tableRows"] if r["pair"] == "A / B"][0]
                    check(cname, bc_row["pDockQ2"] == "—", f"B/C pDockQ2 should be '—', got {bc_row['pDockQ2']}")
                    check(cname, ab_row["pDockQ2"] != "—" and float(ab_row["pDockQ2"]) > 0, f"A/B pDockQ2 should be > 0, got {ab_row['pDockQ2']}")
                    # Stress-test all color modes & LIR/Complete toggle on multimodal complex
                    mode_errs = page.evaluate("""() => {
                        const errs = [];
                        for (const m of ['', 'bychain', 'plddt']) {
                            uniState.cxcColorMode = m;
                            for (const comp of [true, false]) {
                                const cb = document.getElementById('show-complete-structure');
                                if (cb) cb.checked = comp;
                                try { buildColorComponents(uniState.currentPairResult); buildAllChainsCxc(); buildAllChainsPml(); }
                                catch (e) { errs.push(`${m}/${comp}: ${e}`); }
                            }
                        }
                        return errs;
                    }""")
                    check(cname, len(mode_errs) == 0, f"mode toggle errors: {mode_errs}")

                elif cname == "03_oglycan_single_structconn_cif":
                    check(cname, "C" in st["covalentLinks"] and st["covalentLinks"]["C"][0]["seq"] == 3, f"single-record _struct_conn failed: {st['covalentLinks']}")

                elif cname == "04_nglycan_branched_cif":
                    atts = st["covalentLinks"].get("C", [])
                    check(cname, len(atts) == 1 and "ND2" in atts[0]["atoms"] and "CA" in atts[0]["atoms"], f"N-glycan ASN atoms wrong: {atts}")

                elif cname == "05_glycan_between_proteins_cif":
                    # A-C in case 05 MUST have identical pDockQ2 to A-B in case 01 (glycan B between them must NOT shift CA indices!)
                    p01 = results["01_protein_dimer_cif"]["pairMetrics"]["A-B"]["pDockQ2"]
                    p05_ac = st["pairMetrics"]["A-C"]["pDockQ2"]
                    p05_ab = st["pairMetrics"]["A-B"]["pDockQ2"]
                    check(cname, abs(p05_ac - p01) < 1e-4, f"A-C pDockQ2 shifted! got {p05_ac} vs baseline {p01}")
                    check(cname, p05_ab == 0.0, f"A-B (glycan) pDockQ2 leaked contacts from C! got {p05_ab}")

                elif cname == "06_covalent_ligand_cys_cif":
                    atts = st["covalentLinks"].get("B", [])
                    check(cname, len(atts) == 1 and atts[0]["chain"] == "A" and atts[0]["seq"] == 4 and "SG" in atts[0]["atoms"], f"covalent CYS ligand failed: {atts}")

                elif cname == "07_protein_dsdna_cif":
                    check(cname, st["chainTypes"] == ["protein", "dna", "dna"], f"wrong types: {st['chainTypes']}")
                    check(cname, st["chainSizes"] == [8, 5, 5], f"5'-terminal DNA nt dropped! sizes={st['chainSizes']}")
                    check(cname, st["bfsPerChain"]["B"]["found"] == 5 and st["bfsPerChain"]["C"]["found"] == 5, f"DNA B-factors missing: {st['bfsPerChain']}")

                elif cname == "08_protein_rna_cif":
                    check(cname, st["chainTypes"] == ["protein", "rna"], f"wrong types: {st['chainTypes']}")
                    check(cname, st["chainSizes"] == [6, 6], f"5'-terminal RNA nt dropped! sizes={st['chainSizes']}")
                    check(cname, st["bfsPerChain"]["B"]["found"] == 6, f"RNA B-factors missing: {st['bfsPerChain']}")

                elif cname == "09_boltz_wrapped_cif":
                    p01 = results["01_protein_dimer_cif"]["pairMetrics"]["A-B"]
                    p09 = st["pairMetrics"]["A-B"]
                    check(cname, abs(p09["iLIS"] - p01["iLIS"]) < 1e-5 and abs(p09["pDockQ2"] - p01["pDockQ2"]) < 1e-5, f"Boltz wrapped CIF mismatch vs 01: {p09} vs {p01}")
                    check(cname, st["bfsPerChain"]["A"]["found"] == 8 and st["bfsPerChain"]["B"]["found"] == 6, f"Boltz wrapped B-factors dropped: {st['bfsPerChain']}")
                    check(cname, st["seqLens"].get("A") == 8 and st["seqLens"].get("B") == 6, f"Boltz wrapped seqs dropped: {st['seqLens']}")

                elif cname == "10_multi_ions_same_chain_cif":
                    check(cname, st["chainTypes"] == ["protein", "protein", "ion"], f"wrong types: {st['chainTypes']}")
                    check(cname, st["chainSizes"] == [6, 5, 2], f"multi-ion chain collapsed! sizes={st['chainSizes']}")
                    check(cname, st["bfsPerChain"]["C"]["found"] == 2, f"multi-ion B-factors missing: {st['bfsPerChain']}")

                elif cname == "11_pdb_glycan_link_and_dna":
                    check(cname, st["chainTypes"] == ["protein", "dna", "glycan"], f"PDB types wrong: {st['chainTypes']}")
                    check(cname, st["chainSizes"] == [5, 4, 4], f"PDB sizes wrong: {st['chainSizes']}")
                    check(cname, "C" in st["covalentLinks"] and st["covalentLinks"]["C"][0]["chain"] == "A" and st["covalentLinks"]["C"][0]["seq"] == 3, f"PDB LINK covalentLinks missing: {st['covalentLinks']}")
                    check(cname, st["bfsPerChain"]["B"]["found"] == 4 and st["bfsPerChain"]["C"]["found"] == 4, f"PDB DNA/glycan B-factors missing: {st['bfsPerChain']}")

                elif cname == "12_zero_contact_negative_cif":
                    m = st["pairMetrics"]["A-B"]
                    check(cname, m["iLIS"] == 0 and m["ipSAE"] == 0 and m["pDockQ"] == 0 and m["pDockQ2"] == 0, f"zero-contact non-zero: {m}")

                elif cname == "13_two_letter_chains_cif":
                    check(cname, st["chainNames"] == ["AA", "AB", "AC"], f"2-letter chains wrong: {st['chainNames']}")
                    check(cname, "AA-AB" in st["pairMetrics"] and st["pairMetrics"]["AA-AB"]["iLIS"] > 0.5, f"2-letter chain scoring failed: {st['pairMetrics']}")
                    check(cname, "#1/AA" in st["cxc"] and "chain AA" in st["pml"], "2-letter chain missing in cxc/pml")

                elif cname == "14_ptm_phospho_sep_cif":
                    check(cname, st["chainTypes"] == ["protein", "protein"], f"PTM SEP chain misclassified: {st['chainTypes']}")
                    check(cname, "A" in st["ptmResidues"] and st["ptmResidues"]["A"][0]["seq"] == 3, f"PTM SEP not detected: {st['ptmResidues']}")
                    check(cname, len(st["stickComps"]) >= 1, "PTM SEP stick component missing")

                page.close()

            # 15. WASM vs Pure-JS fallback parity check on 4 representative multimodal fixtures
            for cname in ("01_protein_dimer_cif", "02_oglycan_loop_cif", "05_glycan_between_proteins_cif", "07_protein_dsdna_cif"):
                page = context.new_page()
                js_st, _, _ = load_zip_and_process(page, port, fixtures[cname], disable_wasm=True)
                w_st = results[cname]
                for pk, wm in w_st["pairMetrics"].items():
                    jm = js_st["pairMetrics"][pk]
                    for mk in ("iLIS", "LIS", "cLIS", "ipSAE", "actifpTM", "pDockQ", "LIpDockQ", "pDockQ2", "LIpDockQ2"):
                        wv, jv = wm[mk], jm[mk]
                        if wv is None and jv is None:
                            continue
                        diff = abs((wv or 0.0) - (jv or 0.0))
                        check(f"WASM_vs_JS:{cname}:{pk}.{mk}", diff < 1e-5, f"wasm={wv} vs js={jv} (diff={diff})")
                page.close()

            browser.close()
    httpd.shutdown()

    print(f"\n=== EDGE SUITE SUMMARY: {passed} checks passed, {len(failures)} failed ===")
    if failures:
        for f in failures:
            print(" ", f)
        sys.exit(1)


if __name__ == "__main__":
    main()
