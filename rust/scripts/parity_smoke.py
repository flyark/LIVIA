#!/usr/bin/env python3
"""Compare lis-rs native CLI output against python/lis.py across all example prediction bundles."""

import csv
import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

RUST_ROOT = Path(__file__).resolve().parents[1]
REPO_ROOT = RUST_ROOT.parent
EXAMPLES_DIR = REPO_ROOT / "examples"
LIS_PY = REPO_ROOT / "python" / "lis.py"
LIS_RS = RUST_ROOT / "target" / "release" / "lis-rs"
OUT_DIR = RUST_ROOT / "test_out" / "parity_smoke" 


def parse_csv(path: Path):
    with path.open("r", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def compare_tables(name: str, py_rows, rs_rows):
    errors = []
    if len(py_rows) != len(rs_rows):
        errors.append(f"{name}: row count mismatch py={len(py_rows)} vs rs={len(rs_rows)}")
        return errors

    key_cols = ("saved folder", "model_number", "chain_i", "chain_j")
    py_map = {tuple(r.get(k, "") for k in key_cols): r for r in py_rows}
    rs_map = {tuple(r.get(k, "") for k in key_cols): r for r in rs_rows}

    if set(py_map.keys()) != set(rs_map.keys()):
        only_py = list(set(py_map.keys()) - set(rs_map.keys()))[:5]
        only_rs = list(set(rs_map.keys()) - set(py_map.keys()))[:5]
        errors.append(f"{name}: key mismatch only_py={only_py} only_rs={only_rs}")
        return errors

    num_cols = [
        "iLIS", "iLIA", "iLISA", "LIS", "cLIS", "LIA", "cLIA",
        "LIR", "cLIR", "LIR_A", "LIR_B", "cLIR_A", "cLIR_B",
        "ipSAE", "actifpTM", "ipTM", "pTM", "confidence",
        "LIpLDDT", "cLIpLDDT", "LIpLDDT_A", "LIpLDDT_B", "cLIpLDDT_A", "cLIpLDDT_B",
        "pLDDT_A", "pLDDT_B", "pDockQ", "LIpDockQ", "pDockQ2", "LIpDockQ2",
        "len_A", "len_B",
    ]
    str_cols = [
        "platform",
        "LIR_indices_A", "LIR_indices_B",
        "cLIR_indices_A", "cLIR_indices_B",
    ]

    max_diff = 0.0
    for k, pr in py_map.items():
        rr = rs_map[k]
        for col in str_cols:
            if col in pr or col in rr:
                pv = (pr.get(col) or "").strip()
                rv = (rr.get(col) or "").strip()
                if pv != rv:
                    errors.append(f"{name} {k} [{col}]: py={pv!r} vs rs={rv!r}")
        for col in num_cols:
            if col not in pr and col not in rr:
                continue
            ps = (pr.get(col) or "").strip()
            rs_val = (rr.get(col) or "").strip()
            if ps == "" and rs_val == "":
                continue
            if (ps == "") != (rs_val == ""):
                errors.append(f"{name} {k} [{col}]: blank mismatch py={ps!r} vs rs={rs_val!r}")
                continue
            pf = float(ps)
            rf = float(rs_val)
            if math.isnan(pf) and math.isnan(rf):
                continue
            diff = abs(pf - rf)
            if diff > max_diff:
                max_diff = diff
            tol = 2e-3 if "pLDDT" in col else 1.5e-3
            if diff > tol:
                errors.append(f"{name} {k} [{col}]: diff={diff:.6g} (py={ps}, rs={rs_val})")
    return errors, max_diff


def main():
    if OUT_DIR.exists():
        shutil.rmtree(OUT_DIR)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    zips = sorted(
        p for p in EXAMPLES_DIR.glob("*.zip")
        if "interactome" not in p.name and "fly_pathway" not in p.name
    )
    print(f"Testing {len(zips)} prediction bundles in {EXAMPLES_DIR}...")
    total_py_s = 0.0
    total_rs_s = 0.0
    all_ok = True

    for zp in zips:
        py_csv = OUT_DIR / f"{zp.stem}.py.csv"
        rs_csv = OUT_DIR / f"{zp.stem}.rs.csv"

        py_exe = os.environ.get("PYTHON", sys.executable)
        t0 = time.perf_counter()
        p_py = subprocess.run(
            [py_exe, str(LIS_PY), str(zp), "-d", str(OUT_DIR), "-o", py_csv.name],
            capture_output=True,
            text=True,
        )
        dt_py = time.perf_counter() - t0

        t1 = time.perf_counter()
        p_rs = subprocess.run(
            [str(LIS_RS), str(zp), "-d", str(OUT_DIR), "-o", rs_csv.name],
            capture_output=True,
            text=True,
        )
        dt_rs = time.perf_counter() - t1

        if p_py.returncode != 0:
            print(f"  [SKIP non-structure archive] {zp.name}: lis.py rc={p_py.returncode}")
            continue
        if p_rs.returncode != 0:
            print(f"  [FAIL] {zp.name}: lis-rs rc={p_rs.returncode}\n{p_rs.stderr}")
            all_ok = False
            continue

        total_py_s += dt_py
        total_rs_s += dt_rs
        py_rows = parse_csv(py_csv)
        rs_rows = parse_csv(rs_csv)
        res = compare_tables(zp.name, py_rows, rs_rows)
        if isinstance(res, list):
            errs, max_diff = res, 999.0
        else:
            errs, max_diff = res
        speedup = dt_py / max(dt_rs, 1e-6)
        if errs:
            all_ok = False
            print(f"  [DIFF] {zp.name}: {len(errs)} diffs (first 5):")
            for e in errs[:5]:
                print("    ", e)
        else:
            print(
                f"  [PASS] {zp.name:<34s} rows={len(py_rows):3d} "
                f"max_abs_diff={max_diff:.2e}  py={dt_py*1000:6.1f}ms  rs={dt_rs*1000:5.1f}ms ({speedup:5.1f}x)"
            )

    print(
        f"\nSummary: py={total_py_s:.3f}s vs lis-rs={total_rs_s:.3f}s "
        f"({total_py_s / max(total_rs_s, 1e-6):.1f}x overall speedup)"
    )
    sys.exit(0 if all_ok else 1)


if __name__ == "__main__":
    main()
