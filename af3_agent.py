"""
AF3 → ChimeraX Agent
=====================
Process AlphaFold3 zip/folder outputs:
  1. Parse PAE matrix + chain boundaries from full_data JSON
  2. Parse Cβ coordinates from CIF → compute contact map
  3. Calculate LIS, cLIS, iLIS, LIR, cLIR for all chain pairs
  4. Average across 5 models
  5. Generate ChimeraX .cxc script for selected chain pair

Usage:
    python af3_agent.py <path_to_zip_or_folder> [--chain-pair A,B] [--output-dir DIR]
"""

import argparse
import glob
import json
import os
import re
import sys
import zipfile
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.spatial.distance import pdist, squareform


# ── PAE Transform ──────────────────────────────────────────────────────────

def transform_pae_matrix(pae_matrix, pae_cutoff=12):
    """Transform PAE to confidence scores: 0→1 (best), cutoff→0, >cutoff→0."""
    transformed = np.zeros_like(pae_matrix)
    within = pae_matrix < pae_cutoff
    transformed[within] = 1 - (pae_matrix[within] / pae_cutoff)
    return transformed


# ── Contact Map from CIF ──────────────────────────────────────────────────

def parse_cif_cb_coords(cif_text):
    """Extract one Cβ (or Cα for GLY) coordinate per residue from CIF text.

    Priority: CB > CA (for GLY) > P (for nucleic acids).
    Returns list of dicts sorted by chain and resnum.
    """
    residues = {}  # (chain, resnum) -> {x, y, z, has_P}
    for line in cif_text.split('\n'):
        if not line.startswith('ATOM') and not line.startswith('HETATM'):
            continue
        parts = line.split()
        if len(parts) < 14:
            continue

        atom_name = parts[3]
        comp_id = parts[5]
        chain = parts[6]
        resnum = int(parts[8])
        key = (chain, resnum)
        x, y, z = float(parts[10]), float(parts[11]), float(parts[12])

        # CB always wins
        if atom_name == 'CB':
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_P': False}
        # CA for GLY only if no CB yet
        elif atom_name == 'CA' and comp_id == 'GLY' and key not in residues:
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_P': False}
        # P for nucleic acids only if nothing else
        elif atom_name == 'P' and key not in residues:
            residues[key] = {'chain': chain, 'resnum': resnum, 'x': x, 'y': y, 'z': z, 'has_P': True}

    return list(residues.values())


def compute_contact_map(cif_text, distance_threshold=8):
    """Compute NxN binary contact map from CIF Cβ coordinates."""
    coords = parse_cif_cb_coords(cif_text)
    n = len(coords)
    if n == 0:
        return np.zeros((0, 0)), []

    xyz = np.array([[c['x'], c['y'], c['z']] for c in coords])
    has_p = np.array([c['has_P'] for c in coords])

    distances = squareform(pdist(xyz))

    # Adjust threshold for phosphorus-containing residues
    adjusted = np.where(
        has_p[:, np.newaxis] | has_p[np.newaxis, :],
        distances - 4, distances
    )

    contact_map = np.where(adjusted < distance_threshold, 1, 0).astype(np.float32)
    return contact_map, coords


# ── Mean LIS per chain pair ───────────────────────────────────────────────

def calculate_mean_lis(matrix, subunit_sizes):
    """Calculate mean of non-zero values per chain-pair submatrix."""
    cum = np.cumsum(subunit_sizes)
    starts = np.concatenate(([0], cum[:-1]))
    n = len(subunit_sizes)
    result = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            sub = matrix[starts[i]:cum[i], starts[j]:cum[j]]
            nz = sub[sub > 0]
            result[i, j] = nz.mean() if len(nz) > 0 else 0.0
    return result


# ── File Discovery ────────────────────────────────────────────────────────

def discover_files(path):
    """Find AF3 output files from a folder or zip.

    Returns dict of model_index → {full_data, summary, cif} file contents.
    """
    models = {}

    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as zf:
            names = zf.namelist()
            for n in names:
                m = re.search(r'_full_data_(\d+)\.json$', n)
                if m:
                    idx = int(m.group(1))
                    prefix = n.replace(f'_full_data_{idx}.json', '')
                    models[idx] = {
                        'full_data': zf.read(n).decode('utf-8'),
                        'summary': zf.read(f'{prefix}_summary_confidences_{idx}.json').decode('utf-8'),
                        'cif': zf.read(f'{prefix}_model_{idx}.cif').decode('utf-8'),
                        'cif_filename': os.path.basename(f'{prefix}_model_{idx}.cif'),
                    }
    else:
        # Folder path
        matches = glob.glob(os.path.join(path, '*_full_data_0.json'))
        if not matches:
            print(f"[AF3] No full_data files found in {path}")
            return {}
        prefix = matches[0].replace('_full_data_0.json', '')
        for idx in range(5):
            fd = f'{prefix}_full_data_{idx}.json'
            sc = f'{prefix}_summary_confidences_{idx}.json'
            cif = f'{prefix}_model_{idx}.cif'
            if os.path.exists(fd) and os.path.exists(sc) and os.path.exists(cif):
                models[idx] = {
                    'full_data': open(fd).read(),
                    'summary': open(sc).read(),
                    'cif': open(cif).read(),
                    'cif_filename': os.path.basename(cif),
                }
    return models


# ── Per-Model Analysis ────────────────────────────────────────────────────

def analyze_model(full_data_str, summary_str, cif_str, pae_cutoff=12, distance_cutoff=8):
    """Analyze a single AF3 model. Returns metrics for all chain pairs."""
    fd = json.loads(full_data_str)
    sc = json.loads(summary_str)

    # Chain boundaries
    chain_ids = fd['token_chain_ids']
    chain_counts = Counter(chain_ids)
    chain_names = list(chain_counts.keys())
    subunit_sizes = list(chain_counts.values())
    cum = np.cumsum(subunit_sizes)
    starts = np.concatenate(([0], cum[:-1]))

    # PAE
    pae = np.array(fd['pae'], dtype=np.float32)
    pae = np.nan_to_num(pae)
    transformed = transform_pae_matrix(pae, pae_cutoff)
    transformed = np.nan_to_num(transformed)
    lia_map = np.where(transformed > 0, 1, 0)

    # Contact map
    contact_map, coords = compute_contact_map(cif_str, distance_cutoff)
    n_residues = sum(subunit_sizes)
    if contact_map.shape[0] != n_residues:
        print(f"[AF3] WARNING: CIF has {contact_map.shape[0]} residues, PAE has {n_residues}. Padding/truncating.")
        cm = np.zeros((n_residues, n_residues), dtype=np.float32)
        mn = min(contact_map.shape[0], n_residues)
        cm[:mn, :mn] = contact_map[:mn, :mn]
        contact_map = cm

    # Combined map (cLIA)
    combined = np.where((transformed > 0) & (contact_map == 1), transformed, 0)

    # LIS/cLIS matrices
    lis_matrix = calculate_mean_lis(transformed, subunit_sizes)
    lis_matrix = np.nan_to_num(lis_matrix)
    clis_matrix = calculate_mean_lis(combined, subunit_sizes)
    clis_matrix = np.nan_to_num(clis_matrix)

    # ipTM
    iptm_matrix = np.array(sc['chain_pair_iptm'], dtype=float)
    iptm_matrix = np.nan_to_num(iptm_matrix)

    # Per chain pair metrics
    n_chains = len(subunit_sizes)
    pairs = {}
    for i in range(n_chains):
        for j in range(n_chains):
            if i == j:
                continue
            si, ei = starts[i], cum[i]
            sj, ej = starts[j], cum[j]

            # LIR indices
            lia_sub = lia_map[si:ei, sj:ej]
            lir_i = set(int(r + 1) for r in np.unique(np.where(lia_sub > 0)[0]))
            lir_j = set(int(r + 1) for r in np.unique(np.where(lia_sub > 0)[1]))

            # cLIR indices
            combined_sub = combined[si:ei, sj:ej]
            clir_i = set(int(r + 1) for r in np.unique(np.where(combined_sub > 0)[0]))
            clir_j = set(int(r + 1) for r in np.unique(np.where(combined_sub > 0)[1]))

            # LIA/cLIA counts
            lia_count = int(np.count_nonzero(lia_sub))
            clia_count = int(np.count_nonzero(combined_sub))

            lis_val = float(lis_matrix[i, j])
            clis_val = float(clis_matrix[i, j])
            ilis_val = float(np.sqrt(lis_val * clis_val))
            iptm_val = float(iptm_matrix[i, j]) if i < iptm_matrix.shape[0] and j < iptm_matrix.shape[1] else 0.0

            key = (chain_names[i], chain_names[j])
            pairs[key] = {
                'chain_i': chain_names[i],
                'chain_j': chain_names[j],
                'LIS': lis_val,
                'cLIS': clis_val,
                'iLIS': ilis_val,
                'ipTM': iptm_val,
                'LIA': lia_count,
                'cLIA': clia_count,
                'LIR_i': lir_i,
                'LIR_j': lir_j,
                'cLIR_i': clir_i,
                'cLIR_j': clir_j,
                'len_i': subunit_sizes[i],
                'len_j': subunit_sizes[j],
            }

    return pairs, chain_names, subunit_sizes


# ── Average Across Models ─────────────────────────────────────────────────

def average_models(all_model_pairs):
    """Average metrics across models, symmetrize (i,j)/(j,i)."""
    if not all_model_pairs:
        return {}

    # Collect all pair keys
    all_keys = set()
    for mp in all_model_pairs:
        all_keys.update(mp.keys())

    averaged = {}
    for key in all_keys:
        vals = [mp[key] for mp in all_model_pairs if key in mp]
        n = len(vals)
        if n == 0:
            continue
        avg = {
            'chain_i': vals[0]['chain_i'],
            'chain_j': vals[0]['chain_j'],
            'LIS': sum(v['LIS'] for v in vals) / n,
            'cLIS': sum(v['cLIS'] for v in vals) / n,
            'ipTM': sum(v['ipTM'] for v in vals) / n,
            'LIA': sum(v['LIA'] for v in vals) / n,
            'cLIA': sum(v['cLIA'] for v in vals) / n,
            'len_i': vals[0]['len_i'],
            'len_j': vals[0]['len_j'],
            # Union of residue sets across models
            'LIR_i': set().union(*(v['LIR_i'] for v in vals)),
            'LIR_j': set().union(*(v['LIR_j'] for v in vals)),
            'cLIR_i': set().union(*(v['cLIR_i'] for v in vals)),
            'cLIR_j': set().union(*(v['cLIR_j'] for v in vals)),
        }
        avg['iLIS'] = float(np.sqrt(avg['LIS'] * avg['cLIS']))
        averaged[key] = avg

    # Symmetrize: for (A,B) and (B,A), combine into a single entry
    symmetric = {}
    seen = set()
    for key, val in averaged.items():
        ci, cj = key
        canon = tuple(sorted(key))
        if canon in seen:
            continue
        seen.add(canon)

        reverse = (cj, ci)
        if reverse in averaged:
            rv = averaged[reverse]
            sym = {
                'chain_i': ci,
                'chain_j': cj,
                'LIS': (val['LIS'] + rv['LIS']) / 2,
                'cLIS': (val['cLIS'] + rv['cLIS']) / 2,
                'ipTM': (val['ipTM'] + rv['ipTM']) / 2,
                'len_i': val['len_i'],
                'len_j': val['len_j'],
                # LIR for chain i: from (i→j) direction
                # LIR for chain j: from (j→i) direction
                'LIR_i': val['LIR_i'],
                'LIR_j': rv['LIR_i'],
                'cLIR_i': val['cLIR_i'],
                'cLIR_j': rv['cLIR_i'],
            }
            sym['iLIS'] = float(np.sqrt(sym['LIS'] * sym['cLIS']))
        else:
            sym = val

        symmetric[(ci, cj)] = sym

    return symmetric


# ── ChimeraX Script Generation ────────────────────────────────────────────

def _res_spec(positions, chain):
    if not positions:
        return ""
    sorted_pos = sorted(positions)
    ranges = []
    start = end = sorted_pos[0]
    for p in sorted_pos[1:]:
        if p == end + 1:
            end = p
        else:
            ranges.append(f"{start}-{end}" if start != end else str(start))
            start = end = p
    ranges.append(f"{start}-{end}" if start != end else str(start))
    return f"/{chain}:" + ",".join(ranges)


def _fill_gaps(positions, max_gap=20):
    if not positions:
        return set()
    sorted_pos = sorted(positions)
    filled = set(sorted_pos)
    for i in range(len(sorted_pos) - 1):
        gap = sorted_pos[i + 1] - sorted_pos[i] - 1
        if 0 < gap <= max_gap:
            for j in range(sorted_pos[i] + 1, sorted_pos[i + 1]):
                filled.add(j)
    return filled


def _ranges_str(positions):
    if not positions:
        return "none"
    sorted_pos = sorted(positions)
    ranges = []
    start = end = sorted_pos[0]
    for p in sorted_pos[1:]:
        if p <= end + 15:
            end = p
        else:
            ranges.append(f"{start}-{end}")
            start = end = p
    ranges.append(f"{start}-{end}")
    return ", ".join(ranges)


def generate_chimerax_script(
    cif_filename, pair_data,
    colors=None,
):
    """Generate .cxc script for an AF3 chain pair."""
    if colors is None:
        colors = {
            'lir_a': '#b3d4e8', 'clir_a': '#2471A3',
            'lir_b': '#f5cba7', 'clir_b': '#E67E22',
        }

    ci = pair_data['chain_i']
    cj = pair_data['chain_j']
    lir_i = pair_data['LIR_i']
    lir_j = pair_data['LIR_j']
    clir_i = pair_data['cLIR_i']
    clir_j = pair_data['cLIR_j']
    iLIS = pair_data['iLIS']
    ipTM = pair_data['ipTM']
    len_i = pair_data['len_i']
    len_j = pair_data['len_j']

    lir_i_str = _ranges_str(lir_i)
    lir_j_str = _ranges_str(lir_j)

    iLIS_s = f"{iLIS:.3f}"
    ipTM_s = f"{ipTM:.2f}"

    # Title label positions
    part_a = f"Chain {ci} ({lir_i_str})"
    part_b = f"Chain {cj} ({lir_j_str})"
    part_c = f"iLIS: {iLIS_s}  ipTM: {ipTM_s}"
    xpos_b = 0.03 + len(part_a) * 0.009 + 0.015
    xpos_c = xpos_b + len(part_b) * 0.009 + 0.015

    # Fill gaps for continuous cartoon
    lir_i_filled = _fill_gaps(lir_i, max_gap=20)
    lir_j_filled = _fill_gaps(lir_j, max_gap=20)
    lir_i_spec = _res_spec(lir_i_filled, ci)
    lir_j_spec = _res_spec(lir_j_filled, cj)
    clir_i_spec = _res_spec(clir_i, ci)
    clir_j_spec = _res_spec(clir_j, cj)

    lines = [
        f"# ChimeraX visualization: Chain {ci} — Chain {cj} (iLIS={iLIS_s}, ipTM={ipTM_s})",
        f"# Chain /{ci}: {len_i} residues, LIR: {lir_i_str}",
        f"# Chain /{cj}: {len_j} residues, LIR: {lir_j_str}",
        f"# Generated by AF3-ChimeraX Agent",
        f"",
        f"close",
        f"open {cif_filename}",
        f"",
        f"# — Setup: white background, hide everything first —",
        f"set bgColor white",
        f"graphics silhouettes true",
        f"lighting soft",
        f"color #1/{ci} {colors['lir_a']}",
        f"color #1/{cj} {colors['lir_b']}",
        f"hide atoms",
        f"hide cartoons",
        f"",
        f"# — Title —",
        f'2dlabels create title text "{part_a}" xpos 0.03 ypos 0.95 size 16 color {colors["clir_a"]} bold true',
        f'2dlabels create title2 text "{part_b}" xpos {xpos_b:.3f} ypos 0.95 size 16 color {colors["clir_b"]} bold true',
        f'2dlabels create title3 text "{part_c}" xpos {xpos_c:.3f} ypos 0.95 size 16 color black bold true',
        f"",
        f"# — Show ONLY LIR regions as cartoons —",
        f"# LIR chain {ci}: {len(lir_i)} residues",
    ]
    if lir_i_spec:
        lines.append(f"show #1{lir_i_spec} cartoons")
    lines.append("")
    lines.append(f"# LIR chain {cj}: {len(lir_j)} residues")
    if lir_j_spec:
        lines.append(f"show #1{lir_j_spec} cartoons")
    lines.append("")

    lines.append(f"# — Highlight cLIR residues (contact interface) —")
    if clir_i:
        lines.append(f"# cLIR chain {ci}: {len(clir_i)} residues")
        lines.append(f"color #1{clir_i_spec} {colors['clir_a']}")
        lines.append("")
    if clir_j:
        lines.append(f"# cLIR chain {cj}: {len(clir_j)} residues")
        lines.append(f"color #1{clir_j_spec} {colors['clir_b']}")
        lines.append("")

    lines.append("# — Final view —")
    lines.append("view")
    lines.append("lighting soft depthCue true")
    lines.append("")
    lines.append(f"# — Save figure —")
    lines.append(f"save chain{ci}_chain{cj}.png transparentBackground true")
    lines.append("")

    lines.append("# — Color Legend —")
    lines.append(f"# LIR  = Local Interaction Region (PAE <= 12)")
    lines.append(f"# cLIR = contact LIR (PAE <= 12 & C-beta distance <= 8 \u00c5)")
    lines.append(f"#")
    lines.append(f"# LIR cartoons  ({colors['lir_a']}) = Chain {ci} LIR")
    lines.append(f"# cLIR cartoons ({colors['clir_a']}) = Chain {ci} cLIR")
    lines.append(f"# LIR cartoons  ({colors['lir_b']}) = Chain {cj} LIR")
    lines.append(f"# cLIR cartoons ({colors['clir_b']}) = Chain {cj} cLIR")
    lines.append(f"# Non-LIR regions are HIDDEN")

    return "\n".join(lines)


# ── Main Pipeline ──────────────────────────────────────────────────────────

def run(path, chain_pair=None, output_dir=None, pae_cutoff=12, distance_cutoff=8):
    """Run the full AF3 → ChimeraX pipeline."""
    path = str(path)
    print(f"[AF3] Processing: {path}")

    models = discover_files(path)
    if not models:
        print("[AF3] ERROR: No AF3 model files found")
        sys.exit(1)

    print(f"[AF3] Found {len(models)} models")

    # Analyze each model
    all_pairs = []
    chain_names = None
    for idx in sorted(models.keys()):
        m = models[idx]
        print(f"[AF3] Analyzing model {idx}...")
        pairs, cn, sizes = analyze_model(
            m['full_data'], m['summary'], m['cif'],
            pae_cutoff, distance_cutoff
        )
        all_pairs.append(pairs)
        if chain_names is None:
            chain_names = cn

    # Average
    averaged = average_models(all_pairs)

    # Print summary
    print(f"\n[AF3] Chain pair summary (averaged across {len(models)} models):")
    print(f"{'Pair':<10} {'iLIS':>8} {'ipTM':>8} {'LIS':>8} {'cLIS':>8} {'LIR(i)':>8} {'LIR(j)':>8} {'cLIR(i)':>8} {'cLIR(j)':>8}")
    print("-" * 82)
    for key, v in sorted(averaged.items()):
        ci, cj = key
        print(f"{ci}-{cj:<8} {v['iLIS']:8.3f} {v['ipTM']:8.2f} {v['LIS']:8.3f} {v['cLIS']:8.3f} {len(v['LIR_i']):8d} {len(v['LIR_j']):8d} {len(v['cLIR_i']):8d} {len(v['cLIR_j']):8d}")

    # Select chain pair
    if chain_pair:
        ci, cj = chain_pair.split(',')
        key = (ci.strip(), cj.strip())
    else:
        # Pick pair with highest iLIS
        key = max(averaged.keys(), key=lambda k: averaged[k]['iLIS'])

    if key not in averaged:
        print(f"[AF3] Chain pair {key} not found")
        sys.exit(1)

    pair_data = averaged[key]
    ci, cj = key
    print(f"\n[AF3] Selected pair: {ci}-{cj} (iLIS={pair_data['iLIS']:.3f}, ipTM={pair_data['ipTM']:.2f})")

    # Output
    if output_dir is None:
        basename = os.path.basename(path).replace('.zip', '')
        output_dir = f"output/{basename}"
    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Use model 0 CIF
    cif_filename = models[0]['cif_filename']

    # Save CIF
    cif_path = Path(output_dir) / cif_filename
    if not cif_path.exists():
        cif_path.write_text(models[0]['cif'])
        print(f"[AF3] CIF saved: {cif_path}")

    # Generate CXC
    script = generate_chimerax_script(cif_filename, pair_data)
    cxc_path = Path(output_dir) / f"chain{ci}_chain{cj}_interface.cxc"
    cxc_path.write_text(script)
    print(f"[AF3] ChimeraX script saved: {cxc_path}")

    return averaged


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="AF3 → ChimeraX Agent")
    parser.add_argument("path", help="Path to AF3 zip file or extracted folder")
    parser.add_argument("--chain-pair", default=None, help="Chain pair to visualize, e.g. 'A,B'")
    parser.add_argument("--output-dir", default=None, help="Output directory")
    parser.add_argument("--pae-cutoff", type=float, default=12, help="PAE cutoff (default: 12)")
    parser.add_argument("--distance-cutoff", type=float, default=8, help="Cβ distance cutoff (default: 8)")
    args = parser.parse_args()

    run(args.path, chain_pair=args.chain_pair, output_dir=args.output_dir,
        pae_cutoff=args.pae_cutoff, distance_cutoff=args.distance_cutoff)
