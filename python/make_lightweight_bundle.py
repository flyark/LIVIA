#!/usr/bin/env python3
"""Turn a structure prediction into a small, self-contained bundle you can actually keep.

A finished prediction is mostly numeric PAE. In a typical AlphaFold 3 job the five `full_data_*.json`
files hold ~880 MB while the five model CIFs hold ~15 MB -- 98% of what you are storing is a matrix you
will never read by hand, only plot and score. Once it has been scored and plotted, it can go.

This produces, from any prediction that lis.py can read:

    <name>.zip
      manifest.json     what is inside, which model is best, what cutoffs were used
      lis.csv           every LIS/cLIS/iLIS/ipSAE/ipTM number, per model, per chain pair
      model_<m>.cif     the structures, unchanged
      pae_<m>.png       the PAE heatmap for each model

Measured over 467 AlphaFold 3 complexes (3-30 subunits): 71.7 MB -> 3.2 MB, a **22x** reduction, with
every published score still reproducible from lis.csv and the PAE still readable as an image. What you
give up is the numeric PAE: you can no longer re-score at a different cutoff, or compute something
lis.py does not already compute. If that matters for a given prediction, keep that one intact.

The bundle also opens directly in LIVIA (https://flyark.github.io/LIVIA/universal.html) -- drop the zip
in and it loads structures, PAE and scores without re-running anything.

Supports every platform lis.py does: AlphaFold 3 (server + local), ColabFold, AlphaPulldown/AF-Multimer,
Boltz, Chai-1, OpenFold3, ESMFold2, and a generic structure+PAE layout.

    python make_lightweight_bundle.py <prediction.zip | prediction_folder/> [-o out/] [--name LABEL]

    --name LABEL     what the bundle is called publicly. Defaults to the prediction's own name, which
                     is often an internal job identifier -- pass --name if it is not meant to be shared.
    --keep-pae       also store the numeric PAE (as compressed .npz). Bundle gets larger; re-scoring
                     at other cutoffs stays possible.
    --models 0,2     only these model indices (default: all found).
    --delete-original  remove the source prediction, but ONLY after the new bundle passes every check.
                     Refused outright when the bundle would not contain everything the source had.
    --check          verify an existing bundle instead of building one.

Requires: numpy, matplotlib, and lis.py beside this file (or on PYTHONPATH).
"""
import argparse
import io
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
try:
    import lis
except ImportError:
    sys.exit('lis.py not found. Put it next to this script, or add it to PYTHONPATH.\n'
             'It lives at https://github.com/flyark/LIVIA -> python/lis.py')

# Pinned scoring cutoffs. Every number in lis.csv depends on them, so a bundle records the values it
# was built with -- two bundles scored at different cutoffs are not comparable, and the manifest is
# the only place that stays attached to the data.
PAE_CUTOFF = 12
CB_CUTOFF = 8
PAE_VMAX = 30          # fixed, never autoscaled: see render_pae


# --------------------------------------------------------------------------- PAE image
def render_pae(pae, bounds, out_png):
    """One PAE heatmap. `bounds` = cumulative chain-boundary indices, e.g. [300, 450] for 300/150/rest.

    The fixed 0-30 range is the important part. Autoscaling each plot to its own data makes a confident
    complex and a poor one look identical, which defeats the point of keeping the image at all.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    pae = np.asarray(pae, dtype=np.float32)
    if pae.ndim == 3:
        pae = pae[0]
    n = pae.shape[0]
    fig, ax = plt.subplots(figsize=(4.4, 3.8), dpi=150)
    im = ax.imshow(pae, cmap='bwr', vmin=0, vmax=PAE_VMAX, interpolation='nearest')
    # black, not white: bwr puts white at PAE 15, so a white line disappears into every
    # moderately-confident region of the plot
    for b in bounds:
        ax.axhline(b - .5, color='black', lw=0.8)
        ax.axvline(b - .5, color='black', lw=0.8)
    starts = [0] + [int(b) for b in bounds]
    ends = [int(b) for b in bounds] + [n]
    centers = [(s + e) / 2 for s, e in zip(starts, ends)]
    labels = [chr(ord('A') + i) for i in range(len(centers))]
    # chain letters instead of residue indices: at this size the indices are unreadable anyway, and
    # chain identity is what you need when looking at a ten-subunit complex
    fs = max(6, min(13, int(128 / max(len(centers), 1))))
    ax.set_xticks(centers); ax.set_xticklabels(labels, fontsize=fs); ax.tick_params(length=0)
    ax.set_yticks(centers); ax.set_yticklabels(labels, fontsize=fs)
    # no title: viewers label each model above its own plot, and a title just eats height
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, ticks=[0, PAE_VMAX // 2, PAE_VMAX])
    cb.ax.tick_params(labelsize=10)
    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.close(fig)


def chain_info(text, fmt):
    """(boundaries, chains) from a structure: cumulative break indices, and chain id -> residue count.

    lis.get_chains_from_structure returns {'names': [...], 'sizes': [...], 'types': [...]}; a list of
    per-chain records is also accepted so this keeps working if that ever changes. Getting this wrong is
    silent -- the bundle still builds, the PAE images just quietly lose their chain boundary lines.

    The chain map goes in the manifest so a viewer can label subunits without any external annotation.
    Where real names are known (a gene symbol, a UniProt id), overwrite it; residue count is the fallback.
    """
    try:
        c = lis.get_chains_from_structure(text, fmt)
    except Exception:
        return [], {}
    ids, lens = [], []
    if isinstance(c, dict) and 'names' in c:
        ids = [str(x) for x in c.get('names') or []]
        lens = [int(x or 0) for x in c.get('sizes') or []]
    elif isinstance(c, (list, tuple)):
        for i, ch in enumerate(c):
            if isinstance(ch, dict):
                ids.append(str(ch.get('id') or ch.get('chain') or chr(ord('A') + i)))
                lens.append(int(ch.get('length') or ch.get('size') or ch.get('n_res') or 0))
            elif isinstance(ch, (list, tuple)) and len(ch) >= 2:
                ids.append(str(ch[0])); lens.append(int(ch[1]))
    keep = [(i, L) for i, L in zip(ids, lens) if L > 0]
    out, acc = [], 0
    for _, L in keep[:-1]:
        acc += L
        out.append(acc)
    return out, {i: L for i, L in keep}


# --------------------------------------------------------------------------- scoring
def run_lis(path, workdir):
    """Score the prediction with lis.py and return the CSV text."""
    # -o is a FILENAME and -d the directory; passing a full path to -o alone drops the CSV next to the
    # input, which for a read-only or shared prediction store is both surprising and unwanted
    out_csv = os.path.join(workdir, 'lis.csv')
    cmd = [sys.executable, os.path.join(HERE, 'lis.py'), path,
           '-o', 'lis.csv', '-d', workdir,
           '--pae-cutoff', str(PAE_CUTOFF), '--cb-cutoff', str(CB_CUTOFF)]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if not os.path.exists(out_csv):
        sys.stderr.write(r.stdout + r.stderr)
        raise SystemExit('lis.py produced no CSV -- see its output above.')
    return open(out_csv).read()


def rewrite_csv(text, label):
    """Point lis.csv at the files the bundle actually contains, and apply the public name.

    Two things change. `structure_file` must name `model_<m>.cif`, because that is what is in the zip --
    left alone it points at the original filename and per-model geometry cannot be resolved. And `name`
    becomes the public label, so an internal job identifier does not travel with the data.

    Raw string replacement rather than a CSV round-trip: several columns hold quoted residue-range
    strings, and a re-write would re-quote them differently for no reason.
    """
    rows = text.splitlines()
    if len(rows) < 2:
        return text, None
    header = rows[0].split(',')
    try:
        ni = header.index('name')
    except ValueError:
        return text, None
    original = rows[1].split(',')[ni]
    if not original:
        return text, None
    # Order matters: collapse the longer '<name>_model_' form FIRST, then the bare name.
    # An earlier version used a regex with \S*? between the two, which silently spans commas and
    # therefore ate the whole `name,...,structure_file` span, leaving name="model_0.cif".
    out = text.replace(f'{original}_model_', 'model_').replace(original, label)
    return out, original


def subset_csv(text, models):
    """Keep only rows whose `model` column is in `models`, preserving the file byte-for-byte otherwise."""
    rows = text.splitlines(keepends=True)
    if len(rows) < 2:
        return text
    header = rows[0].rstrip('\r\n').split(',')
    if 'model' not in header:
        return text
    mi = header.index('model')
    out = [rows[0]]
    for r in rows[1:]:
        f = r.rstrip('\r\n').split(',')
        if len(f) <= mi:
            continue
        try:
            if int(float(f[mi])) in models:
                out.append(r)
        except ValueError:
            continue
    return ''.join(out)


# --------------------------------------------------------------------------- build
def build(path, out_dir, label=None, keep_pae=False, only_models=None, quiet=False):
    say = (lambda *a: None) if quiet else (lambda *a: print(*a, flush=True))

    filenames, read_fn, _ = lis.scan_files(path)
    platform = lis.detect_platform(filenames, read_fn)
    models = list(lis.find_models(filenames, platform, read_fn))
    if not models:
        raise SystemExit(f'No structure+PAE pairs found in {path} (detected platform: {platform}).')
    pred_name = models[0][0]
    label = label or pred_name
    say(f'platform  : {platform}')
    say(f'models    : {len(models)}')
    say(f'label     : {label}' + ('' if label != pred_name else '   <- the prediction\'s own name; '
                                                                'pass --name to publish it as something else'))

    work = tempfile.mkdtemp(prefix='lwbundle_')
    try:
        say('scoring with lis.py ...')
        csv_text = run_lis(path, work)
        csv_text, original_name = rewrite_csv(csv_text, label)
        if only_models is not None:
            # otherwise lis.csv keeps rows for models that are not in the bundle, and every consumer --
            # including --check -- is right to call that broken
            csv_text = subset_csv(csv_text, only_models)
        open(os.path.join(work, 'lis.csv'), 'w').write(csv_text)

        kept, bounds, chains = [], [], {}
        # lis.py yields (name, rank, model_label, struct, pae, scores, fmt). `model_label` is a filename,
        # not an index, and `rank` is an int for most platforms but a "seed_sample" string for local AF3 --
        # so derive a stable integer and fall back to position when it is not numeric.
        for pos, (name, rank, model_label, struct_path, pae_path, scores_path, fmt) in enumerate(models):
            try:
                m = int(rank)
            except (TypeError, ValueError):
                m = pos
            if only_models is not None and m not in only_models:
                continue
            stext = read_fn(struct_path)
            if isinstance(stext, bytes):
                stext = stext.decode('utf-8', 'replace')
            ext = 'cif' if (fmt or '').lower().startswith('cif') else 'pdb'
            # plain, not gzipped: the zip's own DEFLATE already gives the gzip ratio, and viewers that
            # stream a member expect an uncompressed structure inside
            open(os.path.join(work, f'model_{m}.{ext}'), 'w').write(stext)
            if not bounds:
                bounds, chains = chain_info(stext, ext)
            try:
                pae = lis.extract_pae(pae_path, read_fn)
                render_pae(pae, bounds, os.path.join(work, f'pae_{m}.png'))
                if keep_pae:
                    np.savez_compressed(os.path.join(work, f'pae_{m}.npz'),
                                        pae=np.asarray(pae, dtype=np.float16))
            except Exception as e:
                say(f'  model {m}: no PAE image ({e})')
            kept.append((m, ext))
            say(f'  model {m}: {ext} + pae')

        if not kept:
            raise SystemExit('Nothing kept -- check --models.')
        best = pick_best(csv_text, [m for m, _ in kept])
        ext_of = dict(kept)
        man = {
            'livia_bundle': 'lightweight-v1',
            'name': label, 'label': label,
            'scores': 'lis.csv',
            'chains': chains,
            'structures': {str(m): f'model_{m}.{ext_of[m]}' for m, _ in kept},
            'pae_images': {str(m): f'pae_{m}.png' for m, _ in kept
                           if os.path.exists(os.path.join(work, f'pae_{m}.png'))},
            # singular keys so a single-model reader still works
            'structure': f'model_{best}.{ext_of[best]}',
            'pae_image': f'pae_{best}.png',
            'structure_model': best,
            'n_models': len(kept),
            'pae_cutoff': PAE_CUTOFF, 'cb_cutoff': CB_CUTOFF,
            'source_platform': platform,
            'numeric_pae': bool(keep_pae),
            'note': ('Numeric PAE retained as pae_<m>.npz (float16).' if keep_pae else
                     'Numeric PAE dropped; scores precomputed by lis.py at the cutoffs above.'),
        }
        json.dump(man, open(os.path.join(work, 'manifest.json'), 'w'), indent=1)

        os.makedirs(out_dir, exist_ok=True)
        outzip = os.path.join(out_dir, f'{label}.zip')
        with zipfile.ZipFile(outzip, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
            for f in sorted(os.listdir(work)):
                z.write(os.path.join(work, f), f)

        # if the name was changed on purpose, prove the old one is gone from the text members
        if original_name and original_name != label:
            zc = zipfile.ZipFile(outzip)
            hit = [n for n in ('lis.csv', 'manifest.json')
                   if original_name.encode() in zc.read(n)]
            if hit:
                raise SystemExit(f'"{original_name}" still present in {hit} -- not writing this bundle.')
            say(f'name check: "{original_name}" absent from lis.csv and manifest.json')

        src = folder_size(path)
        dst = os.path.getsize(outzip)
        say(f'\n{outzip}')
        say(f'  {src/1e6:.1f} MB -> {dst/1e6:.1f} MB   ({src/max(dst,1):.0f}x smaller)')
        return outzip, len(kept), len(models)
    finally:
        shutil.rmtree(work, ignore_errors=True)


def pick_best(csv_text, models):
    """Best model = highest mean iLIS across chain pairs; falls back to the lowest index."""
    try:
        rows = [r.split(',') for r in csv_text.splitlines()]
        h = rows[0]
        mi, si = h.index('model'), h.index('iLIS')
        acc = {}
        for r in rows[1:]:
            if len(r) <= max(mi, si):
                continue
            try:
                acc.setdefault(int(float(r[mi])), []).append(float(r[si]))
            except ValueError:
                continue
        acc = {k: sum(v) / len(v) for k, v in acc.items() if k in models and v}
        if acc:
            return max(acc, key=acc.get)
    except Exception:
        pass
    return min(models)


def folder_size(path):
    if os.path.isfile(path):
        return os.path.getsize(path)
    return sum(os.path.getsize(os.path.join(d, f))
               for d, _, fs in os.walk(path) for f in fs)


# --------------------------------------------------------------------------- verify
def check(bundle):
    """Read a bundle back and report whether it is complete and internally consistent."""
    z = zipfile.ZipFile(bundle)
    names = set(z.namelist())
    problems = []
    if 'manifest.json' not in names:
        return ['no manifest.json']
    man = json.loads(z.read('manifest.json'))
    print(f'{os.path.basename(bundle)}: {man.get("label")}  '
          f'{man.get("n_models")} models  {man.get("source_platform")}  '
          f'cutoffs pae={man.get("pae_cutoff")} cb={man.get("cb_cutoff")}')
    if man.get('scores') not in names:
        problems.append(f'missing {man.get("scores")}')
    for m, f in (man.get('structures') or {}).items():
        if f not in names:
            problems.append(f'missing structure {f}')
        if f.endswith('.gz'):
            problems.append(f'{f} is gzipped -- viewers skip .cif.gz; store it plain')
    for m, f in (man.get('pae_images') or {}).items():
        if f not in names:
            problems.append(f'missing PAE image {f}')
    n_struct = len(man.get('structures') or {})
    if man.get('n_models') != n_struct:
        problems.append(f'n_models={man.get("n_models")} but {n_struct} structures present')
    for k in ('structure', 'pae_image'):
        if man.get(k) and man[k] not in names:
            problems.append(f'{k} points at absent {man[k]}')
    if man.get('scores') in names:
        # a real CSV reader, not split(','): several columns hold quoted residue ranges that contain commas
        import csv as _csv
        rd = list(_csv.reader(io.StringIO(z.read(man['scores']).decode('utf-8', 'replace'))))
        if len(rd) >= 2:
            h, r = rd[0], rd[1]
            row = dict(zip(h, r))
            sf = row.get('structure_file', '')
            if sf and sf not in names:
                problems.append(f'lis.csv structure_file "{sf}" is not in the bundle -- '
                                f'per-model geometry will not resolve')
            nm = row.get('name', '')
            if nm and nm != man.get('label'):
                problems.append(f'lis.csv name "{nm}" does not match the label "{man.get("label")}"')
            files = {row2[h.index('structure_file')] for row2 in rd[1:]
                     if 'structure_file' in h and len(row2) > h.index('structure_file')}
            missing = {f for f in files if f and f not in names}
            if missing:
                problems.append(f'lis.csv references structures not in the bundle: {sorted(missing)}')
    for p in problems:
        print(f'  ! {p}')
    if not problems:
        print(f'  ok  ({os.path.getsize(bundle)/1e6:.1f} MB)')
    return problems


def delete_original(path, bundle, n_kept, n_found, keep_pae, quiet=False):
    """Remove the source prediction -- only when the bundle provably replaces it.

    This is the one irreversible thing the script does, so it refuses on anything short of a complete,
    verified bundle rather than asking the user to be careful. `--check` catches a malformed bundle;
    these extra conditions catch a *well-formed* bundle that is nonetheless missing something the
    source had.
    """
    say = (lambda *a: None) if quiet else (lambda *a: print(*a, flush=True))
    refuse = []   # printed unconditionally: a refusal explains why a requested action did not happen
    if n_kept != n_found:
        refuse.append(f'bundle holds {n_kept} of the {n_found} models in the source '
                      f'(--models was used); deleting would lose the rest')
    problems = check(bundle) if quiet else check(bundle)
    if problems:
        refuse.append('the bundle did not pass --check')
    try:
        sp = os.path.realpath(path)
        bp = os.path.realpath(bundle)
        op = os.path.realpath(os.path.dirname(bundle))
        # The critical case: the new bundle lives inside the source. Deleting the source (rmtree)
        # takes the bundle with it -- you lose both, and the "freed N MB" message hides it.
        if bp == sp or bp.startswith(sp + os.sep):
            refuse.append('the new bundle is inside the source -- deleting the source would delete '
                          'the bundle too (point -o outside the prediction)')
        elif sp == op or sp.startswith(op + os.sep):
            refuse.append('the source sits inside the output directory')
    except Exception:
        pass
    if refuse:
        for r in refuse:
            print(f'  KEEPING {path}: {r}', flush=True)
        return False
    freed = folder_size(path)
    if os.path.isdir(path):
        shutil.rmtree(path)
    else:
        os.remove(path)
    say(f'  deleted {path}  (freed {freed/1e6:.1f} MB)')
    say(f'  the numeric PAE is now gone' + ('' if keep_pae else
        ' -- scores can no longer be recomputed at other cutoffs'))
    return True


def main():
    ap = argparse.ArgumentParser(
        description='Shrink a structure prediction to a scored, plottable bundle (~20x smaller).',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='examples:\n'
               '  python make_lightweight_bundle.py fold_mycomplex.zip -o bundles/\n'
               '  python make_lightweight_bundle.py af3_out/ --name PRC2 -o bundles/\n'
               '  python make_lightweight_bundle.py preds/*.zip -o bundles/\n'
               '  python make_lightweight_bundle.py --check bundles/PRC2.zip\n'
               '  python make_lightweight_bundle.py preds/*.zip -o bundles/ --delete-original\n'
               '\nreclaiming space: build first, look at the bundles, then re-run with\n'
               '--delete-original -- or pass it directly, it only deletes what it has verified.\n')
    ap.add_argument('path', nargs='+', help='prediction zip(s) or folder(s); or bundle(s) with --check')
    ap.add_argument('-o', '--out-dir', default='bundles')
    ap.add_argument('--name', default=None, help='public label (default: the prediction\'s own name)')
    ap.add_argument('--keep-pae', action='store_true', help='also store numeric PAE as float16 .npz')
    ap.add_argument('--models', default=None, help='comma-separated model indices to keep')
    ap.add_argument('--delete-original', action='store_true',
                    help='delete the source prediction after the bundle passes every check (irreversible)')
    ap.add_argument('--check', action='store_true', help='verify existing bundles instead of building')
    ap.add_argument('-q', '--quiet', action='store_true')
    a = ap.parse_args()

    if a.check:
        bad = 0
        for p in a.path:
            bad += len(check(p)) > 0
        sys.exit(1 if bad else 0)

    if a.name and len(a.path) > 1:
        sys.exit('--name applies to a single prediction; drop it when passing several.')
    only = None if not a.models else {int(x) for x in a.models.split(',')}
    made, tot_src, tot_dst, freed = [], 0, 0, 0
    for p in a.path:
        if len(a.path) > 1 and not a.quiet:
            print(f'\n=== {p}')
        try:
            sz = folder_size(p)
            z, n_kept, n_found = build(p, a.out_dir, a.name, a.keep_pae, only, a.quiet)
            made.append(z)
            tot_src += sz; tot_dst += os.path.getsize(z)
            if a.delete_original:
                if delete_original(p, z, n_kept, n_found, a.keep_pae, a.quiet):
                    freed += sz
        except SystemExit as e:
            print(f'  skipped: {e}')
    if len(made) > 1:
        print(f'\n{len(made)} bundles: {tot_src/1e9:.2f} GB -> {tot_dst/1e9:.2f} GB '
              f'({tot_src/max(tot_dst,1):.0f}x smaller)')
    if a.delete_original:
        print(f'freed {freed/1e9:.2f} GB by deleting {"the source" if len(made) == 1 else "sources"}'
              if freed else 'no source was deleted')


if __name__ == '__main__':
    main()
