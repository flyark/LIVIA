#!/usr/bin/env python3
"""Canonical LIVIA lightweight-v1 bundle format — manifest builder + validator.

The single place the `livia_bundle: "lightweight-v1"` contract lives in code, shared by both
generators so the format cannot drift between them:

  * LIVIA        python/make_lightweight_bundle.py   (one prediction -> one bundle)
  * lis-toolkit  scripts/screen_site_bundles.py      (a whole screen -> one bundle per pair)

Kept in sync as a twin across those repos (the same convention as lis.py). Spec:
LIVIA/docs/lightweight-bundle.md. Loader: LIVIA/universal.html
(detectLightweightBundle -> processPrecomputedBundle).

Pure stdlib (json / zipfile / csv / io) so it is safe to vendor anywhere.
"""
from __future__ import annotations

import csv as _csv
import io
import json
import zipfile

TAG = "lightweight-v1"
SCORES_NAME = "lis.csv"
MANIFEST_NAME = "manifest.json"
CSV_SIGNATURE = "cLIR_indices_i"   # a column only a lis.py CSV carries — how LIVIA recognises the scores


def build_manifest(*, name, chains, structures, best_model,
                   pae_cutoff, cb_cutoff, source_platform,
                   pae_images=None, structure=None, scores=SCORES_NAME, **extra):
    """Build a lightweight-v1 manifest dict.

    structures : {model_key(str): filename} — one structure file per model.
    best_model : the representative model key; picks the singular `structure` / `pae_image`
                 a single-model reader falls back to (may be None).
    pae_images : {model_key: filename}  — per-model PAE images (LIVIA renders these), or
                 a single filename str  — one image shared across models (ColabFold all-ranks), or
                 None.
    structure  : override the singular structure filename (defaults to structures[best_model]).
    **extra    : any additional keys stored verbatim — e.g. label, n_models, numeric_pae, note.
    """
    structures = {str(k): v for k, v in structures.items()}
    bm = str(best_model)
    man = {
        "livia_bundle": TAG,
        "name": name,
        "scores": scores,
        "chains": dict(chains),
        "structures": structures,
        "structure": structure or structures.get(bm) or (next(iter(structures.values())) if structures else None),
        "structure_model": best_model,
        "pae_cutoff": pae_cutoff,
        "cb_cutoff": cb_cutoff,
        "source_platform": source_platform,
    }
    if isinstance(pae_images, dict):
        pm = {str(k): v for k, v in pae_images.items()}
        man["pae_images"] = pm
        man["pae_image"] = pm.get(bm) or (next(iter(pm.values())) if pm else None)
    elif pae_images:
        man["pae_image"] = pae_images        # single image shared across models
    man.update(extra)
    return man


def validate_bundle(zip_path):
    """Return a list of problems for a lightweight-v1 zip (empty list = valid).

    A superset of both generators' checks, each gated on the keys actually present, so it
    validates single-model (screen-site) and multi-model (LIVIA) bundles alike.
    """
    problems = []
    try:
        z = zipfile.ZipFile(zip_path)
    except Exception as e:                                    # noqa: BLE001
        return [f"cannot open zip: {e!r}"]
    names = set(z.namelist())
    if z.testzip() is not None:
        problems.append("corrupt archive")
    if MANIFEST_NAME not in names:
        return problems + [f"no {MANIFEST_NAME}"]
    try:
        man = json.loads(z.read(MANIFEST_NAME))
    except Exception as e:                                    # noqa: BLE001
        return problems + [f"unreadable {MANIFEST_NAME}: {e!r}"]

    if man.get("livia_bundle") != TAG:
        problems.append(f"bad livia_bundle tag: {man.get('livia_bundle')!r}")

    structs = man.get("structures") or {}
    for f in structs.values():
        if f not in names:
            problems.append(f"structures map broken: {f}")
        if str(f).endswith(".gz"):
            problems.append(f"{f} is gzipped — store the structure plain")
    for f in (man.get("pae_images") or {}).values():
        if f not in names:
            problems.append(f"missing PAE image {f}")
    for k in ("structure", "pae_image"):
        if man.get(k) and man[k] not in names:
            problems.append(f"{k} points at absent {man[k]}")
    if man.get("n_models") is not None and structs and man["n_models"] != len(structs):
        problems.append(f"n_models={man['n_models']} but {len(structs)} structures present")

    scores = man.get("scores", SCORES_NAME)
    if scores not in names:
        problems.append(f"missing scores {scores}")
    else:
        rows = list(_csv.reader(io.StringIO(z.read(scores).decode("utf-8", "replace"))))
        if not rows or CSV_SIGNATURE not in rows[0]:
            problems.append(f"{scores} lacks the lis.py signature column ({CSV_SIGNATURE})")
        if len(rows) < 2:
            problems.append(f"{scores} has no data rows")
        else:
            if man.get("label"):
                nm = dict(zip(rows[0], rows[1])).get("name", "")
                if nm and nm != man["label"]:
                    problems.append(f'{scores} name "{nm}" does not match label "{man["label"]}"')
            if "structure_file" in rows[0]:
                i = rows[0].index("structure_file")
                refs = {r[i] for r in rows[1:] if len(r) > i and r[i]}
                miss = {f for f in refs if f not in names}
                if miss:
                    problems.append(f"{scores} references structures not in the bundle: {sorted(miss)}")
    return problems
