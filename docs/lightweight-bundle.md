# LIVIA lightweight bundle — authoring guide

**Goal:** turn a full AlphaFold3 (or ColabFold / Boltz / Chai / OpenFold3 …) prediction —
which for a higher-order complex can be tens to hundreds of MB, dominated by the numeric
PAE — into a **~1–3 MB** bundle that `universal.html` renders with **no loss of any figure**.

This works because `universal.html` never displays the numeric PAE. It uses PAE at exactly
one point — to compute the LIS metrics and the interface residue sets — and **`lis.py`
already writes both of those to its CSV**. Hand LIVIA that CSV (all models) + one structure —
or one per model — + a PAE *image*, and every card (3D viewer, chord, contact map, score
matrix, sequence viewer) is identical to a full load. The PAE image is a bonus: the full path
shows no PAE at all.

---

## 1. The bundle

**One bundle carries every model.** All N models' scores + interface residues live in the one
`lis.csv` (AF3 ships 5 models → all 5 are in the bundle), so the score matrix, the per-model
selector, and the cross-model average all work exactly like a full load. You choose only how
many **structures** to include — the CSV is unchanged either way:

**(a) One structure — smallest, the usual choice:**

```
<name>.zip
├── manifest.json     # describes the bundle (recommended; see §6)
├── lis.csv           # lis.py output — ALL models' scores + interface residues
├── model.cif         # ONE structure (.pdb/.cif); every model shares it for 3D/chord/contact
└── pae.png           # PAE plot as a raster image (see §4)
```

**(b) One structure per model — full per-model geometry (still ~10× smaller than the full bundle):**

```
<name>.zip
├── manifest.json          # with a `structures` map (+ optional `pae_images` map) — see §5/§6
├── lis.csv                # ALL models' scores + interface residues
├── model_0.cif.gz         # one structure per model; .gz is optional (LIVIA decompresses).
├── model_1.cif.gz         #   The manifest `structures` map ties each file to its model, so
├── model_2.cif.gz         #   switching models in 3D shows that model's real coordinates.
├── model_3.cif.gz
├── model_4.cif.gz
├── pae_0.png              # optional per-model PAE images (`pae_images` map) …
├── …                      #   … or a single shared pae.png for all models
└── pae_4.png
```

Typical sizes: `lis.csv` tens of KB (holds every model), each structure 1–3 MB, `pae.png` ~50 KB.
The numeric PAE — the thing you are dropping — was the multi-tens-of-MB part. Example (PRC2, 5
models): full **17 MB** → one-structure **0.5 MB** → five-structure **1.7 MB**.

---

## 2. Step 1 — run lis.py to get `lis.csv`

`lis.py` lives at `LIVIA/python/lis.py` (twin: `AFM-LIS/lis.py`). Point it at the full
prediction folder or zip:

```bash
python LIVIA/python/lis.py /path/to/full_prediction/ -o lis.csv
# defaults: --pae-cutoff 12 --cb-cutoff 8  ← DO NOT change these
```

**Keep the default cutoffs (12 and 8).** They are the exact defaults `universal.html` uses
(`pae-cutoff` input = 12, `cb-cutoff` input = 8), so your hosted scores equal what an
in-browser full load would compute. Changing them silently makes your bundle disagree with LIVIA.

The CSV has one row **per chain pair per model** (a 4-chain complex → C(4,2)=6 pairs × 5
models = 30 rows). Its header is exactly:

```
name,rank,model,chain_i,chain_j,iLIS,iLIA,iLISA,ipSAE,actifpTM,LIS,cLIS,LIA,cLIA,ipTM,
pLDDT_i,pLDDT_j,pTM,LIR_i,LIR_j,cLIR_i,cLIR_j,LIpLDDT_i,LIpLDDT_j,cLIpLDDT_i,cLIpLDDT_j,
len_i,len_j,LIR_indices_i,LIR_indices_j,cLIR_indices_i,cLIR_indices_j,structure_file
```

The `*_indices_*` columns are the interface residue sets, written as a quoted, bracketed,
compact range string, e.g. `"[1-4,10-11,13]"` means residues 1,2,3,4,10,11,13. These are
1-based residue **numbers as they appear in the structure** — which is why you must ship the
same structure(s) lis.py read (see §9). **Do not edit or reformat the CSV.**

> Keep **all** model rows in the CSV (it is tiny). LIVIA's default view averages the scores
> across models and unions the residue sets — exactly like a full load. If instead you want the
> display to reflect only one model, delete the other models' rows before zipping.

---

## 3. Step 2 — keep the structure(s)

All five models' **scores** always come from the CSV. The only choice here is how much
**geometry** to ship — a size/fidelity trade-off:

**One structure (smallest — the usual choice).** Copy out a single structure — the top-ranked
model (rank 1; for AF3 `*_model_0.cif`), or the highest-`iLISA` model in `lis.csv`. Every
model's scores still show; the 3D/chord/contact views all use this one geometry. Record which
model it is in the manifest (`structure_model`). You may rename it (`model.cif`).

**One structure per model (larger — full per-model geometry).** Copy out every model's
structure so switching models in the 3D viewer shows that model's real coordinates. Still far
smaller than the full bundle (no numeric PAE — e.g. PRC2: 5 structures ≈ 1.7 MB vs 17 MB full).
Tell LIVIA which file is which model in one of two ways:
- **manifest `structures` map (recommended, robust):** `"structures": {"0":"model_0.cif.gz", …}`.
  This works even if you renamed or gzipped the files, and `structure` points at a representative
  (e.g. model 0).
- **matching filenames:** keep the structures' **original names** so they match the CSV
  `structure_file` column, and LIVIA maps them automatically (no `structures` map needed).

Either way:
- `.pdb` and `.cif` both work, and you may **gzip** them (`.cif.gz` / `.pdb.gz`) to save more
  space — LIVIA decompresses them on load.
- **Do not alter coordinates, chain IDs, or residue numbering** — the CSV's residue indices
  are keyed to them (gzipping is fine; it doesn't change the content).
- Chain IDs must match `chain_i`/`chain_j` in the CSV (they will, since lis.py derived the CSV
  from these structures).

---

## 4. Step 3 — render `pae.png`

The full prediction still has the numeric PAE at authoring time — render it once to a PNG,
then drop the numeric file. Use the **same model** whose structure you kept. Any PAE plot
image works (the `bwr` heatmap below matches LIVIA's own PAE maps — blue = low error):

```python
import json, numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

def load_pae(path):
    """Return an (N,N) float PAE matrix from a prediction's PAE file."""
    if path.endswith('.json'):
        d = json.load(open(path))
        if isinstance(d, dict):
            for k in ('pae', 'predicted_aligned_error'):
                if k in d: return np.asarray(d[k], dtype=np.float32)
        if isinstance(d, list):
            if d and isinstance(d[0], dict) and 'predicted_aligned_error' in d[0]:
                return np.asarray(d[0]['predicted_aligned_error'], dtype=np.float32)
            return np.asarray(d, dtype=np.float32)          # bare matrix
    if path.endswith('.npz'):
        z = np.load(path)
        for k in ('pae', 'predicted_aligned_error', *z.files):
            if k in z: return np.asarray(z[k], dtype=np.float32)
    if path.endswith('.npy'):
        return np.asarray(np.load(path), dtype=np.float32)
    raise ValueError(f'Unrecognized PAE file: {path}')

def render_pae_png(pae, out_png, chain_lengths=None, max_pae=30):
    pae = np.asarray(pae, dtype=np.float32)
    if pae.ndim == 3: pae = pae[0]
    fig, ax = plt.subplots(figsize=(5, 5), dpi=100)
    im = ax.imshow(pae, cmap='bwr', vmin=0, vmax=max_pae,
                   interpolation='nearest', origin='upper')          # blue = low PAE (confident) → red = high; matches LIVIA's PAE maps
    if chain_lengths:                                                # optional: white chain-boundary lines
        b = 0
        for L in chain_lengths[:-1]:
            b += L
            ax.axhline(b - 0.5, color='white', lw=0.8)
            ax.axvline(b - 0.5, color='white', lw=0.8)
    ax.set_xlabel('Scored residue'); ax.set_ylabel('Aligned residue')
    cb = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cb.set_label('Predicted aligned error (Å)')
    fig.tight_layout(); fig.savefig(out_png, dpi=100, bbox_inches='tight'); plt.close(fig)

# --- usage ---
pae = load_pae('/path/to/full_prediction/..._full_data_0.json')   # AF3 example
render_pae_png(pae, 'pae.png', chain_lengths=[300, 150, 90, 120]) # chain_lengths optional
```

`chain_lengths` (residues per chain, in order) just adds boundary lines; you can read them
from the CSV `len_i`/`len_j` columns or count residues per chain in the structure. If your
platform already produced a PAE PNG (e.g. ColabFold `*_pae.png`), just use that.

**Per-model PAE images (optional).** Render one PNG per model (`pae_0.png … pae_4.png`) and list
them in the manifest as `"pae_images": {"0":"pae_0.png", …}`; the PAE plot then switches with the
selected model. Otherwise a single representative `pae_image` is shown for all models.

---

## 5. Step 4 — write `manifest.json`

Recommended — it makes the bundle self-describing and lets LIVIA label chains with gene
names. (LIVIA also **auto-detects** a manifest-less bundle: it finds the CSV by its header
signature, the structure by extension, and the PNG by name. The manifest is the robust path.)

```json
{
  "livia_bundle": "lightweight-v1",
  "name": "CG4679-mRpS9-mRpS29-mRpS31",
  "scores": "lis.csv",
  "structure": "model.pdb",
  "pae_image": "pae.png",
  "structure_model": 0,
  "chains": { "A": "CG4679", "B": "mRpS9", "C": "mRpS29", "D": "mRpS31" },
  "pae_cutoff": 12,
  "cb_cutoff": 8,
  "source_platform": "alphafold3",
  "note": "Numeric PAE dropped to reduce size; scores precomputed by lis.py."
}
```

---

## 6. manifest.json fields

| Field | Required | Meaning |
|---|---|---|
| `livia_bundle` | **yes** | Format tag. Use the literal `"lightweight-v1"`. This is how LIVIA positively identifies the bundle. |
| `name` | recommended | Display label (used as the title; equivalent to `&name=` on the URL). |
| `scores` | recommended | Filename of the lis.py CSV inside the zip (default guess: any `*.csv` with the signature header). |
| `structure` | recommended | Filename of the representative `.pdb`/`.cif[.gz]` — used for the chain layout and as the fallback geometry. |
| `structures` | optional | Map of `model` → structure filename, e.g. `{"0":"model_0.cif.gz", …}`. **The robust way to give per-model geometry** — LIVIA maps each model to its own file by this (preferred over the CSV `structure_file` column, which may not match if you renamed/gzipped the files). |
| `pae_image` | optional | Filename of the representative PAE image. Omit if you shipped none — the PAE card just hides. |
| `pae_images` | optional | Map of `model` → PAE image filename, e.g. `{"0":"pae_0.png", …}`. When present, the PAE plot switches with the selected model. |
| `structure_model` | optional | Which model the representative `structure`/`pae_image` correspond to (the `model` value from the CSV). |
| `chains` | optional | Map of chain ID → gene/protein label, shown in the figures. |
| `pae_cutoff`, `cb_cutoff` | optional | Provenance. Should be `12` and `8` unless you deliberately changed them. |
| `source_platform` | optional | Provenance (`alphafold3`, `colabfold`, …). |
| `note` | optional | Free text. |

---

## 7. How LIVIA loads it

Host the zip anywhere with CORS (GitHub Pages, OSF, S3 with CORS, …) and open:

```
https://flyark.github.io/LIVIA/universal.html?data=<url-encoded-zip-url>&name=<label>
```

`name` on the URL is optional if the manifest has one. GitHub Pages and AFDB send CORS
headers and load directly; hosts that don't (e.g. OSF) go through LIVIA's allow-listed proxy
automatically — if you use a **new** host, it must be added to the proxy allow-list
(`LIVIA/tools/cloudflare-worker/worker.js`) or the fetch is refused.

---

## 8. Validation checklist (run before hosting)

- [ ] `lis.csv` first line is the exact header in §2 (the `cLIR_indices_i` column is present).
- [ ] `lis.csv` keeps **all** models' rows (the `model` column spans every model, e.g. 0–4);
      `chain_i`/`chain_j` cover every pair of the complex.
- [ ] The structure parses and its chain IDs match the CSV's `chain_i`/`chain_j`.
- [ ] The structure's residue numbering is **unchanged** from what lis.py read.
- [ ] *(per-model structures only)* every shipped `.cif/.pdb` filename matches a `structure_file`
      value in the CSV, so LIVIA maps each model to its own geometry.
- [ ] `pae.png` opens as an image (if you shipped one).
- [ ] `manifest.json` is valid JSON with `"livia_bundle": "lightweight-v1"` and correct filenames.
- [ ] Total zip is small (target a few MB — if it's tens of MB, a numeric PAE probably slipped in).
- [ ] Smoke test: serve the zip on `http://localhost:8000/` and open
      `universal.html?data=http://localhost:8000/<name>.zip` — all cards render, the model
      selector lists all N models, scores match a full load, and the PAE image shows.

---

## 9. Fidelity notes & gotchas

- **Lossless for every figure.** Scores and interface residues come from the CSV; geometry
  (3D, chord, contact-map lines within 8 Å, sequence highlights) comes from the one structure.
  Nothing LIVIA renders needs the numeric PAE.
- **Geometry follows how many structures you ship.** All models' *scores* always appear. Ship
  **one** structure and every model's 3D/chord/contact uses it (the "Average" view unions the
  interface residues across models — same as a full load). Ship **one per model** and each
  model's views use its own coordinates.
- **Trust model.** Scores are taken from the CSV, not recomputed — the bundle is only as
  correct as the CSV. Generate it with unmodified `lis.py` and default cutoffs.
- **Don't hand-edit the CSV** (especially the residue-range strings) — the loader parses them
  strictly.
- **Residue numbering must match** between CSV and structure. Ship the *same* structure lis.py
  processed; don't renumber or re-relax it.

---

## 10. Ready-made prompt for another session

Fill in the placeholders. The prompt keeps **all models' scores** either way; the bracketed
clause in step 2 is where you pick single-structure (default) vs. one-structure-per-model.

> You have a full <PLATFORM> prediction of a multi-chain complex at `<PATH>` (it has <N> models —
> AF3 = 5). Produce a **LIVIA lightweight bundle** following `LIVIA/docs/lightweight-bundle.md`:
>
> 1. Run `LIVIA/python/lis.py <PATH> -o lis.csv` with default cutoffs. The CSV **must keep all
>    <N> models' rows** — do not filter it (that's what makes the bundle multi-model).
> 2. Copy out the structure(s), unaltered:
>    - **default (smallest):** just the top-ranked model as `model.cif`; **or**
>    - **per-model geometry:** **all <N>** model structures, keeping their **original filenames**
>      so they match the CSV `structure_file` column.
> 3. Render one model's numeric PAE to `pae.png` with the §4 snippet.
> 4. Write `manifest.json` (§5) with a `chains` gene-label map (point `structure` at model 0).
> 5. Zip everything as `<name>.zip` and run the §8 checklist.
>
> Never put a numeric PAE file (`*_full_data_*.json`, `*.npz`, `*.npy`) in the zip. Report the
> final zip path, its size, and how many models/structures it contains.
