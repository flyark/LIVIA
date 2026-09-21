# Hosting bundles on GitHub Pages

**Goal:** publish a catalog of lightweight bundles (see
[`lightweight-bundle.md`](lightweight-bundle.md)) or cLIP CSV+FASTA bundles as plain files on
GitHub Pages, and deep-link straight into `universal.html` or `clip.html` for each one.

This is the simple, default hosting path — one file per bundle, served over a static site with
no packing and no special extraction code. Reach for the
[Zenodo packed-bundle guide](hosting-zenodo-packed-bundles.md) instead only when you specifically
need Zenodo's permanent DOI and are willing to deal with its 100-files-per-record quota.

---

## 1. Why GitHub Pages works fine at scale

GitHub Pages has no meaningful per-repo file-count limit for this use case. Two production
catalogs prove it:

- **[`FlyPredictome-clip-data`](https://github.com/flyark/FlyPredictome-clip-data)** — 5,059
  per-gene cLIP bundles (`<FBgn>.zip`, one per FlyPredictome network node) plus an
  `available.json` index. 5,187 files total.
- **[`Visco-ID-interactome-screen-structures`](https://github.com/flyark/Visco-ID-interactome-screen-structures)**
  — 115 lightweight structure bundles, 187 MB, still served this way for its cLIP data even
  after the structure bundles themselves moved to Zenodo (see that guide for why).

Git's own soft limits (a single file capped around 100 MB, a repo getting unwieldy in the
low-GB range) are the real ceiling, not GitHub Pages itself. A lightweight bundle is typically
1–3 MB (see `lightweight-bundle.md` §1), so thousands of them fit comfortably.

## 2. Layout

One file per bundle, named by a stable identifier — not an internal job name (see
`lightweight-bundle.md` §11.1 for why). Two conventions already in use:

**By gene/protein ID** (`FlyPredictome-clip-data`):

```
FBgn0000008.zip
FBgn0000014.zip
...
available.json          # { "FBgn0000008": "symbol", ... } — what exists, and its display name
```

**By subdomain, no index file** (`Visco-ID-interactome-screen-structures/bundles/<screen>/`):

```
bundles/TOMM20/TOM20_HUMAN___ACPM_HUMAN.zip
bundles/TOMM20/TOM20_HUMAN___ARFG1_HUMAN.zip
...
```

Here the catalog lives in the *consuming* repo instead (`Visco-ID-interactome-screen`'s
`data/interactions.csv` has a `bundle` column naming the file) — an index file is convenient
when a page needs to know what's available without fetching a second repo's file listing, but
isn't required if the catalog already lives elsewhere.

## 3. Deep-link URL patterns

**Lightweight structure bundle → `universal.html`:**

```
https://flyark.github.io/LIVIA/universal.html?data=<bundle-url>&name=<label>
```

`universal.html` fetches `data` directly (same-origin or CORS-enabled — GitHub Pages sends
permissive CORS headers by default), falling back to a Cloudflare Worker proxy
(`livia-proxy.flyark.workers.dev`) if the direct fetch fails. See `universal.html`'s
`?data=` handler (search for `URL deep-link` in the file) — no `entry=` parameter here, that's
only for the Zenodo packed-catalog variant.

**cLIP CSV+FASTA bundle → `clip.html`:**

```
https://flyark.github.io/LIVIA/clip.html?data=<bundle-url>&gene=<symbol>&species=<taxid>&ilis=<cutoff>
```

- `gene` picks which protein in the bundle is the bait (falls back to whichever protein
  appears in the most rows if omitted or not found).
- `species` pins the organism for symbol resolution — needed when a symbol is ambiguous across
  species (e.g. a human ortholog vs. the fly gene). Inferred as fly (`7227`) automatically if
  the URL contains an `FBgn` ID or `flypredictome`, otherwise explicit.
- `ilis` pre-sets the iLIS cutoff filter.

Handled by `clip.html`'s `dataParam()` function (search `Direct data hand-off` in the file) — a
plain `fetch()` of the whole bundle, no partial extraction. Fine here because cLIP bundles
(CSV + FASTA) are small; a `TOMM20.zip` cLIP bundle is ~1.4 MB uncompressed for ~2,700 rows.

## 4. Publishing checklist

1. Build each bundle (`make_lightweight_bundle.py` for structures; see
   `lightweight-bundle.md` §12 — or the CSV+FASTA convention above for cLIP).
2. Name files by a stable public identifier, never an internal job name.
3. Commit and push to a dedicated GitHub Pages repo (`git add`, avoid `git add -A` — check
   `git status` for anything that shouldn't ship).
4. Enable GitHub Pages on that repo (Settings → Pages → deploy from the branch root) if not
   already on.
5. Test one deep link end to end before publishing the catalog page that links to all of them.
6. If cataloging with an index file, regenerate it as bundles are added — `available.json`'s
   shape (`{id: displayName}`) is a reasonable default; keep it small since some consumers may
   fetch the whole thing on page load.
