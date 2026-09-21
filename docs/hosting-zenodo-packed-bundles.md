# Hosting bundles on Zenodo (packed catalog + RemoteZip)

**Goal:** archive a catalog of lightweight bundles (see
[`lightweight-bundle.md`](lightweight-bundle.md)) on Zenodo for a permanent DOI, without
running into Zenodo's 100-files-per-record limit — by packing every bundle into **one** zip
and having LIVIA pull out just the one requested bundle via HTTP Range requests.

For the common case — no DOI requirement, just serving files — use the simpler
[GitHub Pages guide](hosting-github-pages.md) instead. GitHub Pages has no comparable
file-count limit; Zenodo does, which is the entire reason this pattern exists.

Reference implementation: [`Visco-ID-interactome-screen`](https://github.com/flyark/Visco-ID-interactome-screen),
DOI [10.5281/zenodo.22821570](https://doi.org/10.5281/zenodo.22821570) — ~115 structure bundles
packed into one `viscoid-interactome-bundles.zip`. (Its cLIP bundles stayed on GitHub Pages —
smaller, no DOI need, see that repo's `Visco-ID-interactome-screen-structures/clip/`.)

---

## 1. Why packing is necessary

Zenodo caps a single deposit at **100 files**. A screen with hundreds of interacting pairs
blows past that immediately if each bundle is uploaded as its own file. Packing every bundle
into one zip sidesteps the limit — Zenodo only ever sees one file — while a normal unzip tool
still opens it as an ordinary archive for anyone who downloads it directly.

The tradeoff: LIVIA can no longer just `fetch()` one bundle's URL, because there isn't one —
every bundle is now an *entry inside* a single large archive. That's what
[`js/remotezip.js`](../js/remotezip.js) solves.

## 2. How RemoteZip works

`RemoteZip.remoteZipExtract(url, entryName)` fetches exactly one named entry out of a remote
zip using only `Range:` HTTP requests — never downloading the archive as a whole:

1. `HEAD`/range-probe the URL to get the total file size.
2. Range-fetch the last ~64 KiB (grown and retried if needed) to find the End Of Central
   Directory record and parse the central directory — the list of every entry's name, size,
   compression method, and byte offset.
3. Range-fetch just that one entry's local file header (30 bytes) to find where its data
   actually starts, then range-fetch the entry's compressed bytes.
4. Inflate in-browser via the native `DecompressionStream('deflate-raw')` if the entry was
   DEFLATE-compressed; pass through unchanged if STORED (uncompressed).

Total network cost for one bundle: a handful of small range requests, independent of how many
other bundles or how large the whole archive is. **Requires the host to honor `Range:` requests
with a real `206 Partial Content` response** — confirmed working against both Zenodo's
`/api/records/<id>/files/<name>/content` endpoint and GitHub Pages.

## 3. Pack the bundles

Standard zip tools work — RemoteZip supports both STORED (method 0) and DEFLATE (method 8)
compression, so either `zip` or Python's `zipfile` is fine. From a directory of already-built
lightweight bundles:

```bash
cd bundles/
zip -r -X ../catalog-bundles.zip *.zip   # -X strips extra metadata zip doesn't need
```

or in Python, if you're assembling the catalog programmatically:

```python
import zipfile
with zipfile.ZipFile('catalog-bundles.zip', 'w', zipfile.ZIP_DEFLATED) as z:
    for bundle_path in sorted(bundle_dir.glob('*.zip')):
        z.write(bundle_path, arcname=bundle_path.name)   # arcname = the `entry=` value later
```

The `arcname` (the filename *inside* the packed zip, not its path on disk) is exactly what you
pass as `entry=` in the deep-link URL — keep it a flat filename, no directories, matching
whatever your CSV/catalog already calls the bundle (e.g. `interactions.csv`'s `bundle` column
in the Visco-ID screen).

Nested zips are fine — a lightweight bundle is already a zip, and it just becomes one STORED or
DEFLATED entry inside the outer catalog zip. RemoteZip extracts the entry's raw bytes and hands
them to LIVIA as-is; LIVIA then unzips *that* normally, same as any other bundle.

## 4. Upload to Zenodo

1. New upload → add files → the **one** packed zip (not the individual bundles).
2. Fill in metadata, publish → note the resulting **record ID** (the numeric ID in the
   record's URL) and the exact **filename** you uploaded it as.
3. Construct the direct-content URL:
   ```
   https://zenodo.org/api/records/<record-id>/files/<packed-zip-filename>/content
   ```
   This is the URL RemoteZip range-fetches against — confirm it in a browser first (it should
   download the full zip if opened directly; Range support is what makes LIVIA's partial fetch
   work, not this URL shape itself).

## 5. Deep-link URL pattern

```
https://flyark.github.io/LIVIA/universal.html?data=<zenodo-content-url>&entry=<bundle-filename>&name=<label>
```

- `data` — the packed zip's Zenodo content URL (same for every row/link in your catalog).
- `entry` — the one bundle's filename *inside* the packed zip (the `arcname` from §3).
- `name` — display label for the loading status message.

`universal.html`'s handler (search `Packed-catalog variant` in the file) tries RemoteZip's
Range extraction first; if the host doesn't support Range, or anything else fails, it falls
back to downloading the whole archive once via `fetch()` (proxied through
`livia-proxy.flyark.workers.dev` if direct fetch is blocked) and pulling the entry out locally
with JSZip. Slower on failure, but never broken outright.

This URL scheme is identical whether the catalog has 10 bundles or 10,000 — only `entry=`
changes per row, so a table like `Visco-ID-interactome-screen`'s just needs one `data=` constant
and a `bundle` column to build every link.

## 6. Publishing checklist

1. Build every bundle (`lightweight-bundle.md` §12).
2. Pack them into one zip, `entry=`-friendly flat filenames (§3).
3. Upload the one packed zip to Zenodo, publish, get the DOI + content URL (§4).
4. Test one `?data=&entry=` deep link end to end — confirm it actually range-fetches (check
   the browser's network tab for `206` responses, not a full-file `200` download) before
   publishing the catalog page.
5. Keep the packed zip's Zenodo upload as the version of record — a new Zenodo *version* gets a
   new record ID, so an existing catalog's `data=` links break unless updated. Uploading a
   correction as a new version, not editing files in place, is how Zenodo enforces its
   permanence guarantee — plan the catalog page's update process around that.
