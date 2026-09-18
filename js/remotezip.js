// RemoteZip — extract ONE named entry out of a large remote ZIP using only HTTP Range requests,
// without downloading the whole archive. Lets a static deep-link catalog pack hundreds of lightweight
// bundles into a single uploaded file (e.g. to stay under a host's per-record file-count quota, such
// as Zenodo's 100-files-per-record limit) while still streaming just one bundle per click.
//
// Requires the host to honor `Range:` GET requests with a real 206 response (confirmed for both
// Zenodo's `/api/records/<id>/files/<name>/content` endpoint and GitHub Pages). Verified byte-identical
// against the original standalone files, end to end through LIVIA's normal load path.
//
// ZIP layout (from the end of the file):
//   [ ... file data ... ] [ central directory records ] [ End Of Central Directory (EOCD), 22 bytes ]
// EOCD gives the central directory's offset + size + entry count. Each central directory record gives
// a filename, compressed size, compression method, and the offset of that entry's LOCAL file header
// (which itself precedes the entry's actual data by a further 30+ bytes, variable on the local header's
// own filename/extra-field lengths — hence the second small fetch in extractEntry()).

const EOCD_SIG = 0x06054b50;
const CD_SIG = 0x02014b50;
const LFH_SIG = 0x04034b50;

async function rzRangeFetch(url, start, end) {
  // end is inclusive, matching HTTP Range semantics
  const resp = await fetch(url, { headers: { Range: `bytes=${start}-${end}` } });
  if (resp.status !== 206) throw new Error(`expected 206 Partial Content, got ${resp.status} (server may not support Range)`);
  return await resp.arrayBuffer();
}

async function rzGetRemoteSize(url) {
  const resp = await fetch(url, { headers: { Range: 'bytes=0-0' } });
  const cr = resp.headers.get('Content-Range'); // may be unexposed by CORS on some hosts — don't rely on it
  await resp.arrayBuffer();
  if (cr) { const m = /\/(\d+)$/.exec(cr); if (m) return parseInt(m[1], 10); }
  const head = await fetch(url, { method: 'HEAD' });
  const cl = head.headers.get('Content-Length');
  if (cl) return parseInt(cl, 10);
  throw new Error('could not determine remote file size');
}

async function rzFetchCentralDirectory(url, totalSize) {
  // 64 KiB is generous for a zip comment + EOCD + a few hundred entries' worth of central-directory
  // records; grow the window and retry if the EOCD signature isn't found (e.g. a large zip comment).
  let windowSize = 65536;
  for (let attempt = 0; attempt < 3; attempt++) {
    const start = Math.max(0, totalSize - windowSize);
    const buf = await rzRangeFetch(url, start, totalSize - 1);
    const dv = new DataView(buf);
    for (let i = buf.byteLength - 22; i >= 0; i--) {
      if (dv.getUint32(i, true) === EOCD_SIG) {
        const entryCount = dv.getUint16(i + 10, true);
        const cdSize = dv.getUint32(i + 12, true);
        const cdOffset = dv.getUint32(i + 16, true);
        const cdStart = cdOffset;
        const cdEnd = cdOffset + cdSize;
        let cdBuf;
        if (cdStart >= start) {
          cdBuf = buf.slice(cdStart - start, cdEnd - start);   // central directory already inside this window
        } else {
          cdBuf = await rzRangeFetch(url, cdStart, cdEnd - 1); // window was too small — one more fetch
        }
        return rzParseCentralDirectory(cdBuf, entryCount);
      }
    }
    windowSize *= 4; // EOCD not found in this window (unusually large comment) — widen and retry
  }
  throw new Error('EOCD signature not found — not a valid zip, or comment too large');
}

function rzParseCentralDirectory(buf, entryCount) {
  const dv = new DataView(buf);
  const entries = [];
  let p = 0;
  for (let i = 0; i < entryCount; i++) {
    if (dv.getUint32(p, true) !== CD_SIG) throw new Error(`bad central directory record signature at entry ${i}`);
    const method = dv.getUint16(p + 10, true);
    const compSize = dv.getUint32(p + 20, true);
    const uncompSize = dv.getUint32(p + 24, true);
    const nameLen = dv.getUint16(p + 28, true);
    const extraLen = dv.getUint16(p + 30, true);
    const commentLen = dv.getUint16(p + 32, true);
    const localHeaderOffset = dv.getUint32(p + 42, true);
    const nameBytes = new Uint8Array(buf, p + 46, nameLen);
    const name = new TextDecoder().decode(nameBytes);
    entries.push({ name, method, compSize, uncompSize, localHeaderOffset });
    p += 46 + nameLen + extraLen + commentLen;
  }
  return entries;
}

async function rzExtractEntry(url, entry) {
  // Local file header repeats name/extra fields (sometimes with different extra-field content than
  // the central directory copy), so its total size must be read from the local header itself before
  // we know where the actual entry data starts. Fetch the fixed 30-byte header first.
  const lfhBuf = await rzRangeFetch(url, entry.localHeaderOffset, entry.localHeaderOffset + 29);
  const dv = new DataView(lfhBuf);
  if (dv.getUint32(0, true) !== LFH_SIG) throw new Error('bad local file header signature');
  const nameLen = dv.getUint16(26, true);
  const extraLen = dv.getUint16(28, true);
  const dataStart = entry.localHeaderOffset + 30 + nameLen + extraLen;
  const dataEnd = dataStart + entry.compSize - 1;
  const raw = await rzRangeFetch(url, dataStart, dataEnd);
  if (entry.method === 0) return raw; // STORED — already the original bytes, no inflate needed
  if (entry.method === 8) {
    // DEFLATE — browsers can inflate raw deflate streams natively via DecompressionStream.
    const ds = new DecompressionStream('deflate-raw');
    const stream = new Blob([raw]).stream().pipeThrough(ds);
    return await new Response(stream).arrayBuffer();
  }
  throw new Error(`unsupported compression method ${entry.method}`);
}

// Convenience: fetch just the named entry's bytes out of a remote zip in one call. Throws if the host
// doesn't support Range, or the entry isn't found — callers should fall back to a full download + local
// unzip (e.g. via JSZip, already used elsewhere in LIVIA) on failure rather than surface this directly.
async function remoteZipExtract(url, entryName) {
  const size = await rzGetRemoteSize(url);
  const entries = await rzFetchCentralDirectory(url, size);
  const entry = entries.find(e => e.name === entryName);
  if (!entry) throw new Error(`entry "${entryName}" not found in remote zip (${entries.length} entries present)`);
  return await rzExtractEntry(url, entry);
}

window.RemoteZip = { getRemoteSize: rzGetRemoteSize, fetchCentralDirectory: rzFetchCentralDirectory, extractEntry: rzExtractEntry, remoteZipExtract };
