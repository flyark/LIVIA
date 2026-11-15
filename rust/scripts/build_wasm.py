#!/usr/bin/env python3
"""Bundle livia_wasm.wasm into a self-contained universal JavaScript wrapper (livia-wasm.js)."""

import base64
import gzip
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
WASM_PATH = ROOT / "target" / "wasm32-unknown-unknown" / "wasm-release" / "livia_wasm.wasm"
OUT_JS = ROOT / "dist" / "livia-wasm.js"

JS_TEMPLATE = r"""/* LIVIA WebAssembly Compute Engine (livia-wasm.js)
 * Self-contained universal module for main thread, Web Workers, file://, offline reports, and strict CSP.
 * Automatically falls back to pure JavaScript when WebAssembly is unavailable or blocked.
 */
(function (root) {
  'use strict';

  const WASM_B64 = '__WASM_B64__';

  let exportsObj = null;
  let initPromise = null;
  let initFailed = false;

  const decodeBase64 = function (b64) {
    if (typeof atob === 'function') {
      const bin = atob(b64);
      const len = bin.length;
      const bytes = new Uint8Array(len);
      for (let i = 0; i < len; i++) bytes[i] = bin.charCodeAt(i);
      return bytes;
    }
    if (typeof Buffer !== 'undefined') {
      const buf = Buffer.from(b64, 'base64');
      return new Uint8Array(buf.buffer, buf.byteOffset, buf.byteLength);
    }
    const lut = new Int16Array(128).fill(-1);
    const alpha = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/';
    for (let i = 0; i < 64; i++) lut[alpha.charCodeAt(i)] = i;
    let cleanLen = b64.length;
    while (cleanLen > 0 && b64.charCodeAt(cleanLen - 1) === 61) cleanLen--;
    const outLen = Math.floor((cleanLen * 3) / 4);
    const out = new Uint8Array(outLen);
    let p = 0;
    for (let i = 0; i < cleanLen; i += 4) {
      const a = lut[b64.charCodeAt(i)];
      const b = i + 1 < cleanLen ? lut[b64.charCodeAt(i + 1)] : 0;
      const c = i + 2 < cleanLen ? lut[b64.charCodeAt(i + 2)] : 0;
      const d = i + 3 < cleanLen ? lut[b64.charCodeAt(i + 3)] : 0;
      const trip = (a << 18) | (b << 12) | (c << 6) | d;
      if (p < outLen) out[p++] = (trip >>> 16) & 255;
      if (p < outLen) out[p++] = (trip >>> 8) & 255;
      if (p < outLen) out[p++] = trip & 255;
    }
    return out;
  };

  let cachedBytes = null;
  const getBytes = function () {
    if (!cachedBytes) cachedBytes = decodeBase64(WASM_B64);
    return cachedBytes;
  };

  const trySyncInit = function () {
    if (exportsObj) return true;
    if (initFailed) return false;
    if (typeof WebAssembly !== 'object' || WebAssembly === null) {
      initFailed = true;
      return false;
    }
    const isWorker = typeof WorkerGlobalScope !== 'undefined' && root instanceof WorkerGlobalScope;
    const isNode = typeof process !== 'undefined' && !!(process.versions && process.versions.node);
    if (!isWorker && !isNode) return false;
    try {
      const mod = new WebAssembly.Module(getBytes());
      const inst = new WebAssembly.Instance(mod, {});
      exportsObj = inst.exports;
      cachedBytes = null;
      return true;
    } catch (_err) {
      return false;
    }
  };

  const ready = function () {
    if (exportsObj) return Promise.resolve(true);
    if (initFailed) return Promise.resolve(false);
    if (trySyncInit()) return Promise.resolve(true);
    if (initPromise) return initPromise;
    if (typeof WebAssembly !== 'object' || WebAssembly === null || typeof WebAssembly.instantiate !== 'function') {
      initFailed = true;
      return Promise.resolve(false);
    }
    initPromise = WebAssembly.instantiate(getBytes(), {})
      .then(function (res) {
        const inst = res.instance || res;
        exportsObj = inst.exports;
        cachedBytes = null;
        return true;
      })
      .catch(function () {
        initFailed = true;
        return false;
      });
    return initPromise;
  };

  trySyncInit();
  if (!exportsObj && !initFailed) {
    ready();
  }

  const isReady = function () {
    return exportsObj !== null;
  };

  const allocF64 = function (arr) {
    const byteLen = arr.length * 8;
    const ptr = exportsObj.wasm_alloc(byteLen);
    const view = new Float64Array(exportsObj.memory.buffer, ptr, arr.length);
    view.set(arr);
    return { ptr: ptr, byteLen: byteLen };
  };

  const allocF32 = function (arr) {
    const byteLen = arr.length * 4;
    const ptr = exportsObj.wasm_alloc(byteLen);
    const view = new Float32Array(exportsObj.memory.buffer, ptr, arr.length);
    view.set(arr);
    return { ptr: ptr, byteLen: byteLen };
  };

  const allocU32 = function (arr) {
    const byteLen = arr.length * 4;
    const ptr = exportsObj.wasm_alloc(byteLen);
    const view = new Uint32Array(exportsObj.memory.buffer, ptr, arr.length);
    view.set(arr);
    return { ptr: ptr, byteLen: byteLen };
  };

  const freeBuf = function (b) {
    if (b && b.ptr) exportsObj.wasm_free(b.ptr, b.byteLen);
  };

  const readOutU32 = function () {
    const ptr = exportsObj.wasm_out_u32_ptr();
    const len = exportsObj.wasm_out_u32_len();
    return new Uint32Array(exportsObj.memory.buffer, ptr, len).slice();
  };

  const readOutF64 = function () {
    const ptr = exportsObj.wasm_out_f64_ptr();
    const len = exportsObj.wasm_out_f64_len();
    return new Float64Array(exportsObj.memory.buffer, ptr, len).slice();
  };

  const readOutF32 = function () {
    const ptr = exportsObj.wasm_out_f32_ptr();
    const len = exportsObj.wasm_out_f32_len();
    return new Float32Array(exportsObj.memory.buffer, ptr, len).slice();
  };

  const makePrng = function (seed) {
    let a = (seed >>> 0);
    return function () {
      a = (a + 0x6D2B79F5) | 0;
      let t = Math.imul(a ^ (a >>> 15), 1 | a);
      t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
  };

  const shuffledOrder = function (n, rng) {
    const ord = new Array(n);
    for (let i = 0; i < n; i++) ord[i] = i;
    for (let i = n - 1; i > 0; i--) {
      const j = Math.floor(rng() * (i + 1));
      const tmp = ord[i];
      ord[i] = ord[j];
      ord[j] = tmp;
    }
    return ord;
  };

  const jsLeidenSingle = function (n, simpleEdges, degW, m2, gamma, useCpm, rng) {
    let curN = n;
    let curEdges = simpleEdges.slice();
    let curSizes = new Float64Array(n).fill(1);
    let pPart = new Array(n);
    let members = new Array(n);
    for (let i = 0; i < n; i++) {
      pPart[i] = i;
      members[i] = [i];
    }
    if (!useCpm && m2 <= 0) {
      return { part: pPart, score: 0 };
    }
    for (let level = 0; level < 20; level++) {
      const adj = new Array(curN);
      for (let i = 0; i < curN; i++) adj[i] = [];
      const kW = new Float64Array(curN);
      for (let e = 0; e < curEdges.length; e++) {
        const u = curEdges[e][0], v = curEdges[e][1], w = curEdges[e][2];
        if (u === v) {
          kW[u] += 2 * w;
        } else {
          adj[u].push([v, w]);
          adj[v].push([u, w]);
          kW[u] += w;
          kW[v] += w;
        }
      }
      const maxC = Math.max.apply(null, pPart) + 1;
      const totW = new Float64Array(Math.max(curN, maxC));
      const totSz = new Float64Array(Math.max(curN, maxC));
      for (let i = 0; i < curN; i++) {
        totW[pPart[i]] += kW[i];
        totSz[pPart[i]] += curSizes[i];
      }
      let moved = true;
      let rounds = 0;
      while (moved && rounds < 50) {
        moved = false;
        rounds++;
        const ord = shuffledOrder(curN, rng);
        for (let oi = 0; oi < curN; oi++) {
          const i = ord[oi];
          const ci = pPart[i];
          const wt = new Map();
          const nb = adj[i];
          for (let ni = 0; ni < nb.length; ni++) {
            const cj = pPart[nb[ni][0]];
            wt.set(cj, (wt.get(cj) || 0) + nb[ni][1]);
          }
          totW[ci] -= kW[i];
          totSz[ci] -= curSizes[i];
          const stayPen = useCpm ? gamma * curSizes[i] * totSz[ci] : (gamma * kW[i] * totW[ci]) / m2;
          let bestC = ci;
          let bestGain = (wt.get(ci) || 0) - stayPen;
          for (const entry of wt.entries()) {
            const cj = entry[0], wTo = entry[1];
            const pen = useCpm ? gamma * curSizes[i] * totSz[cj] : (gamma * kW[i] * totW[cj]) / m2;
            const g = wTo - pen;
            if (g > bestGain + 1e-12) {
              bestGain = g;
              bestC = cj;
            }
          }
          totW[bestC] += kW[i];
          totSz[bestC] += curSizes[i];
          if (bestC !== ci) {
            pPart[i] = bestC;
            moved = true;
          }
        }
      }
      const commMap = new Map();
      for (let i = 0; i < curN; i++) {
        const c = pPart[i];
        if (!commMap.has(c)) commMap.set(c, []);
        commMap.get(c).push(i);
      }
      if (commMap.size === curN) break;
      const nextN = commMap.size;
      const cIds = new Map();
      let idxCounter = 0;
      for (const c of commMap.keys()) cIds.set(c, idxCounter++);
      const nextMembers = new Array(nextN);
      const nextSizes = new Float64Array(nextN);
      for (let i = 0; i < nextN; i++) nextMembers[i] = [];
      for (let i = 0; i < curN; i++) {
        const ni = cIds.get(pPart[i]);
        nextMembers[ni].push.apply(nextMembers[ni], members[i]);
        nextSizes[ni] += curSizes[i];
      }
      const edgeAcc = new Map();
      for (let e = 0; e < curEdges.length; e++) {
        const a = cIds.get(pPart[curEdges[e][0]]);
        const b = cIds.get(pPart[curEdges[e][1]]);
        const lo = a < b ? a : b;
        const hi = a < b ? b : a;
        const key = lo * nextN + hi;
        edgeAcc.set(key, (edgeAcc.get(key) || 0) + curEdges[e][2]);
      }
      curEdges = [];
      for (const entry of edgeAcc.entries()) {
        const key = entry[0];
        curEdges.push([Math.floor(key / nextN), key % nextN, entry[1]]);
      }
      curN = nextN;
      curSizes = nextSizes;
      members = nextMembers;
      pPart = new Array(curN);
      for (let i = 0; i < curN; i++) pPart[i] = i;
    }
    const flatPart = new Array(n).fill(0);
    const groups = new Map();
    for (let i = 0; i < curN; i++) {
      const c = pPart[i];
      if (!groups.has(c)) groups.set(c, []);
      groups.get(c).push.apply(groups.get(c), members[i]);
    }
    let cid = 0;
    for (const grp of groups.values()) {
      for (let j = 0; j < grp.length; j++) flatPart[grp[j]] = cid;
      cid++;
    }
    let internalW = 0;
    for (let e = 0; e < simpleEdges.length; e++) {
      if (flatPart[simpleEdges[e][0]] === flatPart[simpleEdges[e][1]]) {
        internalW += simpleEdges[e][2];
      }
    }
    let score = 0;
    if (useCpm) {
      const sz = new Float64Array(cid);
      for (let i = 0; i < n; i++) sz[flatPart[i]] += 1;
      let pen = 0;
      for (let c = 0; c < cid; c++) pen += gamma * (sz[c] * (sz[c] - 1)) * 0.5;
      score = internalW - pen;
    } else if (m2 > 0) {
      const tw = new Float64Array(cid);
      for (let i = 0; i < n; i++) tw[flatPart[i]] += degW[i];
      let pen = 0;
      for (let c = 0; c < cid; c++) pen += gamma * (tw[c] / m2) * (tw[c] / m2);
      score = (2 * internalW) / m2 - pen;
    }
    return { part: flatPart, score: score };
  };

  const jsLeidenNetworkFallback = function (n, edges, gamma, useCpm, nIterations, seed) {
    if (n === 0) return { membership: [], degree: [] };
    const pairMap = new Map();
    const degSet = new Array(n);
    for (let i = 0; i < n; i++) degSet[i] = new Set();
    for (let i = 0; i < edges.length; i++) {
      const u = edges[i][0] | 0, v = edges[i][1] | 0, w = +edges[i][2];
      if (u === v || u < 0 || v < 0 || u >= n || v >= n) continue;
      degSet[u].add(v);
      degSet[v].add(u);
      const lo = u < v ? u : v;
      const hi = u < v ? v : u;
      const key = lo * n + hi;
      pairMap.set(key, (pairMap.get(key) || 0) + w);
    }
    const degree = degSet.map(function (s) { return s.size; });
    const simpleEdges = [];
    const degW = new Float64Array(n);
    let m2 = 0;
    for (const entry of pairMap.entries()) {
      const u = Math.floor(entry[0] / n), v = entry[0] % n, w = entry[1];
      simpleEdges.push([u, v, w]);
      degW[u] += w;
      degW[v] += w;
      m2 += 2 * w;
    }
    const rng = makePrng((seed >>> 0) || 42);
    const runs = Math.max(1, nIterations | 0);
    let bestPart = null;
    let bestScore = -Infinity;
    for (let r = 0; r < runs; r++) {
      const res = jsLeidenSingle(n, simpleEdges, degW, m2, gamma, useCpm, rng);
      if (res.score > bestScore + 1e-12 || bestPart === null) {
        bestScore = res.score;
        bestPart = res.part;
      }
    }
    const groups = new Map();
    for (let i = 0; i < n; i++) {
      const c = bestPart[i];
      if (!groups.has(c)) groups.set(c, []);
      groups.get(c).push(i);
    }
    const sorted = Array.from(groups.values()).sort(function (a, b) {
      return b.length - a.length || a[0] - b[0];
    });
    const membership = new Array(n).fill(0);
    for (let rank = 0; rank < sorted.length; rank++) {
      const g = sorted[rank];
      for (let j = 0; j < g.length; j++) membership[g[j]] = rank;
    }
    return { membership: membership, degree: degree };
  };

  const jsGraphLayoutFallback = function (method, n, edges, niter, seed) {
    if (n === 0) return [];
    if (n === 1) return [[0, 0]];
    const rng = makePrng((seed >>> 0) || 42);
    const px = new Float64Array(n);
    const py = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      px[i] = rng() - 0.5;
      py[i] = rng() - 0.5;
    }
    if (method === 'kk') {
      const adj = new Array(n);
      for (let i = 0; i < n; i++) adj[i] = [];
      for (let e = 0; e < edges.length; e++) {
        const u = edges[e][0] | 0, v = edges[e][1] | 0;
        if (u !== v && u < n && v < n) {
          adj[u].push(v);
          adj[v].push(u);
        }
      }
      const dm = new Int16Array(n * n).fill(-1);
      let dmax = 1;
      for (let s = 0; s < n; s++) {
        const base = s * n;
        dm[base + s] = 0;
        const q = [s];
        for (let h = 0; h < q.length; h++) {
          const u = q[h];
          const nb = adj[u];
          for (let k = 0; k < nb.length; k++) {
            const v = nb[k];
            if (dm[base + v] < 0) {
              dm[base + v] = dm[base + u] + 1;
              if (dm[base + v] > dmax) dmax = dm[base + v];
              q.push(v);
            }
          }
        }
      }
      for (let i = 0; i < n; i++) {
        px[i] *= dmax;
        py[i] *= dmax;
      }
      const rounds = Math.max(50, Math.min(300, niter || 200));
      for (let it = 0; it < rounds; it++) {
        for (let i = 0; i < n; i++) {
          let sx = 0, sy = 0, sw = 0;
          const base = i * n;
          for (let j = 0; j < n; j++) {
            if (i === j) continue;
            const dij = dm[base + j] < 0 ? dmax + 1 : dm[base + j];
            const w = 1 / (dij * dij);
            const dx = px[i] - px[j], dy = py[i] - py[j];
            const l = Math.max(1e-6, Math.hypot(dx, dy));
            sx += w * (px[j] + (dij * dx) / l);
            sy += w * (py[j] + (dij * dy) / l);
            sw += w;
          }
          if (sw > 0) {
            px[i] = sx / sw;
            py[i] = sy / sw;
          }
        }
      }
    } else {
      const iters = Math.max(50, niter || 500);
      const kIdeal = Math.sqrt(1 / Math.max(1, n));
      let wmax = 1e-9;
      for (let e = 0; e < edges.length; e++) {
        if (edges[e][2] > wmax) wmax = edges[e][2];
      }
      const dxArr = new Float64Array(n);
      const dyArr = new Float64Array(n);
      for (let it = 0; it < iters; it++) {
        const temp = 0.12 * (1 - it / iters) + 0.0015;
        dxArr.fill(0);
        dyArr.fill(0);
        for (let i = 0; i < n; i++) {
          for (let j = i + 1; j < n; j++) {
            const dx = px[i] - px[j], dy = py[i] - py[j];
            const d2 = Math.max(1e-8, dx * dx + dy * dy);
            const f = (kIdeal * kIdeal) / d2;
            dxArr[i] += dx * f;
            dyArr[i] += dy * f;
            dxArr[j] -= dx * f;
            dyArr[j] -= dy * f;
          }
        }
        for (let e = 0; e < edges.length; e++) {
          const u = edges[e][0] | 0, v = edges[e][1] | 0, w = +edges[e][2];
          if (u === v || u >= n || v >= n) continue;
          const dx = px[u] - px[v], dy = py[u] - py[v];
          const d = Math.max(1e-4, Math.hypot(dx, dy));
          const f = (d / kIdeal) * (0.3 + (0.7 * w) / wmax);
          dxArr[u] -= dx * f;
          dyArr[u] -= dy * f;
          dxArr[v] += dx * f;
          dyArr[v] += dy * f;
        }
        for (let i = 0; i < n; i++) {
          const l = Math.max(1e-9, Math.hypot(dxArr[i], dyArr[i]));
          const step = Math.min(l, temp) / l;
          px[i] = (px[i] + dxArr[i] * step) * 0.99;
          py[i] = (py[i] + dyArr[i] * step) * 0.99;
        }
      }
    }
    const out = new Array(n);
    for (let i = 0; i < n; i++) out[i] = [px[i], py[i]];
    return out;
  };

  const packEdges = function (edges) {
    const m = edges.length;
    const flat = new Float64Array(m * 3);
    for (let i = 0; i < m; i++) {
      const e = edges[i];
      flat[i * 3] = e[0];
      flat[i * 3 + 1] = e[1];
      flat[i * 3 + 2] = e.length > 2 ? e[2] : 1;
    }
    return flat;
  };

  const leidenNetwork = function (n, edges, gamma, objective, nIterations, seed) {
    const useCpm = objective === 'CPM' || objective === 'cpm';
    if (!exportsObj) {
      return jsLeidenNetworkFallback(n, edges, gamma, useCpm, nIterations || 10, seed || 42);
    }
    const flat = packEdges(edges);
    const eBuf = allocF64(flat);
    try {
      exportsObj.wasm_leiden_network(
        n >>> 0,
        eBuf.ptr,
        edges.length >>> 0,
        +gamma,
        useCpm ? 1 : 0,
        (nIterations || 10) >>> 0,
        (seed >>> 0)
      );
      const u32 = readOutU32();
      return {
        membership: Array.from(u32.subarray(0, n)),
        degree: Array.from(u32.subarray(n, 2 * n)),
      };
    } finally {
      freeBuf(eBuf);
    }
  };

  const leidenAtlas = function (n, edges, gamma, seed) {
    if (!exportsObj) return null;
    const flat = packEdges(edges);
    const eBuf = allocF64(flat);
    try {
      exportsObj.wasm_leiden_atlas(
        n >>> 0,
        eBuf.ptr,
        edges.length >>> 0,
        gamma == null ? 1.0 : +gamma,
        (seed || 0) >>> 0
      );
      return Array.from(readOutU32());
    } finally {
      freeBuf(eBuf);
    }
  };

  const steinerAdd = function (nNodes, edges, terms, prizesArr, cap, directed) {
    if (!exportsObj) return null;
    const isFlat = edges instanceof Float64Array;
    const flatE = isFlat ? edges : packEdges(edges);
    const nEdges = isFlat ? ((edges.length / 3) | 0) : edges.length;
    const termsU32 = Uint32Array.from(terms);
    const prizesF64 = prizesArr && prizesArr.length ? Float64Array.from(prizesArr) : new Float64Array(0);
    const eBuf = allocF64(flatE);
    const tBuf = allocU32(termsU32);
    const pBuf = prizesF64.length ? allocF64(prizesF64) : { ptr: 0, byteLen: 0 };
    try {
      exportsObj.wasm_steiner_add(
        nNodes >>> 0,
        eBuf.ptr,
        nEdges >>> 0,
        tBuf.ptr,
        termsU32.length >>> 0,
        pBuf.ptr,
        prizesF64.length >>> 0,
        cap >>> 0,
        directed ? 1 : 0
      );
      const f64 = readOutF64();
      const u32 = readOutU32();
      return {
        joined: f64[0] | 0,
        cost: f64[1],
        nodes: Array.from(u32),
      };
    } finally {
      freeBuf(eBuf);
      freeBuf(tBuf);
      freeBuf(pBuf);
    }
  };

  const graphLayout = function (method, n, edges, niter, seed) {
    if (!exportsObj) {
      return jsGraphLayoutFallback(method, n, edges, niter || 500, seed || 42);
    }
    const code = method === 'kk' ? 1 : (method === 'drl' ? 2 : 0);
    const flat = packEdges(edges);
    const eBuf = allocF64(flat);
    try {
      exportsObj.wasm_graph_layout(
        code >>> 0,
        n >>> 0,
        eBuf.ptr,
        edges.length >>> 0,
        (niter || 500) >>> 0,
        (seed >>> 0)
      );
      const f64 = readOutF64();
      const out = new Array(n);
      for (let i = 0; i < n; i++) {
        out[i] = [f64[i * 2], f64[i * 2 + 1]];
      }
      return out;
    } finally {
      freeBuf(eBuf);
    }
  };

  const packFingerprintsCsr = function (fingerprints) {
    const n = fingerprints.length;
    const offsets = new Uint32Array(n + 1);
    let total = 0;
    const lists = new Array(n);
    for (let i = 0; i < n; i++) {
      const fp = fingerprints[i];
      const arr = Array.isArray(fp) ? fp : Array.from(fp);
      lists[i] = arr;
      total += arr.length;
      offsets[i + 1] = total;
    }
    const flat = new Uint32Array(total);
    let pos = 0;
    for (let i = 0; i < n; i++) {
      const arr = lists[i];
      for (let j = 0; j < arr.length; j++) {
        flat[pos++] = arr[j] >>> 0;
      }
    }
    return { offsets: offsets, flat: flat };
  };

  const unpackLinkage = function (f64, n) {
    if (n <= 1) return [];
    const rows = n - 1;
    const Z = new Array(rows);
    for (let i = 0; i < rows; i++) {
      const b = i * 4;
      Z[i] = [f64[b], f64[b + 1], f64[b + 2], f64[b + 3]];
    }
    return Z;
  };

  const clipCosineDistUpper = function (fingerprints) {
    if (!exportsObj) return null;
    const csr = packFingerprintsCsr(fingerprints);
    const oBuf = allocU32(csr.offsets);
    const fBuf = allocU32(csr.flat);
    try {
      exportsObj.wasm_clip_cosine_dist(oBuf.ptr, csr.offsets.length >>> 0, fBuf.ptr, csr.flat.length >>> 0);
      return readOutF64();
    } finally {
      freeBuf(oBuf);
      freeBuf(fBuf);
    }
  };

  const clipFromFingerprints = function (fingerprints, maxK, dropThreshold) {
    if (!exportsObj) return null;
    const n = fingerprints.length;
    if (n === 0) {
      return { optimalK: 0, labels: [], rawLabels: [], order: [], linkage: [] };
    }
    if (n === 1) {
      return { optimalK: 1, labels: [1], rawLabels: [1], order: [0], linkage: [] };
    }
    const csr = packFingerprintsCsr(fingerprints);
    const oBuf = allocU32(csr.offsets);
    const fBuf = allocU32(csr.flat);
    try {
      exportsObj.wasm_clip_full(
        oBuf.ptr,
        csr.offsets.length >>> 0,
        fBuf.ptr,
        csr.flat.length >>> 0,
        (maxK || 30) >>> 0,
        dropThreshold == null ? 0.05 : +dropThreshold
      );
      const u32 = readOutU32();
      const f64 = readOutF64();
      const optimalK = u32[0];
      const labels = Array.from(u32.subarray(1, 1 + n));
      const rawLabels = Array.from(u32.subarray(1 + n, 1 + 2 * n));
      const order = Array.from(u32.subarray(1 + 2 * n, 1 + 3 * n));
      const linkage = unpackLinkage(f64, n);
      exportsObj.wasm_clear_outputs();
      return {
        optimalK: optimalK,
        labels: labels,
        rawLabels: rawLabels,
        order: order,
        linkage: linkage,
      };
    } finally {
      freeBuf(oBuf);
      freeBuf(fBuf);
    }
  };

  const clipFromDistMatrix = function (D, maxK, dropThreshold) {
    if (!exportsObj) return null;
    const n = D.length;
    if (n <= 1) {
      return {
        optimalK: n,
        labels: n === 1 ? [1] : [],
        rawLabels: n === 1 ? [1] : [],
        order: n === 1 ? [0] : [],
        linkage: [],
      };
    }
    const upperLen = (n * (n - 1)) / 2;
    const upper = new Float64Array(upperLen);
    let p = 0;
    for (let i = 0; i < n; i++) {
      const row = D[i];
      for (let j = i + 1; j < n; j++) {
        upper[p++] = row[j];
      }
    }
    const uBuf = allocF64(upper);
    try {
      exportsObj.wasm_clip_from_dist(
        uBuf.ptr,
        upperLen >>> 0,
        n >>> 0,
        (maxK || 10) >>> 0,
        dropThreshold == null ? 0.05 : +dropThreshold
      );
      const u32 = readOutU32();
      const f64 = readOutF64();
      const optimalK = u32[0];
      const labels = Array.from(u32.subarray(1, 1 + n));
      const rawLabels = Array.from(u32.subarray(1 + n, 1 + 2 * n));
      const order = Array.from(u32.subarray(1 + 2 * n, 1 + 3 * n));
      const linkage = unpackLinkage(f64, n);
      exportsObj.wasm_clear_outputs();
      return {
        optimalK: optimalK,
        labels: labels,
        rawLabels: rawLabels,
        order: order,
        linkage: linkage,
      };
    } finally {
      freeBuf(uBuf);
    }
  };

  const paeChainPairIptm = function (pae, starts, ends, nc) {
    if (!exportsObj) return null;
    const dim = pae.length;
    const flatPae = new Float64Array(dim * dim);
    for (let i = 0; i < dim; i++) {
      const row = pae[i];
      const rlen = row ? Math.min(dim, row.length) : 0;
      const base = i * dim;
      for (let j = 0; j < rlen; j++) flatPae[base + j] = row[j];
    }
    const stU32 = Uint32Array.from(starts);
    const enU32 = Uint32Array.from(ends);
    const pBuf = allocF64(flatPae);
    const sBuf = allocU32(stU32);
    const eBuf = allocU32(enU32);
    try {
      exportsObj.wasm_pae_chain_pair_iptm(pBuf.ptr, dim >>> 0, sBuf.ptr, eBuf.ptr, nc >>> 0);
      const f64 = readOutF64();
      const out = new Array(nc);
      for (let i = 0; i < nc; i++) {
        out[i] = Array.from(f64.subarray(i * nc, (i + 1) * nc));
      }
      return out;
    } finally {
      freeBuf(pBuf);
      freeBuf(sBuf);
      freeBuf(eBuf);
    }
  };

  const packCoords = function (coords) {
    const n = coords ? coords.length : 0;
    const flat = new Float64Array(n * 4);
    for (let i = 0; i < n; i++) {
      const c = coords[i];
      const b = i * 4;
      if (c) {
        flat[b] = c.x;
        flat[b + 1] = c.y;
        flat[b + 2] = c.z;
        flat[b + 3] = c.hasP ? 1 : 0;
      } else {
        flat[b] = NaN;
        flat[b + 1] = NaN;
        flat[b + 2] = NaN;
        flat[b + 3] = 0;
      }
    }
    return flat;
  };

  const computeModelPairs = function (opts) {
    if (!exportsObj) return null;
    const pae = opts.pae;
    const dim = pae.length;
    const nTotal = opts.nTotal;
    const chainNames = opts.chainNames;
    const nc = chainNames.length;
    const starts = opts.starts;
    const sizes = opts.sizes;
    const paeCutoff = opts.paeCutoff;
    const cbCutoff = opts.cbCutoff;

    const flatPae = new Float64Array(dim * dim);
    for (let i = 0; i < dim; i++) {
      const row = pae[i];
      const rlen = row ? Math.min(dim, row.length) : 0;
      const base = i * dim;
      for (let j = 0; j < rlen; j++) flatPae[base + j] = row[j];
    }

    const cbFlat = packCoords(opts.coords);
    const caFlat = packCoords(opts.caCoords);

    let totalChainRes = 0;
    for (let i = 0; i < nc; i++) totalChainRes += sizes[i];
    const plddtFlat = new Float64Array(totalChainRes).fill(NaN);
    if (opts.structText && opts.bfsAll && opts.bfsAll.size > 0) {
      let cur = 0;
      for (let c = 0; c < nc; c++) {
        const cname = chainNames[c];
        const sz = sizes[c];
        for (let r = 1; r <= sz; r++) {
          const v = opts.bfsAll.get(cname + ':' + r);
          if (v != null && !isNaN(v)) {
            plddtFlat[cur + r - 1] = v;
          }
        }
        cur += sz;
      }
    }

    const pairIptmFlat = new Float64Array(nc * nc);
    const iptmMat = opts.iptmMatrix;
    const gIptm = opts.globalIptm || 0;
    for (let i = 0; i < nc; i++) {
      for (let j = 0; j < nc; j++) {
        const cell = iptmMat && iptmMat[i] && iptmMat[i][j] != null ? iptmMat[i][j] : gIptm;
        pairIptmFlat[i * nc + j] = cell;
      }
    }

    const stU32 = Uint32Array.from(starts);
    const szU32 = Uint32Array.from(sizes);

    const pBuf = allocF64(flatPae);
    const cbBuf = allocF64(cbFlat);
    const caBuf = allocF64(caFlat);
    const bfBuf = allocF64(plddtFlat);
    const sBuf = allocU32(stU32);
    const zBuf = allocU32(szU32);
    const ipBuf = allocF64(pairIptmFlat);

    try {
      const nPairs = exportsObj.wasm_compute_model_pairs(
        pBuf.ptr,
        dim >>> 0,
        nTotal >>> 0,
        cbBuf.ptr,
        (opts.coords ? opts.coords.length : 0) >>> 0,
        caBuf.ptr,
        (opts.caCoords ? opts.caCoords.length : 0) >>> 0,
        bfBuf.ptr,
        totalChainRes >>> 0,
        sBuf.ptr,
        zBuf.ptr,
        nc >>> 0,
        ipBuf.ptr,
        +paeCutoff,
        +cbCutoff
      );

      const f32View = new Float32Array(
        exportsObj.memory.buffer,
        exportsObj.wasm_out_f32_ptr(),
        exportsObj.wasm_out_f32_len()
      );
      const mapCells = nTotal * nTotal;
      const lisMap = new Array(nTotal);
      const clisMap = new Array(nTotal);
      for (let r = 0; r < nTotal; r++) {
        const rOff = r * nTotal;
        lisMap[r] = f32View.slice(rOff, rOff + nTotal);
        clisMap[r] = f32View.slice(mapCells + rOff, mapCells + rOff + nTotal);
      }

      const f64 = readOutF64();
      const u32 = readOutU32();
      const optVal = function (v) { return isNaN(v) ? null : v; };
      const makeSet = function (off, len) {
        const s = new Set();
        for (let k = 0; k < len; k++) s.add(u32[off + k]);
        return s;
      };

      const pairs = {};
      for (let p = 0; p < nPairs; p++) {
        const b = p * 30;
        const ciIdx = f64[b] | 0;
        const cjIdx = f64[b + 1] | 0;
        const ci = chainNames[ciIdx];
        const cj = chainNames[cjIdx];
        const hb = p * 8;
        const lirI = makeSet(u32[hb], u32[hb + 1]);
        const lirJ = makeSet(u32[hb + 2], u32[hb + 3]);
        const clirI = makeSet(u32[hb + 4], u32[hb + 5]);
        const clirJ = makeSet(u32[hb + 6], u32[hb + 7]);

        pairs[ci + ',' + cj] = {
          ci: ci,
          cj: cj,
          LIS: f64[b + 2],
          cLIS: f64[b + 3],
          iLIS: f64[b + 4],
          LIA: f64[b + 5],
          cLIA: f64[b + 6],
          iLIA: f64[b + 7],
          iLISA: f64[b + 8],
          ipTM: f64[b + 9],
          ipSAE: f64[b + 10],
          actifpTM: f64[b + 11],
          LIpLDDT_i: f64[b + 12],
          LIpLDDT_j: f64[b + 13],
          cLIpLDDT_i: f64[b + 14],
          cLIpLDDT_j: f64[b + 15],
          geom_plddt_i: f64[b + 16],
          geom_plddt_j: f64[b + 17],
          pDockQ: optVal(f64[b + 18]),
          LIpDockQ: optVal(f64[b + 19]),
          pDockQ2_i: optVal(f64[b + 20]),
          pDockQ2_j: optVal(f64[b + 21]),
          LIpDockQ2_i: optVal(f64[b + 22]),
          LIpDockQ2_j: optVal(f64[b + 23]),
          geomNContacts: f64[b + 24],
          LIpDockQN: f64[b + 25],
          geomIfISize: f64[b + 26],
          geomIfJSize: f64[b + 27],
          lirI: lirI,
          lirJ: lirJ,
          clirI: clirI,
          clirJ: clirJ,
          lenI: f64[b + 28],
          lenJ: f64[b + 29],
        };
      }

      if (nTotal > 800) {
        exportsObj.wasm_clear_outputs();
      }
      return { pairs: pairs, lisMap: lisMap, clisMap: clisMap };
    } finally {
      freeBuf(pBuf);
      freeBuf(cbBuf);
      freeBuf(caBuf);
      freeBuf(bfBuf);
      freeBuf(sBuf);
      freeBuf(zBuf);
      freeBuf(ipBuf);
    }
  };

  root.LiviaWasm = {
    version: '1.0.0',
    ready: ready,
    isReady: isReady,
    leidenNetwork: leidenNetwork,
    leidenAtlas: leidenAtlas,
    steinerAdd: steinerAdd,
    graphLayout: graphLayout,
    clipCosineDistUpper: clipCosineDistUpper,
    clipFromFingerprints: clipFromFingerprints,
    clipFromDistMatrix: clipFromDistMatrix,
    paeChainPairIptm: paeChainPairIptm,
    computeModelPairs: computeModelPairs,
  };
})(typeof globalThis !== 'undefined' ? globalThis : self);
"""


def main():
    wasm_bytes = WASM_PATH.read_bytes()
    b64 = base64.b64encode(wasm_bytes).decode("ascii")
    OUT_JS.parent.mkdir(parents=True, exist_ok=True)
    js_code = JS_TEMPLATE.replace("__WASM_B64__", b64)
    OUT_JS.write_text(js_code, encoding="utf-8")
    gz_size = len(gzip.compress(js_code.encode("utf-8"), compresslevel=9))
    print(
        f"Wrote {OUT_JS}: raw={len(js_code)/1024:.1f} KB, "
        f"wasm={len(wasm_bytes)/1024:.1f} KB, gzip={gz_size/1024:.1f} KB"
    )


if __name__ == "__main__":
    main()
