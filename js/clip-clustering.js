/*
 * cLIP clustering — contact-fingerprint hierarchical clustering.
 * Contact-fingerprint hierarchical clustering: cosine distance + average
 * linkage (UPGMA) + silhouette-based optimal-k, matching scipy/sklearn.
 *
 * Input: fingerprints = array of arrays of 1-indexed residue positions
 * (the target-chain cLIR_indice for each partner). protein_len is unused
 * for cosine on binary vectors (kept for API parity).
 */
(function (root, factory) {
  if (typeof module !== 'undefined' && module.exports) module.exports = factory();
  else root.CLIPCluster = factory();
})(typeof self !== 'undefined' ? self : this, function () {
  'use strict';

  // Cosine distance between binary vectors, computed from residue-index sets:
  //   cos_sim = |A∩B| / (sqrt|A| * sqrt|B|)   (since ||u|| = sqrt(#ones))
  //   cos_dist = 1 - cos_sim
  function cosineDistMatrix(fingerprints) {
    const n = fingerprints.length;
    const sets = fingerprints.map((f) => new Set(f));
    const norm = fingerprints.map((f) => Math.sqrt(f.length));
    const D = Array.from({ length: n }, () => new Float64Array(n));
    for (let i = 0; i < n; i++) {
      for (let j = i + 1; j < n; j++) {
        let inter = 0;
        const a = sets[i];
        const bf = fingerprints[j];
        for (let t = 0; t < bf.length; t++) if (a.has(bf[t])) inter++;
        let d = norm[i] === 0 || norm[j] === 0 ? 1 : 1 - inter / (norm[i] * norm[j]);
        if (d < 0) d = 0;                // clamp tiny negatives
        D[i][j] = d;
        D[j][i] = d;
      }
    }
    return D;
  }

  // Average-linkage (UPGMA) agglomeration via Lance-Williams update.
  // Returns Z: [[idA, idB, dist, size], ...] with merged ids n, n+1, ...
  function averageLinkage(D0, n) {
    const size = {};
    const D = {};
    const active = new Set();
    for (let i = 0; i < n; i++) { size[i] = 1; D[i] = {}; active.add(i); }
    for (let i = 0; i < n; i++)
      for (let j = i + 1; j < n; j++) { D[i][j] = D0[i][j]; D[j][i] = D0[i][j]; }

    const Z = [];
    let nextId = n;
    while (active.size > 1) {
      const ids = [...active].sort((a, b) => a - b); // deterministic tie-break
      let best = Infinity, bi = -1, bj = -1;
      for (let x = 0; x < ids.length; x++) {
        const i = ids[x];
        const Di = D[i];
        for (let y = x + 1; y < ids.length; y++) {
          const j = ids[y];
          const d = Di[j];
          if (d < best - 1e-12) { best = d; bi = i; bj = j; }
        }
      }
      const ni = size[bi], nj = size[bj];
      const id = nextId++;
      size[id] = ni + nj;
      D[id] = {};
      active.delete(bi);
      active.delete(bj);
      for (const k of active) {
        const dn = (ni * D[bi][k] + nj * D[bj][k]) / (ni + nj);
        D[id][k] = dn;
        D[k][id] = dn;
      }
      active.add(id);
      Z.push([bi, bj, best, ni + nj]);
    }
    return Z;
  }

  // Flatten to k clusters by undoing the last k-1 merges (UPGMA is monotonic,
  // so cutting after the first n-k merges yields exactly k groups). 1-indexed labels.
  function fcluster(Z, n, k) {
    const uf = {};
    const find = (x) => { let r = x; while (uf[r] !== undefined) r = uf[r]; return r; };
    const nMerges = n - k;
    for (let m = 0; m < nMerges; m++) {
      const a = Z[m][0], b = Z[m][1];
      const ra = find(a), rb = find(b);
      if (ra !== rb) uf[ra] = rb;
      uf[n + m] = find(b);            // merged-cluster id → current root
    }
    const labels = new Array(n);
    const rootLabel = {};
    let next = 1;
    for (let i = 0; i < n; i++) {
      const r = find(i);
      if (rootLabel[r] === undefined) rootLabel[r] = next++;
      labels[i] = rootLabel[r];
    }
    return labels;
  }

  // Mean silhouette (cosine distances = precomputed D). Singleton sample → 0.
  function silhouette(D, labels) {
    const n = labels.length;
    const clusters = {};
    for (let i = 0; i < n; i++) (clusters[labels[i]] || (clusters[labels[i]] = [])).push(i);
    const cids = Object.keys(clusters);
    if (cids.length < 2) return -1;
    let total = 0;
    for (let i = 0; i < n; i++) {
      const own = clusters[labels[i]];
      let a = 0;
      if (own.length > 1) {
        let s = 0;
        for (const j of own) if (j !== i) s += D[i][j];
        a = s / (own.length - 1);
      }
      let b = Infinity;
      for (const cid of cids) {
        if (+cid === labels[i]) continue;
        const grp = clusters[cid];
        let s = 0;
        for (const j of grp) s += D[i][j];
        const mean = s / grp.length;
        if (mean < b) b = mean;
      }
      total += own.length === 1 ? 0 : (b - a) / Math.max(a, b);
    }
    return total / n;
  }

  // Dendrogram leaf order (left→right) from the merge tree, for heatmap row ordering.
  function leafOrder(Z, n) {
    if (n <= 1) return n ? [0] : [];
    const child = {};
    for (let m = 0; m < Z.length; m++) child[n + m] = [Z[m][0], Z[m][1]];
    const order = [];
    const stack = [n + Z.length - 1];
    while (stack.length) {
      const id = stack.pop();
      if (id < n) { order.push(id); continue; }
      const c = child[id];
      stack.push(c[1], c[0]); // push right then left → left popped first
    }
    return order;
  }

  function findOptimalK(Z, D, n, maxK, dropThreshold) {
    maxK = maxK || 30;
    dropThreshold = dropThreshold == null ? 0.05 : dropThreshold;
    if (n < 3) return Math.min(2, n);
    maxK = Math.min(maxK, n - 1);
    if (maxK < 2) return 2;
    // Phase 1: silhouette for every k, tracking the argmax.
    const kScores = [];
    let bestScore = -1, optimalK = 2;
    for (let k = 2; k <= maxK; k++) {
      const labels = fcluster(Z, n, k);
      if (new Set(labels).size < 2) continue;
      const score = silhouette(D, labels);
      kScores.push([k, score]);
      if (score > bestScore) { bestScore = score; optimalK = k; }
    }
    // Phase 2: back off to the k just before the first significant consecutive drop ("elbow").
    for (let i = 1; i < kScores.length; i++) {
      if (kScores[i - 1][1] - kScores[i][1] > dropThreshold) { optimalK = kScores[i - 1][0]; break; }
    }
    return optimalK;
  }

  // Full pipeline: fingerprints (residue-index arrays) → {optimalK, labels, Z, D}.
  function clusterInteractions(fingerprints, opts) {
    opts = opts || {};
    const n = fingerprints.length;
    if (n < 2) return { optimalK: n, labels: n ? [1] : [], Z: [], D: [] };
    const D = cosineDistMatrix(fingerprints);
    const Z = averageLinkage(D, n);
    const optimalK = findOptimalK(Z, D, n, opts.maxK || 30, opts.dropThreshold == null ? 0.05 : opts.dropThreshold);
    const labels = fcluster(Z, n, optimalK);
    return { optimalK, labels, Z, D };
  }

  return { cosineDistMatrix, averageLinkage, fcluster, silhouette, findOptimalK, clusterInteractions, leafOrder };
});
