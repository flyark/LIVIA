//! High-speed cLIP contact-fingerprint clustering:
//! - Hardware/Wasm `u64::count_ones()` (`popcnt`) binary cosine distance
//! - Packed upper-triangle `f64` matrix (`N*(N-1)/2`, halving memory for N=10,000+)
//! - Exact O(N^2) nearest-neighbor-cached Lance-Williams UPGMA matching `clipfast.js`
//!   and `js/clip-clustering.js` (`EPS = 1e-12` tie-breaking)
//! - Incremental O(N^2 + K^2 * N) silhouette scan for `find_optimal_k`
//! - `fcluster`, `leaf_order`, and size-sorted cluster remapping

/// Linkage row: `(cluster_a, cluster_b, distance, merged_size)`.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LinkageRow {
    pub c1: u32,
    pub c2: u32,
    pub dist: f64,
    pub count: u32,
}

/// Complete output of the cLIP clustering pipeline.
#[derive(Clone, Debug)]
pub struct ClipResult {
    pub optimal_k: usize,
    /// Size-remapped 1-based cluster labels (`1` = largest cluster; ties: lower raw label first).
    pub labels: Vec<u32>,
    /// Raw 1-based `fcluster` labels before size remapping.
    pub raw_labels: Vec<u32>,
    /// Leaf order (`0..n-1`) from in-order dendrogram traversal.
    pub order: Vec<u32>,
    /// `(n - 1)` linkage rows.
    pub linkage: Vec<LinkageRow>,
}

/// Compact upper-triangle symmetric distance matrix (`i < j` stored at `tri_idx(i, j, n)`).
pub struct UpperDistMatrix {
    pub n: usize,
    pub data: Vec<f64>,
}

impl UpperDistMatrix {
    #[inline]
    pub fn with_size(n: usize) -> Self {
        let len = n.saturating_mul(n.saturating_sub(1)) / 2;
        Self {
            n,
            data: vec![0.0_f64; len],
        }
    }

    #[inline]
    pub fn idx(i: usize, j: usize, n: usize) -> usize {
        debug_assert!(i < j && j < n);
        // Row offset for i: i * (2*n - i - 1) / 2 + (j - i - 1)
        i * (2 * n - i - 1) / 2 + (j - i - 1)
    }

    #[inline]
    pub fn get(&self, i: usize, j: usize) -> f64 {
        if i == j {
            0.0
        } else if i < j {
            self.data[Self::idx(i, j, self.n)]
        } else {
            self.data[Self::idx(j, i, self.n)]
        }
    }

    #[inline]
    pub fn set_upper(&mut self, i: usize, j: usize, v: f64) {
        let pos = Self::idx(i, j, self.n);
        self.data[pos] = v;
    }

    pub fn to_full_vec(&self) -> Vec<f64> {
        let mut out = vec![0.0_f64; self.n * self.n];
        for i in 0..self.n {
            for j in (i + 1)..self.n {
                let d = self.data[Self::idx(i, j, self.n)];
                out[i * self.n + j] = d;
                out[j * self.n + i] = d;
            }
        }
        out
    }
}

/// Pack variable-length 1-based (or 0-based non-negative) integer residue fingerprints
/// into `u64` bitset rows `[n * words]`.
/// `offsets` has length `n + 1`, `flat_residues` holds concatenated residue numbers.
pub fn cosine_dist_from_csr(offsets: &[u32], flat_residues: &[u32]) -> UpperDistMatrix {
    let n = offsets.len().saturating_sub(1);
    let mut dmat = UpperDistMatrix::with_size(n);
    if n <= 1 {
        return dmat;
    }
    let mut max_r = 0u32;
    for &r in flat_residues {
        if r > max_r {
            max_r = r;
        }
    }
    let words = ((max_r as usize) / 64) + 1;
    let mut bits = vec![0u64; n * words];
    let mut pop = vec![0u32; n];

    for i in 0..n {
        let start = offsets[i] as usize;
        let end = offsets[i + 1] as usize;
        let row = &mut bits[i * words..(i + 1) * words];
        for &r in &flat_residues[start..end] {
            let idx = r as usize;
            let w = idx >> 6;
            let mask = 1u64 << (idx & 63);
            if (row[w] & mask) == 0 {
                row[w] |= mask;
                pop[i] += 1;
            }
        }
    }

    for i in 0..n {
        let pi = pop[i];
        let row_i = &bits[i * words..(i + 1) * words];
        for j in (i + 1)..n {
            let pj = pop[j];
            let d = if pi == 0 || pj == 0 {
                1.0
            } else {
                let row_j = &bits[j * words..(j + 1) * words];
                let mut dot = 0u32;
                for w in 0..words {
                    dot += (row_i[w] & row_j[w]).count_ones();
                }
                1.0 - (dot as f64) / ((pi as f64) * (pj as f64)).sqrt()
            };
            dmat.set_upper(i, j, d);
        }
    }
    dmat
}

/// Nearest-neighbor-cached Lance-Williams UPGMA average linkage matching `clipfast.js`
/// and `js/clip-clustering.js` **100% identically** (`EPS = 1e-12` tie-breaking:
/// smallest distance, then lowest cluster id, then lowest partner cluster id).
pub fn average_linkage_upgma(dmat: &UpperDistMatrix) -> Vec<LinkageRow> {
    let n = dmat.n;
    if n < 2 {
        return Vec::new();
    }
    // Working upper-triangle mutable copy (plus active diagonal doesn't need storage).
    // Because merged clusters reuse slot `i` (`i < j` or we store `w(a, b)` symmetrically via helper),
    // let's use a compact upper-triangle working buffer `work`: `get(a, b)` reads `work[idx(min(a,b), max(a,b))]`!
    // Wait: `work` is ONLY 400 MB at N = 10,000 (or we can mutate `UpperDistMatrix` in place if owned!).
    average_linkage_in_place(&mut UpperDistMatrix {
        n,
        data: dmat.data.clone(),
    })
}

pub fn average_linkage_in_place(w: &mut UpperDistMatrix) -> Vec<LinkageRow> {
    let n = w.n;
    if n < 2 {
        return Vec::new();
    }
    const EPS: f64 = 1e-12;

    let mut id: Vec<u32> = (0..n as u32).collect();
    let mut size: Vec<u32> = vec![1; n];
    let mut act: Vec<usize> = (0..n).collect();
    let mut nn: Vec<usize> = vec![0; n];
    let mut nd: Vec<f64> = vec![f64::INFINITY; n];

    let pick_nn = |i: usize, act: &[usize], id: &[u32], w: &UpperDistMatrix| -> (usize, f64) {
        let mut best_d = f64::INFINITY;
        let mut best_j = 0usize;
        let mut best_id = u32::MAX;
        for &j in act {
            if j == i {
                continue;
            }
            let d = w.get(i, j);
            let jid = id[j];
            if d < best_d - EPS || ((d - best_d).abs() <= EPS && jid < best_id) {
                best_d = d;
                best_j = j;
                best_id = jid;
            }
        }
        (best_j, best_d)
    };

    for i in 0..n {
        let (bj, bd) = pick_nn(i, &act, &id, w);
        nn[i] = bj;
        nd[i] = bd;
    }

    let mut z = Vec::with_capacity(n - 1);
    for step in 0..(n - 1) {
        let mut best_d = f64::INFINITY;
        let mut bi = 0usize;
        let mut bj = 0usize;
        let mut bi_id = u32::MAX;
        let mut bj_id = u32::MAX;

        for &a in &act {
            let d = nd[a];
            let aid = id[a];
            let bid = id[nn[a]];
            if d < best_d - EPS
                || ((d - best_d).abs() <= EPS
                    && (aid < bi_id || (aid == bi_id && bid < bj_id)))
            {
                best_d = d;
                bi = a;
                bj = nn[a];
                bi_id = aid;
                bj_id = bid;
            }
        }

        let (c1, c2) = if id[bi] < id[bj] {
            (id[bi], id[bj])
        } else {
            (id[bj], id[bi])
        };
        let si = size[bi] as f64;
        let sj = size[bj] as f64;
        let tot = size[bi] + size[bj];
        z.push(LinkageRow {
            c1,
            c2,
            dist: best_d,
            count: tot,
        });

        // Remove `bj` from `act` via swap_remove or order-preserving remove (order of `act` doesn't
        // affect tie-breaks because tie-breaks explicitly compare `id[a]` and `id[j]`!).
        if let Some(pos) = act.iter().position(|&x| x == bj) {
            act.swap_remove(pos);
        }

        let inv_tot = 1.0 / (tot as f64);
        for &k in &act {
            if k == bi {
                continue;
            }
            let new_d = (w.get(bi, k) * si + w.get(bj, k) * sj) * inv_tot;
            if bi < k {
                w.set_upper(bi, k, new_d);
            } else {
                w.set_upper(k, bi, new_d);
            }
        }
        id[bi] = (n + step) as u32;
        size[bi] = tot;

        if act.len() > 1 {
            let (nbj, nbd) = pick_nn(bi, &act, &id, w);
            nn[bi] = nbj;
            nd[bi] = nbd;
            for &k in &act {
                if k == bi {
                    continue;
                }
                if nn[k] == bi || nn[k] == bj {
                    let (nkj, nkd) = pick_nn(k, &act, &id, w);
                    nn[k] = nkj;
                    nd[k] = nkd;
                } else {
                    let d = w.get(k, bi);
                    if d < nd[k] - EPS || ((d - nd[k]).abs() <= EPS && id[bi] < id[nn[k]]) {
                        nn[k] = bi;
                        nd[k] = d;
                    }
                }
            }
        }
    }

    z
}

/// Cut dendrogram `z` into `k` flat clusters (1-based labels `1..=actual_k` in leaf encounter order),
/// matching `CLIPCluster.fcluster(Z, n, k)` identically.
pub fn fcluster(z: &[LinkageRow], n: usize, k: usize) -> Vec<u32> {
    if n == 0 {
        return Vec::new();
    }
    let k = k.max(1);
    if k >= n {
        return (1..=n as u32).collect();
    }
    let total = 2 * n - 1;
    let mut parent: Vec<usize> = (0..total).collect();

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    let num_merges = n.saturating_sub(k).min(z.len());
    for s in 0..num_merges {
        let new_id = n + s;
        let ra = find(&mut parent, z[s].c1 as usize);
        let rb = find(&mut parent, z[s].c2 as usize);
        parent[ra] = new_id;
        parent[rb] = new_id;
    }

    let mut label_map = vec![0u32; total];
    let mut next_label = 1u32;
    let mut out = Vec::with_capacity(n);
    for i in 0..n {
        let r = find(&mut parent, i);
        if label_map[r] == 0 {
            label_map[r] = next_label;
            next_label += 1;
        }
        out.push(label_map[r]);
    }
    out
}

/// Exact mean silhouette score matching `CLIPCluster.silhouette(D, labels)`.
pub fn silhouette(dmat: &UpperDistMatrix, labels: &[u32]) -> f64 {
    let n = labels.len();
    if n < 2 {
        return 0.0;
    }
    let mut max_lbl = 0usize;
    for &l in labels {
        if (l as usize) > max_lbl {
            max_lbl = l as usize;
        }
    }
    // Map arbitrary label integers to compact `0..K-1`.
    let mut slot_of = vec![usize::MAX; max_lbl + 1];
    let mut counts: Vec<usize> = Vec::new();
    let mut cl = Vec::with_capacity(n);
    for &l in labels {
        let li = l as usize;
        if slot_of[li] == usize::MAX {
            slot_of[li] = counts.len();
            counts.push(0);
        }
        let s = slot_of[li];
        counts[s] += 1;
        cl.push(s);
    }
    let k = counts.len();
    if k < 2 {
        return 0.0;
    }

    let mut sums = vec![0.0_f64; n * k];
    for i in 0..n {
        let ci = cl[i];
        for j in (i + 1)..n {
            let d = dmat.get(i, j);
            sums[i * k + cl[j]] += d;
            sums[j * k + ci] += d;
        }
    }

    let mut total = 0.0_f64;
    for i in 0..n {
        let own = cl[i];
        let own_cnt = counts[own];
        if own_cnt <= 1 {
            continue;
        }
        let a = sums[i * k + own] / ((own_cnt - 1) as f64);
        let mut b = f64::INFINITY;
        let row = &sums[i * k..(i + 1) * k];
        for c in 0..k {
            if c == own || counts[c] == 0 {
                continue;
            }
            let mean_d = row[c] / (counts[c] as f64);
            if mean_d < b {
                b = mean_d;
            }
        }
        let denom = a.max(b);
        if denom > 0.0 {
            total += (b - a) / denom;
        }
    }
    total / (n as f64)
}

/// Incremental `O(N^2 + K^2 * N)` optimal-K finder matching
/// `CLIPCluster.findOptimalK(Z, D, n, maxK = 30, dropThreshold = 0.05)`.
/// Scans `k = 2..=min(max_k, n - 1)` using a single triangle pass over `dmat`!
pub fn find_optimal_k(
    z: &[LinkageRow],
    dmat: &UpperDistMatrix,
    max_k: usize,
    drop_threshold: f64,
) -> usize {
    let n = dmat.n;
    if n <= 2 {
        return n.min(2);
    }
    let cap = max_k.min(n - 1);
    if cap < 2 {
        return 2;
    }

    // Step 1: Build Union-Find up to `n - cap` merges (yielding `cap` clusters).
    let total_nodes = 2 * n - 1;
    let mut parent: Vec<usize> = (0..total_nodes).collect();
    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    let base_merges = n - cap;
    for s in 0..base_merges {
        let new_id = n + s;
        let ra = find(&mut parent, z[s].c1 as usize);
        let rb = find(&mut parent, z[s].c2 as usize);
        parent[ra] = new_id;
        parent[rb] = new_id;
    }

    // Assign compact slot `0..cap-1` to each active root at `k = cap`.
    let mut root_to_slot = vec![usize::MAX; total_nodes];
    let mut next_slot = 0usize;
    let mut cl = vec![0usize; n];
    let mut counts = vec![0usize; cap];
    for i in 0..n {
        let r = find(&mut parent, i);
        if root_to_slot[r] == usize::MAX {
            root_to_slot[r] = next_slot;
            next_slot += 1;
        }
        let sl = root_to_slot[r];
        cl[i] = sl;
        counts[sl] += 1;
    }
    let active_slots = next_slot;

    // Step 2: Single upper-triangle pass over `dmat` (`O(N^2 / 2)` reads total!).
    let mut sums = vec![0.0_f64; n * active_slots];
    for i in 0..n {
        let ci = cl[i];
        for j in (i + 1)..n {
            let d = dmat.get(i, j);
            sums[i * active_slots + cl[j]] += d;
            sums[j * active_slots + ci] += d;
        }
    }

    // Step 3: Evaluate silhouette at `k = cap`, then merge one step at a time down to `k = 2`.
    let eval_sil = |cl: &[usize], counts: &[usize], sums: &[f64]| -> Option<f64> {
        let active_cnt = counts.iter().filter(|&&c| c > 0).count();
        if active_cnt < 2 {
            return None;
        }
        let mut total = 0.0_f64;
        for i in 0..n {
            let own = cl[i];
            let own_c = counts[own];
            if own_c <= 1 {
                continue;
            }
            let row = &sums[i * active_slots..(i + 1) * active_slots];
            let a = row[own] / ((own_c - 1) as f64);
            let mut b = f64::INFINITY;
            for c in 0..active_slots {
                let cc = counts[c];
                if c == own || cc == 0 {
                    continue;
                }
                let mean_d = row[c] / (cc as f64);
                if mean_d < b {
                    b = mean_d;
                }
            }
            let denom = a.max(b);
            if denom > 0.0 {
                total += (b - a) / denom;
            }
        }
        Some(total / (n as f64))
    };

    let mut scores_by_k = vec![None; cap + 1];
    scores_by_k[cap] = eval_sil(&cl, &counts, &sums);

    for k_cur in (3..=cap).rev() {
        let merge_idx = n - k_cur;
        let ra = find(&mut parent, z[merge_idx].c1 as usize);
        let rb = find(&mut parent, z[merge_idx].c2 as usize);
        let new_id = n + merge_idx;
        let sa = root_to_slot[ra];
        let sb = root_to_slot[rb];
        parent[ra] = new_id;
        parent[rb] = new_id;
        root_to_slot[new_id] = sa;

        if sa != sb {
            for i in 0..n {
                let base = i * active_slots;
                sums[base + sa] += sums[base + sb];
                if cl[i] == sb {
                    cl[i] = sa;
                }
            }
            counts[sa] += counts[sb];
            counts[sb] = 0;
        }
        scores_by_k[k_cur - 1] = eval_sil(&cl, &counts, &sums);
    }

    // Phase 1 + Phase 2 matching `CLIPCluster.findOptimalK` in `js/clip-clustering.js`.
    let mut best_k = 2usize;
    let mut best_score = -1.0_f64;
    let mut k_scores: Vec<(usize, f64)> = Vec::with_capacity(cap);
    for k in 2..=cap {
        if let Some(s) = scores_by_k[k] {
            k_scores.push((k, s));
            if s > best_score {
                best_score = s;
                best_k = k;
            }
        }
    }

    for idx in 1..k_scores.len() {
        if k_scores[idx - 1].1 - k_scores[idx].1 > drop_threshold {
            best_k = k_scores[idx - 1].0;
            break;
        }
    }
    best_k
}

/// Iterative in-order dendrogram leaf order matching `CLIPCluster.leafOrder(Z, n)`.
pub fn leaf_order(z: &[LinkageRow], n: usize) -> Vec<u32> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 || z.is_empty() {
        return (0..n as u32).collect();
    }
    let mut order = Vec::with_capacity(n);
    let mut stack = vec![(2 * n - 2) as u32];
    while let Some(node) = stack.pop() {
        if (node as usize) < n {
            order.push(node);
        } else {
            let row = &z[(node as usize) - n];
            stack.push(row.c2);
            stack.push(row.c1);
        }
    }
    order
}

/// Renumber 1-based cluster labels by descending size (`1` = largest; ties: lower raw label first),
/// matching `clip.html` and `clipworker.js`.
pub fn remap_labels_by_size(raw: &[u32]) -> Vec<u32> {
    if raw.is_empty() {
        return Vec::new();
    }
    let max_l = raw.iter().copied().max().unwrap_or(0) as usize;
    let mut cnt = vec![0usize; max_l + 1];
    for &l in raw {
        cnt[l as usize] += 1;
    }
    let mut present: Vec<usize> = (1..=max_l).filter(|&l| cnt[l] > 0).collect();
    present.sort_by(|&x, &y| cnt[y].cmp(&cnt[x]).then_with(|| x.cmp(&y)));
    let mut remap = vec![0u32; max_l + 1];
    for (idx, old) in present.into_iter().enumerate() {
        remap[old] = (idx + 1) as u32;
    }
    raw.iter().map(|&l| remap[l as usize]).collect()
}

/// End-to-end cLIP pipeline from CSR fingerprints (`offsets`, `flat_residues`) in a single memory space
/// (never materializes an `N x N` JS array on the heap!).
pub fn run_clip_from_csr(
    offsets: &[u32],
    flat_residues: &[u32],
    max_k: usize,
    drop_threshold: f64,
) -> ClipResult {
    let n = offsets.len().saturating_sub(1);
    if n == 0 {
        return ClipResult {
            optimal_k: 0,
            labels: Vec::new(),
            raw_labels: Vec::new(),
            order: Vec::new(),
            linkage: Vec::new(),
        };
    }
    if n == 1 {
        return ClipResult {
            optimal_k: 1,
            labels: vec![1],
            raw_labels: vec![1],
            order: vec![0],
            linkage: Vec::new(),
        };
    }
    let mut dmat = cosine_dist_from_csr(offsets, flat_residues);
    // Notice: `find_optimal_k` needs original `dmat` while `average_linkage_in_place` mutates its working matrix.
    // Cloning `UpperDistMatrix` temporarily uses `N*(N-1)/2 * 8` bytes during linkage, then drops immediately.
    let linkage = average_linkage_upgma(&dmat);
    let optimal_k = find_optimal_k(&linkage, &dmat, max_k, drop_threshold);
    // Free `dmat` early before building outputs.
    dmat.data = Vec::new();
    let raw_labels = fcluster(&linkage, n, optimal_k);
    let labels = remap_labels_by_size(&raw_labels);
    let order = leaf_order(&linkage, n);
    ClipResult {
        optimal_k,
        labels,
        raw_labels,
        order,
        linkage,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_incremental_silhouette_matches_exact() {
        // Synthetic 60-prediction protein with 4 distinct binding patches + noisy overlaps.
        let mut offsets = vec![0u32];
        let mut flat = Vec::new();
        for i in 0..60u32 {
            let base = (i % 4) * 25 + 1;
            for r in base..(base + 10 + (i % 5)) {
                flat.push(r);
            }
            if i % 7 == 0 {
                flat.push(15);
                flat.push(40);
            }
            offsets.push(flat.len() as u32);
        }
        let dmat = cosine_dist_from_csr(&offsets, &flat);
        let z = average_linkage_upgma(&dmat);
        let fast_k = find_optimal_k(&z, &dmat, 20, 0.05);

        let mut brute_best_k = 2;
        let mut brute_best_s = -1.0;
        let mut k_scores = Vec::new();
        for k in 2..=20 {
            let lbl = fcluster(&z, 60, k);
            let s = silhouette(&dmat, &lbl);
            k_scores.push((k, s));
            if s > brute_best_s {
                brute_best_s = s;
                brute_best_k = k;
            }
        }
        for idx in 1..k_scores.len() {
            if k_scores[idx - 1].1 - k_scores[idx].1 > 0.05 {
                brute_best_k = k_scores[idx - 1].0;
                break;
            }
        }
        assert_eq!(fast_k, brute_best_k);
        assert_eq!(fast_k, 4);
    }
}
