//! Deterministic community detection (Leiden with Modularity and CPM objectives),
//! prize-collecting / Takahashi-Matsuyama Steiner connector, and force-directed
//! graph layouts (Fruchterman-Reingold, Kamada-Kawai, and DrL/OpenOrd).

use std::cmp::Ordering;
use std::collections::{HashMap, VecDeque};

/// Exact Mulberry32 PRNG matching LIVIA Atlas `seededRandom(seed)` bit-for-bit.
#[derive(Clone, Copy)]
pub struct Mulberry32 {
    state: u32,
}

impl Mulberry32 {
    pub fn new(seed: u32) -> Self {
        Self { state: seed }
    }

    #[inline]
    pub fn next_f64(&mut self) -> f64 {
        self.state = self.state.wrapping_add(0x6D2B79F5);
        let mut t = (self.state ^ (self.state >> 15)).wrapping_mul(1 | self.state);
        t = (t.wrapping_add((t ^ (t >> 7)).wrapping_mul(61 | t))) ^ t;
        let out = t ^ (t >> 14);
        (out as f64) / 4294967296.0
    }
}

/// Fisher-Yates shuffle matching `app.js` `shuffled(n, rand)`.
/// When `rng` is `None`, returns `[0, 1, ..., n - 1]`.
pub fn shuffled(n: usize, rng: Option<&mut Mulberry32>) -> Vec<usize> {
    let mut o: Vec<usize> = (0..n).collect();
    if let Some(r) = rng {
        for i in (1..n).rev() {
            let j = (r.next_f64() * ((i + 1) as f64)).floor() as usize;
            let j = j.min(i);
            o.swap(i, j);
        }
    }
    o
}

/// Exact bit-for-bit match to `livia-atlas/app.js` `leiden(n, edges, gamma, rand)`.
/// `seed == 0` uses fixed index order (`rand = null`); `seed > 0` uses `Mulberry32(seed)`.
pub fn leiden_atlas(n: usize, edges: &[(usize, usize, f64)], gamma: f64, seed: u32) -> Vec<u32> {
    if n == 0 {
        return Vec::new();
    }
    let mut rng_opt = if seed > 0 { Some(Mulberry32::new(seed)) } else { None };
    let mut cur_n = n;
    let mut cur_edges: Vec<(usize, usize, f64)> = edges
        .iter()
        .copied()
        .filter(|&(u, v, _)| u < n && v < n)
        .collect();
    let mut p: Vec<usize> = (0..n).collect();
    let mut member: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();

    for _level in 0..20 {
        // Preserve insertion order in adjacency lists just like JS `Map`.
        let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); cur_n];
        let mut adj_pos: Vec<HashMap<usize, usize>> = vec![HashMap::new(); cur_n];
        let mut k = vec![0.0_f64; cur_n];
        let mut m2 = 0.0_f64;

        for &(u, v, w) in &cur_edges {
            if let Some(&idx) = adj_pos[u].get(&v) {
                adj[u][idx].1 += w;
            } else {
                adj_pos[u].insert(v, adj[u].len());
                adj[u].push((v, w));
            }
            if u != v {
                if let Some(&idx) = adj_pos[v].get(&u) {
                    adj[v][idx].1 += w;
                } else {
                    adj_pos[v].insert(u, adj[v].len());
                    adj[v].push((u, w));
                }
            } else {
                let idx = adj_pos[u][&u];
                adj[u][idx].1 += w;
            }
            k[u] += w;
            k[v] += w;
            m2 += 2.0 * w;
        }
        if m2 == 0.0 {
            break;
        }

        let mut tot: HashMap<usize, f64> = HashMap::new();
        for i in 0..cur_n {
            *tot.entry(p[i]).or_insert(0.0) += k[i];
        }

        let ord = shuffled(cur_n, rng_opt.as_mut());
        let mut moved = true;
        let mut rounds = 0usize;
        while moved && rounds < 50 {
            rounds += 1;
            moved = false;
            for &i in &ord {
                let ci = p[i];
                let mut wt_keys: Vec<usize> = Vec::new();
                let mut wt_map: HashMap<usize, usize> = HashMap::new();
                let mut wt_vals: Vec<f64> = Vec::new();
                for &(j, w) in &adj[i] {
                    if j != i {
                        let cj = p[j];
                        if let Some(&pos) = wt_map.get(&cj) {
                            wt_vals[pos] += w;
                        } else {
                            wt_map.insert(cj, wt_keys.len());
                            wt_keys.push(cj);
                            wt_vals.push(w);
                        }
                    }
                }
                let tot_ci_after = tot.get(&ci).copied().unwrap_or(0.0) - k[i];
                tot.insert(ci, tot_ci_after);
                let wt_ci = wt_map.get(&ci).map(|&pos| wt_vals[pos]).unwrap_or(0.0);
                let mut best = ci;
                let mut gain = wt_ci - gamma * tot_ci_after * k[i] / m2;
                for (pos, &cj) in wt_keys.iter().enumerate() {
                    let w = wt_vals[pos];
                    let tot_cj = tot.get(&cj).copied().unwrap_or(0.0);
                    let g = w - gamma * tot_cj * k[i] / m2;
                    if g > gain + 1e-12 {
                        gain = g;
                        best = cj;
                    }
                }
                *tot.entry(best).or_insert(0.0) += k[i];
                if best != ci {
                    p[i] = best;
                    moved = true;
                }
            }
        }

        // Preserve community encounter order like JS `Map`.
        let mut comm_order: Vec<usize> = Vec::new();
        let mut comms: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..cur_n {
            comms.entry(p[i]).or_insert_with(|| {
                comm_order.push(p[i]);
                Vec::new()
            }).push(i);
        }
        if comm_order.len() == cur_n {
            break;
        }

        let mut r: Vec<usize> = (0..cur_n).collect();
        let mut rtot = k.clone();
        let mut rsize = vec![1usize; cur_n];
        let mut s_out = vec![0.0_f64; cur_n];
        let mut in_c = vec![false; cur_n];

        for cid in &comm_order {
            let nodes = &comms[cid];
            let mut tot_c = 0.0_f64;
            for &i in nodes {
                in_c[i] = true;
                tot_c += k[i];
            }
            let mut w_c = vec![0.0_f64; cur_n];
            for &i in nodes {
                let mut s = 0.0_f64;
                for &(j, w) in &adj[i] {
                    if j != i && in_c[j] {
                        s += w;
                    }
                }
                w_c[i] = s;
                s_out[i] = s;
            }
            let visit: Vec<usize> = if let Some(ref mut rng) = rng_opt {
                shuffled(nodes.len(), Some(rng)).into_iter().map(|x| nodes[x]).collect()
            } else {
                nodes.clone()
            };
            for v in visit {
                if rsize[r[v]] != 1 || r[v] != v {
                    continue;
                }
                if w_c[v] < gamma * k[v] * (tot_c - k[v]) / m2 {
                    continue;
                }
                let mut ws_keys: Vec<usize> = Vec::new();
                let mut ws_map: HashMap<usize, usize> = HashMap::new();
                let mut ws_vals: Vec<f64> = Vec::new();
                for &(j, w) in &adj[v] {
                    if j != v && in_c[j] {
                        let rj = r[j];
                        if let Some(&pos) = ws_map.get(&rj) {
                            ws_vals[pos] += w;
                        } else {
                            ws_map.insert(rj, ws_keys.len());
                            ws_keys.push(rj);
                            ws_vals.push(w);
                        }
                    }
                }
                let mut best = v;
                let mut gain = 0.0_f64;
                for (pos, &s) in ws_keys.iter().enumerate() {
                    let w = ws_vals[pos];
                    if s_out[s] < gamma * rtot[s] * (tot_c - rtot[s]) / m2 {
                        continue;
                    }
                    let g = w - gamma * k[v] * rtot[s] / m2;
                    if g > gain + 1e-12 {
                        gain = g;
                        best = s;
                    }
                }
                if best != v {
                    let w = ws_vals[ws_map[&best]];
                    r[v] = best;
                    rtot[best] += k[v];
                    rsize[best] += 1;
                    rsize[v] = 0;
                    s_out[best] = s_out[best] + w_c[v] - 2.0 * w;
                }
            }
            for &i in nodes {
                in_c[i] = false;
            }
        }

        let mut ids: HashMap<usize, usize> = HashMap::new();
        for &x in &r {
            let next_id = ids.len();
            ids.entry(x).or_insert(next_id);
        }
        let n_sub = ids.len();
        if n_sub == cur_n {
            break;
        }

        let mut next_member: Vec<Vec<usize>> = vec![Vec::new(); n_sub];
        let mut next_p = vec![0usize; n_sub];
        let mut pid: HashMap<usize, usize> = HashMap::new();
        for i in 0..cur_n {
            let a = ids[&r[i]];
            next_member[a].extend_from_slice(&member[i]);
            let next_pid = pid.len();
            let pv = *pid.entry(p[i]).or_insert(next_pid);
            next_p[a] = pv;
        }

        let mut agg_order: Vec<(usize, usize)> = Vec::new();
        let mut agg_map: HashMap<(usize, usize), usize> = HashMap::new();
        let mut agg_w: Vec<f64> = Vec::new();
        for &(u, v, w) in &cur_edges {
            let a = ids[&r[u]];
            let b = ids[&r[v]];
            let key = if a < b { (a, b) } else { (b, a) };
            if let Some(&pos) = agg_map.get(&key) {
                agg_w[pos] += w;
            } else {
                agg_map.insert(key, agg_order.len());
                agg_order.push(key);
                agg_w.push(w);
            }
        }
        cur_edges = agg_order
            .into_iter()
            .zip(agg_w)
            .map(|((a, b), w)| (a, b, w))
            .collect();
        cur_n = n_sub;
        member = next_member;
        p = next_p;
    }

    let mut group_order: Vec<usize> = Vec::new();
    let mut groups: HashMap<usize, Vec<usize>> = HashMap::new();
    for i in 0..cur_n {
        groups
            .entry(p[i])
            .or_insert_with(|| {
                group_order.push(p[i]);
                Vec::new()
            })
            .extend_from_slice(&member[i]);
    }
    let mut group_list: Vec<Vec<usize>> = group_order
        .into_iter()
        .map(|gid| groups.remove(&gid).unwrap())
        .collect();
    // Stable sort by descending size (matches JS Array.prototype.sort `(a, b) => b.length - a.length`).
    group_list.sort_by(|a, b| b.len().cmp(&a.len()));

    let mut out = vec![0u32; n];
    for (rank, g) in group_list.iter().enumerate() {
        for &x in g {
            out[x] = rank as u32;
        }
    }
    out
}

/// Full multi-iteration Leiden community detection for `LIVIA/network.html`
/// (replaces Pyodide + `python-igraph` `community_leiden`).
///
/// Supports both:
/// - `use_cpm = false` (`"modularity"`): Reichardt-Bornholdt modularity with resolution `gamma`.
/// - `use_cpm = true` (`"CPM"`): Constant Potts Model `sum_c (e_c - gamma * n_c * (n_c - 1) / 2)`.
///
/// Runs `n_iterations` outer Leiden passes (starting each iteration from the previous partition
/// or exploring seeded multi-start moves) and guarantees internally connected communities.
/// Returns `(membership, degree)` with communities numbered `0..K-1` largest-first.
pub fn leiden_network(
    n: usize,
    edges: &[(usize, usize, f64)],
    gamma: f64,
    use_cpm: bool,
    n_iterations: usize,
    seed: u32,
) -> (Vec<u32>, Vec<u32>) {
    if n == 0 {
        return (Vec::new(), Vec::new());
    }
    // Simple unweighted degree of unique incident edges (matches `igraph.Graph.degree()`).
    let mut degree = vec![0u32; n];
    let mut valid_edges: Vec<(usize, usize, f64)> = Vec::with_capacity(edges.len());
    for &(u, v, w) in edges {
        if u < n && v < n {
            degree[u] += 1;
            degree[v] += 1;
            valid_edges.push((u, v, w));
        }
    }
    if valid_edges.is_empty() {
        return ((0..n as u32).collect(), degree);
    }

    let mut rng = Mulberry32::new(if seed == 0 { 1 } else { seed });
    let iters = n_iterations.max(1);

    let mut part: Vec<usize> = (0..n).collect();
    let mut best_part = part.clone();
    let mut best_q = f64::NEG_INFINITY;

    for outer in 0..iters {
        // Multi-level Leiden pass starting from `part` on iteration > 0 (or fresh singletons on restart).
        if outer > 0 && outer % 4 == 0 {
            // Periodic seeded restart escapes shallow local optima while keeping best_part.
            part = (0..n).collect();
        }
        part = leiden_one_pass(n, &valid_edges, &part, gamma, use_cpm, &mut rng);
        part = split_disconnected(n, &valid_edges, &part);
        let q = quality_score(n, &valid_edges, &part, gamma, use_cpm);
        if q > best_q + 1e-12 {
            best_q = q;
            best_part = part.clone();
        } else if (q - best_q).abs() <= 1e-12 {
            // Converged on current trajectory.
            part = best_part.clone();
        }
    }

    best_part = split_disconnected(n, &valid_edges, &best_part);

    // Remap communities by descending size (ties: first appearance order, matching Python Counter.most_common()).
    let mut first_seen: Vec<usize> = Vec::new();
    let mut counts: HashMap<usize, usize> = HashMap::new();
    for &c in &best_part {
        counts
            .entry(c)
            .and_modify(|cnt| *cnt += 1)
            .or_insert_with(|| {
                first_seen.push(c);
                1
            });
    }
    first_seen.sort_by(|a, b| counts[b].cmp(&counts[a]));
    let mut remap: HashMap<usize, u32> = HashMap::with_capacity(first_seen.len());
    for (idx, c) in first_seen.into_iter().enumerate() {
        remap.insert(c, idx as u32);
    }
    let membership: Vec<u32> = best_part.iter().map(|c| remap[c]).collect();
    (membership, degree)
}

fn quality_score(
    n: usize,
    edges: &[(usize, usize, f64)],
    part: &[usize],
    gamma: f64,
    use_cpm: bool,
) -> f64 {
    let mut k = vec![0.0_f64; n];
    let mut m2 = 0.0_f64;
    let mut in_w: HashMap<usize, f64> = HashMap::new();
    for &(u, v, w) in edges {
        k[u] += w;
        k[v] += w;
        m2 += 2.0 * w;
        if part[u] == part[v] {
            *in_w.entry(part[u]).or_insert(0.0) += 2.0 * w;
        }
    }
    if m2 <= 0.0 {
        return 0.0;
    }
    if use_cpm {
        let mut sz: HashMap<usize, f64> = HashMap::new();
        for &c in part {
            *sz.entry(c).or_insert(0.0) += 1.0;
        }
        let mut q = 0.0;
        for (c, s) in sz {
            let ew = in_w.get(&c).copied().unwrap_or(0.0) * 0.5;
            q += ew - gamma * s * (s - 1.0) * 0.5;
        }
        q
    } else {
        let mut tot: HashMap<usize, f64> = HashMap::new();
        for i in 0..n {
            *tot.entry(part[i]).or_insert(0.0) += k[i];
        }
        let mut q = 0.0;
        for (c, t) in tot {
            let ew = in_w.get(&c).copied().unwrap_or(0.0);
            let norm = t / m2;
            q += ew / m2 - gamma * norm * norm;
        }
        q
    }
}

fn split_disconnected(n: usize, edges: &[(usize, usize, f64)], part: &[usize]) -> Vec<usize> {
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(u, v, _) in edges {
        if u != v && part[u] == part[v] {
            adj[u].push(v);
            adj[v].push(u);
        }
    }
    let mut out = vec![usize::MAX; n];
    let mut next_cid = 0usize;
    let mut q = VecDeque::new();
    for s in 0..n {
        if out[s] != usize::MAX {
            continue;
        }
        out[s] = next_cid;
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            for &v in &adj[u] {
                if out[v] == usize::MAX {
                    out[v] = next_cid;
                    q.push_back(v);
                }
            }
        }
        next_cid += 1;
    }
    out
}

fn leiden_one_pass(
    n: usize,
    edges: &[(usize, usize, f64)],
    init_part: &[usize],
    gamma: f64,
    use_cpm: bool,
    rng: &mut Mulberry32,
) -> Vec<usize> {
    let mut cur_n = n;
    let mut cur_edges = edges.to_vec();
    let mut p = init_part.to_vec();
    let mut node_size = vec![1.0_f64; n];
    let mut member: Vec<Vec<usize>> = (0..n).map(|i| vec![i]).collect();

    for _level in 0..20 {
        let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); cur_n];
        let mut adj_pos: Vec<HashMap<usize, usize>> = vec![HashMap::new(); cur_n];
        let mut k = vec![0.0_f64; cur_n];
        let mut m2 = 0.0_f64;

        for &(u, v, w) in &cur_edges {
            if u != v {
                if let Some(&idx) = adj_pos[u].get(&v) {
                    adj[u][idx].1 += w;
                } else {
                    adj_pos[u].insert(v, adj[u].len());
                    adj[u].push((v, w));
                }
                if let Some(&idx) = adj_pos[v].get(&u) {
                    adj[v][idx].1 += w;
                } else {
                    adj_pos[v].insert(u, adj[v].len());
                    adj[v].push((u, w));
                }
            }
            k[u] += w;
            k[v] += w;
            m2 += 2.0 * w;
        }
        if m2 <= 0.0 {
            break;
        }

        // Penalty weight of each node: degree k[i] for modularity, node_size[i] for CPM;
        // scale = gamma / m2 for modularity, gamma for CPM.
        let scale = if use_cpm { gamma } else { gamma / m2 };
        let nw: Vec<f64> = if use_cpm { node_size.clone() } else { k.clone() };

        let mut tot: HashMap<usize, f64> = HashMap::new();
        for i in 0..cur_n {
            *tot.entry(p[i]).or_insert(0.0) += nw[i];
        }

        // Fast local move queue (Traag et al. 2019 queue-based local move).
        let ord = shuffled(cur_n, Some(rng));
        let mut in_queue = vec![true; cur_n];
        let mut queue: VecDeque<usize> = ord.into_iter().collect();
        let mut steps = 0usize;
        let max_steps = cur_n.saturating_mul(60).max(100);

        while let Some(i) = queue.pop_front() {
            in_queue[i] = false;
            steps += 1;
            if steps > max_steps {
                break;
            }
            let ci = p[i];
            let mut wt_keys: Vec<usize> = Vec::new();
            let mut wt_map: HashMap<usize, usize> = HashMap::new();
            let mut wt_vals: Vec<f64> = Vec::new();
            for &(j, w) in &adj[i] {
                let cj = p[j];
                if let Some(&pos) = wt_map.get(&cj) {
                    wt_vals[pos] += w;
                } else {
                    wt_map.insert(cj, wt_keys.len());
                    wt_keys.push(cj);
                    wt_vals.push(w);
                }
            }
            let tot_ci_after = tot.get(&ci).copied().unwrap_or(0.0) - nw[i];
            tot.insert(ci, tot_ci_after);
            let wt_ci = wt_map.get(&ci).map(|&pos| wt_vals[pos]).unwrap_or(0.0);

            // Also allow moving to an empty singleton community (gain = 0.0).
            let stay_gain = wt_ci - scale * tot_ci_after * nw[i];
            let mut best = ci;
            let mut best_gain = stay_gain;
            let mut use_empty = false;
            if 0.0 > best_gain + 1e-12 && tot_ci_after > 0.0 {
                best_gain = 0.0;
                use_empty = true;
            }

            for (pos, &cj) in wt_keys.iter().enumerate() {
                let w = wt_vals[pos];
                let tot_cj = tot.get(&cj).copied().unwrap_or(0.0);
                let g = w - scale * tot_cj * nw[i];
                if g > best_gain + 1e-12 {
                    best_gain = g;
                    best = cj;
                    use_empty = false;
                }
            }

            if use_empty {
                let new_c = tot.keys().copied().max().unwrap_or(0) + 1;
                best = new_c;
            }
            *tot.entry(best).or_insert(0.0) += nw[i];
            if best != ci {
                p[i] = best;
                for &(j, _) in &adj[i] {
                    if p[j] != best && !in_queue[j] {
                        in_queue[j] = true;
                        queue.push_back(j);
                    }
                }
            }
        }

        // Refinement phase: split into well-connected subcommunities within each community of `p`.
        let mut comms: HashMap<usize, Vec<usize>> = HashMap::new();
        for i in 0..cur_n {
            comms.entry(p[i]).or_default().push(i);
        }
        if comms.len() == cur_n {
            break;
        }

        let mut r: Vec<usize> = (0..cur_n).collect();
        let mut rtot = nw.clone();
        let mut rsize = vec![1usize; cur_n];
        let mut s_out = vec![0.0_f64; cur_n];
        let mut in_c = vec![false; cur_n];

        let mut comm_keys: Vec<usize> = comms.keys().copied().collect();
        comm_keys.sort_unstable();
        for cid in comm_keys {
            let nodes = &comms[&cid];
            let mut tot_c = 0.0_f64;
            for &i in nodes {
                in_c[i] = true;
                tot_c += nw[i];
            }
            let mut w_c = vec![0.0_f64; cur_n];
            for &i in nodes {
                let mut s = 0.0_f64;
                for &(j, w) in &adj[i] {
                    if in_c[j] {
                        s += w;
                    }
                }
                w_c[i] = s;
                s_out[i] = s;
            }
            let visit: Vec<usize> = shuffled(nodes.len(), Some(rng))
                .into_iter()
                .map(|idx| nodes[idx])
                .collect();
            for v in visit {
                if rsize[r[v]] != 1 || r[v] != v {
                    continue;
                }
                if w_c[v] < scale * nw[v] * (tot_c - nw[v]) {
                    continue;
                }
                let mut ws_keys: Vec<usize> = Vec::new();
                let mut ws_map: HashMap<usize, usize> = HashMap::new();
                let mut ws_vals: Vec<f64> = Vec::new();
                for &(j, w) in &adj[v] {
                    if in_c[j] {
                        let rj = r[j];
                        if let Some(&pos) = ws_map.get(&rj) {
                            ws_vals[pos] += w;
                        } else {
                            ws_map.insert(rj, ws_keys.len());
                            ws_keys.push(rj);
                            ws_vals.push(w);
                        }
                    }
                }
                let mut best = v;
                let mut best_g = 0.0_f64;
                for (pos, &s) in ws_keys.iter().enumerate() {
                    let w = ws_vals[pos];
                    if s_out[s] < scale * rtot[s] * (tot_c - rtot[s]) {
                        continue;
                    }
                    let g = w - scale * nw[v] * rtot[s];
                    if g > best_g + 1e-12 {
                        best_g = g;
                        best = s;
                    }
                }
                if best != v {
                    let w = ws_vals[ws_map[&best]];
                    r[v] = best;
                    rtot[best] += nw[v];
                    rsize[best] += 1;
                    rsize[v] = 0;
                    s_out[best] = s_out[best] + w_c[v] - 2.0 * w;
                }
            }
            for &i in nodes {
                in_c[i] = false;
            }
        }

        let mut ids: HashMap<usize, usize> = HashMap::new();
        for &x in &r {
            let next_id = ids.len();
            ids.entry(x).or_insert(next_id);
        }
        let n_sub = ids.len();
        if n_sub == cur_n {
            break;
        }

        let mut next_member = vec![Vec::new(); n_sub];
        let mut next_size = vec![0.0_f64; n_sub];
        let mut next_p = vec![0usize; n_sub];
        let mut pid: HashMap<usize, usize> = HashMap::new();
        for i in 0..cur_n {
            let a = ids[&r[i]];
            next_member[a].extend_from_slice(&member[i]);
            next_size[a] += node_size[i];
            let next_pid = pid.len();
            next_p[a] = *pid.entry(p[i]).or_insert(next_pid);
        }

        let mut agg_map: HashMap<(usize, usize), f64> = HashMap::new();
        for &(u, v, w) in &cur_edges {
            let a = ids[&r[u]];
            let b = ids[&r[v]];
            let key = if a < b { (a, b) } else { (b, a) };
            *agg_map.entry(key).or_insert(0.0) += w;
        }
        let mut agg_edges: Vec<(usize, usize, f64)> =
            agg_map.into_iter().map(|((a, b), w)| (a, b, w)).collect();
        agg_edges.sort_unstable_by(|x, y| x.0.cmp(&y.0).then(x.1.cmp(&y.1)));

        cur_edges = agg_edges;
        cur_n = n_sub;
        member = next_member;
        node_size = next_size;
        p = next_p;
    }

    let mut out = vec![0usize; n];
    for i in 0..cur_n {
        let cid = p[i];
        for &orig in &member[i] {
            out[orig] = cid;
        }
    }
    out
}

/// Binary min-heap matching `livia-atlas/app.js` `class MinHeap` sifting order bit-for-bit.
struct JsMinHeap {
    a: Vec<(f64, usize)>,
}

impl JsMinHeap {
    #[inline]
    fn new() -> Self {
        Self { a: Vec::new() }
    }

    #[inline]
    fn push(&mut self, k: f64, v: usize) {
        self.a.push((k, v));
        let mut i = self.a.len() - 1;
        while i > 0 {
            let p = (i - 1) >> 1;
            if self.a[p].0 <= self.a[i].0 {
                break;
            }
            self.a.swap(p, i);
            i = p;
        }
    }

    #[inline]
    fn pop(&mut self) -> Option<(f64, usize)> {
        if self.a.is_empty() {
            return None;
        }
        let top = self.a[0];
        let last = self.a.pop().unwrap();
        if !self.a.is_empty() {
            self.a[0] = last;
            let len = self.a.len();
            let mut i = 0;
            loop {
                let l = 2 * i + 1;
                let r = l + 1;
                let mut m = i;
                if l < len && self.a[l].0 < self.a[m].0 {
                    m = l;
                }
                if r < len && self.a[r].0 < self.a[m].0 {
                    m = r;
                }
                if m == i {
                    break;
                }
                self.a.swap(m, i);
                i = m;
            }
        }
        Some(top)
    }
}

/// Steiner / prize-collecting minimal-cost connector matching `livia-atlas/app.js` `steinerAdd`.
/// `edges` contains passing undirected edges `(u, v, best_ilis)`.
/// `prizes` is either empty (unweighted Takahashi-Matsuyama) or length `terminals.len()` parallel to `terminals`.
/// Returns `(joined_count, total_cost, nodes_in_tree_insertion_order)`.
pub fn steiner_add(
    n_nodes: usize,
    edges: &[(usize, usize, f64)],
    terminals: &[usize],
    prizes: &[f64],
    cap: usize,
) -> (usize, f64, Vec<usize>) {
    steiner_add_opts(n_nodes, edges, terminals, prizes, cap, false)
}

/// Steiner / prize-collecting minimal-cost connector with optional pre-expanded directed CSR edges.
/// When `directed == true`, each `(u, v, best)` appends only `(v, best)` to `adj[u]`, preserving exact CSR order.
pub fn steiner_add_opts(
    n_nodes: usize,
    edges: &[(usize, usize, f64)],
    terminals: &[usize],
    prizes: &[f64],
    cap: usize,
    directed: bool,
) -> (usize, f64, Vec<usize>) {
    if terminals.is_empty() || n_nodes == 0 {
        return (0, 0.0, Vec::new());
    }
    let has_prize = prizes.len() == terminals.len() && !prizes.is_empty();
    let mut adj: Vec<Vec<(usize, f64)>> = vec![Vec::new(); n_nodes];
    if directed {
        for &(u, v, best) in edges {
            if u < n_nodes && v < n_nodes && u != v {
                adj[u].push((v, best));
            }
        }
    } else {
        for &(u, v, best) in edges {
            if u < n_nodes && v < n_nodes && u != v {
                adj[u].push((v, best));
                adj[v].push((u, best));
            }
        }
    }

    let mut is_q = vec![false; n_nodes];
    let mut is_term = vec![false; n_nodes];
    let mut prize_of = vec![0.0_f64; n_nodes];
    let mut unique_terms: Vec<usize> = Vec::new();
    for (idx, &t) in terminals.iter().enumerate() {
        if t < n_nodes && !is_q[t] {
            is_q[t] = true;
            is_term[t] = true;
            if has_prize {
                prize_of[t] = prizes[idx];
            }
            unique_terms.push(t);
        }
    }
    if unique_terms.is_empty() {
        return (0, 0.0, Vec::new());
    }

    let qn = |i: usize| -> usize {
        adj[i].iter().filter(|&&(j, _)| j != i && is_q[j]).count()
    };

    unique_terms.sort_by(|&a, &b| {
        let cmp = if has_prize {
            prize_of[b].partial_cmp(&prize_of[a]).unwrap_or(Ordering::Equal)
        } else {
            qn(b).cmp(&qn(a))
        };
        cmp.then_with(|| a.cmp(&b))
    });

    let start = unique_terms[0];
    let mut in_t = vec![false; n_nodes];
    let mut tree_nodes: Vec<usize> = vec![start];
    in_t[start] = true;
    is_term[start] = false;
    let mut terms_left = unique_terms.len() - 1;

    let mut total_cost = 0.0_f64;
    let mut joined = 1usize;

    let mut dist = vec![f64::INFINITY; n_nodes];
    let mut hop = vec![0u8; n_nodes];
    let mut prev = vec![usize::MAX; n_nodes];
    let mut touched: Vec<usize> = Vec::new();

    while terms_left > 0 && tree_nodes.len() < cap {
        for &u in &touched {
            dist[u] = f64::INFINITY;
            hop[u] = 0;
            prev[u] = usize::MAX;
        }
        touched.clear();

        let bound = if has_prize {
            unique_terms
                .iter()
                .filter(|&&t| is_term[t])
                .map(|&t| prize_of[t])
                .fold(f64::NEG_INFINITY, f64::max)
        } else {
            f64::INFINITY
        };

        let mut heap = JsMinHeap::new();
        for &t in &tree_nodes {
            dist[t] = 0.0;
            hop[t] = 0;
            touched.push(t);
            heap.push(0.0, t);
        }

        let mut found: Vec<(usize, f64)> = Vec::new();
        let mut found_pos: HashMap<usize, usize> = HashMap::new();
        while let Some((d, u)) = heap.pop() {
            if d > dist[u] {
                continue;
            }
            if d >= bound {
                break;
            }
            if is_term[u] {
                if let Some(&pos) = found_pos.get(&u) {
                    found[pos].1 = d;
                } else {
                    found_pos.insert(u, found.len());
                    found.push((u, d));
                }
                if !has_prize {
                    break;
                }
                continue;
            }
            if hop[u] >= 4 {
                continue;
            }
            for &(v, best) in &adj[u] {
                if v == u || in_t[v] || (is_q[v] && !is_term[v]) {
                    continue;
                }
                let nd = d + 1.05 - best.min(1.0);
                if nd < dist[v] {
                    if dist[v].is_infinite() {
                        touched.push(v);
                    }
                    dist[v] = nd;
                    hop[v] = hop[u] + 1;
                    prev[v] = u;
                    heap.push(nd, v);
                }
            }
        }

        let mut pick: Option<(usize, f64)> = None;
        let mut best_gain = f64::NEG_INFINITY;
        for &(t, d) in &found {
            let g = if has_prize { prize_of[t] - d } else { -d };
            let better = if g > best_gain + 1e-12 {
                true
            } else if (g - best_gain).abs() <= 1e-12 {
                pick.map_or(true, |(pt, _)| t < pt)
            } else {
                false
            };
            if better {
                best_gain = g;
                pick = Some((t, d));
            }
        }

        let Some((chosen, chosen_d)) = pick else { break };
        if has_prize && best_gain <= 0.0 {
            break;
        }

        let mut cur = chosen;
        while cur != usize::MAX && !in_t[cur] {
            in_t[cur] = true;
            tree_nodes.push(cur);
            cur = prev[cur];
        }
        total_cost += chosen_d;
        is_term[chosen] = false;
        terms_left -= 1;
        joined += 1;
    }

    (joined, total_cost, tree_nodes)
}

/// Fruchterman-Reingold weighted force-directed layout (`method = "fr"`).
/// Uses exact all-pairs repulsion up to `n = 600` and spatial grid Barnes-Hut-like cell
/// aggregation above `600` so 5,000+ node interactomes finish in milliseconds.
pub fn layout_fr(
    n: usize,
    edges: &[(usize, usize, f64)],
    niter: usize,
    seed: u32,
) -> Vec<(f64, f64)> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![(0.0, 0.0)];
    }
    let iters = niter.max(10).min(5000);
    let mut rng = Mulberry32::new(if seed == 0 { 1 } else { seed });
    let width = (n as f64).sqrt();
    let k = (width * width / (n as f64)).sqrt();
    let k2 = k * k;
    let wmax = edges
        .iter()
        .map(|e| e.2)
        .fold(1e-9_f64, f64::max);

    let mut px: Vec<f64> = (0..n).map(|_| (rng.next_f64() - 0.5) * width).collect();
    let mut py: Vec<f64> = (0..n).map(|_| (rng.next_f64() - 0.5) * width).collect();
    let mut dx = vec![0.0_f64; n];
    let mut dy = vec![0.0_f64; n];

    let temp0 = width * 0.12;
    for it in 0..iters {
        let temp = temp0 * (1.0 - (it as f64) / (iters as f64)) + 1e-3 * k;
        dx.fill(0.0);
        dy.fill(0.0);

        if n <= 600 {
            for i in 0..n {
                let xi = px[i];
                let yi = py[i];
                for j in (i + 1)..n {
                    let rx = xi - px[j];
                    let ry = yi - py[j];
                    let d2 = (rx * rx + ry * ry).max(1e-6);
                    let f = k2 / d2;
                    dx[i] += rx * f;
                    dy[i] += ry * f;
                    dx[j] -= rx * f;
                    dy[j] -= ry * f;
                }
            }
        } else {
            grid_repulsion(&px, &py, &mut dx, &mut dy, k2, 24);
        }

        for &(u, v, w) in edges {
            if u >= n || v >= n || u == v {
                continue;
            }
            let rx = px[u] - px[v];
            let ry = py[u] - py[v];
            let d = (rx * rx + ry * ry).sqrt().max(1e-4);
            let wf = 0.3 + 0.7 * (w / wmax);
            let f = (d / k) * wf;
            let fx = rx * f;
            let fy = ry * f;
            dx[u] -= fx;
            dy[u] -= fy;
            dx[v] += fx;
            dy[v] += fy;
        }

        for i in 0..n {
            let len = (dx[i] * dx[i] + dy[i] * dy[i]).sqrt().max(1e-9);
            let step = len.min(temp) / len;
            px[i] += dx[i] * step;
            py[i] += dy[i] * step;
            // Gentle centering gravity keeps disconnected islands from drifting infinitely far apart.
            px[i] *= 0.992;
            py[i] *= 0.992;
        }
    }

    px.into_iter().zip(py).collect()
}

fn grid_repulsion(
    px: &[f64],
    py: &[f64],
    dx: &mut [f64],
    dy: &mut [f64],
    k2: f64,
    g: usize,
) {
    let n = px.len();
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;
    for i in 0..n {
        min_x = min_x.min(px[i]);
        max_x = max_x.max(px[i]);
        min_y = min_y.min(py[i]);
        max_y = max_y.max(py[i]);
    }
    let span_x = (max_x - min_x).max(1e-3);
    let span_y = (max_y - min_y).max(1e-3);
    let cell_w = span_x / (g as f64);
    let cell_h = span_y / (g as f64);

    let cells_n = g * g;
    let mut cell_nodes: Vec<Vec<usize>> = vec![Vec::new(); cells_n];
    let mut cell_mass = vec![0.0_f64; cells_n];
    let mut cell_cx = vec![0.0_f64; cells_n];
    let mut cell_cy = vec![0.0_f64; cells_n];
    let mut node_cell = vec![0usize; n];

    for i in 0..n {
        let cx = (((px[i] - min_x) / cell_w).floor() as isize).clamp(0, (g - 1) as isize) as usize;
        let cy = (((py[i] - min_y) / cell_h).floor() as isize).clamp(0, (g - 1) as isize) as usize;
        let c = cy * g + cx;
        node_cell[i] = c;
        cell_nodes[c].push(i);
        cell_mass[c] += 1.0;
        cell_cx[c] += px[i];
        cell_cy[c] += py[i];
    }
    for c in 0..cells_n {
        if cell_mass[c] > 0.0 {
            cell_cx[c] /= cell_mass[c];
            cell_cy[c] /= cell_mass[c];
        }
    }

    for i in 0..n {
        let ci = node_cell[i];
        let cix = (ci % g) as isize;
        let ciy = (ci / g) as isize;
        let xi = px[i];
        let yi = py[i];
        for c in 0..cells_n {
            if cell_mass[c] == 0.0 {
                continue;
            }
            let cx = (c % g) as isize;
            let cy = (c / g) as isize;
            if (cx - cix).abs() <= 1 && (cy - ciy).abs() <= 1 {
                for &j in &cell_nodes[c] {
                    if j != i {
                        let rx = xi - px[j];
                        let ry = yi - py[j];
                        let d2 = (rx * rx + ry * ry).max(1e-6);
                        let f = k2 / d2;
                        dx[i] += rx * f;
                        dy[i] += ry * f;
                    }
                }
            } else {
                let rx = xi - cell_cx[c];
                let ry = yi - cell_cy[c];
                let d2 = (rx * rx + ry * ry).max(1e-4);
                let f = cell_mass[c] * k2 / d2;
                dx[i] += rx * f;
                dy[i] += ry * f;
            }
        }
    }
}

/// Kamada-Kawai layout via Gansner-Koren-North stress majorization (`method = "kk"`).
/// For large graphs (`n > 700`), uses pivot stress majorization so shortest-path memory
/// and runtime remain linear-logarithmic instead of stalling the browser.
pub fn layout_kk(
    n: usize,
    edges: &[(usize, usize, f64)],
    niter: usize,
    seed: u32,
) -> Vec<(f64, f64)> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![(0.0, 0.0)];
    }
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for &(u, v, _) in edges {
        if u < n && v < n && u != v {
            adj[u].push(v);
            adj[v].push(u);
        }
    }

    let mut rng = Mulberry32::new(if seed == 0 { 1 } else { seed });

    if n > 700 {
        // Landmark / pivot stress majorization + local edge refinement for huge graphs.
        let n_pivots = 96.min(n);
        let step = (n / n_pivots).max(1);
        let pivots: Vec<usize> = (0..n).step_by(step).take(n_pivots).collect();
        let mut pdist: Vec<Vec<i32>> = Vec::with_capacity(pivots.len());
        let mut dmax = 1i32;
        let mut q = VecDeque::with_capacity(n);
        for &p in &pivots {
            let mut d = vec![-1i32; n];
            d[p] = 0;
            q.push_back(p);
            while let Some(u) = q.pop_front() {
                for &v in &adj[u] {
                    if d[v] < 0 {
                        d[v] = d[u] + 1;
                        dmax = dmax.max(d[v]);
                        q.push_back(v);
                    }
                }
            }
            pdist.push(d);
        }
        let scale = dmax as f64;
        let mut px: Vec<f64> = (0..n).map(|_| (rng.next_f64() - 0.5) * scale).collect();
        let mut py: Vec<f64> = (0..n).map(|_| (rng.next_f64() - 0.5) * scale).collect();
        let iters = niter.min(180).max(30);
        for _ in 0..iters {
            for i in 0..n {
                let mut sx = 0.0_f64;
                let mut sy = 0.0_f64;
                let mut sw = 0.0_f64;
                for (pi_idx, &p) in pivots.iter().enumerate() {
                    if p == i {
                        continue;
                    }
                    let raw = pdist[pi_idx][i];
                    let dij = (if raw < 0 { dmax + 1 } else { raw }) as f64;
                    let w = 1.0 / (dij * dij);
                    let rx = px[i] - px[p];
                    let ry = py[i] - py[p];
                    let l = (rx * rx + ry * ry).sqrt().max(1e-6);
                    sx += w * (px[p] + dij * rx / l);
                    sy += w * (py[p] + dij * ry / l);
                    sw += w;
                }
                for &nb in &adj[i] {
                    let rx = px[i] - px[nb];
                    let ry = py[i] - py[nb];
                    let l = (rx * rx + ry * ry).sqrt().max(1e-6);
                    sx += 2.0 * (px[nb] + rx / l);
                    sy += 2.0 * (py[nb] + ry / l);
                    sw += 2.0;
                }
                if sw > 0.0 {
                    px[i] = sx / sw;
                    py[i] = sy / sw;
                }
            }
        }
        return px.into_iter().zip(py).collect();
    }

    let mut dm = vec![-1i16; n * n];
    let mut dmax = 1i16;
    let mut q = VecDeque::with_capacity(n);
    for s in 0..n {
        let row = &mut dm[s * n..(s + 1) * n];
        row[s] = 0;
        q.clear();
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            let du = row[u];
            for &v in &adj[u] {
                if row[v] < 0 {
                    row[v] = du + 1;
                    if row[v] > dmax {
                        dmax = row[v];
                    }
                    q.push_back(v);
                }
            }
        }
    }

    let scale = dmax as f64;
    let mut px: Vec<f64> = (0..n).map(|_| (rng.next_f64() - 0.5) * scale).collect();
    let mut py: Vec<f64> = (0..n).map(|_| (rng.next_f64() - 0.5) * scale).collect();
    let iters = niter.min(300).max(20);
    let fallback_d = (dmax + 1) as f64;

    for _ in 0..iters {
        for i in 0..n {
            let mut sx = 0.0_f64;
            let mut sy = 0.0_f64;
            let mut sw = 0.0_f64;
            let row = &dm[i * n..(i + 1) * n];
            let xi = px[i];
            let yi = py[i];
            for j in 0..n {
                if j == i {
                    continue;
                }
                let raw = row[j];
                let dij = if raw < 0 { fallback_d } else { raw as f64 };
                let w = 1.0 / (dij * dij);
                let rx = xi - px[j];
                let ry = yi - py[j];
                let l = (rx * rx + ry * ry).sqrt().max(1e-6);
                sx += w * (px[j] + dij * rx / l);
                sy += w * (py[j] + dij * ry / l);
                sw += w;
            }
            if sw > 0.0 {
                px[i] = sx / sw;
                py[i] = sy / sw;
            }
        }
    }
    px.into_iter().zip(py).collect()
}

/// DrL / OpenOrd multi-stage annealing layout (`method = "drl"`),
/// specifically tailored to highlight modular PPI clusters cleanly.
/// Stages (Martin et al. 2011): Liquid -> Expansion -> Cooldown -> Crunch -> Simmer.
pub fn layout_drl(
    n: usize,
    edges: &[(usize, usize, f64)],
    niter: usize,
    seed: u32,
) -> Vec<(f64, f64)> {
    if n == 0 {
        return Vec::new();
    }
    if n == 1 {
        return vec![(0.0, 0.0)];
    }
    let mut rng = Mulberry32::new(if seed == 0 { 1 } else { seed });
    let total_iters = niter.max(100).min(2000);
    let width = (n as f64).sqrt() * 2.0;
    let k = (width * width / (n as f64)).sqrt();
    let k2 = k * k;
    let wmax = edges.iter().map(|e| e.2).fold(1e-9_f64, f64::max);

    let mut px: Vec<f64> = (0..n).map(|_| (rng.next_f64() - 0.5) * 0.2 * width).collect();
    let mut py: Vec<f64> = (0..n).map(|_| (rng.next_f64() - 0.5) * 0.2 * width).collect();
    let mut dx = vec![0.0_f64; n];
    let mut dy = vec![0.0_f64; n];

    // (fraction_of_iters, temp_factor, attraction_power, repulsion_scale, damping)
    let stages: [(f64, f64, f64, f64, f64); 5] = [
        (0.25, 0.30, 2.0, 1.0, 0.995), // Liquid
        (0.25, 0.25, 1.2, 2.2, 0.997), // Expansion (separates modules)
        (0.25, 0.12, 1.5, 1.2, 0.995), // Cooldown
        (0.10, 0.04, 1.8, 0.7, 0.992), // Crunch
        (0.15, 0.015, 1.0, 1.0, 0.995), // Simmer
    ];

    for &(frac, tfac, attr_p, rep_s, damp) in &stages {
        let stage_iters = ((total_iters as f64) * frac).round().max(5.0) as usize;
        for it in 0..stage_iters {
            let progress = (it as f64) / (stage_iters as f64);
            let temp = width * tfac * (1.0 - 0.5 * progress);
            dx.fill(0.0);
            dy.fill(0.0);

            if n <= 500 {
                for i in 0..n {
                    let xi = px[i];
                    let yi = py[i];
                    for j in (i + 1)..n {
                        let rx = xi - px[j];
                        let ry = yi - py[j];
                        let d2 = (rx * rx + ry * ry).max(1e-5);
                        let f = rep_s * k2 / d2;
                        dx[i] += rx * f;
                        dy[i] += ry * f;
                        dx[j] -= rx * f;
                        dy[j] -= ry * f;
                    }
                }
            } else {
                grid_repulsion(&px, &py, &mut dx, &mut dy, rep_s * k2, 28);
            }

            for &(u, v, w) in edges {
                if u >= n || v >= n || u == v {
                    continue;
                }
                let rx = px[u] - px[v];
                let ry = py[u] - py[v];
                let d = (rx * rx + ry * ry).sqrt().max(1e-4);
                let wf = (w / wmax).powf(1.25).max(0.15);
                let f = (d / k).powf(attr_p) * wf;
                dx[u] -= rx * f;
                dy[u] -= ry * f;
                dx[v] += rx * f;
                dy[v] += ry * f;
            }

            for i in 0..n {
                let len = (dx[i] * dx[i] + dy[i] * dy[i]).sqrt().max(1e-9);
                let step = len.min(temp) / len;
                px[i] = (px[i] + dx[i] * step) * damp;
                py[i] = (py[i] + dy[i] * step) * damp;
            }
        }
    }

    px.into_iter().zip(py).collect()
}
