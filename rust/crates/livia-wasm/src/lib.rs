//! Zero-overhead TypedArray WebAssembly ABI over `livia-core`.
//!
//! All input allocations created via `wasm_alloc` are 8-byte (`u64`) aligned so JavaScript
//! `Float64Array`, `Float32Array`, and `Uint32Array` views over `memory.buffer` never hit
//! alignment faults on any browser engine (V8, JavaScriptCore, SpiderMonkey).

use std::cell::UnsafeCell;

use livia_core::clip::{
    average_linkage_upgma, cosine_dist_from_csr, fcluster, find_optimal_k, leaf_order,
    remap_labels_by_size, run_clip_from_csr, UpperDistMatrix,
};
use livia_core::graph::{
    layout_drl, layout_fr, layout_kk, leiden_atlas, leiden_network, steiner_add_opts,
};
use livia_core::metrics::{calc_pae_chain_pair_iptm, compute_directed_pairs_kernel};
use livia_core::structure::ResidueCoord;

struct ScratchBuffers {
    out_f64: Vec<f64>,
    out_f32: Vec<f32>,
    out_u32: Vec<u32>,
}

struct SyncScratch(UnsafeCell<ScratchBuffers>);

// SAFETY: WebAssembly instances execute on a single thread (either main thread or a dedicated Web Worker).
unsafe impl Sync for SyncScratch {}

static SCRATCH: SyncScratch = SyncScratch(UnsafeCell::new(ScratchBuffers {
    out_f64: Vec::new(),
    out_f32: Vec::new(),
    out_u32: Vec::new(),
}));

#[inline]
fn scratch_mut() -> &'static mut ScratchBuffers {
    unsafe { &mut *SCRATCH.0.get() }
}

/// Allocate `byte_len` bytes aligned to 8 bytes (`u64`).
#[no_mangle]
pub extern "C" fn wasm_alloc(byte_len: usize) -> *mut u8 {
    let words = byte_len.div_ceil(8).max(1);
    let mut buf = vec![0u64; words];
    let ptr = buf.as_mut_ptr() as *mut u8;
    core::mem::forget(buf);
    ptr
}

/// Free a buffer previously allocated by `wasm_alloc(byte_len)`.
#[no_mangle]
pub extern "C" fn wasm_free(ptr: *mut u8, byte_len: usize) {
    if ptr.is_null() {
        return;
    }
    let words = byte_len.div_ceil(8).max(1);
    unsafe {
        let _ = Vec::from_raw_parts(ptr as *mut u64, words, words);
    }
}

/// Release retained capacity in internal output buffers after large runs.
#[no_mangle]
pub extern "C" fn wasm_clear_outputs() {
    let s = scratch_mut();
    s.out_f64 = Vec::new();
    s.out_f32 = Vec::new();
    s.out_u32 = Vec::new();
}

#[no_mangle]
pub extern "C" fn wasm_out_f64_ptr() -> *const f64 {
    scratch_mut().out_f64.as_ptr()
}

#[no_mangle]
pub extern "C" fn wasm_out_f64_len() -> usize {
    scratch_mut().out_f64.len()
}

#[no_mangle]
pub extern "C" fn wasm_out_f32_ptr() -> *const f32 {
    scratch_mut().out_f32.as_ptr()
}

#[no_mangle]
pub extern "C" fn wasm_out_f32_len() -> usize {
    scratch_mut().out_f32.len()
}

#[no_mangle]
pub extern "C" fn wasm_out_u32_ptr() -> *const u32 {
    scratch_mut().out_u32.as_ptr()
}

#[no_mangle]
pub extern "C" fn wasm_out_u32_len() -> usize {
    scratch_mut().out_u32.len()
}

#[inline]
unsafe fn slice_or_empty<'a, T>(ptr: *const T, len: usize) -> &'a [T] {
    if len == 0 || ptr.is_null() {
        &[]
    } else {
        core::slice::from_raw_parts(ptr, len)
    }
}

fn decode_edges(flat: &[f64]) -> Vec<(usize, usize, f64)> {
    let m = flat.len() / 3;
    let mut out = Vec::with_capacity(m);
    for idx in 0..m {
        let u = flat[idx * 3] as usize;
        let v = flat[idx * 3 + 1] as usize;
        let w = flat[idx * 3 + 2];
        out.push((u, v, w));
    }
    out
}

/// Multi-iteration Leiden community detection (`modularity` or `CPM`) for `LIVIA/network.html`.
/// Writes `[membership[0..n], degree[0..n]]` into `out_u32`.
#[no_mangle]
pub extern "C" fn wasm_leiden_network(
    n: u32,
    edges_ptr: *const f64,
    n_edges: u32,
    gamma: f64,
    use_cpm: u32,
    n_iterations: u32,
    seed: u32,
) -> u32 {
    let flat = unsafe { slice_or_empty(edges_ptr, (n_edges as usize) * 3) };
    let edges = decode_edges(flat);
    let (membership, degree) = leiden_network(
        n as usize,
        &edges,
        gamma,
        use_cpm != 0,
        n_iterations as usize,
        seed,
    );
    let s = scratch_mut();
    s.out_u32.clear();
    s.out_u32.reserve(membership.len() + degree.len());
    s.out_u32.extend_from_slice(&membership);
    s.out_u32.extend_from_slice(&degree);
    n
}

/// Exact deterministic Leiden matching `livia-atlas/app.js` `leiden(n, edges, gamma, rand)`.
/// `seed == 0` uses fixed order (`rand = null`); `seed > 0` uses `Mulberry32(seed)`.
/// Writes `membership[0..n]` into `out_u32`.
#[no_mangle]
pub extern "C" fn wasm_leiden_atlas(
    n: u32,
    edges_ptr: *const f64,
    n_edges: u32,
    gamma: f64,
    seed: u32,
) -> u32 {
    let flat = unsafe { slice_or_empty(edges_ptr, (n_edges as usize) * 3) };
    let edges = decode_edges(flat);
    let membership = leiden_atlas(n as usize, &edges, gamma, seed);
    let s = scratch_mut();
    s.out_u32.clear();
    s.out_u32.extend_from_slice(&membership);
    n
}

/// Exact Steiner / prize-collecting minimal connector matching `livia-atlas/app.js` `steinerAdd`.
/// Writes `[joined as f64, total_cost]` into `out_f64` and `tree_nodes[0..len]` into `out_u32`.
#[no_mangle]
pub extern "C" fn wasm_steiner_add(
    n_nodes: u32,
    edges_ptr: *const f64,
    n_edges: u32,
    terms_ptr: *const u32,
    n_terms: u32,
    prizes_ptr: *const f64,
    n_prizes: u32,
    cap: u32,
    directed: u32,
) -> u32 {
    let flat_edges = unsafe { slice_or_empty(edges_ptr, (n_edges as usize) * 3) };
    let terms_u32 = unsafe { slice_or_empty(terms_ptr, n_terms as usize) };
    let prizes = unsafe { slice_or_empty(prizes_ptr, n_prizes as usize) };
    let edges = decode_edges(flat_edges);
    let terms: Vec<usize> = terms_u32.iter().map(|&t| t as usize).collect();

    let (joined, cost, nodes) = steiner_add_opts(
        n_nodes as usize,
        &edges,
        &terms,
        prizes,
        cap as usize,
        directed != 0,
    );
    let s = scratch_mut();
    s.out_f64.clear();
    s.out_f64.push(joined as f64);
    s.out_f64.push(cost);
    s.out_u32.clear();
    s.out_u32.extend(nodes.into_iter().map(|x| x as u32));
    s.out_u32.len() as u32
}

/// Graph layout (`0 = fr`, `1 = kk`, `2 = drl`).
/// Writes `[x0, y0, x1, y1, ...]` (`2 * n` f64s) into `out_f64`.
#[no_mangle]
pub extern "C" fn wasm_graph_layout(
    method: u32,
    n: u32,
    edges_ptr: *const f64,
    n_edges: u32,
    niter: u32,
    seed: u32,
) -> u32 {
    let flat = unsafe { slice_or_empty(edges_ptr, (n_edges as usize) * 3) };
    let edges = decode_edges(flat);
    let coords = match method {
        1 => layout_kk(n as usize, &edges, niter as usize, seed),
        2 => layout_drl(n as usize, &edges, niter as usize, seed),
        _ => layout_fr(n as usize, &edges, niter as usize, seed),
    };
    let s = scratch_mut();
    s.out_f64.clear();
    s.out_f64.reserve(coords.len() * 2);
    for (x, y) in coords {
        s.out_f64.push(x);
        s.out_f64.push(y);
    }
    n
}

/// Compute only the upper-triangle bitset cosine distance matrix (`n * (n - 1) / 2` f64s in `out_f64`).
#[no_mangle]
pub extern "C" fn wasm_clip_cosine_dist(
    offsets_ptr: *const u32,
    offsets_len: u32,
    flat_ptr: *const u32,
    flat_len: u32,
) -> u32 {
    let offsets = unsafe { slice_or_empty(offsets_ptr, offsets_len as usize) };
    let flat = unsafe { slice_or_empty(flat_ptr, flat_len as usize) };
    let dmat = cosine_dist_from_csr(offsets, flat);
    let n = dmat.n as u32;
    let s = scratch_mut();
    s.out_f64 = dmat.data;
    n
}

/// Full end-to-end cLIP clustering from CSR fingerprints (`offsets`, `flat_residues`).
/// Writes `[optimal_k, labels[0..n], raw_labels[0..n], order[0..n]]` into `out_u32`
/// and flattened linkage `Z` (`4 * (n - 1)` f64s) into `out_f64`.
#[no_mangle]
pub extern "C" fn wasm_clip_full(
    offsets_ptr: *const u32,
    offsets_len: u32,
    flat_ptr: *const u32,
    flat_len: u32,
    max_k: u32,
    drop_threshold: f64,
) -> u32 {
    let offsets = unsafe { slice_or_empty(offsets_ptr, offsets_len as usize) };
    let flat = unsafe { slice_or_empty(flat_ptr, flat_len as usize) };
    let res = run_clip_from_csr(offsets, flat, max_k as usize, drop_threshold);

    let s = scratch_mut();
    s.out_u32.clear();
    s.out_u32.reserve(1 + res.labels.len() * 3);
    s.out_u32.push(res.optimal_k as u32);
    s.out_u32.extend_from_slice(&res.labels);
    s.out_u32.extend_from_slice(&res.raw_labels);
    s.out_u32.extend_from_slice(&res.order);

    s.out_f64.clear();
    s.out_f64.reserve(res.linkage.len() * 4);
    for row in &res.linkage {
        s.out_f64.push(row.c1 as f64);
        s.out_f64.push(row.c2 as f64);
        s.out_f64.push(row.dist);
        s.out_f64.push(row.count as f64);
    }
    res.optimal_k as u32
}

/// Run UPGMA linkage + incremental optimal-K + labels + leafOrder from an upper-triangle distance matrix.
/// Writes `[optimal_k, labels[0..n], raw_labels[0..n], order[0..n]]` into `out_u32`
/// and flattened linkage `Z` (`4 * (n - 1)` f64s) into `out_f64`.
#[no_mangle]
pub extern "C" fn wasm_clip_from_dist(
    upper_ptr: *const f64,
    upper_len: u32,
    n: u32,
    max_k: u32,
    drop_threshold: f64,
) -> u32 {
    let n_usize = n as usize;
    let slice = unsafe { slice_or_empty(upper_ptr, upper_len as usize) };
    let dmat = UpperDistMatrix {
        n: n_usize,
        data: slice.to_vec(),
    };
    let linkage = average_linkage_upgma(&dmat);
    let optimal_k = find_optimal_k(&linkage, &dmat, max_k as usize, drop_threshold);
    let raw_labels = fcluster(&linkage, n_usize, optimal_k);
    let labels = remap_labels_by_size(&raw_labels);
    let order = leaf_order(&linkage, n_usize);

    let s = scratch_mut();
    s.out_u32.clear();
    s.out_u32.reserve(1 + n_usize * 3);
    s.out_u32.push(optimal_k as u32);
    s.out_u32.extend_from_slice(&labels);
    s.out_u32.extend_from_slice(&raw_labels);
    s.out_u32.extend_from_slice(&order);

    s.out_f64.clear();
    s.out_f64.reserve(linkage.len() * 4);
    for row in &linkage {
        s.out_f64.push(row.c1 as f64);
        s.out_f64.push(row.c2 as f64);
        s.out_f64.push(row.dist);
        s.out_f64.push(row.count as f64);
    }
    optimal_k as u32
}

/// PAE-derived per-chain-pair ipTM fallback (`nc * nc` f64s written to `out_f64`).
#[no_mangle]
pub extern "C" fn wasm_pae_chain_pair_iptm(
    pae_ptr: *const f64,
    pae_dim: u32,
    starts_ptr: *const u32,
    ends_ptr: *const u32,
    nc: u32,
) -> u32 {
    let dim = pae_dim as usize;
    let pae = unsafe { slice_or_empty(pae_ptr, dim * dim) };
    let st_u32 = unsafe { slice_or_empty(starts_ptr, nc as usize) };
    let en_u32 = unsafe { slice_or_empty(ends_ptr, nc as usize) };
    let starts: Vec<usize> = st_u32.iter().map(|&x| x as usize).collect();
    let ends: Vec<usize> = en_u32.iter().map(|&x| x as usize).collect();

    let mat = calc_pae_chain_pair_iptm(pae, dim, &starts, &ends);
    let s = scratch_mut();
    s.out_f64.clear();
    s.out_f64.reserve((nc as usize) * (nc as usize));
    for row in mat {
        s.out_f64.extend_from_slice(&row);
    }
    nc
}

fn unpack_coords(flat: &[f64]) -> Vec<ResidueCoord> {
    let n = flat.len() / 4;
    let mut out = Vec::with_capacity(n);
    for idx in 0..n {
        out.push(ResidueCoord {
            chain: String::new(),
            resnum: (idx + 1) as i32,
            x: flat[idx * 4],
            y: flat[idx * 4 + 1],
            z: flat[idx * 4 + 2],
            has_p: flat[idx * 4 + 3] != 0.0,
        });
    }
    out
}

/// Compute all directed `(i, j)` pair metrics + `lis_map` + `clis_map` for `universal.html`.
///
/// Outputs:
/// - `out_f32`: `[lis_map (n_total * n_total), clis_map (n_total * n_total)]`
/// - `out_f64`: `30 * n_pairs` numbers per directed pair `(i, j)`
/// - `out_u32`: `8 * n_pairs` `(offset, len)` descriptors followed by 1-based residue index arrays
/// Returns `n_pairs`.
#[no_mangle]
pub extern "C" fn wasm_compute_model_pairs(
    pae_ptr: *const f64,
    pae_dim: u32,
    n_total: u32,
    cb_ptr: *const f64,
    n_cb: u32,
    ca_ptr: *const f64,
    n_ca: u32,
    plddt_ptr: *const f64,
    n_plddt: u32,
    starts_ptr: *const u32,
    sizes_ptr: *const u32,
    nc: u32,
    pair_iptm_ptr: *const f64,
    pae_cutoff: f64,
    cb_cutoff: f64,
) -> u32 {
    let dim = pae_dim as usize;
    let nt = n_total as usize;
    let nc_usize = nc as usize;

    let pae = unsafe { slice_or_empty(pae_ptr, dim * dim) };
    let cb_flat = unsafe { slice_or_empty(cb_ptr, (n_cb as usize) * 4) };
    let ca_flat = unsafe { slice_or_empty(ca_ptr, (n_ca as usize) * 4) };
    let plddt_flat = unsafe { slice_or_empty(plddt_ptr, n_plddt as usize) };
    let starts_u32 = unsafe { slice_or_empty(starts_ptr, nc_usize) };
    let sizes_u32 = unsafe { slice_or_empty(sizes_ptr, nc_usize) };
    let pair_iptm = unsafe { slice_or_empty(pair_iptm_ptr, nc_usize * nc_usize) };

    let cb_coords = unpack_coords(cb_flat);
    let ca_coords = unpack_coords(ca_flat);
    let starts: Vec<usize> = starts_u32.iter().map(|&x| x as usize).collect();
    let sizes: Vec<usize> = sizes_u32.iter().map(|&x| x as usize).collect();

    let mut chain_plddts: Vec<Vec<f64>> = Vec::with_capacity(nc_usize);
    let mut cur = 0usize;
    for &sz in &sizes {
        let end = (cur + sz).min(plddt_flat.len());
        let mut v = Vec::with_capacity(sz);
        if cur < end {
            v.extend_from_slice(&plddt_flat[cur..end]);
        }
        while v.len() < sz {
            v.push(f64::NAN);
        }
        chain_plddts.push(v);
        cur += sz;
    }

    let iptm_lookup = |i: usize, j: usize| -> f64 {
        pair_iptm.get(i * nc_usize + j).copied().unwrap_or(0.0)
    };

    let (directed, lis_map, clis_map) = compute_directed_pairs_kernel(
        pae,
        dim,
        nt,
        &cb_coords,
        &ca_coords,
        &chain_plddts,
        &starts,
        &sizes,
        pae_cutoff,
        cb_cutoff,
        true,
        &iptm_lookup,
    );

    let n_pairs = directed.len();
    let s = scratch_mut();

    s.out_f32.clear();
    s.out_f32.reserve(lis_map.len() + clis_map.len());
    s.out_f32.extend_from_slice(&lis_map);
    s.out_f32.extend_from_slice(&clis_map);

    s.out_f64.clear();
    s.out_f64.reserve(n_pairs * 30);

    let header_words = n_pairs * 8;
    let total_idx_words: usize = directed
        .iter()
        .map(|d| d.lir_i.len() + d.lir_j.len() + d.clir_i.len() + d.clir_j.len())
        .sum();
    s.out_u32.clear();
    s.out_u32.resize(header_words, 0);
    s.out_u32.reserve(total_idx_words);

    let opt_or_nan = |v: Option<f64>| v.unwrap_or(f64::NAN);

    for (p_idx, d) in directed.iter().enumerate() {
        s.out_f64.extend_from_slice(&[
            d.i as f64,
            d.j as f64,
            d.lis,
            d.clis,
            d.ilis,
            d.lia as f64,
            d.clia as f64,
            d.ilia,
            d.ilisa,
            d.iptm,
            d.ipsae,
            d.actifptm,
            d.liplddt_i,
            d.liplddt_j,
            d.cliplddt_i,
            d.cliplddt_j,
            d.geom_plddt_i,
            d.geom_plddt_j,
            opt_or_nan(d.pdockq),
            opt_or_nan(d.lipdockq),
            opt_or_nan(d.pdockq2_i),
            opt_or_nan(d.pdockq2_j),
            opt_or_nan(d.lipdockq2_i),
            opt_or_nan(d.lipdockq2_j),
            d.geom_n_contacts as f64,
            d.lipdockq_n as f64,
            d.geom_if_i_count as f64,
            d.geom_if_j_count as f64,
            d.len_i as f64,
            d.len_j as f64,
        ]);

        let hbase = p_idx * 8;
        let slices = [&d.lir_i, &d.lir_j, &d.clir_i, &d.clir_j];
        for (k, arr) in slices.iter().enumerate() {
            let off = s.out_u32.len() as u32;
            let len = arr.len() as u32;
            s.out_u32[hbase + k * 2] = off;
            s.out_u32[hbase + k * 2 + 1] = len;
            s.out_u32.extend_from_slice(arr);
        }
    }

    n_pairs as u32
}
