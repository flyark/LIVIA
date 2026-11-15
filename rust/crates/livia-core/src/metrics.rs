//! Single-source-of-truth LIS / cLIS / iLIS / ipSAE / actifpTM / pDockQ / LIpDockQ /
//! pDockQ2 / LIpDockQ2 metric engine matching `python/lis.py` and `universal.html`
//! to floating-point precision.

use std::collections::HashMap;

use crate::structure::{
    format_indices, get_chains_from_structure, parse_bfactors_per_residue,
    parse_structure_ca_coords, parse_structure_coords, ResidueCoord,
};

pub const CSV_HEADER: &str = "name,rank,model,chain_i,chain_j,iLIS,iLIA,iLISA,ipSAE,actifpTM,LIS,cLIS,LIA,cLIA,ipTM,pLDDT_i,pLDDT_j,pLDDT,pTM,LIR_i,LIR_j,cLIR_i,cLIR_j,LIpLDDT_i,LIpLDDT_j,LIpLDDT,cLIpLDDT_i,cLIpLDDT_j,cLIpLDDT,pDockQ,LIpDockQ,pDockQ2_i,pDockQ2_j,pDockQ2,LIpDockQ2_i,LIpDockQ2_j,LIpDockQ2,len_i,len_j,LIR_indices_i,LIR_indices_j,cLIR_indices_i,cLIR_indices_j,structure_file";

/// Directed `(i, j)` chain-pair result (before symmetrization), used directly by
/// `universal.html` (`LiviaWasm.computeModelPairs`) AND symmetrized for `lis-rs` CSV output.
#[derive(Clone, Debug)]
pub struct DirectedPairMetrics {
    pub i: usize,
    pub j: usize,
    pub lis: f64,
    pub clis: f64,
    pub ilis: f64,
    pub lia: u32,
    pub clia: u32,
    pub ilia: f64,
    pub ilisa: f64,
    pub iptm: f64,
    pub ipsae: f64,
    pub actifptm: f64,
    pub liplddt_i: f64,
    pub liplddt_j: f64,
    pub cliplddt_i: f64,
    pub cliplddt_j: f64,
    pub geom_plddt_i: f64,
    pub geom_plddt_j: f64,
    pub pdockq: Option<f64>,
    pub lipdockq: Option<f64>,
    pub pdockq2_i: Option<f64>,
    pub pdockq2_j: Option<f64>,
    pub lipdockq2_i: Option<f64>,
    pub lipdockq2_j: Option<f64>,
    pub geom_n_contacts: u32,
    pub lipdockq_n: u32,
    pub lir_i: Vec<u32>,
    pub lir_j: Vec<u32>,
    pub clir_i: Vec<u32>,
    pub clir_j: Vec<u32>,
    pub geom_if_i_count: u32,
    pub geom_if_j_count: u32,
    pub len_i: usize,
    pub len_j: usize,
}

/// Symmetrized chain-pair metrics matching one `lis.py` CSV row.
#[derive(Clone, Debug)]
pub struct SymmetricPairMetrics {
    pub ci: String,
    pub cj: String,
    pub lis: f64,
    pub clis: f64,
    pub ilis: f64,
    pub lia: u32,
    pub clia: u32,
    pub ilia: f64,
    pub ilisa: f64,
    pub iptm: f64,
    pub ipsae: f64,
    pub actifptm: f64,
    pub plddt_i: Option<f64>,
    pub plddt_j: Option<f64>,
    pub ptm: Option<f64>,
    pub liplddt_i: f64,
    pub liplddt_j: f64,
    pub cliplddt_i: f64,
    pub cliplddt_j: f64,
    pub geom_plddt_i: f64,
    pub geom_plddt_j: f64,
    pub geom_if_i_len: usize,
    pub geom_if_j_len: usize,
    pub geom_n_contacts: u32,
    pub lipdockq_n: u32,
    pub pdockq2_i: Option<f64>,
    pub pdockq2_j: Option<f64>,
    pub lipdockq2_i: Option<f64>,
    pub lipdockq2_j: Option<f64>,
    pub lir_i: Vec<u32>,
    pub lir_j: Vec<u32>,
    pub clir_i: Vec<u32>,
    pub clir_j: Vec<u32>,
    pub len_i: usize,
    pub len_j: usize,
}

#[derive(Clone, Debug, Default)]
pub struct ConfidenceInputs {
    pub ptm: Option<f64>,
    pub iptm: Option<f64>,
    pub global_plddt: Option<f64>,
    /// Optional `nc x nc` per-chain-pair ipTM matrix.
    pub chain_pair_iptm: Option<Vec<Vec<Option<f64>>>>,
    /// Optional `(chain_a, chain_b)` named ipTM lookup (OpenFold3 format).
    pub named_pair_iptm: Option<HashMap<(String, String), f64>>,
}

/// Compute `NxN` binary contact map (`1` where `dist - (has_p ? 4.0 : 0.0) <= threshold`).
/// `NaN` coordinates naturally yield `false` (`NaN <= threshold` is false in IEEE-754).
pub fn compute_contact_map(coords: &[ResidueCoord], threshold: f64) -> Vec<u8> {
    let n = coords.len();
    let mut contact = vec![0u8; n * n];
    for i in 0..n {
        let ci = &coords[i];
        if ci.x.is_nan() {
            continue;
        }
        // Diagonal: distance 0 (-4 if has_p) <= threshold when threshold >= 0.
        if (if ci.has_p { -4.0 } else { 0.0 }) <= threshold {
            contact[i * n + i] = 1;
        }
        for j in (i + 1)..n {
            let cj = &coords[j];
            if cj.x.is_nan() {
                continue;
            }
            let dx = ci.x - cj.x;
            let dy = ci.y - cj.y;
            let dz = ci.z - cj.z;
            let mut dist = (dx * dx + dy * dy + dz * dz).sqrt();
            if ci.has_p || cj.has_p {
                dist -= 4.0;
            }
            if dist <= threshold {
                contact[i * n + j] = 1;
                contact[j * n + i] = 1;
            }
        }
    }
    contact
}

/// Scalar trait for PAE buffers (`f32` in compact batch CLI, `f64` in Wasm/JS exact parity).
pub trait PaeScalar: Copy {
    fn to_f64(self) -> f64;
}
impl PaeScalar for f32 {
    #[inline]
    fn to_f64(self) -> f64 {
        self as f64
    }
}
impl PaeScalar for f64 {
    #[inline]
    fn to_f64(self) -> f64 {
        self
    }
}

/// PAE-based per-chain-pair ipTM fallback (`calc_pae_chain_pair_iptm` in `lis.py`).
pub fn calc_pae_chain_pair_iptm<P: PaeScalar>(
    pae: &[P],
    pae_dim: usize,
    starts: &[usize],
    ends: &[usize],
) -> Vec<Vec<f64>> {
    let nc = starts.len();
    let mut out = vec![vec![0.0_f64; nc]; nc];
    for a in 0..nc {
        let sa = starts[a].min(pae_dim);
        let ea = ends[a].min(pae_dim);
        if ea <= sa {
            continue;
        }
        let na = ea - sa;
        for b in 0..nc {
            let sb = starts[b].min(pae_dim);
            let eb = ends[b].min(pae_dim);
            if eb <= sb {
                continue;
            }
            let nb = eb - sb;
            let n_for_d0 = if a == b { na.max(19) } else { (na + nb).max(19) };
            let d0 = 1.24 * ((n_for_d0 as f64 - 15.0).cbrt()) - 1.8;
            if d0 <= 0.0 {
                continue;
            }
            let d0sq = d0 * d0;
            let mut best = 0.0_f64;
            if a == b {
                for i in sa..ea {
                    let row = &pae[i * pae_dim..i * pae_dim + ea];
                    let mut sum = 0.0_f64;
                    for &vp in &row[sa..ea] {
                        let v = vp.to_f64();
                        sum += 1.0 / (1.0 + (v * v) / d0sq);
                    }
                    let avg = sum / (na as f64);
                    if avg > best {
                        best = avg;
                    }
                }
            } else {
                for i in sa..ea {
                    let row = &pae[i * pae_dim..i * pae_dim + eb];
                    let mut sum = 0.0_f64;
                    for &vp in &row[sb..eb] {
                        let v = vp.to_f64();
                        sum += 1.0 / (1.0 + (v * v) / d0sq);
                    }
                    let avg = sum / (nb as f64);
                    if avg > best {
                        best = avg;
                    }
                }
                for i in sb..eb {
                    let row = &pae[i * pae_dim..i * pae_dim + ea];
                    let mut sum = 0.0_f64;
                    for &vp in &row[sa..ea] {
                        let v = vp.to_f64();
                        sum += 1.0 / (1.0 + (v * v) / d0sq);
                    }
                    let avg = sum / (na as f64);
                    if avg > best {
                        best = avg;
                    }
                }
            }
            out[a][b] = best;
        }
    }
    out
}

/// Dunbrack `ipSAE` (`d0res` method, strict `PAE < cutoff`, `cutoff = 10.0`).
pub fn calc_ipsae<P: PaeScalar>(
    pae: &[P],
    pae_dim: usize,
    si: usize,
    ei: usize,
    sj: usize,
    ej: usize,
    pae_cutoff: f64,
) -> f64 {
    let ei = ei.min(pae_dim);
    let ej = ej.min(pae_dim);
    if ei <= si || ej <= sj {
        return 0.0;
    }

    let one_dir = |r_start: usize, r_end: usize, c_start: usize, c_end: usize| -> f64 {
        let mut best = 0.0_f64;
        for r in r_start..r_end {
            let row = &pae[r * pae_dim..];
            let mut good = 0usize;
            for c in c_start..c_end {
                if row[c].to_f64() < pae_cutoff {
                    good += 1;
                }
            }
            if good > 0 {
                let l = good.max(27) as f64;
                let d0 = (1.24 * (l - 15.0).cbrt() - 1.8).max(1.0);
                let d0sq = d0 * d0;
                let mut sum = 0.0_f64;
                for c in c_start..c_end {
                    let v = row[c].to_f64();
                    if v < pae_cutoff {
                        sum += 1.0 / (1.0 + (v * v) / d0sq);
                    }
                }
                let s = sum / (good as f64);
                if s > best {
                    best = s;
                }
            }
        }
        best
    };

    one_dir(si, ei, sj, ej).max(one_dir(sj, ej, si, ei))
}

/// Approximate `actifpTM` (Varga & Ovchinnikov 2025) from PAE + binary `Cb` contact map.
pub fn calc_actifptm<P: PaeScalar>(
    pae: &[P],
    pae_dim: usize,
    contact: &[u8],
    n_coords: usize,
    n_use: usize,
    si: usize,
    ei: usize,
    sj: usize,
    ej: usize,
) -> f64 {
    let clipped = pae_dim.max(19) as f64;
    let d0 = 1.24 * (clipped - 15.0).cbrt() - 1.8;
    let d0sq = d0 * d0;

    let one_dir = |r_start: usize, r_end: usize, s_start: usize, s_end: usize| -> f64 {
        let ri_end = r_end.min(n_use).min(pae_dim);
        let si_end = s_end.min(n_use).min(pae_dim);
        if ri_end <= r_start || si_end <= s_start {
            return 0.0;
        }
        let mut best = 0.0_f64;
        for r in r_start..ri_end {
            let crow = &contact[r * n_coords..];
            let mut wsum = 0usize;
            for s in s_start..si_end {
                if crow[s] != 0 {
                    wsum += 1;
                }
            }
            if wsum == 0 {
                continue;
            }
            let prow = &pae[r * pae_dim..];
            let mut tm_sum = 0.0_f64;
            for s in s_start..si_end {
                if crow[s] != 0 {
                    let p = prow[s].to_f64();
                    tm_sum += 1.0 / (1.0 + (p * p) / d0sq);
                }
            }
            let score = tm_sum / (wsum as f64);
            if score > best {
                best = score;
            }
        }
        best
    };

    one_dir(si, ei, sj, ej).max(one_dir(sj, ej, si, ei))
}

/// Bryant et al. 2022 `pDockQ` / `LIpDockQ` formula (`0.0` when `n_if_contacts == 0`).
#[inline]
pub fn calc_pdockq(avg_if_plddt: Option<f64>, n_if_contacts: u32) -> Option<f64> {
    if n_if_contacts == 0 {
        return Some(0.0);
    }
    let avg = avg_if_plddt?;
    if avg.is_nan() {
        return None;
    }
    let x = avg * (n_if_contacts as f64).log10();
    Some(0.724 / (1.0 + (-0.052 * (x - 152.611)).exp()) + 0.018)
}

/// Zhu et al. 2023 `pDockQ2` / `LIpDockQ2` per-chain logistic transform from separate sums.
#[inline]
pub fn pdockq2_from_sums(decay_sum: f64, plddt_sum: f64, n: usize, empty_chain: bool) -> Option<f64> {
    if empty_chain {
        return None;
    }
    if n == 0 {
        return Some(0.0);
    }
    let x = (decay_sum / (n as f64)) * (plddt_sum / (n as f64));
    Some(1.31034849 / (1.0 + (-0.07471577 * (x - 84.73262390)).exp()) + 0.00501886)
}

#[inline]
pub fn wmean_plddt(vi: Option<f64>, vj: Option<f64>, mut wi: usize, mut wj: usize) -> Option<f64> {
    let ok_i = vi.filter(|v| !v.is_nan());
    let ok_j = vj.filter(|v| !v.is_nan());
    if ok_i.is_none() && ok_j.is_none() {
        return None;
    }
    let a = if let Some(v) = ok_i {
        v
    } else {
        wi = 0;
        0.0
    };
    let b = if let Some(v) = ok_j {
        v
    } else {
        wj = 0;
        0.0
    };
    let w = wi + wj;
    if w == 0 {
        return None;
    }
    Some(((wi as f64) * a + (wj as f64) * b) / (w as f64))
}

#[inline]
pub fn mean2_nullable(a: Option<f64>, b: Option<f64>) -> Option<f64> {
    match (a, b) {
        (None, None) => None,
        (Some(v), None) => Some(v),
        (None, Some(v)) => Some(v),
        (Some(va), Some(vb)) => Some((va + vb) * 0.5),
    }
}

/// Core flat-array kernel shared by BOTH `livia-wasm` (browser) and `lis-rs` (native CLI).
/// Takes pre-flattened arrays so Wasm callers can pass TypedArrays with zero string/JSON overhead:
/// - `pae`: `pae_dim * pae_dim` row-major `f32`
/// - `cb_coords`: `ResidueCoord` slice (Cb / Gly-CA / nucleic-P / HETATM)
/// - `ca_coords`: `ResidueCoord` slice (CA / nucleic-P for `pDockQ2`)
/// - `chain_plddts`: length `nc` array of `Vec<f64>` (`f64::NAN` when a residue lacks a B-factor)
/// - `starts`, `sizes`: length `nc` chain start offsets and lengths
pub fn compute_directed_pairs_kernel<P: PaeScalar>(
    pae: &[P],
    pae_dim: usize,
    n_total: usize,
    cb_coords: &[ResidueCoord],
    ca_coords: &[ResidueCoord],
    chain_plddts: &[Vec<f64>],
    starts: &[usize],
    sizes: &[usize],
    pae_cutoff: f64,
    cb_cutoff: f64,
    f32_map_accum: bool,
    iptm_lookup: &dyn Fn(usize, usize) -> f64,
) -> (Vec<DirectedPairMetrics>, Vec<f32>, Vec<f32>) {
    let n_eff = n_total.min(pae_dim);
    let contact = compute_contact_map(cb_coords, cb_cutoff);
    let n_coords = cb_coords.len();
    let n_use = n_eff.min(n_coords);

    let ca_contact = compute_contact_map(ca_coords, 8.0);
    let n_ca = ca_coords.len();
    let n_use_ca = n_eff.min(n_ca);

    // Build transformed LIS map (`n_total x n_total`, matching `universal.html` display maps & `lis.py` `transformed`).
    let mut lis_map = vec![0.0_f32; n_total * n_total];
    let mut clis_map = vec![0.0_f32; n_total * n_total];
    let mut within = vec![0u8; n_total * n_total];

    for i in 0..n_eff {
        let prow = &pae[i * pae_dim..i * pae_dim + n_eff];
        let lrow_off = i * n_total;
        for j in 0..n_eff {
            let v = prow[j].to_f64();
            if v <= pae_cutoff {
                let tv = (1.0 - v / pae_cutoff) as f32;
                lis_map[lrow_off + j] = tv;
                within[lrow_off + j] = 1;
                if tv > 0.0 && i < n_use && j < n_use && contact[i * n_coords + j] != 0 {
                    clis_map[lrow_off + j] = tv;
                }
            }
        }
    }

    let nc = starts.len();
    let mut out = Vec::with_capacity(nc * nc.saturating_sub(1));

    for i in 0..nc {
        for j in 0..nc {
            if i == j {
                continue;
            }
            let si = starts[i];
            let ei = (starts[i] + sizes[i]).min(n_total);
            let sj = starts[j];
            let ej = (starts[j] + sizes[j]).min(n_total);

            let len_ri = ei.saturating_sub(si);
            let len_rj = ej.saturating_sub(sj);

            let mut lir_i_mask = vec![0u8; len_ri];
            let mut lir_j_mask = vec![0u8; len_rj];
            let mut clir_i_mask = vec![0u8; len_ri];
            let mut clir_j_mask = vec![0u8; len_rj];
            let mut geom_i_mask = vec![0u8; len_ri];
            let mut geom_j_mask = vec![0u8; len_rj];

            let mut lis_sum = 0.0_f64;
            let mut lis_count_avg = 0u32;
            let mut clis_sum = 0.0_f64;
            let mut clis_count_avg = 0u32;
            let mut lis_count_ab = 0u32;
            let mut lis_count_ba = 0u32;
            let mut clis_count_ab = 0u32;
            let mut clis_count_ba = 0u32;
            let mut geom_n_contacts = 0u32;
            let mut lipdockq_n = 0u32;

            let mut pdq2_i_decay = 0.0_f64;
            let mut pdq2_i_plddt = 0.0_f64;
            let mut pdq2_i_n = 0usize;
            let mut pdq2_j_decay = 0.0_f64;
            let mut pdq2_j_plddt = 0.0_f64;
            let mut pdq2_j_n = 0usize;

            let mut cpdq2_i_decay = 0.0_f64;
            let mut cpdq2_i_plddt = 0.0_f64;
            let mut cpdq2_i_n = 0usize;
            let mut cpdq2_j_decay = 0.0_f64;
            let mut cpdq2_j_plddt = 0.0_f64;
            let mut cpdq2_j_n = 0usize;

            let plddt_vec_i = &chain_plddts[i];
            let plddt_vec_j = &chain_plddts[j];

            for ri in si..ei {
                let local_i = ri - si;
                let pi = plddt_vec_i.get(local_i).copied().unwrap_or(f64::NAN);
                let pi_ok = !pi.is_nan();
                for rj in sj..ej {
                    let local_j = rj - sj;
                    let in_contact =
                        ri < n_use && rj < n_use && contact[ri * n_coords + rj] != 0;

                    if within[ri * n_total + rj] != 0 {
                        let tv64 = if f32_map_accum {
                            lis_map[ri * n_total + rj] as f64
                        } else {
                            1.0 - pae[ri * pae_dim + rj].to_f64() / pae_cutoff
                        };
                        lis_sum += tv64;
                        lis_count_avg += 1;
                        if in_contact {
                            clis_sum += tv64;
                            clis_count_avg += 1;
                        }
                    }

                    if in_contact {
                        geom_i_mask[local_i] = 1;
                        geom_j_mask[local_j] = 1;
                        geom_n_contacts += 1;
                    }

                    if ri < pae_dim && rj < pae_dim {
                        let pae_ab = pae[ri * pae_dim + rj].to_f64();
                        let pae_ba = pae[rj * pae_dim + ri].to_f64();
                        let ab_ok = pae_ab <= pae_cutoff;
                        let ba_ok = pae_ba <= pae_cutoff;

                        if ab_ok {
                            lis_count_ab += 1;
                        }
                        if ba_ok {
                            lis_count_ba += 1;
                        }
                        if in_contact {
                            if ab_ok {
                                clis_count_ab += 1;
                            }
                            if ba_ok {
                                clis_count_ba += 1;
                            }
                        }

                        let pj = plddt_vec_j.get(local_j).copied().unwrap_or(f64::NAN);
                        let pj_ok = !pj.is_nan();

                        if ab_ok || ba_ok {
                            lir_i_mask[local_i] = 1;
                            lir_j_mask[local_j] = 1;
                            if in_contact {
                                clir_i_mask[local_i] = 1;
                                clir_j_mask[local_j] = 1;
                                lipdockq_n += 1;
                                if pi_ok {
                                    cpdq2_i_decay += 1.0 / (1.0 + (pae_ab / 10.0).powi(2));
                                    cpdq2_i_plddt += pi;
                                    cpdq2_i_n += 1;
                                }
                                if pj_ok {
                                    cpdq2_j_decay += 1.0 / (1.0 + (pae_ba / 10.0).powi(2));
                                    cpdq2_j_plddt += pj;
                                    cpdq2_j_n += 1;
                                }
                            }
                        }

                        if ri < n_use_ca && rj < n_use_ca && ca_contact[ri * n_ca + rj] != 0 {
                            if pi_ok {
                                pdq2_i_decay += 1.0 / (1.0 + (pae_ab / 10.0).powi(2));
                                pdq2_i_plddt += pi;
                                pdq2_i_n += 1;
                            }
                            if pj_ok {
                                pdq2_j_decay += 1.0 / (1.0 + (pae_ba / 10.0).powi(2));
                                pdq2_j_plddt += pj;
                                pdq2_j_n += 1;
                            }
                        }
                    }
                }
            }

            let c_ei = ei.min(n_use);
            let c_ej = ej.min(n_use);
            let empty_cb = c_ei <= si || c_ej <= sj;
            let ca_ei = ei.min(n_use_ca);
            let ca_ej = ej.min(n_use_ca);
            let empty_ca = ca_ei <= si || ca_ej <= sj;

            let pdockq2_i = pdockq2_from_sums(pdq2_i_decay, pdq2_i_plddt, pdq2_i_n, empty_ca);
            let pdockq2_j = pdockq2_from_sums(pdq2_j_decay, pdq2_j_plddt, pdq2_j_n, empty_ca);
            let lipdockq2_i =
                pdockq2_from_sums(cpdq2_i_decay, cpdq2_i_plddt, cpdq2_i_n, empty_cb);
            let lipdockq2_j =
                pdockq2_from_sums(cpdq2_j_decay, cpdq2_j_plddt, cpdq2_j_n, empty_cb);

            let lis = if lis_count_avg > 0 {
                lis_sum / (lis_count_avg as f64)
            } else {
                0.0
            };
            let clis = if clis_count_avg > 0 {
                clis_sum / (clis_count_avg as f64)
            } else {
                0.0
            };
            let ilis = (lis * clis).sqrt();
            let lia = lis_count_ab + lis_count_ba;
            let clia = clis_count_ab + clis_count_ba;
            let ilia = ((lia as f64) * (clia as f64)).sqrt();
            let ilisa = ilis * ilia;

            let iptm = iptm_lookup(i, j);
            let ipsae = calc_ipsae(pae, pae_dim, si, ei, sj, ej, 10.0);
            let actifptm =
                calc_actifptm(pae, pae_dim, &contact, n_coords, n_use, si, ei, sj, ej);

            let mask_to_indices = |m: &[u8]| -> Vec<u32> {
                m.iter()
                    .enumerate()
                    .filter_map(|(idx, &b)| if b != 0 { Some((idx + 1) as u32) } else { None })
                    .collect()
            };
            let avg_plddt_for = |indices: &[u32], pvec: &[f64]| -> f64 {
                let mut s = 0.0_f64;
                let mut c = 0usize;
                for &r1 in indices {
                    if let Some(&v) = pvec.get((r1 - 1) as usize) {
                        if !v.is_nan() {
                            s += v;
                            c += 1;
                        }
                    }
                }
                if c > 0 {
                    s / (c as f64)
                } else {
                    f64::NAN
                }
            };

            let lir_i = mask_to_indices(&lir_i_mask);
            let lir_j = mask_to_indices(&lir_j_mask);
            let clir_i = mask_to_indices(&clir_i_mask);
            let clir_j = mask_to_indices(&clir_j_mask);
            let geom_i = mask_to_indices(&geom_i_mask);
            let geom_j = mask_to_indices(&geom_j_mask);

            let liplddt_i = avg_plddt_for(&lir_i, plddt_vec_i);
            let liplddt_j = avg_plddt_for(&lir_j, plddt_vec_j);
            let cliplddt_i = avg_plddt_for(&clir_i, plddt_vec_i);
            let cliplddt_j = avg_plddt_for(&clir_j, plddt_vec_j);
            let geom_plddt_i = avg_plddt_for(&geom_i, plddt_vec_i);
            let geom_plddt_j = avg_plddt_for(&geom_j, plddt_vec_j);

            let geom_pair = wmean_plddt(
                Some(geom_plddt_i),
                Some(geom_plddt_j),
                geom_i.len(),
                geom_j.len(),
            );
            let pdockq = calc_pdockq(geom_pair, geom_n_contacts);

            let clir_pair = wmean_plddt(
                Some(cliplddt_i),
                Some(cliplddt_j),
                clir_i.len(),
                clir_j.len(),
            );
            let lipdockq = calc_pdockq(clir_pair, lipdockq_n);

            out.push(DirectedPairMetrics {
                i,
                j,
                lis,
                clis,
                ilis,
                lia,
                clia,
                ilia,
                ilisa,
                iptm,
                ipsae,
                actifptm,
                liplddt_i,
                liplddt_j,
                cliplddt_i,
                cliplddt_j,
                geom_plddt_i,
                geom_plddt_j,
                pdockq,
                lipdockq,
                pdockq2_i,
                pdockq2_j,
                lipdockq2_i,
                lipdockq2_j,
                geom_n_contacts,
                lipdockq_n,
                lir_i,
                lir_j,
                clir_i,
                clir_j,
                geom_if_i_count: geom_i.len() as u32,
                geom_if_j_count: geom_j.len() as u32,
                len_i: sizes[i],
                len_j: sizes[j],
            });
        }
    }

    (out, lis_map, clis_map)
}

/// Full single-model analysis matching `lis.py` `analyze_single_model`.
pub fn analyze_single_model(
    struct_text: &str,
    is_pdb: bool,
    pae: &[f32],
    pae_dim: usize,
    token_chain_ids: Option<&[String]>,
    scores: &ConfidenceInputs,
    pae_cutoff: f64,
    cb_cutoff: f64,
) -> Vec<SymmetricPairMetrics> {
    let chain_info = get_chains_from_structure(struct_text, is_pdb);
    let (chain_names, sizes, mut n_total) = if let Some(tcids) = token_chain_ids {
        let mut names: Vec<String> = Vec::new();
        let mut counts: HashMap<String, usize> = HashMap::new();
        for c in tcids {
            if let Some(cnt) = counts.get_mut(c) {
                *cnt += 1;
            } else {
                names.push(c.clone());
                counts.insert(c.clone(), 1);
            }
        }
        let sz: Vec<usize> = names.iter().map(|n| counts[n]).collect();
        (names, sz, tcids.len())
    } else {
        let nt: usize = chain_info.sizes.iter().sum();
        (chain_info.names, chain_info.sizes, nt)
    };

    if pae_dim != n_total {
        n_total = pae_dim;
    }

    let nc = chain_names.len();
    if nc < 2 {
        return Vec::new();
    }

    let mut starts = Vec::with_capacity(nc);
    let mut ends = Vec::with_capacity(nc);
    let mut acc = 0usize;
    for &s in &sizes {
        starts.push(acc);
        acc += s;
        ends.push(acc);
    }

    let cb_coords = parse_structure_coords(struct_text, is_pdb);
    let ca_coords = parse_structure_ca_coords(struct_text, is_pdb);
    let bfs = parse_bfactors_per_residue(struct_text, is_pdb);

    // Per-chain mean pLDDT and per-residue indexed arrays.
    let mut chain_mean_plddt: HashMap<String, f64> = HashMap::new();
    let mut chain_vals_map: HashMap<&str, (f64, usize)> = HashMap::new();
    for (k, &val) in &bfs {
        if let Some((ch, _)) = k.split_once(':') {
            let entry = chain_vals_map.entry(ch).or_insert((0.0, 0));
            entry.0 += val;
            entry.1 += 1;
        }
    }
    for (ch, (s, c)) in chain_vals_map {
        if c > 0 {
            chain_mean_plddt.insert(ch.to_string(), s / (c as f64));
        }
    }

    let mut chain_plddts: Vec<Vec<f64>> = Vec::with_capacity(nc);
    for (idx, cname) in chain_names.iter().enumerate() {
        let len = sizes[idx];
        let mut v = vec![f64::NAN; len];
        for p in 0..len {
            if let Some(&bf) = bfs.get(&format!("{cname}:{}", p + 1)) {
                v[p] = bf;
            }
        }
        chain_plddts.push(v);
    }

    let fallback_iptm = if scores.chain_pair_iptm.is_none()
        && scores.named_pair_iptm.is_none()
        && nc >= 3
    {
        Some(calc_pae_chain_pair_iptm(pae, pae_dim, &starts, &ends))
    } else {
        None
    };
    let global_iptm = scores.iptm.unwrap_or(0.0);

    let iptm_lookup = |i: usize, j: usize| -> f64 {
        if let Some(ref named) = scores.named_pair_iptm {
            if let Some(&v) = named.get(&(chain_names[i].clone(), chain_names[j].clone())) {
                return v;
            }
        }
        if let Some(ref mat) = scores.chain_pair_iptm {
            if let Some(Some(v)) = mat.get(i).and_then(|r| r.get(j)).copied() {
                return v;
            }
        }
        if let Some(ref est) = fallback_iptm {
            return est[i][j];
        }
        global_iptm
    };

    let (directed, _, _) = compute_directed_pairs_kernel(
        pae,
        pae_dim,
        n_total,
        &cb_coords,
        &ca_coords,
        &chain_plddts,
        &starts,
        &sizes,
        pae_cutoff,
        cb_cutoff,
        false,
        &iptm_lookup,
    );

    let mut by_ij: HashMap<(usize, usize), DirectedPairMetrics> = HashMap::new();
    for d in directed {
        by_ij.insert((d.i, d.j), d);
    }

    let mut seen: std::collections::HashSet<(String, String)> = std::collections::HashSet::new();
    let mut keyed: Vec<(String, SymmetricPairMetrics)> = Vec::new();

    for i in 0..nc {
        for j in 0..nc {
            if i == j {
                continue;
            }
            let ci = &chain_names[i];
            let cj = &chain_names[j];
            let canon = if ci <= cj {
                (ci.clone(), cj.clone())
            } else {
                (cj.clone(), ci.clone())
            };
            if !seen.insert(canon) {
                continue;
            }
            let val = &by_ij[&(i, j)];
            let rv = &by_ij[&(j, i)];

            let lis = (val.lis + rv.lis) * 0.5;
            let clis = (val.clis + rv.clis) * 0.5;
            let ilis = (lis * clis).sqrt();
            let lia = val.lia;
            let clia = val.clia;
            let ilia = ((lia as f64) * (clia as f64)).sqrt();
            let ilisa = ilis * ilia;
            let iptm = (val.iptm + rv.iptm) * 0.5;
            let ipsae = val.ipsae.max(rv.ipsae);
            let actifptm = val.actifptm.max(rv.actifptm);

            let union_sorted = |a: &[u32], b: &[u32]| -> Vec<u32> {
                let mut u = Vec::with_capacity(a.len() + b.len());
                u.extend_from_slice(a);
                u.extend_from_slice(b);
                u.sort_unstable();
                u.dedup();
                u
            };

            let lir_i = union_sorted(&val.lir_i, &rv.lir_j);
            let lir_j = union_sorted(&val.lir_j, &rv.lir_i);
            let clir_i = union_sorted(&val.clir_i, &rv.clir_j);
            let clir_j = union_sorted(&val.clir_j, &rv.clir_i);

            let plddt_i = chain_mean_plddt.get(ci).copied().or(scores.global_plddt);
            let plddt_j = chain_mean_plddt.get(cj).copied().or(scores.global_plddt);

            keyed.push((
                format!("{ci}-{cj}"),
                SymmetricPairMetrics {
                    ci: ci.clone(),
                    cj: cj.clone(),
                    lis,
                    clis,
                    ilis,
                    lia,
                    clia,
                    ilia,
                    ilisa,
                    iptm,
                    ipsae,
                    actifptm,
                    plddt_i,
                    plddt_j,
                    ptm: scores.ptm,
                    liplddt_i: val.liplddt_i,
                    liplddt_j: rv.liplddt_i,
                    cliplddt_i: val.cliplddt_i,
                    cliplddt_j: rv.cliplddt_i,
                    geom_plddt_i: val.geom_plddt_i,
                    geom_plddt_j: rv.geom_plddt_i,
                    geom_if_i_len: val.geom_if_i_count as usize,
                    geom_if_j_len: val.geom_if_j_count as usize,
                    geom_n_contacts: val.geom_n_contacts,
                    lipdockq_n: val.lipdockq_n,
                    pdockq2_i: val.pdockq2_i,
                    pdockq2_j: rv.pdockq2_i,
                    lipdockq2_i: val.lipdockq2_i,
                    lipdockq2_j: rv.lipdockq2_i,
                    lir_i,
                    lir_j,
                    clir_i,
                    clir_j,
                    len_i: val.len_i,
                    len_j: val.len_j,
                },
            ));
        }
    }

    keyed.sort_by(|a, b| a.0.cmp(&b.0));
    keyed.into_iter().map(|(_, v)| v).collect()
}

pub fn extract_model_num(struct_filename: &str) -> String {
    for prefix in ["model_idx_", "model_", "sample_"] {
        if let Some(pos) = struct_filename.find(prefix) {
            let rest = &struct_filename[pos + prefix.len()..];
            let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
            if !digits.is_empty() {
                // Note: check `model_` before `model_idx_` to mirror `lis.py`'s exact regex order!
                if prefix == "model_idx_" && struct_filename.contains("model_") {
                    // handled below
                } else {
                    return digits;
                }
            }
        }
    }
    for prefix in ["model_", "sample_", "model_idx_"] {
        let mut slice = struct_filename;
        while let Some(pos) = slice.find(prefix) {
            let after = &slice[pos + prefix.len()..];
            let digits: String = after.chars().take_while(|c| c.is_ascii_digit()).collect();
            if !digits.is_empty() {
                return digits;
            }
            slice = after;
        }
    }
    String::new()
}

pub fn format_csv_row(
    name: &str,
    rank: &str,
    struct_file: &str,
    pair: &SymmetricPairMetrics,
) -> String {
    let fmt_plddt = |v: Option<f64>| -> String {
        match v {
            Some(x) if !x.is_nan() => format!("{x:.1}"),
            _ => String::new(),
        }
    };
    let fmt_f4 = |v: Option<f64>| -> String {
        match v {
            Some(x) if !x.is_nan() => format!("{x:.4}"),
            _ => String::new(),
        }
    };

    let model_num = extract_model_num(struct_file);
    let plddt_pair = wmean_plddt(pair.plddt_i, pair.plddt_j, pair.len_i, pair.len_j);
    let liplddt_pair = wmean_plddt(
        Some(pair.liplddt_i),
        Some(pair.liplddt_j),
        pair.lir_i.len(),
        pair.lir_j.len(),
    );
    let cliplddt_pair = wmean_plddt(
        Some(pair.cliplddt_i),
        Some(pair.cliplddt_j),
        pair.clir_i.len(),
        pair.clir_j.len(),
    );
    let geom_pair = wmean_plddt(
        Some(pair.geom_plddt_i),
        Some(pair.geom_plddt_j),
        pair.geom_if_i_len,
        pair.geom_if_j_len,
    );

    let pdockq_orig = calc_pdockq(geom_pair, pair.geom_n_contacts);
    let lipdockq = calc_pdockq(cliplddt_pair, pair.lipdockq_n);
    let pdockq2_pair = mean2_nullable(pair.pdockq2_i, pair.pdockq2_j);
    let lipdockq2_pair = mean2_nullable(pair.lipdockq2_i, pair.lipdockq2_j);

    let cols = [
        name.to_string(),
        rank.to_string(),
        model_num,
        pair.ci.clone(),
        pair.cj.clone(),
        format!("{:.4}", pair.ilis),
        format!("{:.1}", pair.ilia),
        format!("{:.1}", pair.ilisa),
        format!("{:.4}", pair.ipsae),
        format!("{:.4}", pair.actifptm),
        format!("{:.4}", pair.lis),
        format!("{:.4}", pair.clis),
        format!("{:.1}", pair.lia as f64),
        format!("{:.1}", pair.clia as f64),
        format!("{:.4}", pair.iptm),
        fmt_plddt(pair.plddt_i),
        fmt_plddt(pair.plddt_j),
        fmt_plddt(plddt_pair),
        pair.ptm
            .filter(|v| !v.is_nan())
            .map(|v| format!("{v:.3}"))
            .unwrap_or_default(),
        pair.lir_i.len().to_string(),
        pair.lir_j.len().to_string(),
        pair.clir_i.len().to_string(),
        pair.clir_j.len().to_string(),
        fmt_plddt(Some(pair.liplddt_i)),
        fmt_plddt(Some(pair.liplddt_j)),
        fmt_plddt(liplddt_pair),
        fmt_plddt(Some(pair.cliplddt_i)),
        fmt_plddt(Some(pair.cliplddt_j)),
        fmt_plddt(cliplddt_pair),
        fmt_f4(pdockq_orig),
        fmt_f4(lipdockq),
        fmt_f4(pair.pdockq2_i),
        fmt_f4(pair.pdockq2_j),
        fmt_f4(pdockq2_pair),
        fmt_f4(pair.lipdockq2_i),
        fmt_f4(pair.lipdockq2_j),
        fmt_f4(lipdockq2_pair),
        pair.len_i.to_string(),
        pair.len_j.to_string(),
        format_indices(&pair.lir_i),
        format_indices(&pair.lir_j),
        format_indices(&pair.clir_i),
        format_indices(&pair.clir_j),
        struct_file.to_string(),
    ];
    cols.join(",")
}
