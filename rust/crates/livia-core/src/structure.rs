//! PDB and mmCIF dynamic-header coordinate, chain-boundary, and B-factor parsers
//! matching `python/lis.py` and `universal.html`.

use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct ResidueCoord {
    pub chain: String,
    pub resnum: i32,
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub has_p: bool,
}

#[derive(Clone, Debug)]
pub struct ChainInfo {
    pub names: Vec<String>,
    pub sizes: Vec<usize>,
    pub types: Vec<String>,
}

#[inline]
pub fn is_ion_name(comp_id: &str) -> bool {
    matches!(
        comp_id,
        "ZN" | "CA" | "MG" | "MN" | "FE" | "CU" | "NA" | "K" | "CL" | "NI" | "CO" | "CD"
    )
}

#[inline]
pub fn is_mod_aa(comp_id: &str) -> bool {
    matches!(
        comp_id,
        "SEP" | "TPO" | "PTR" | "MSE" | "CSO" | "HYP" | "MLY" | "M3L" | "ALY" | "PCA" | "KCX"
            | "LLP" | "CME" | "CSD" | "OCS" | "SAH" | "SAM" | "FME" | "ACE" | "NME"
    )
}

struct OrderedResidues {
    keys: Vec<String>,
    index: HashMap<String, usize>,
    vals: Vec<ResidueCoord>,
}

impl OrderedResidues {
    fn new() -> Self {
        Self {
            keys: Vec::new(),
            index: HashMap::new(),
            vals: Vec::new(),
        }
    }

    fn contains(&self, key: &str) -> bool {
        self.index.contains_key(key)
    }

    fn set(&mut self, key: String, val: ResidueCoord) {
        if let Some(&idx) = self.index.get(&key) {
            self.vals[idx] = val;
        } else {
            let idx = self.vals.len();
            self.index.insert(key.clone(), idx);
            self.keys.push(key);
            self.vals.push(val);
        }
    }
}

/// Fill missing residue numbers inside each polymer chain's `[min, max]` range with `NaN`
/// placeholders, matching `lis.py` `_fill_residue_gaps` so ColabFold `'X'` residue gaps
/// stay aligned with the PAE matrix.
fn fill_residue_gaps(ord: OrderedResidues) -> Vec<ResidueCoord> {
    let mut chain_order: Vec<String> = Vec::new();
    let mut polymer: HashMap<String, HashMap<i32, ResidueCoord>> = HashMap::new();
    let mut het: Vec<ResidueCoord> = Vec::new();

    for (key, r) in ord.keys.into_iter().zip(ord.vals) {
        if key.starts_with("het:") {
            het.push(r);
        } else {
            if !polymer.contains_key(&r.chain) {
                chain_order.push(r.chain.clone());
                polymer.insert(r.chain.clone(), HashMap::new());
            }
            polymer.get_mut(&r.chain).unwrap().insert(r.resnum, r);
        }
    }

    let mut out = Vec::new();
    for chain in chain_order {
        let byres = &polymer[&chain];
        let lo = *byres.keys().min().unwrap();
        let hi = *byres.keys().max().unwrap();
        for n in lo..=hi {
            if let Some(r) = byres.get(&n) {
                out.push(r.clone());
            } else {
                out.push(ResidueCoord {
                    chain: chain.clone(),
                    resnum: n,
                    x: f64::NAN,
                    y: f64::NAN,
                    z: f64::NAN,
                    has_p: false,
                });
            }
        }
    }
    out.extend(het);
    out
}

pub fn parse_pdb_coords(pdb_text: &str) -> Vec<ResidueCoord> {
    let mut residues = OrderedResidues::new();
    let mut het_idx = 0usize;

    for line in pdb_text.lines() {
        let is_atom = line.starts_with("ATOM");
        let is_het = !is_atom && line.starts_with("HETATM");
        if (!is_atom && !is_het) || line.len() < 54 {
            continue;
        }
        let Some(atom_s) = line.get(12..16) else { continue };
        let Some(comp_s) = line.get(17..20) else { continue };
        let Some(chain_s) = line.get(21..22) else { continue };
        let Some(xs) = line.get(30..38) else { continue };
        let Some(ys) = line.get(38..46) else { continue };
        let Some(zs) = line.get(46..54) else { continue };

        let atom_name = atom_s.trim();
        let comp_id = comp_s.trim();
        let chain_trim = chain_s.trim();
        let chain = if chain_trim.is_empty() { "A" } else { chain_trim };

        let (Ok(x), Ok(y), Ok(z)) = (
            xs.trim().parse::<f64>(),
            ys.trim().parse::<f64>(),
            zs.trim().parse::<f64>(),
        ) else {
            continue;
        };

        if is_het {
            if is_ion_name(comp_id) {
                let k = format!("{chain}:1");
                if !residues.contains(&k) {
                    residues.set(
                        k,
                        ResidueCoord {
                            chain: chain.to_string(),
                            resnum: 1,
                            x,
                            y,
                            z,
                            has_p: false,
                        },
                    );
                }
            } else {
                let k = format!("het:{chain}:{het_idx}");
                het_idx += 1;
                residues.set(
                    k,
                    ResidueCoord {
                        chain: chain.to_string(),
                        resnum: 1,
                        x,
                        y,
                        z,
                        has_p: atom_name.contains('P'),
                    },
                );
            }
            continue;
        }

        let Some(res_s) = line.get(22..26) else { continue };
        let Ok(resnum) = res_s.trim().parse::<i32>() else { continue };
        let key = format!("{chain}:{resnum}");

        if atom_name == "CB" {
            residues.set(
                key,
                ResidueCoord {
                    chain: chain.to_string(),
                    resnum,
                    x,
                    y,
                    z,
                    has_p: false,
                },
            );
        } else if atom_name == "CA" && comp_id == "GLY" && !residues.contains(&key) {
            residues.set(
                key,
                ResidueCoord {
                    chain: chain.to_string(),
                    resnum,
                    x,
                    y,
                    z,
                    has_p: false,
                },
            );
        } else if atom_name == "P" && !residues.contains(&key) {
            residues.set(
                key,
                ResidueCoord {
                    chain: chain.to_string(),
                    resnum,
                    x,
                    y,
                    z,
                    has_p: true,
                },
            );
        }
    }

    fill_residue_gaps(residues)
}

#[derive(Default)]
struct CifCols {
    group_pdb: Option<usize>,
    atom_id: Option<usize>,
    comp_id: Option<usize>,
    asym_id: Option<usize>,
    seq_id: Option<usize>,
    x: Option<usize>,
    y: Option<usize>,
    z: Option<usize>,
    bfac: Option<usize>,
}

impl CifCols {
    fn from_names(names: &[&str]) -> Self {
        let mut c = Self::default();
        for (i, &n) in names.iter().enumerate() {
            match n {
                "group_PDB" => c.group_pdb = Some(i),
                "label_atom_id" => c.atom_id = Some(i),
                "label_comp_id" => c.comp_id = Some(i),
                "label_asym_id" => c.asym_id = Some(i),
                "label_seq_id" => c.seq_id = Some(i),
                "Cartn_x" => c.x = Some(i),
                "Cartn_y" => c.y = Some(i),
                "Cartn_z" => c.z = Some(i),
                "B_iso_or_equiv" => c.bfac = Some(i),
                _ => {}
            }
        }
        c
    }
}

#[inline]
fn get_tok<'a>(parts: &[&'a str], idx: Option<usize>) -> &'a str {
    idx.and_then(|i| parts.get(i).copied()).unwrap_or("")
}

pub fn parse_cif_coords(cif_text: &str) -> Vec<ResidueCoord> {
    let mut residues = OrderedResidues::new();
    let mut het_idx = 0usize;
    let mut in_atom_site = false;
    let mut col_names: Vec<&str> = Vec::new();
    let mut cols = CifCols::default();
    let mut cols_ready = false;

    for raw_line in cif_text.lines() {
        if let Some(rest) = raw_line.strip_prefix("_atom_site.") {
            if !in_atom_site {
                in_atom_site = true;
                col_names.clear();
                cols_ready = false;
            }
            col_names.push(rest.trim());
            continue;
        }
        if in_atom_site {
            let trimmed = raw_line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            if trimmed.starts_with("loop_") || trimmed.starts_with('_') {
                in_atom_site = false;
                continue;
            }
            if !cols_ready {
                cols = CifCols::from_names(&col_names);
                cols_ready = true;
            }
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.is_empty() || (parts[0] != "ATOM" && parts[0] != "HETATM") {
                continue;
            }
            let group_pdb = get_tok(&parts, cols.group_pdb);
            if group_pdb != "ATOM" && group_pdb != "HETATM" {
                continue;
            }
            let atom_name = get_tok(&parts, cols.atom_id);
            let comp_id = get_tok(&parts, cols.comp_id);
            let chain = get_tok(&parts, cols.asym_id);
            let res_seq = get_tok(&parts, cols.seq_id);

            let (Ok(x), Ok(y), Ok(z)) = (
                get_tok(&parts, cols.x).parse::<f64>(),
                get_tok(&parts, cols.y).parse::<f64>(),
                get_tok(&parts, cols.z).parse::<f64>(),
            ) else {
                continue;
            };

            if group_pdb == "HETATM" {
                if is_ion_name(comp_id) {
                    let k = format!("{chain}:1");
                    if !residues.contains(&k) {
                        residues.set(
                            k,
                            ResidueCoord {
                                chain: chain.to_string(),
                                resnum: 1,
                                x,
                                y,
                                z,
                                has_p: false,
                            },
                        );
                    }
                } else {
                    let k = format!("het:{chain}:{het_idx}");
                    het_idx += 1;
                    residues.set(
                        k,
                        ResidueCoord {
                            chain: chain.to_string(),
                            resnum: 1,
                            x,
                            y,
                            z,
                            has_p: atom_name.contains('P'),
                        },
                    );
                }
                continue;
            }

            let Ok(resnum) = res_seq.parse::<i32>() else { continue };
            let key = format!("{chain}:{resnum}");
            if atom_name == "CB" {
                residues.set(
                    key,
                    ResidueCoord {
                        chain: chain.to_string(),
                        resnum,
                        x,
                        y,
                        z,
                        has_p: false,
                    },
                );
            } else if atom_name == "CA" && comp_id == "GLY" && !residues.contains(&key) {
                residues.set(
                    key,
                    ResidueCoord {
                        chain: chain.to_string(),
                        resnum,
                        x,
                        y,
                        z,
                        has_p: false,
                    },
                );
            } else if atom_name == "P" && !residues.contains(&key) {
                residues.set(
                    key,
                    ResidueCoord {
                        chain: chain.to_string(),
                        resnum,
                        x,
                        y,
                        z,
                        has_p: true,
                    },
                );
            }
        }
    }
    residues.vals
}

pub fn parse_structure_coords(text: &str, is_pdb: bool) -> Vec<ResidueCoord> {
    if is_pdb {
        parse_pdb_coords(text)
    } else {
        parse_cif_coords(text)
    }
}

pub fn parse_pdb_ca_coords(pdb_text: &str) -> Vec<ResidueCoord> {
    let mut residues = OrderedResidues::new();
    for line in pdb_text.lines() {
        if (!line.starts_with("ATOM") && !line.starts_with("HETATM")) || line.len() < 54 {
            continue;
        }
        let Some(atom_s) = line.get(12..16) else { continue };
        let Some(chain_s) = line.get(21..22) else { continue };
        let Some(res_s) = line.get(22..26) else { continue };
        let Some(xs) = line.get(30..38) else { continue };
        let Some(ys) = line.get(38..46) else { continue };
        let Some(zs) = line.get(46..54) else { continue };

        let atom_name = atom_s.trim();
        if atom_name != "CA" && atom_name != "P" {
            continue;
        }
        let chain_trim = chain_s.trim();
        let chain = if chain_trim.is_empty() { "A" } else { chain_trim };
        let Ok(resnum) = res_s.trim().parse::<i32>() else { continue };
        let (Ok(x), Ok(y), Ok(z)) = (
            xs.trim().parse::<f64>(),
            ys.trim().parse::<f64>(),
            zs.trim().parse::<f64>(),
        ) else {
            continue;
        };
        let key = format!("{chain}:{resnum}");
        if atom_name == "CA" {
            residues.set(
                key,
                ResidueCoord {
                    chain: chain.to_string(),
                    resnum,
                    x,
                    y,
                    z,
                    has_p: false,
                },
            );
        } else if atom_name == "P" && !residues.contains(&key) {
            residues.set(
                key,
                ResidueCoord {
                    chain: chain.to_string(),
                    resnum,
                    x,
                    y,
                    z,
                    has_p: true,
                },
            );
        }
    }
    residues.vals
}

pub fn parse_cif_ca_coords(cif_text: &str) -> Vec<ResidueCoord> {
    let mut residues = OrderedResidues::new();
    let mut in_atom_site = false;
    let mut col_names: Vec<&str> = Vec::new();
    let mut cols = CifCols::default();
    let mut cols_ready = false;

    for raw_line in cif_text.lines() {
        if let Some(rest) = raw_line.strip_prefix("_atom_site.") {
            if !in_atom_site {
                in_atom_site = true;
                col_names.clear();
                cols_ready = false;
            }
            col_names.push(rest.trim());
            continue;
        }
        if in_atom_site {
            let trimmed = raw_line.trim();
            if trimmed.is_empty() || trimmed.starts_with('#') {
                continue;
            }
            if trimmed.starts_with("loop_") || trimmed.starts_with('_') {
                in_atom_site = false;
                continue;
            }
            if !cols_ready {
                cols = CifCols::from_names(&col_names);
                cols_ready = true;
            }
            let parts: Vec<&str> = trimmed.split_whitespace().collect();
            if parts.is_empty() || (parts[0] != "ATOM" && parts[0] != "HETATM") {
                continue;
            }
            let group_pdb = get_tok(&parts, cols.group_pdb);
            if group_pdb != "ATOM" && group_pdb != "HETATM" {
                continue;
            }
            let atom_name = get_tok(&parts, cols.atom_id);
            if atom_name != "CA" && atom_name != "P" {
                continue;
            }
            let chain = get_tok(&parts, cols.asym_id);
            let Ok(resnum) = get_tok(&parts, cols.seq_id).parse::<i32>() else { continue };
            let (Ok(x), Ok(y), Ok(z)) = (
                get_tok(&parts, cols.x).parse::<f64>(),
                get_tok(&parts, cols.y).parse::<f64>(),
                get_tok(&parts, cols.z).parse::<f64>(),
            ) else {
                continue;
            };
            let key = format!("{chain}:{resnum}");
            if atom_name == "CA" {
                residues.set(
                    key,
                    ResidueCoord {
                        chain: chain.to_string(),
                        resnum,
                        x,
                        y,
                        z,
                        has_p: false,
                    },
                );
            } else if !residues.contains(&key) {
                residues.set(
                    key,
                    ResidueCoord {
                        chain: chain.to_string(),
                        resnum,
                        x,
                        y,
                        z,
                        has_p: true,
                    },
                );
            }
        }
    }
    residues.vals
}

pub fn parse_structure_ca_coords(text: &str, is_pdb: bool) -> Vec<ResidueCoord> {
    if is_pdb {
        parse_pdb_ca_coords(text)
    } else {
        parse_cif_ca_coords(text)
    }
}

pub fn get_chains_from_structure(text: &str, is_pdb: bool) -> ChainInfo {
    if is_pdb {
        let coords = parse_pdb_coords(text);
        let mut names: Vec<String> = Vec::new();
        let mut counts: HashMap<String, usize> = HashMap::new();
        for c in coords {
            if let Some(cnt) = counts.get_mut(&c.chain) {
                *cnt += 1;
            } else {
                names.push(c.chain.clone());
                counts.insert(c.chain, 1);
            }
        }
        let sizes = names.iter().map(|n| counts[n]).collect();
        let types = vec!["protein".to_string(); names.len()];
        ChainInfo { names, sizes, types }
    } else {
        let mut names: Vec<String> = Vec::new();
        let mut counts: HashMap<String, usize> = HashMap::new();
        let mut types_map: HashMap<String, String> = HashMap::new();
        let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
        let mut in_atom_site = false;
        let mut col_names: Vec<&str> = Vec::new();
        let mut cols = CifCols::default();
        let mut cols_ready = false;

        for raw_line in text.lines() {
            if let Some(rest) = raw_line.strip_prefix("_atom_site.") {
                if !in_atom_site {
                    in_atom_site = true;
                    col_names.clear();
                    cols_ready = false;
                }
                col_names.push(rest.trim());
                continue;
            }
            if in_atom_site {
                let trimmed = raw_line.trim();
                if trimmed.is_empty() || trimmed.starts_with('#') {
                    continue;
                }
                if trimmed.starts_with("loop_") || trimmed.starts_with('_') {
                    in_atom_site = false;
                    continue;
                }
                if !cols_ready {
                    cols = CifCols::from_names(&col_names);
                    cols_ready = true;
                }
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.is_empty() || (parts[0] != "ATOM" && parts[0] != "HETATM") {
                    continue;
                }
                let group_pdb = get_tok(&parts, cols.group_pdb);
                let atom_name = get_tok(&parts, cols.atom_id);
                let comp_id = get_tok(&parts, cols.comp_id);
                let chain = get_tok(&parts, cols.asym_id);
                let res_seq = get_tok(&parts, cols.seq_id);

                let mut counted = false;
                if group_pdb == "ATOM" && (atom_name == "CA" || atom_name == "P") {
                    let rkey = format!("{chain}:{res_seq}");
                    if seen.insert(rkey) {
                        counted = true;
                        if atom_name == "P" {
                            types_map.insert(
                                chain.to_string(),
                                if comp_id.starts_with('D') { "dna" } else { "rna" }.to_string(),
                            );
                        } else {
                            types_map
                                .entry(chain.to_string())
                                .or_insert_with(|| "protein".to_string());
                        }
                    }
                } else if group_pdb == "HETATM" && is_ion_name(comp_id) {
                    let rkey = format!("{chain}:ion");
                    if seen.insert(rkey) {
                        counted = true;
                        types_map.insert(chain.to_string(), "ion".to_string());
                    }
                } else if group_pdb == "HETATM" {
                    counted = true;
                    if atom_name == "CA" || is_mod_aa(comp_id) {
                        types_map.insert(chain.to_string(), "protein".to_string());
                    } else {
                        types_map
                            .entry(chain.to_string())
                            .or_insert_with(|| "glycan".to_string());
                    }
                }

                if counted {
                    if let Some(cnt) = counts.get_mut(chain) {
                        *cnt += 1;
                    } else {
                        names.push(chain.to_string());
                        counts.insert(chain.to_string(), 1);
                    }
                }
            }
        }
        let sizes = names.iter().map(|n| counts[n]).collect();
        let types = names
            .iter()
            .map(|n| types_map.get(n).cloned().unwrap_or_else(|| "protein".to_string()))
            .collect();
        ChainInfo { names, sizes, types }
    }
}

pub fn parse_bfactors_per_residue(text: &str, is_pdb: bool) -> HashMap<String, f64> {
    let mut bfactors: HashMap<String, f64> = HashMap::new();
    let mut max_bf = f64::NEG_INFINITY;

    if is_pdb {
        for line in text.lines() {
            if !line.starts_with("ATOM") || line.len() < 66 {
                continue;
            }
            let Some(atom_s) = line.get(12..16) else { continue };
            if atom_s.trim() != "CA" {
                continue;
            }
            let Some(chain_s) = line.get(21..22) else { continue };
            let Some(res_s) = line.get(22..26) else { continue };
            let Some(bf_s) = line.get(60..66) else { continue };
            let chain_trim = chain_s.trim();
            let chain = if chain_trim.is_empty() { "A" } else { chain_trim };
            let (Ok(rn), Ok(bf)) = (res_s.trim().parse::<i32>(), bf_s.trim().parse::<f64>())
            else {
                continue;
            };
            if bf > max_bf {
                max_bf = bf;
            }
            bfactors.insert(format!("{chain}:{rn}"), bf);
        }
    } else {
        let mut in_atom_site = false;
        let mut col_names: Vec<&str> = Vec::new();
        let mut cols = CifCols::default();
        let mut cols_ready = false;

        for raw_line in text.lines() {
            if let Some(rest) = raw_line.strip_prefix("_atom_site.") {
                if !in_atom_site {
                    in_atom_site = true;
                    col_names.clear();
                    cols_ready = false;
                }
                col_names.push(rest.trim());
                continue;
            }
            if in_atom_site {
                let trimmed = raw_line.trim();
                if trimmed.is_empty() || trimmed.starts_with('#') {
                    continue;
                }
                if trimmed.starts_with("loop_") || trimmed.starts_with('_') {
                    in_atom_site = false;
                    continue;
                }
                if !cols_ready {
                    cols = CifCols::from_names(&col_names);
                    cols_ready = true;
                }
                let parts: Vec<&str> = trimmed.split_whitespace().collect();
                if parts.is_empty() || parts[0] != "ATOM" {
                    continue;
                }
                if get_tok(&parts, cols.group_pdb) != "ATOM"
                    || get_tok(&parts, cols.atom_id) != "CA"
                {
                    continue;
                }
                let chain = get_tok(&parts, cols.asym_id);
                let (Ok(rn), Ok(bf)) = (
                    get_tok(&parts, cols.seq_id).parse::<i32>(),
                    get_tok(&parts, cols.bfac).parse::<f64>(),
                ) else {
                    continue;
                };
                if bf > max_bf {
                    max_bf = bf;
                }
                bfactors.insert(format!("{chain}:{rn}"), bf);
            }
        }
    }

    if !bfactors.is_empty() && max_bf > 0.0 && max_bf <= 1.0 {
        for v in bfactors.values_mut() {
            *v *= 100.0;
        }
    }
    bfactors
}

/// Format sorted 1-based residue indices as `lis.py`'s `"[1-4,10-11,13]"` CSV cell.
pub fn format_indices(sorted_pos: &[u32]) -> String {
    if sorted_pos.is_empty() {
        return String::new();
    }
    let mut ranges: Vec<String> = Vec::new();
    let mut start = sorted_pos[0];
    let mut end = sorted_pos[0];
    for &p in &sorted_pos[1..] {
        if p == end + 1 {
            end = p;
        } else {
            if start != end {
                ranges.push(format!("{start}-{end}"));
            } else {
                ranges.push(start.to_string());
            }
            start = p;
            end = p;
        }
    }
    if start != end {
        ranges.push(format!("{start}-{end}"));
    } else {
        ranges.push(start.to_string());
    }
    format!("\"[{}]\"", ranges.join(","))
}
