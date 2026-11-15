//! `lis-rs`: Zero-Python, multi-threaded native CLI for LIVIA Local Interaction Score analysis.
//!
//! Produces 44-column CSV output matching `python/lis.py` across AlphaFold 3, Protenix-v2,
//! ColabFold, AlphaPulldown (AF-Multimer), Boltz-1/2, Chai-1, OpenFold 3, ESMFold 2,
//! generic PDB/CIF+JSON/NPZ folders, nested `.zip` archives, and `.tar.{gz,xz,zst}` bundles.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::env;
use std::fs::{self, File};
use std::io::{BufRead, BufReader, BufWriter, Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

use flate2::read::GzDecoder;
use livia_core::metrics::{analyze_single_model, format_csv_row, ConfidenceInputs, CSV_HEADER};
use rayon::prelude::*;
use regex::Regex;
use serde_json::Value;

#[derive(Clone)]
enum VirtualSource {
    Disk(PathBuf),
    Bytes(Arc<Vec<u8>>),
}

struct VirtualFs {
    /// Ordered relative POSIX paths (`a/b/c.cif`).
    filenames: Vec<String>,
    entries: HashMap<String, VirtualSource>,
}

impl VirtualFs {
    fn read_bytes(&self, vpath: &str) -> Option<Vec<u8>> {
        let src = self.entries.get(vpath)?;
        let raw = match src {
            VirtualSource::Disk(p) => fs::read(p).ok()?,
            VirtualSource::Bytes(b) => (**b).clone(),
        };
        let lower = vpath.to_ascii_lowercase();
        if lower.ends_with(".gz") && !lower.ends_with(".tar.gz") {
            let mut dec = GzDecoder::new(Cursor::new(raw));
            let mut out = Vec::new();
            dec.read_to_end(&mut out).ok()?;
            Some(out)
        } else if lower.ends_with(".xz") && !lower.ends_with(".tar.xz") {
            let mut out = Vec::new();
            lzma_rs::xz_decompress(&mut Cursor::new(raw), &mut out).ok()?;
            Some(out)
        } else {
            Some(raw)
        }
    }

    fn read_text(&self, vpath: &str) -> Option<String> {
        let bytes = self.read_bytes(vpath)?;
        String::from_utf8(bytes).ok()
    }
}

fn scan_directory_recursive(root: &Path, cur: &Path, filenames: &mut Vec<String>, entries: &mut HashMap<String, VirtualSource>) {
    let Ok(rd) = fs::read_dir(cur) else { return };
    let mut dir_entries: Vec<_> = rd.filter_map(|e| e.ok()).collect();
    dir_entries.sort_by_key(|e| e.file_name());
    for entry in dir_entries {
        let p = entry.path();
        let Ok(ft) = entry.file_type() else { continue };
        if ft.is_dir() {
            scan_directory_recursive(root, &p, filenames, entries);
        } else if ft.is_file() {
            if let Ok(rel) = p.strip_prefix(root) {
                let vpath = rel
                    .components()
                    .map(|c| c.as_os_str().to_string_lossy().into_owned())
                    .collect::<Vec<_>>()
                    .join("/");
                filenames.push(vpath.clone());
                entries.insert(vpath, VirtualSource::Disk(p));
            }
        }
    }
}

fn ingest_zip_bytes(
    prefix: &str,
    data: &[u8],
    filenames: &mut Vec<String>,
    entries: &mut HashMap<String, VirtualSource>,
    depth: usize,
) -> Result<(), String> {
    let mut archive = zip::ZipArchive::new(Cursor::new(data)).map_err(|e| e.to_string())?;
    for i in 0..archive.len() {
        let mut file = archive.by_index(i).map_err(|e| e.to_string())?;
        if file.is_dir() {
            continue;
        }
        let raw_name = file.name().trim_start_matches('/').to_string();
        if raw_name.is_empty() {
            continue;
        }
        let mut buf = Vec::with_capacity(file.size() as usize);
        file.read_to_end(&mut buf).map_err(|e| e.to_string())?;
        let full_vpath = if prefix.is_empty() {
            raw_name.clone()
        } else {
            format!("{prefix}/{raw_name}")
        };
        let lower = raw_name.to_ascii_lowercase();
        // Expand nested `.zip` archives (not `.npz`!).
        if depth < 3 && lower.ends_with(".zip") {
            let sub_prefix = full_vpath.trim_end_matches(".zip").trim_end_matches(".ZIP");
            let _ = ingest_zip_bytes(sub_prefix, &buf, filenames, entries, depth + 1);
        } else {
            filenames.push(full_vpath.clone());
            entries.insert(full_vpath, VirtualSource::Bytes(Arc::new(buf)));
        }
    }
    Ok(())
}

fn ingest_tar_reader<R: Read>(
    reader: R,
    filenames: &mut Vec<String>,
    entries: &mut HashMap<String, VirtualSource>,
) -> Result<(), String> {
    let mut archive = tar::Archive::new(reader);
    let iter = archive.entries().map_err(|e| e.to_string())?;
    for item in iter {
        let mut entry = item.map_err(|e| e.to_string())?;
        if !entry.header().entry_type().is_file() {
            continue;
        }
        let path_str = entry
            .path()
            .map_err(|e| e.to_string())?
            .to_string_lossy()
            .trim_start_matches("./")
            .trim_start_matches('/')
            .to_string();
        if path_str.is_empty() {
            continue;
        }
        let mut buf = Vec::new();
        entry.read_to_end(&mut buf).map_err(|e| e.to_string())?;
        if path_str.to_ascii_lowercase().ends_with(".zip") {
            let sub_prefix = path_str.trim_end_matches(".zip").trim_end_matches(".ZIP");
            let _ = ingest_zip_bytes(sub_prefix, &buf, filenames, entries, 1);
        } else {
            filenames.push(path_str.clone());
            entries.insert(path_str, VirtualSource::Bytes(Arc::new(buf)));
        }
    }
    Ok(())
}

fn scan_files(path: &Path) -> Result<VirtualFs, String> {
    let mut filenames = Vec::new();
    let mut entries = HashMap::new();

    if path.is_dir() {
        scan_directory_recursive(path, path, &mut filenames, &mut entries);
        return Ok(VirtualFs { filenames, entries });
    }

    let name_lower = path
        .file_name()
        .map(|s| s.to_string_lossy().to_ascii_lowercase())
        .unwrap_or_default();

    if name_lower.ends_with(".zip") {
        let bytes = fs::read(path).map_err(|e| e.to_string())?;
        ingest_zip_bytes("", &bytes, &mut filenames, &mut entries, 0)?;
    } else if name_lower.ends_with(".tar.gz") || name_lower.ends_with(".tgz") {
        let f = File::open(path).map_err(|e| e.to_string())?;
        ingest_tar_reader(GzDecoder::new(BufReader::new(f)), &mut filenames, &mut entries)?;
    } else if name_lower.ends_with(".tar.xz") || name_lower.ends_with(".txz") {
        let f = File::open(path).map_err(|e| e.to_string())?;
        let mut tar_bytes = Vec::new();
        lzma_rs::xz_decompress(&mut BufReader::new(f), &mut tar_bytes).map_err(|e| e.to_string())?;
        ingest_tar_reader(Cursor::new(tar_bytes), &mut filenames, &mut entries)?;
    } else if name_lower.ends_with(".tar.zst") || name_lower.ends_with(".tar.zstd") {
        let f = File::open(path).map_err(|e| e.to_string())?;
        let dec = ruzstd::StreamingDecoder::new(BufReader::new(f)).map_err(|e| e.to_string())?;
        ingest_tar_reader(dec, &mut filenames, &mut entries)?;
    } else if name_lower.ends_with(".tar") {
        let f = File::open(path).map_err(|e| e.to_string())?;
        ingest_tar_reader(BufReader::new(f), &mut filenames, &mut entries)?;
    } else {
        return Err(format!("Unsupported input path: {}", path.display()));
    }

    Ok(VirtualFs { filenames, entries })
}

fn posix_basename(p: &str) -> &str {
    p.rsplit('/').next().unwrap_or(p)
}

fn posix_dirname(p: &str) -> &str {
    match p.rfind('/') {
        Some(idx) => &p[..idx],
        None => "",
    }
}

fn get_toplevel_name(p: &str) -> Option<String> {
    let mut parts = p.split('/');
    let first = parts.next()?;
    let _second = parts.next()?;
    if !first.is_empty() && !first.starts_with("__") {
        Some(first.to_string())
    } else {
        None
    }
}

#[derive(Clone, Debug)]
struct ModelSpec {
    name: String,
    rank: String,
    model_label: String,
    struct_path: String,
    pae_path: Option<String>,
    scores_path: Option<String>,
    is_pdb: bool,
    pae_key: Option<String>,
}

fn detect_platform(filenames: &[String], vfs: &VirtualFs) -> String {
    let basenames: Vec<&str> = filenames.iter().map(|f| posix_basename(f)).collect();

    let rx_tam_m = Regex::new(r"^result_sample_\d+_model\.pdb$").unwrap();
    let rx_tam_c = Regex::new(r"^result_sample_\d+_confidences\.json$").unwrap();
    if basenames.iter().any(|b| rx_tam_m.is_match(b)) && basenames.iter().any(|b| rx_tam_c.is_match(b)) {
        for f in filenames {
            if posix_basename(f) == "experiment_config.json" {
                if let Some(txt) = vfs.read_text(f) {
                    if let Ok(v) = serde_json::from_str::<Value>(&txt) {
                        if v.get("inference_ckpt_path")
                            .and_then(|s| s.as_str())
                            .is_some_and(|s| s.contains("of3"))
                        {
                            return "openfold3".to_string();
                        }
                    }
                }
            }
        }
        return "alphafold3".to_string();
    }

    let rx_af3_m = Regex::new(r"_model_\d+\.cif$").unwrap();
    let rx_af3_s = Regex::new(r"_summary_confidences_\d+\.json$").unwrap();
    let rx_af3_f = Regex::new(r"_full_data_\d+\.json$").unwrap();
    if basenames.iter().any(|b| rx_af3_m.is_match(b))
        && basenames.iter().any(|b| rx_af3_s.is_match(b))
        && basenames.iter().any(|b| rx_af3_f.is_match(b))
    {
        return "alphafold3".to_string();
    }

    let rx_pr_m = Regex::new(r"_sample_\d+\.cif$").unwrap();
    let rx_pr_s = Regex::new(r"_summary_confidence_sample_\d+\.json$").unwrap();
    let rx_pr_f = Regex::new(r"_full_data_sample_\d+\.json$").unwrap();
    if basenames.iter().any(|b| rx_pr_m.is_match(b))
        && basenames.iter().any(|b| rx_pr_s.is_match(b))
        && basenames.iter().any(|b| rx_pr_f.is_match(b))
    {
        return "alphafold3".to_string();
    }

    let rx_srv_m = Regex::new(r"seed-\d+_sample-\d+/model\.cif$").unwrap();
    let rx_srv_c = Regex::new(r"seed-\d+_sample-\d+/confidences\.json$").unwrap();
    if filenames.iter().any(|f| rx_srv_m.is_match(f)) && filenames.iter().any(|f| rx_srv_c.is_match(f)) {
        return "alphafold3".to_string();
    }

    let rx_cf_m = Regex::new(r"_unrelaxed_rank_\d+").unwrap();
    let rx_cf_s = Regex::new(r"_scores_rank_\d+").unwrap();
    if basenames
        .iter()
        .any(|b| rx_cf_m.is_match(b) && (b.ends_with(".pdb") || b.ends_with(".cif")))
        && basenames.iter().any(|b| rx_cf_s.is_match(b) && b.ends_with(".json"))
    {
        return "colabfold".to_string();
    }

    let rx_ap_p = Regex::new(r"^pae_model_\d+_.+\.json$").unwrap();
    let rx_ap_c = Regex::new(r"^confidence_model_\d+_.+\.json$").unwrap();
    if basenames.iter().any(|b| rx_ap_p.is_match(b)) && basenames.iter().any(|b| rx_ap_c.is_match(b)) {
        return "alphapulldown".to_string();
    }

    if basenames.iter().any(|b| b.starts_with("confidence_") && b.ends_with(".json"))
        && basenames.iter().any(|b| b.ends_with(".cif") || b.ends_with(".pdb"))
    {
        return "boltz".to_string();
    }

    let rx_chai_m = Regex::new(r"^pred\.(rank_\d+|model_idx_\d+)\.cif$").unwrap();
    let rx_chai_s = Regex::new(r"^scores\.(rank_\d+|model_idx_\d+)\.(json|npz)$").unwrap();
    if basenames.iter().any(|b| rx_chai_m.is_match(b))
        || (basenames.iter().any(|b| rx_chai_s.is_match(b)) && basenames.iter().any(|b| b.starts_with("pred.")))
    {
        return "chai".to_string();
    }

    if basenames.iter().any(|b| b.ends_with("_model.pdb"))
        && basenames.iter().any(|b| b.ends_with("_pae.json"))
    {
        return "esmfold2".to_string();
    }

    let mut by_dir: HashMap<&str, Vec<&str>> = HashMap::new();
    for f in filenames {
        by_dir.entry(posix_dirname(f)).or_default().push(posix_basename(f));
    }
    let struct_exts = [".pdb", ".cif", ".pdb.gz", ".cif.gz", ".pdb.xz", ".cif.xz"];
    let pae_exts = [".npy", ".npz"];
    for names in by_dir.values() {
        let has_s = names.iter().any(|n| {
            let nl = n.to_ascii_lowercase();
            struct_exts.iter().any(|e| nl.ends_with(e))
        });
        let has_p = names.iter().any(|n| {
            let nl = n.to_ascii_lowercase();
            pae_exts.iter().any(|e| nl.ends_with(e))
        });
        if has_s && has_p {
            return "esmfold2_native".to_string();
        }
    }

    if basenames
        .iter()
        .any(|b| b.ends_with("_confidences.json") && !b.ends_with("_summary_confidences.json"))
        && basenames.iter().any(|b| b.ends_with("_summary_confidences.json"))
        && basenames.iter().any(|b| b.ends_with("_model.cif"))
    {
        return "alphafold3".to_string();
    }

    "generic".to_string()
}

fn find_models(raw_filenames: &[String], platform: &str, vfs: &VirtualFs) -> Vec<ModelSpec> {
    let filenames: Vec<String> = raw_filenames
        .iter()
        .filter(|f| !f.contains("__MACOSX") && !posix_basename(f).starts_with("._"))
        .cloned()
        .collect();
    let mut basenames_map: HashMap<String, String> = HashMap::new();
    for f in &filenames {
        basenames_map.insert(posix_basename(f).to_string(), f.clone());
    }

    match platform {
        "colabfold" => find_colabfold(&filenames),
        "alphafold3" | "openfold3" => {
            let v = find_af3(&filenames, &basenames_map);
            if !v.is_empty() {
                v
            } else {
                find_af3_local(&filenames)
            }
        }
        "alphapulldown" => find_alphapulldown(&filenames, vfs),
        "esmfold2" => find_esmfold2(&filenames, &basenames_map),
        "esmfold2_native" => find_esmfold2_native(&filenames),
        "boltz" => find_boltz(&filenames),
        "chai" => find_chai(&filenames, &basenames_map),
        _ => find_generic(&filenames, vfs),
    }
}

fn find_colabfold(filenames: &[String]) -> Vec<ModelSpec> {
    let rx = Regex::new(r"^(.+)_unrelaxed_rank_(\d+)(.*)\.(pdb|cif)$").unwrap();
    let mut out = Vec::new();
    for name in filenames {
        let base = posix_basename(name);
        let Some(caps) = rx.captures(base) else { continue };
        let prefix = &caps[1];
        let rank_num: usize = caps[2].parse().unwrap_or(0);
        let rank = rank_num.to_string();
        let is_pdb = &caps[4] == "pdb";
        let padded = format!("{rank_num:03}");
        let scores_prefix = format!("{prefix}_scores_");

        let mut scores_path = None;
        for fn_cand in filenames {
            let fb = posix_basename(fn_cand);
            if fb.starts_with(&scores_prefix) && fb.contains(&format!("rank_{padded}")) && fb.ends_with(".json") {
                scores_path = Some(fn_cand.clone());
                break;
            }
        }
        if scores_path.is_none() {
            for fn_cand in filenames {
                let fb = posix_basename(fn_cand);
                if fb.starts_with(&scores_prefix) && fb.contains(&format!("rank_{rank}")) && fb.ends_with(".json") {
                    scores_path = Some(fn_cand.clone());
                    break;
                }
            }
        }
        if scores_path.is_none() {
            for fn_cand in filenames {
                let fb = posix_basename(fn_cand);
                if fb.starts_with(&scores_prefix) && fb.ends_with(".json") {
                    scores_path = Some(fn_cand.clone());
                    break;
                }
            }
        }

        out.push(ModelSpec {
            name: prefix.to_string(),
            rank,
            model_label: base.to_string(),
            struct_path: name.clone(),
            pae_path: scores_path.clone(),
            scores_path,
            is_pdb,
            pae_key: None,
        });
    }
    out
}

fn find_af3(filenames: &[String], basenames_map: &HashMap<String, String>) -> Vec<ModelSpec> {
    let mut out = Vec::new();
    let rx_std = Regex::new(r"^(.+)_model_(\d+)\.cif$").unwrap();
    for name in filenames {
        let base = posix_basename(name);
        if let Some(c) = rx_std.captures(base) {
            let prefix = &c[1];
            let idx = &c[2];
            let fd = basenames_map.get(&format!("{prefix}_full_data_{idx}.json"));
            let sm = basenames_map.get(&format!("{prefix}_summary_confidences_{idx}.json"));
            if let (Some(fd), Some(sm)) = (fd, sm) {
                out.push(ModelSpec {
                    name: prefix.to_string(),
                    rank: idx.to_string(),
                    model_label: base.to_string(),
                    struct_path: name.clone(),
                    pae_path: Some(fd.clone()),
                    scores_path: Some(sm.clone()),
                    is_pdb: false,
                    pae_key: None,
                });
            }
        }
    }

    let rx_prot = Regex::new(r"^(.+)_sample_(\d+)\.cif$").unwrap();
    for name in filenames {
        let base = posix_basename(name);
        if let Some(c) = rx_prot.captures(base) {
            let prefix = &c[1];
            let idx = &c[2];
            let fd = basenames_map.get(&format!("{prefix}_full_data_sample_{idx}.json"));
            let sm = basenames_map.get(&format!("{prefix}_summary_confidence_sample_{idx}.json"));
            if let (Some(fd), Some(sm)) = (fd, sm) {
                out.push(ModelSpec {
                    name: prefix.to_string(),
                    rank: idx.to_string(),
                    model_label: base.to_string(),
                    struct_path: name.clone(),
                    pae_path: Some(fd.clone()),
                    scores_path: Some(sm.clone()),
                    is_pdb: false,
                    pae_key: None,
                });
            }
        }
    }

    let rx_tam = Regex::new(r"^result_sample_(\d+)_model\.pdb$").unwrap();
    for name in filenames {
        let base = posix_basename(name);
        if let Some(c) = rx_tam.captures(base) {
            let idx = &c[1];
            let conf = basenames_map.get(&format!("result_sample_{idx}_confidences.json")).cloned();
            let agg = basenames_map
                .get(&format!("result_sample_{idx}_confidences_aggregated.json"))
                .cloned();
            let pred_name = get_toplevel_name(name).unwrap_or_else(|| "prediction".to_string());
            out.push(ModelSpec {
                name: pred_name,
                rank: idx.to_string(),
                model_label: base.to_string(),
                struct_path: name.clone(),
                pae_path: conf.clone(),
                scores_path: agg.or(conf),
                is_pdb: true,
                pae_key: None,
            });
        }
    }

    let rx_srv = Regex::new(r"^seed-(\d+)_sample-(\d+)").unwrap();
    let fset: HashSet<&str> = filenames.iter().map(|s| s.as_str()).collect();
    for name in filenames {
        let base = posix_basename(name);
        if base != "model.cif" {
            continue;
        }
        let dirpath = posix_dirname(name);
        let sample_dir = posix_basename(dirpath);
        let Some(c) = rx_srv.captures(sample_dir) else { continue };
        let pred_dir = posix_dirname(dirpath);
        let mut pred_name = posix_basename(pred_dir).to_string();
        if pred_name.is_empty() || pred_name == "AF3_outputs" {
            pred_name = "af3_prediction".to_string();
        }
        let conf_cand = if dirpath.is_empty() {
            "confidences.json".to_string()
        } else {
            format!("{dirpath}/confidences.json")
        };
        if !fset.contains(conf_cand.as_str()) {
            continue;
        }
        let summ_cand = if dirpath.is_empty() {
            "summary_confidences.json".to_string()
        } else {
            format!("{dirpath}/summary_confidences.json")
        };
        let summ = if fset.contains(summ_cand.as_str()) {
            summ_cand
        } else {
            conf_cand.clone()
        };
        out.push(ModelSpec {
            name: pred_name,
            rank: format!("{}_{}", &c[1], &c[2]),
            model_label: base.to_string(),
            struct_path: name.clone(),
            pae_path: Some(conf_cand),
            scores_path: Some(summ),
            is_pdb: false,
            pae_key: None,
        });
    }

    out
}

fn find_af3_local(filenames: &[String]) -> Vec<ModelSpec> {
    let fset: HashSet<&str> = filenames.iter().map(|s| s.as_str()).collect();
    let rx_ss = Regex::new(r"seed-(\d+)_sample-(\d+)").unwrap();
    let rx_strip = Regex::new(r"_seed-\d+_sample-\d+$").unwrap();
    let mut seed_models = Vec::new();
    let mut top_models = Vec::new();

    for name in filenames {
        let base = posix_basename(name);
        let Some(prefix) = base.strip_suffix("_model.cif") else { continue };
        let dirpath = posix_dirname(name);
        let join = |f: &str| {
            if dirpath.is_empty() {
                f.to_string()
            } else {
                format!("{dirpath}/{f}")
            }
        };
        let conf = join(&format!("{prefix}_confidences.json"));
        if !fset.contains(conf.as_str()) {
            continue;
        }
        let summ_cand = join(&format!("{prefix}_summary_confidences.json"));
        let summ = if fset.contains(summ_cand.as_str()) {
            summ_cand
        } else {
            conf.clone()
        };
        let stripped = rx_strip.replace(prefix, "").to_string();
        let pred = if !stripped.is_empty() {
            stripped
        } else {
            get_toplevel_name(name).unwrap_or_else(|| "af3_prediction".to_string())
        };
        if let Some(c) = rx_ss.captures(prefix).or_else(|| rx_ss.captures(dirpath)) {
            seed_models.push(ModelSpec {
                name: pred,
                rank: format!("{}_{}", &c[1], &c[2]),
                model_label: base.to_string(),
                struct_path: name.clone(),
                pae_path: Some(conf),
                scores_path: Some(summ),
                is_pdb: false,
                pae_key: None,
            });
        } else {
            top_models.push(ModelSpec {
                name: pred,
                rank: "0".to_string(),
                model_label: base.to_string(),
                struct_path: name.clone(),
                pae_path: Some(conf),
                scores_path: Some(summ),
                is_pdb: false,
                pae_key: None,
            });
        }
    }
    if !seed_models.is_empty() {
        seed_models
    } else {
        top_models
    }
}

fn find_alphapulldown(filenames: &[String], vfs: &VirtualFs) -> Vec<ModelSpec> {
    let mut dir_order: Vec<String> = Vec::new();
    let mut by_dir: HashMap<String, Vec<String>> = HashMap::new();
    for f in filenames {
        let d = posix_dirname(f).to_string();
        by_dir
            .entry(d.clone())
            .or_insert_with(|| {
                dir_order.push(d);
                Vec::new()
            })
            .push(f.clone());
    }

    let rx_pae = Regex::new(r"^pae_(model_\d+.*?)\.json$").unwrap();
    let rx_conf = Regex::new(r"^confidence_(model_\d+.*?)\.json$").unwrap();
    let rx_unrel = Regex::new(r"^unrelaxed_(model_\d+.*?)\.(pdb|cif)$").unwrap();
    let rx_rel = Regex::new(r"^relaxed_(model_\d+.*?)\.(pdb|cif)$").unwrap();
    let rx_ranked = Regex::new(r"^ranked_(\d+)\.(pdb|cif)$").unwrap();

    let mut out = Vec::new();
    for d in dir_order {
        let files = &by_dir[&d];
        let mut local: HashMap<String, String> = HashMap::new();
        for f in files {
            local.insert(posix_basename(f).to_string(), f.clone());
        }
        let mut pae_by_id: HashMap<String, String> = HashMap::new();
        let mut conf_by_id: HashMap<String, String> = HashMap::new();
        let mut struct_by_id: HashMap<String, String> = HashMap::new();
        let mut ranked_by_rank: HashMap<usize, String> = HashMap::new();

        for (b, f) in &local {
            if let Some(c) = rx_pae.captures(b) {
                pae_by_id.insert(c[1].to_string(), f.clone());
            }
            if let Some(c) = rx_conf.captures(b) {
                conf_by_id.insert(c[1].to_string(), f.clone());
            }
            if let Some(c) = rx_unrel.captures(b) {
                struct_by_id.entry(c[1].to_string()).or_insert_with(|| f.clone());
            }
            if let Some(c) = rx_ranked.captures(b) {
                if let Ok(rk) = c[1].parse::<usize>() {
                    ranked_by_rank.insert(rk, f.clone());
                }
            }
        }
        for (b, f) in &local {
            if let Some(c) = rx_rel.captures(b) {
                struct_by_id.insert(c[1].to_string(), f.clone());
            }
        }
        if pae_by_id.is_empty() || (struct_by_id.is_empty() && ranked_by_rank.is_empty()) {
            continue;
        }

        let mut rank_order: Vec<String> = Vec::new();
        if let Some(rd_path) = local.get("ranking_debug.json") {
            if let Some(txt) = vfs.read_text(rd_path) {
                if let Ok(v) = serde_json::from_str::<Value>(&txt) {
                    if let Some(arr) = v.get("order").and_then(|x| x.as_array()) {
                        rank_order = arr
                            .iter()
                            .filter_map(|x| x.as_str().map(|s| s.to_string()))
                            .collect();
                    }
                }
            }
        }

        let pred_name = {
            let bn = posix_basename(&d);
            if bn.is_empty() {
                "prediction".to_string()
            } else {
                bn.to_string()
            }
        };

        if !rank_order.is_empty() && !struct_by_id.is_empty() {
            for (rank_idx, model_id) in rank_order.iter().enumerate() {
                let (Some(sp), Some(pp)) = (struct_by_id.get(model_id), pae_by_id.get(model_id)) else {
                    continue;
                };
                let cp = conf_by_id.get(model_id).unwrap_or(pp);
                out.push(ModelSpec {
                    name: pred_name.clone(),
                    rank: rank_idx.to_string(),
                    model_label: posix_basename(sp).to_string(),
                    struct_path: sp.clone(),
                    pae_path: Some(pp.clone()),
                    scores_path: Some(cp.clone()),
                    is_pdb: !sp.to_ascii_lowercase().ends_with(".cif"),
                    pae_key: None,
                });
            }
        } else if !struct_by_id.is_empty() {
            let mut ids: Vec<String> = struct_by_id
                .keys()
                .filter(|k| pae_by_id.contains_key(*k))
                .cloned()
                .collect();
            ids.sort();
            for (rank_idx, model_id) in ids.iter().enumerate() {
                let sp = &struct_by_id[model_id];
                let pp = &pae_by_id[model_id];
                let cp = conf_by_id.get(model_id).unwrap_or(pp);
                out.push(ModelSpec {
                    name: pred_name.clone(),
                    rank: rank_idx.to_string(),
                    model_label: posix_basename(sp).to_string(),
                    struct_path: sp.clone(),
                    pae_path: Some(pp.clone()),
                    scores_path: Some(cp.clone()),
                    is_pdb: !sp.to_ascii_lowercase().ends_with(".cif"),
                    pae_key: None,
                });
            }
        } else if !rank_order.is_empty() {
            for (rank_idx, model_id) in rank_order.iter().enumerate() {
                let (Some(sp), Some(pp)) = (ranked_by_rank.get(&rank_idx), pae_by_id.get(model_id)) else {
                    continue;
                };
                let cp = conf_by_id.get(model_id).unwrap_or(pp);
                out.push(ModelSpec {
                    name: pred_name.clone(),
                    rank: rank_idx.to_string(),
                    model_label: model_id.clone(),
                    struct_path: sp.clone(),
                    pae_path: Some(pp.clone()),
                    scores_path: Some(cp.clone()),
                    is_pdb: !sp.to_ascii_lowercase().ends_with(".cif"),
                    pae_key: None,
                });
            }
        }
    }
    out
}

fn find_boltz(filenames: &[String]) -> Vec<ModelSpec> {
    let mut dir_groups: BTreeMap<String, Vec<String>> = BTreeMap::new();
    for f in filenames {
        dir_groups
            .entry(posix_dirname(f).to_string())
            .or_default()
            .push(f.clone());
    }
    if dir_groups.len() <= 1 {
        dir_groups.clear();
        dir_groups.insert(String::new(), filenames.to_vec());
    }

    let rx_mod = Regex::new(r"model_(\d+)").unwrap();
    let mut out = Vec::new();
    for (dirpath, dir_files) in dir_groups {
        let mut struct_files: Vec<String> = dir_files
            .iter()
            .filter(|f| f.ends_with(".cif") || f.ends_with(".pdb"))
            .cloned()
            .collect();
        struct_files.sort_by(|a, b| posix_basename(a).cmp(posix_basename(b)));
        if struct_files.is_empty() {
            continue;
        }
        let conf_files: Vec<String> = dir_files
            .iter()
            .filter(|f| posix_basename(f).starts_with("confidence") && f.ends_with(".json"))
            .cloned()
            .collect();
        let pae_files: Vec<String> = dir_files
            .iter()
            .filter(|f| posix_basename(f).starts_with("pae") && f.ends_with(".npz"))
            .cloned()
            .collect();

        let mut pred_name = if dirpath.is_empty() {
            String::new()
        } else {
            posix_basename(&dirpath).to_string()
        };
        if pred_name.is_empty()
            || matches!(
                pred_name.as_str(),
                "Boltz1_outputs" | "boltz_outputs" | "predictions" | "result" | "results"
            )
        {
            pred_name = get_toplevel_name(&struct_files[0]).unwrap_or_else(|| "boltz".to_string());
        }

        for (i, sf) in struct_files.iter().enumerate() {
            let sb = posix_basename(sf);
            let model_idx = rx_mod
                .captures(sb)
                .map(|c| c[1].to_string())
                .unwrap_or_else(|| i.to_string());
            let needle = format!("model_{model_idx}");
            let conf_path = conf_files
                .iter()
                .find(|cf| posix_basename(cf).contains(&needle))
                .or_else(|| conf_files.get(i))
                .cloned();
            let pae_path = pae_files
                .iter()
                .find(|pf| posix_basename(pf).contains(&needle))
                .or_else(|| pae_files.get(i))
                .cloned();

            out.push(ModelSpec {
                name: pred_name.clone(),
                rank: model_idx,
                model_label: sb.to_string(),
                struct_path: sf.clone(),
                pae_path: pae_path.or_else(|| conf_path.clone()),
                scores_path: conf_path,
                is_pdb: sf.ends_with(".pdb"),
                pae_key: None,
            });
        }
    }
    out
}

fn find_chai(filenames: &[String], basenames_map: &HashMap<String, String>) -> Vec<ModelSpec> {
    let rx = Regex::new(r"^pred\.(?:rank_(\d+)|model_idx_(\d+))\.cif$").unwrap();
    let mut out = Vec::new();
    for name in filenames {
        let base = posix_basename(name);
        let Some(c) = rx.captures(base) else { continue };
        let idx = c
            .get(1)
            .or_else(|| c.get(2))
            .map(|m| m.as_str())
            .unwrap_or("0");
        let mut score_path = basenames_map
            .get(&format!("scores.rank_{idx}.json"))
            .or_else(|| basenames_map.get(&format!("scores.model_idx_{idx}.json")))
            .cloned();
        if score_path.is_none() {
            score_path = filenames
                .iter()
                .find(|f| posix_basename(f).contains("scores") && f.ends_with(".json"))
                .cloned();
        }
        let pae_path = basenames_map
            .get(&format!("pae.rank_{idx}.npy"))
            .or_else(|| basenames_map.get(&format!("pae.rank_{idx}.npz")))
            .or_else(|| basenames_map.get(&format!("pae.model_idx_{idx}.npy")))
            .or_else(|| basenames_map.get(&format!("pae.model_idx_{idx}.npz")))
            .cloned();
        let pred_name = get_toplevel_name(name).unwrap_or_else(|| "prediction".to_string());
        out.push(ModelSpec {
            name: pred_name,
            rank: idx.to_string(),
            model_label: base.to_string(),
            struct_path: name.clone(),
            pae_path: pae_path.or_else(|| score_path.clone()),
            scores_path: score_path,
            is_pdb: false,
            pae_key: None,
        });
    }
    out
}

fn find_esmfold2(filenames: &[String], basenames_map: &HashMap<String, String>) -> Vec<ModelSpec> {
    let mut out = Vec::new();
    for name in filenames {
        let base = posix_basename(name);
        let Some(prefix) = base.strip_suffix("_model.pdb") else { continue };
        let Some(pae_path) = basenames_map.get(&format!("{prefix}_pae.json")) else { continue };
        out.push(ModelSpec {
            name: prefix.to_string(),
            rank: "0".to_string(),
            model_label: base.to_string(),
            struct_path: name.clone(),
            pae_path: Some(pae_path.clone()),
            scores_path: Some(pae_path.clone()),
            is_pdb: true,
            pae_key: None,
        });
    }
    out
}

fn find_esmfold2_native(filenames: &[String]) -> Vec<ModelSpec> {
    let mut dir_order: Vec<String> = Vec::new();
    let mut by_dir: HashMap<String, Vec<String>> = HashMap::new();
    for f in filenames {
        let d = posix_dirname(f).to_string();
        by_dir
            .entry(d.clone())
            .or_insert_with(|| {
                dir_order.push(d);
                Vec::new()
            })
            .push(f.clone());
    }

    let rx_s = Regex::new(r"(?i)\.(pdb|cif)(\.gz|\.xz)?$").unwrap();
    let rx_p = Regex::new(r"(?i)\.(npy|npz)$").unwrap();
    let rx_strip = Regex::new(r"(?i)\.(pdb|cif|npy|npz|json)(\.gz|\.xz)?$").unwrap();

    let stem_of = |p: &str| rx_strip.replace(posix_basename(p), "").to_string();
    let lcp = |strs: &[String]| -> String {
        let (Some(s1), Some(s2)) = (strs.iter().min(), strs.iter().max()) else {
            return String::new();
        };
        s1.chars()
            .zip(s2.chars())
            .take_while(|(a, b)| a == b)
            .map(|(a, _)| a)
            .collect()
    };
    let ids_by_lcp = |files: &[String]| -> HashMap<String, String> {
        if files.len() == 1 {
            let mut m = HashMap::new();
            m.insert(String::new(), files[0].clone());
            return m;
        }
        let stems: Vec<String> = files.iter().map(|f| stem_of(f)).collect();
        let pref = lcp(&stems);
        let mut out = HashMap::new();
        for (f, st) in files.iter().zip(stems) {
            let sid = st[pref.len()..].trim_start_matches(['_', '-', '.']).to_string();
            out.insert(sid, f.clone());
        }
        out
    };

    let mut out = Vec::new();
    for d in dir_order {
        let files = &by_dir[&d];
        let mut struct_files: Vec<String> = files
            .iter()
            .filter(|f| rx_s.is_match(posix_basename(f)))
            .cloned()
            .collect();
        let mut pae_files: Vec<String> = files
            .iter()
            .filter(|f| rx_p.is_match(posix_basename(f)))
            .cloned()
            .collect();
        let json_files: Vec<String> = files
            .iter()
            .filter(|f| posix_basename(f).to_ascii_lowercase().ends_with(".json"))
            .cloned()
            .collect();
        if struct_files.is_empty() || pae_files.is_empty() {
            continue;
        }
        struct_files.sort_by_key(|p| {
            let b = posix_basename(p).to_ascii_lowercase();
            usize::from(b.ends_with(".gz") || b.ends_with(".xz"))
        });
        pae_files.sort_by_key(|p| usize::from(!p.to_ascii_lowercase().ends_with(".npy")));

        let pred_name = {
            let bn = posix_basename(&d);
            if bn.is_empty() {
                "prediction".to_string()
            } else {
                bn.to_string()
            }
        };

        let emit_single = |sp: &String, pp: &String| -> ModelSpec {
            let jp = json_files.first().unwrap_or(pp);
            let is_pdb = !posix_basename(sp).to_ascii_lowercase().contains(".cif");
            ModelSpec {
                name: pred_name.clone(),
                rank: "0".to_string(),
                model_label: posix_basename(sp).to_string(),
                struct_path: sp.clone(),
                pae_path: Some(pp.clone()),
                scores_path: Some(jp.clone()),
                is_pdb,
                pae_key: None,
            }
        };

        if struct_files.len() == 1 || pae_files.len() == 1 {
            out.push(emit_single(&struct_files[0], &pae_files[0]));
            continue;
        }

        let s_ids = ids_by_lcp(&struct_files);
        let p_ids = ids_by_lcp(&pae_files);
        let j_ids = if json_files.is_empty() {
            HashMap::new()
        } else {
            ids_by_lcp(&json_files)
        };
        let mut common: Vec<String> = s_ids
            .keys()
            .filter(|k| p_ids.contains_key(*k))
            .cloned()
            .collect();
        if common.is_empty() {
            out.push(emit_single(&struct_files[0], &pae_files[0]));
            continue;
        }
        common.sort();
        for sid in common {
            let sp = &s_ids[&sid];
            let pp = &p_ids[&sid];
            let jp = j_ids
                .get(&sid)
                .or_else(|| json_files.first())
                .unwrap_or(pp);
            let is_pdb = !posix_basename(sp).to_ascii_lowercase().contains(".cif");
            out.push(ModelSpec {
                name: pred_name.clone(),
                rank: if sid.is_empty() { "0".to_string() } else { sid },
                model_label: posix_basename(sp).to_string(),
                struct_path: sp.clone(),
                pae_path: Some(pp.clone()),
                scores_path: Some(jp.clone()),
                is_pdb,
                pae_key: None,
            });
        }
    }
    out
}

fn find_generic(filenames: &[String], vfs: &VirtualFs) -> Vec<ModelSpec> {
    let mut struct_files: Vec<String> = filenames
        .iter()
        .filter(|f| f.ends_with(".cif") || f.ends_with(".pdb"))
        .cloned()
        .collect();
    struct_files.sort_by(|a, b| posix_basename(a).cmp(posix_basename(b)));
    let all_json: Vec<String> = filenames
        .iter()
        .filter(|f| f.ends_with(".json"))
        .cloned()
        .collect();
    let npz_files: Vec<String> = filenames
        .iter()
        .filter(|f| f.ends_with(".npz"))
        .cloned()
        .collect();

    let is_non_pae = |f: &str| -> bool {
        let b = posix_basename(f).to_ascii_lowercase();
        [
            "_data.json",
            "ranking_scores",
            "ranking_debug",
            "config",
            "_summary_confidences.json",
            "summary_confidence",
            "metrics.json",
        ]
        .iter()
        .any(|k| b.contains(k))
    };
    let filtered_json: Vec<String> = all_json.iter().filter(|f| !is_non_pae(f)).cloned().collect();
    let json_files = if filtered_json.is_empty() {
        all_json
    } else {
        filtered_json
    };
    let rx_digits = Regex::new(r"(\d+)").unwrap();

    let mut out = Vec::new();
    for (i, sf) in struct_files.iter().enumerate() {
        let sb = posix_basename(sf);
        let is_pdb = sf.ends_with(".pdb");
        let name_prefix = sb.trim_end_matches(".cif").trim_end_matches(".pdb");
        let idx = rx_digits.captures(name_prefix).map(|c| c[1].to_string());

        let mut pae_source = json_files
            .iter()
            .find(|jf| posix_basename(jf).contains(name_prefix))
            .cloned();
        if pae_source.is_none() {
            if let Some(ref idx_s) = idx {
                pae_source = json_files
                    .iter()
                    .find(|jf| {
                        let fb = posix_basename(jf);
                        fb.contains(idx_s.as_str()) && !fb.contains("aggregated")
                    })
                    .cloned();
            }
        }
        if pae_source.is_none() {
            pae_source = json_files.get(i).or_else(|| json_files.first()).cloned();
        }
        if pae_source.is_none() && !npz_files.is_empty() {
            pae_source = npz_files.get(i).or_else(|| npz_files.first()).cloned();
        }
        if let Some(ref chosen) = pae_source {
            if extract_pae(chosen, None, vfs).is_none() {
                for cf in json_files.iter().chain(npz_files.iter()) {
                    if cf != chosen && extract_pae(cf, None, vfs).is_some() {
                        pae_source = Some(cf.clone());
                        break;
                    }
                }
            }
        }
        out.push(ModelSpec {
            name: "generic".to_string(),
            rank: i.to_string(),
            model_label: sb.to_string(),
            struct_path: sf.clone(),
            pae_path: pae_source.clone(),
            scores_path: pae_source,
            is_pdb,
            pae_key: None,
        });
    }
    out
}

fn f16_to_f32(bits: u16) -> f32 {
    let sign = ((bits >> 15) & 1) as u32;
    let exp = ((bits >> 10) & 0x1f) as u32;
    let mant = (bits & 0x3ff) as u32;
    let f_bits = if exp == 0 {
        if mant == 0 {
            sign << 31
        } else {
            let mut m = mant;
            let mut e = 127u32 - 14;
            while (m & 0x400) == 0 {
                m <<= 1;
                e -= 1;
            }
            m &= 0x3ff;
            (sign << 31) | (e << 23) | (m << 13)
        }
    } else if exp == 31 {
        (sign << 31) | 0x7f80_0000 | (mant << 13)
    } else {
        (sign << 31) | ((exp + 127 - 15) << 23) | (mant << 13)
    };
    f32::from_bits(f_bits)
}

fn parse_npy_matrix(bytes: &[u8]) -> Option<(Vec<f32>, usize)> {
    if bytes.len() < 10 || &bytes[0..6] != b"\x93NUMPY" {
        return None;
    }
    let major = bytes[6];
    let (hdr_len, data_off) = if major == 1 {
        let len = u16::from_le_bytes([bytes[8], bytes[9]]) as usize;
        (len, 10 + len)
    } else {
        if bytes.len() < 12 {
            return None;
        }
        let len = u32::from_le_bytes([bytes[8], bytes[9], bytes[10], bytes[11]]) as usize;
        (len, 12 + len)
    };
    if bytes.len() < data_off {
        return None;
    }
    let header = std::str::from_utf8(&bytes[data_off - hdr_len..data_off]).ok()?;

    let rx_descr = Regex::new(r"'descr'\s*:\s*'([^']+)'").ok()?;
    let descr = rx_descr.captures(header)?.get(1)?.as_str();
    let rx_fort = Regex::new(r"'fortran_order'\s*:\s*(True|False)").ok()?;
    let fortran = rx_fort
        .captures(header)
        .and_then(|c| c.get(1))
        .is_some_and(|m| m.as_str() == "True");
    let rx_shape = Regex::new(r"'shape'\s*:\s*\(([^)]*)\)").ok()?;
    let shape_str = rx_shape.captures(header)?.get(1)?.as_str();
    let dims: Vec<usize> = shape_str
        .split(',')
        .map(|s| s.trim())
        .filter(|s| !s.is_empty())
        .filter_map(|s| s.parse::<usize>().ok())
        .collect();

    let (rows, cols) = match dims.as_slice() {
        [r, c] => (*r, *c),
        [_, r, c] => (*r, *c),
        _ => return None,
    };
    if rows == 0 || rows != cols {
        return None;
    }
    let count = rows * cols;
    let payload = &bytes[data_off..];
    let big_endian = descr.starts_with('>');
    let kind_size = descr.trim_start_matches(['<', '>', '|', '=']);

    let mut flat = Vec::with_capacity(count);
    match kind_size {
        "f2" => {
            if payload.len() < count * 2 {
                return None;
            }
            for i in 0..count {
                let b = [payload[i * 2], payload[i * 2 + 1]];
                let bits = if big_endian { u16::from_be_bytes(b) } else { u16::from_le_bytes(b) };
                let v = f16_to_f32(bits);
                flat.push(if v.is_nan() { 31.0 } else { v });
            }
        }
        "f4" => {
            if payload.len() < count * 4 {
                return None;
            }
            for i in 0..count {
                let b = [
                    payload[i * 4],
                    payload[i * 4 + 1],
                    payload[i * 4 + 2],
                    payload[i * 4 + 3],
                ];
                let v = if big_endian { f32::from_be_bytes(b) } else { f32::from_le_bytes(b) };
                flat.push(if v.is_nan() { 31.0 } else { v });
            }
        }
        "f8" => {
            if payload.len() < count * 8 {
                return None;
            }
            for i in 0..count {
                let mut b = [0u8; 8];
                b.copy_from_slice(&payload[i * 8..i * 8 + 8]);
                let v = (if big_endian { f64::from_be_bytes(b) } else { f64::from_le_bytes(b) }) as f32;
                flat.push(if v.is_nan() { 31.0 } else { v });
            }
        }
        _ => return None,
    }

    if fortran {
        let mut tr = vec![0.0_f32; count];
        for r in 0..rows {
            for c in 0..cols {
                tr[r * cols + c] = flat[c * rows + r];
            }
        }
        flat = tr;
    }

    Some((flat, rows))
}

fn parse_npz_matrix(bytes: &[u8]) -> Option<(Vec<f32>, usize)> {
    let mut za = zip::ZipArchive::new(Cursor::new(bytes)).ok()?;
    let mut names: Vec<String> = (0..za.len())
        .filter_map(|i| za.by_index(i).ok().map(|f| f.name().to_string()))
        .collect();
    // Prefer `pae.npy` first just like `lis.py`.
    names.sort_by_key(|n| usize::from(n != "pae.npy" && n != "pae"));
    for name in names {
        if let Ok(mut f) = za.by_name(&name) {
            let mut buf = Vec::with_capacity(f.size() as usize);
            if f.read_to_end(&mut buf).is_ok() {
                if let Some(res) = parse_npy_matrix(&buf) {
                    return Some(res);
                }
            }
        }
    }
    None
}

fn json_val_to_square_pae(v: &Value) -> Option<(Vec<f32>, usize)> {
    let arr = v.as_array()?;
    if arr.is_empty() {
        return None;
    }
    if let Some(first_row) = arr[0].as_array() {
        let n = arr.len();
        if first_row.len() != n {
            return None;
        }
        let mut flat = Vec::with_capacity(n * n);
        for row in arr {
            let r = row.as_array()?;
            if r.len() != n {
                return None;
            }
            for cell in r {
                let x = cell.as_f64().unwrap_or(31.0) as f32;
                flat.push(if x.is_nan() { 31.0 } else { x });
            }
        }
        Some((flat, n))
    } else {
        let len = arr.len();
        let n = (len as f64).sqrt().round() as usize;
        if n * n != len || n == 0 {
            return None;
        }
        let mut flat = Vec::with_capacity(len);
        for cell in arr {
            let x = cell.as_f64().unwrap_or(31.0) as f32;
            flat.push(if x.is_nan() { 31.0 } else { x });
        }
        Some((flat, n))
    }
}

fn extract_pae(pae_source: &str, pae_key: Option<&str>, vfs: &VirtualFs) -> Option<(Vec<f32>, usize)> {
    let base = posix_basename(pae_source).to_ascii_lowercase();
    if base.ends_with(".npy") {
        let bytes = vfs.read_bytes(pae_source)?;
        return parse_npy_matrix(&bytes);
    }
    if base.ends_with(".npz") {
        let bytes = vfs.read_bytes(pae_source)?;
        return parse_npz_matrix(&bytes);
    }
    let txt = vfs.read_text(pae_source)?;
    let data: Value = serde_json::from_str(&txt).ok()?;

    if let Some(k) = pae_key {
        if let Some(v) = data.get(k) {
            if let Some(res) = json_val_to_square_pae(v) {
                return Some(res);
            }
        }
    }

    if let Some(arr) = data.as_array() {
        if let Some(first) = arr.first() {
            if let Some(v) = first.get("predicted_aligned_error") {
                if let Some(res) = json_val_to_square_pae(v) {
                    return Some(res);
                }
            }
        }
    }

    for key in [
        "pae",
        "predicted_aligned_error",
        "pae_matrix",
        "predicted_aligned_error_matrix",
        "pde",
        "token_pair_pae",
    ] {
        if let Some(v) = data.get(key) {
            if let Some(res) = json_val_to_square_pae(v) {
                return Some(res);
            }
        }
    }
    None
}

fn unwrap_scalar(v: &Value) -> Option<f64> {
    if let Some(x) = v.as_f64() {
        return Some(x);
    }
    if let Some(arr) = v.as_array() {
        if arr.len() == 1 {
            return arr[0].as_f64();
        }
    }
    None
}

fn parse_matrix_2d(v: &Value) -> Option<Vec<Vec<Option<f64>>>> {
    let outer = v.as_array()?;
    if outer.is_empty() {
        return None;
    }
    // Chai-1 sometimes nests `[[[...]]]`.
    let rows_val = if outer[0]
        .as_array()
        .and_then(|r| r.first())
        .is_some_and(|c| c.is_array())
    {
        outer[0].as_array()?
    } else {
        outer
    };
    let mut mat = Vec::with_capacity(rows_val.len());
    for row in rows_val {
        let r = row.as_array()?;
        mat.push(r.iter().map(|c| c.as_f64()).collect());
    }
    Some(mat)
}

fn extract_confidence_scores(scores_path: &str, vfs: &VirtualFs) -> ConfidenceInputs {
    let mut out = ConfidenceInputs::default();
    let Some(txt) = vfs.read_text(scores_path) else {
        return out;
    };
    let Ok(data) = serde_json::from_str::<Value>(&txt) else {
        return out;
    };
    let Some(obj) = data.as_object() else {
        return out;
    };

    if let Some(v) = obj.get("ptm").and_then(unwrap_scalar) {
        out.ptm = Some(v);
    }
    if let Some(v) = obj.get("iptm").and_then(unwrap_scalar) {
        out.iptm = Some(v);
    }
    if out.iptm.is_none() {
        if let Some(v) = obj.get("interface_ptm").and_then(unwrap_scalar) {
            out.iptm = Some(v);
        }
    }
    if let Some(v) = obj.get("ptm_score").and_then(unwrap_scalar) {
        out.ptm = Some(v);
    }
    if let Some(v) = obj.get("iptm_score").and_then(unwrap_scalar) {
        out.iptm = Some(v);
    }

    if let Some(v) = obj.get("chain_pair_iptm").and_then(parse_matrix_2d) {
        out.chain_pair_iptm = Some(v);
    }
    if out.chain_pair_iptm.is_none() {
        if let Some(v) = obj.get("per_chain_pair_iptm").and_then(parse_matrix_2d) {
            out.chain_pair_iptm = Some(v);
        }
    }
    if let Some(raw) = obj.get("pair_chains_iptm") {
        if let Some(mat) = parse_matrix_2d(raw) {
            out.chain_pair_iptm = Some(mat);
        } else if let Some(dict) = raw.as_object() {
            let mut row_keys: Vec<&String> = dict.keys().collect();
            row_keys.sort_by_key(|k| k.parse::<i64>().unwrap_or(0));
            let mut mat = Vec::with_capacity(row_keys.len());
            for rk in row_keys {
                if let Some(robj) = dict[rk].as_object() {
                    let mut col_keys: Vec<&String> = robj.keys().collect();
                    col_keys.sort_by_key(|k| k.parse::<i64>().unwrap_or(0));
                    mat.push(col_keys.into_iter().map(|ck| robj[ck].as_f64()).collect());
                }
            }
            if !mat.is_empty() {
                out.chain_pair_iptm = Some(mat);
            }
        }
    }

    // OpenFold3 / Tamarind `iptm_by_asym_id_pair`: numeric `"(1, 2)"` OR chain-letter `"(A, B)"`.
    if out.chain_pair_iptm.is_none() && out.named_pair_iptm.is_none() {
        if let Some(raw) = obj.get("iptm_by_asym_id_pair").and_then(|v| v.as_object()) {
            let rx_num = Regex::new(r"^\((\d+),\s*(\d+)\)$").unwrap();
            let rx_sym = Regex::new(r"^\(([^,]+),\s*([^)]+)\)$").unwrap();
            let mut num_ids: HashSet<usize> = HashSet::new();
            for k in raw.keys() {
                if let Some(c) = rx_num.captures(k) {
                    if let (Ok(a), Ok(b)) = (c[1].parse::<usize>(), c[2].parse::<usize>()) {
                        num_ids.insert(a);
                        num_ids.insert(b);
                    }
                }
            }
            if !num_ids.is_empty() {
                let mut ids: Vec<usize> = num_ids.into_iter().collect();
                ids.sort_unstable();
                let mut mat = Vec::with_capacity(ids.len());
                for &i in &ids {
                    let mut row = Vec::with_capacity(ids.len());
                    for &j in &ids {
                        let v = raw
                            .get(&format!("({i}, {j})"))
                            .or_else(|| raw.get(&format!("({j}, {i})")))
                            .and_then(|x| x.as_f64())
                            .or(Some(0.0));
                        row.push(v);
                    }
                    mat.push(row);
                }
                out.chain_pair_iptm = Some(mat);
            } else {
                let mut named = HashMap::new();
                for (k, val) in raw {
                    if let (Some(c), Some(fv)) = (rx_sym.captures(k), val.as_f64()) {
                        let a = c[1].trim().to_string();
                        let b = c[2].trim().to_string();
                        named.insert((a.clone(), b.clone()), fv);
                        named.insert((b, a), fv);
                    }
                }
                if !named.is_empty() {
                    out.named_pair_iptm = Some(named);
                }
            }
        }
    }

    let mean_arr = |key: &str| -> Option<f64> {
        let arr = obj.get(key)?.as_array()?;
        if arr.is_empty() {
            return None;
        }
        let mut s = 0.0;
        let mut n = 0usize;
        for v in arr {
            if let Some(x) = v.as_f64() {
                s += x;
                n += 1;
            }
        }
        (n > 0).then(|| s / (n as f64))
    };

    if let Some(m) = mean_arr("atom_plddts") {
        out.global_plddt = Some(m);
    }
    if out.global_plddt.is_none() {
        out.global_plddt = mean_arr("plddt");
    }
    if out.global_plddt.is_none() {
        out.global_plddt = mean_arr("confidenceScore");
    }
    if out.global_plddt.is_none() {
        out.global_plddt = obj.get("avg_plddt").and_then(|v| v.as_f64());
    }
    if out.global_plddt.is_none() {
        out.global_plddt = obj.get("complex_plddt").and_then(|v| v.as_f64());
    }

    out
}

fn extract_alphapulldown_scalars(
    struct_path: &str,
    model_label: &str,
    vfs: &VirtualFs,
) -> (Option<f64>, Option<f64>) {
    let dir = posix_dirname(struct_path);
    let rd_path = if dir.is_empty() {
        "ranking_debug.json".to_string()
    } else {
        format!("{dir}/ranking_debug.json")
    };
    let Some(txt) = vfs.read_text(&rd_path) else {
        return (None, None);
    };
    let Ok(rd) = serde_json::from_str::<Value>(&txt) else {
        return (None, None);
    };
    let rx = Regex::new(r"model_\d+(?:_[A-Za-z0-9]+)*").unwrap();
    let Some(m) = rx.find(model_label) else {
        return (None, None);
    };
    let model_id = m.as_str();
    let iptm = rd
        .get("iptm")
        .and_then(|d| d.get(model_id))
        .and_then(|v| v.as_f64());
    let mut ptm = None;
    if let Some(ip) = iptm {
        if let Some(comp) = rd
            .get("iptm+ptm")
            .and_then(|d| d.get(model_id))
            .and_then(|v| v.as_f64())
        {
            let p = 5.0 * comp - 4.0 * ip;
            if (-0.005..=1.005).contains(&p) {
                ptm = Some(p.clamp(0.0, 1.0));
            }
        }
    }
    (iptm, ptm)
}

fn extract_token_chain_ids(pae_path: &str, vfs: &VirtualFs) -> Option<Vec<String>> {
    if !pae_path.to_ascii_lowercase().ends_with(".json") {
        return None;
    }
    let txt = vfs.read_text(pae_path)?;
    let v: Value = serde_json::from_str(&txt).ok()?;
    let arr = v.get("token_chain_ids")?.as_array()?;
    let mut out = Vec::with_capacity(arr.len());
    for item in arr {
        if let Some(s) = item.as_str() {
            out.push(s.to_string());
        } else {
            return None;
        }
    }
    (!out.is_empty()).then_some(out)
}

fn process_one_model(
    spec: &ModelSpec,
    detected: &str,
    pae_cutoff: f64,
    cb_cutoff: f64,
    vfs: &VirtualFs,
) -> Result<Vec<String>, String> {
    let struct_text = vfs
        .read_text(&spec.struct_path)
        .ok_or_else(|| format!("structure file unreadable: {}", spec.struct_path))?;
    let pae_path = spec
        .pae_path
        .as_deref()
        .ok_or_else(|| "missing PAE path".to_string())?;
    let (pae, pae_dim) = extract_pae(pae_path, spec.pae_key.as_deref(), vfs)
        .ok_or_else(|| format!("PAE not found or unreadable: {pae_path}"))?;

    let mut scores = spec
        .scores_path
        .as_deref()
        .map(|sp| extract_confidence_scores(sp, vfs))
        .unwrap_or_default();

    if spec.scores_path.as_deref() != Some(pae_path) {
        let fallback = extract_confidence_scores(pae_path, vfs);
        if scores.ptm.is_none() {
            scores.ptm = fallback.ptm;
        }
        if scores.iptm.is_none() {
            scores.iptm = fallback.iptm;
        }
        if scores.global_plddt.is_none() {
            scores.global_plddt = fallback.global_plddt;
        }
        if scores.chain_pair_iptm.is_none() {
            scores.chain_pair_iptm = fallback.chain_pair_iptm;
        }
        if scores.named_pair_iptm.is_none() {
            scores.named_pair_iptm = fallback.named_pair_iptm;
        }
    }

    if detected == "alphapulldown" {
        let (ap_iptm, ap_ptm) =
            extract_alphapulldown_scalars(&spec.struct_path, &spec.model_label, vfs);
        if scores.iptm.is_none() {
            scores.iptm = ap_iptm;
        }
        if scores.ptm.is_none() {
            scores.ptm = ap_ptm;
        }
    }

    let token_chain_ids = if matches!(detected, "alphafold3" | "openfold3") {
        extract_token_chain_ids(pae_path, vfs)
    } else {
        None
    };

    let pairs = analyze_single_model(
        &struct_text,
        spec.is_pdb,
        &pae,
        pae_dim,
        token_chain_ids.as_deref(),
        &scores,
        pae_cutoff,
        cb_cutoff,
    );

    Ok(pairs
        .iter()
        .map(|p| format_csv_row(&spec.name, &spec.rank, &spec.model_label, p))
        .collect())
}

fn glob_match_simple(pattern: &str, text: &str) -> bool {
    let mut rx = String::from("^");
    for ch in pattern.chars() {
        match ch {
            '*' => rx.push_str(".*"),
            '?' => rx.push('.'),
            '.' | '+' | '(' | ')' | '[' | ']' | '{' | '}' | '^' | '$' | '|' | '\\' => {
                rx.push('\\');
                rx.push(ch);
            }
            _ => rx.push(ch),
        }
    }
    rx.push('$');
    Regex::new(&rx).is_ok_and(|r| r.is_match(text))
}

fn load_manifest_models(
    manifest_val: &Value,
    manifest_dir: &str,
    filenames: &[String],
) -> Option<Vec<ModelSpec>> {
    let obj = manifest_val.as_object()?;
    if !obj.contains_key("models") && !obj.contains_key("pae") {
        return None;
    }
    let pae_key = obj
        .get("pae_key")
        .and_then(|v| v.as_str())
        .map(|s| s.to_string());

    let resolve = |pat: Option<&str>| -> Vec<String> {
        let Some(p) = pat else { return Vec::new() };
        let mut cands = Vec::new();
        if !manifest_dir.is_empty() {
            cands.push(format!("{manifest_dir}/{p}"));
        }
        cands.push(p.to_string());
        for c in cands {
            if filenames.iter().any(|f| f == &c) {
                return vec![c];
            }
        }
        let bn = posix_basename(p);
        filenames
            .iter()
            .filter(|f| glob_match_simple(bn, posix_basename(f)))
            .cloned()
            .collect()
    };

    let mut out = Vec::new();
    if let Some(models_arr) = obj.get("models").and_then(|v| v.as_array()) {
        for (i, entry) in models_arr.iter().enumerate() {
            let s = resolve(entry.get("structure").and_then(|v| v.as_str()));
            let p = resolve(
                entry
                    .get("pae")
                    .or_else(|| entry.get("structure"))
                    .and_then(|v| v.as_str()),
            );
            if s.is_empty() || p.is_empty() {
                continue;
            }
            let sm = resolve(entry.get("summary").and_then(|v| v.as_str()));
            let name = entry
                .get("name")
                .and_then(|v| v.as_str())
                .map(|s| s.to_string())
                .unwrap_or_else(|| format!("model_{i}"));
            let rank = entry
                .get("rank")
                .map(|v| {
                    v.as_str()
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| v.to_string())
                })
                .unwrap_or_else(|| i.to_string());
            let sp = &s[0];
            out.push(ModelSpec {
                name,
                rank,
                model_label: posix_basename(sp).to_string(),
                struct_path: sp.clone(),
                pae_path: Some(p[0].clone()),
                scores_path: sm.first().cloned().or_else(|| Some(p[0].clone())),
                is_pdb: !sp.ends_with(".cif"),
                pae_key: pae_key.clone(),
            });
        }
        return Some(out);
    }

    let mut structs = resolve(obj.get("structure").and_then(|v| v.as_str()));
    if structs.is_empty() {
        structs = filenames
            .iter()
            .filter(|f| f.ends_with(".cif") || f.ends_with(".pdb"))
            .cloned()
            .collect();
    }
    structs.sort_by(|a, b| posix_basename(a).cmp(posix_basename(b)));
    let pae_all = resolve(obj.get("pae").and_then(|v| v.as_str()));
    let summ_all = resolve(obj.get("summary").and_then(|v| v.as_str()));

    let sibling = |files: &[String], sdir: &str, sname: &str, i: usize| -> Option<String> {
        if files.is_empty() {
            return None;
        }
        let same: Vec<&String> = files.iter().filter(|f| posix_dirname(f) == sdir).collect();
        for f in &same {
            if posix_basename(f).contains(sname) {
                return Some((*f).clone());
            }
        }
        for f in files {
            if posix_basename(f).contains(sname) {
                return Some(f.clone());
            }
        }
        if !same.is_empty() {
            Some(same.get(i).copied().unwrap_or(same[0]).clone())
        } else {
            Some(files.get(i).unwrap_or(&files[0]).clone())
        }
    };

    for (i, sp) in structs.iter().enumerate() {
        let sdir = posix_dirname(sp);
        let sname = posix_basename(sp)
            .trim_end_matches(".cif")
            .trim_end_matches(".pdb");
        let pp = sibling(&pae_all, sdir, sname, i);
        let sm = sibling(&summ_all, sdir, sname, i).or_else(|| pp.clone());
        out.push(ModelSpec {
            name: sname.to_string(),
            rank: i.to_string(),
            model_label: posix_basename(sp).to_string(),
            struct_path: sp.clone(),
            pae_path: pp,
            scores_path: sm,
            is_pdb: !sp.ends_with(".cif"),
            pae_key: pae_key.clone(),
        });
    }
    Some(out)
}

fn split_first_two_csv_cols(line: &str) -> (&str, &str) {
    let mut it = line.splitn(3, ',');
    let c0 = it.next().unwrap_or("");
    let c1 = it.next().unwrap_or("");
    (c0, c1)
}

fn sort_csv_file(path: &Path) -> Result<(), String> {
    let f = File::open(path).map_err(|e| e.to_string())?;
    let mut lines = BufReader::new(f).lines();
    let Some(Ok(header)) = lines.next() else {
        return Ok(());
    };
    let mut rows: Vec<String> = lines
        .filter_map(|l| l.ok())
        .filter(|l| !l.trim().is_empty())
        .collect();
    if rows.is_empty() {
        return Ok(());
    }
    rows.sort_by(|a, b| {
        let (na, ra) = split_first_two_csv_cols(a);
        let (nb, rb) = split_first_two_csv_cols(b);
        let za = if ra.is_empty() {
            String::new()
        } else {
            format!("{ra:0>5}")
        };
        let zb = if rb.is_empty() {
            String::new()
        } else {
            format!("{rb:0>5}")
        };
        na.cmp(nb).then_with(|| za.cmp(&zb))
    });
    let out = File::create(path).map_err(|e| e.to_string())?;
    let mut w = BufWriter::new(out);
    writeln!(w, "{header}").map_err(|e| e.to_string())?;
    for r in rows {
        writeln!(w, "{r}").map_err(|e| e.to_string())?;
    }
    Ok(())
}

fn print_usage() {
    eprintln!(
        "Usage: lis-rs <path> [-o OUTPUT] [-d OUTPUT_DIR] [--pae-cutoff 12] [--cb-cutoff 8] \
         [--platform PLATFORM] [--manifest FILE] [-w WORKERS] [--no-skip-existing] [-v]"
    );
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    if args.is_empty() || args.iter().any(|a| a == "-h" || a == "--help") {
        print_usage();
        if args.is_empty() {
            std::process::exit(1);
        }
        return;
    }

    let mut input_path: Option<String> = None;
    let mut output: Option<String> = None;
    let mut output_dir: Option<String> = None;
    let mut pae_cutoff = 12.0_f64;
    let mut cb_cutoff = 8.0_f64;
    let mut platform_arg: Option<String> = None;
    let mut manifest_arg: Option<String> = None;
    let mut workers: usize = 0; // 0 = auto all cores!
    let mut skip_existing = true;
    let mut verbose = false;

    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-o" | "--output" => {
                i += 1;
                output = args.get(i).cloned();
            }
            "-d" | "--output-dir" => {
                i += 1;
                output_dir = args.get(i).cloned();
            }
            "--pae-cutoff" => {
                i += 1;
                pae_cutoff = args
                    .get(i)
                    .and_then(|s| s.parse().ok())
                    .unwrap_or(12.0);
            }
            "--cb-cutoff" => {
                i += 1;
                cb_cutoff = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(8.0);
            }
            "--platform" => {
                i += 1;
                platform_arg = args.get(i).cloned();
            }
            "--manifest" => {
                i += 1;
                manifest_arg = args.get(i).cloned();
            }
            "-w" | "--workers" => {
                i += 1;
                workers = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(1);
            }
            "--no-skip-existing" => {
                skip_existing = false;
            }
            "-v" | "--verbose" => {
                verbose = true;
            }
            "--allow-pickle" => {}
            other if !other.starts_with('-') && input_path.is_none() => {
                input_path = Some(other.to_string());
            }
            other => {
                eprintln!("[LIS] Unknown argument: {other}");
                print_usage();
                std::process::exit(1);
            }
        }
        i += 1;
    }

    let Some(raw_path) = input_path else {
        print_usage();
        std::process::exit(1);
    };
    let path = PathBuf::from(&raw_path);
    if !path.exists() {
        eprintln!("[LIS] ERROR: Path not found: {}", path.display());
        std::process::exit(1);
    }

    if workers > 0 {
        let _ = rayon::ThreadPoolBuilder::new().num_threads(workers).build_global();
    }

    let t_start = Instant::now();
    println!("[LIS] Scanning files from: {}", path.display());
    let vfs = match scan_files(&path) {
        Ok(v) => v,
        Err(e) => {
            eprintln!("[LIS] ERROR: {e}");
            std::process::exit(1);
        }
    };
    println!("[LIS] Found {} files", vfs.filenames.len());
    if vfs.filenames.is_empty() {
        eprintln!("[LIS] ERROR: No files found");
        std::process::exit(1);
    }

    // Optional manifest check.
    let mut manifest_pair: Option<(Value, String)> = None;
    if let Some(ref mp) = manifest_arg {
        if let Ok(txt) = fs::read_to_string(mp) {
            if let Ok(val) = serde_json::from_str::<Value>(&txt) {
                let mdir = Path::new(mp)
                    .parent()
                    .map(|p| p.to_string_lossy().into_owned())
                    .unwrap_or_default();
                manifest_pair = Some((val, mdir));
            }
        }
    } else {
        for f in &vfs.filenames {
            let bn = posix_basename(f).to_ascii_lowercase();
            if bn == "lis.json" || bn == "lis_manifest.json" {
                if let Some(txt) = vfs.read_text(f) {
                    if let Ok(val) = serde_json::from_str::<Value>(&txt) {
                        manifest_pair = Some((val, posix_dirname(f).to_string()));
                        break;
                    }
                }
            }
        }
    }

    let (detected, mut models) = if let Some((ref mval, ref mdir)) = manifest_pair {
        if let Some(m_models) = load_manifest_models(mval, mdir, &vfs.filenames) {
            let det = mval
                .get("platform")
                .and_then(|v| v.as_str())
                .unwrap_or("manifest")
                .to_string();
            println!("[LIS] Platform: {det} (from lis.json manifest)");
            (det, m_models)
        } else {
            let forced = platform_arg.or_else(|| {
                mval.get("platform")
                    .and_then(|v| v.as_str())
                    .map(|s| s.to_string())
            });
            let det = forced.unwrap_or_else(|| detect_platform(&vfs.filenames, &vfs));
            println!("[LIS] Platform: {det}");
            let ms = find_models(&vfs.filenames, &det, &vfs);
            (det, ms)
        }
    } else {
        let det = platform_arg.unwrap_or_else(|| detect_platform(&vfs.filenames, &vfs));
        println!("[LIS] Platform: {det}");
        let ms = find_models(&vfs.filenames, &det, &vfs);
        (det, ms)
    };

    models.retain(|m| !m.name.starts_with("._"));
    if models.is_empty() {
        eprintln!("[LIS] ERROR: No prediction models found");
        std::process::exit(1);
    }
    println!("[LIS] {} model(s) found", models.len());

    // Resolve output CSV path identically to `lis.py` (`-o` overrides filename; if `-o` has a parent dir and `-d` is unset, respect `-o` directly when it contains a separator).
    let out_path: PathBuf = if let (Some(ref out_s), None) = (&output, &output_dir) {
        let op = PathBuf::from(out_s);
        if op.parent().is_some_and(|p| !p.as_os_str().is_empty()) {
            if let Some(parent) = op.parent() {
                let _ = fs::create_dir_all(parent);
            }
            op
        } else {
            let dir = if path.is_dir() {
                path.clone()
            } else {
                path.canonicalize()
                    .ok()
                    .and_then(|p| p.parent().map(|x| x.to_path_buf()))
                    .unwrap_or_else(|| PathBuf::from("."))
            };
            let _ = fs::create_dir_all(&dir);
            dir.join(op.file_name().unwrap())
        }
    } else {
        let rx_archive = Regex::new(r"\.(zip|tar\.gz|tar\.xz|tar\.zstd|tar\.zst|tgz|txz|tar)$").unwrap();
        let default_name = {
            let bn = path
                .file_name()
                .map(|s| s.to_string_lossy().into_owned())
                .unwrap_or_else(|| "output".to_string());
            let stem = rx_archive.replace(&bn, "");
            format!("{stem}_lis_analysis.csv")
        };
        let fname = output
            .as_deref()
            .and_then(|s| Path::new(s).file_name())
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or(default_name);
        let dir = if let Some(ref d) = output_dir {
            PathBuf::from(d)
        } else if path.is_dir() {
            path.clone()
        } else {
            path.canonicalize()
                .ok()
                .and_then(|p| p.parent().map(|x| x.to_path_buf()))
                .unwrap_or_else(|| PathBuf::from("."))
        };
        let _ = fs::create_dir_all(&dir);
        dir.join(fname)
    };

    let mut existing_keys: HashSet<(String, String)> = HashSet::new();
    if skip_existing && out_path.exists() {
        if let Ok(f) = File::open(&out_path) {
            for line in BufReader::new(f).lines().map_while(Result::ok) {
                let (c0, c1) = split_first_two_csv_cols(&line);
                if !c0.is_empty() && c0 != "name" {
                    existing_keys.insert((c0.to_string(), c1.to_string()));
                }
            }
        }
    }

    let models_todo: Vec<&ModelSpec> = models
        .iter()
        .filter(|m| !existing_keys.contains(&(m.name.clone(), m.rank.clone())))
        .collect();
    let total_skipped = models.len() - models_todo.len();
    if models_todo.is_empty() {
        println!("[LIS] All {total_skipped} model(s) already processed.");
        let _ = sort_csv_file(&out_path);
        return;
    }

    if !out_path.exists() || existing_keys.is_empty() {
        let mut f = File::create(&out_path).expect("cannot create output CSV");
        writeln!(f, "{CSV_HEADER}").expect("cannot write CSV header");
    }

    let results: Vec<(&ModelSpec, Result<Vec<String>, String>)> = models_todo
        .par_iter()
        .map(|spec| (*spec, process_one_model(spec, &detected, pae_cutoff, cb_cutoff, &vfs)))
        .collect();

    let mut total_done = 0usize;
    let mut total_failed = 0usize;
    {
        let f = fs::OpenOptions::new()
            .append(true)
            .open(&out_path)
            .expect("cannot append output CSV");
        let mut w = BufWriter::new(f);
        for (spec, res) in results {
            match res {
                Ok(rows) => {
                    total_done += 1;
                    for r in rows {
                        let _ = writeln!(w, "{r}");
                    }
                }
                Err(err) => {
                    total_failed += 1;
                    if verbose || total_failed <= 20 {
                        eprintln!("[LIS] FAIL {} rank={}: {err}", spec.name, spec.rank);
                    }
                }
            }
        }
    }

    let _ = sort_csv_file(&out_path);
    let elapsed = t_start.elapsed().as_secs_f64();
    println!(
        "[LIS] {total_done} done{}{ }in {elapsed:.2}s — {}",
        if total_skipped > 0 {
            format!(", {total_skipped} skipped")
        } else {
            String::new()
        },
        if total_failed > 0 {
            format!(", {total_failed} failed ")
        } else {
            " ".to_string()
        },
        out_path.display()
    );
}
