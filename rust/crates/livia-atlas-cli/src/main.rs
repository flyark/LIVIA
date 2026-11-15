//! `livia-atlas-cli`: High-throughput, low-RAM pure-Rust index & shard builder for LIVIA Atlas.
//!
//! Subcommands:
//!   tested-index <species_dir> <pairs.tsv|csv|.gz|- ...>
//!   paralog-index <species_dir> <cluster.tsv> <identity_note>
//!   extract-set <dataset.zip> [set_id]
//!   offsets-shards <site_repo>
//!   compact-index <site_repo> [--check]
//!   edges-shards <site_repo> [--check]
//!   dedup-stream <inputs...> -o <out.tsv> [--min-ilis 0.0]

use std::collections::{BTreeSet, HashMap, HashSet};
use std::env;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Cursor, Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

use flate2::read::GzDecoder;
use md5::Md5;
use regex::Regex;
use serde_json::{json, Map, Value};
use sha1::{Digest, Sha1};

fn sha1_hex12(keys: &[String]) -> String {
    let mut hasher = Sha1::new();
    hasher.update(keys.join("\n").as_bytes());
    let digest = hasher.finalize();
    let mut s = String::with_capacity(12);
    for b in &digest[..6] {
        s.push_str(&format!("{b:02x}"));
    }
    s
}

fn md5_hex(bytes: &[u8]) -> String {
    let mut hasher = Md5::new();
    hasher.update(bytes);
    let digest = hasher.finalize();
    let mut s = String::with_capacity(32);
    for b in digest {
        s.push_str(&format!("{b:02x}"));
    }
    s
}

/// Format JSON with 1-space indent matching Python `json.dump(..., indent=1)` (`": "` after colon, no trailing newline unless requested).
fn to_json_indent1(val: &Value) -> String {
    let mut buf = Vec::new();
    let formatter = serde_json::ser::PrettyFormatter::with_indent(b" ");
    let mut ser = serde_json::Serializer::with_formatter(&mut buf, formatter);
    serde::Serialize::serialize(val, &mut ser).expect("json serialize");
    String::from_utf8(buf).expect("utf8 json")
}

fn load_species_keys(sp_dir: &Path) -> Result<(Vec<String>, HashMap<String, usize>), String> {
    let p_path = sp_dir.join("proteins.json");
    let txt = fs::read_to_string(&p_path).map_err(|e| format!("{}: {e}", p_path.display()))?;
    let v: Value = serde_json::from_str(&txt).map_err(|e| e.to_string())?;
    let cols = v
        .get("columns")
        .and_then(|c| c.as_array())
        .ok_or("proteins.json missing columns")?;
    let kc = cols
        .iter()
        .position(|c| c.as_str() == Some("key"))
        .ok_or("proteins.json missing 'key' column")?;
    let rows = v
        .get("rows")
        .and_then(|r| r.as_array())
        .ok_or("proteins.json missing rows")?;

    let mut keys = Vec::with_capacity(rows.len());
    let mut row_map = HashMap::with_capacity(rows.len());
    for (idx, r) in rows.iter().enumerate() {
        let k = r
            .as_array()
            .and_then(|a| a.get(kc))
            .and_then(|x| x.as_str())
            .ok_or("invalid key cell")?
            .to_string();
        row_map.insert(k.clone(), idx);
        keys.push(k);
    }
    Ok((keys, row_map))
}

fn open_line_reader(path_str: &str) -> Result<Box<dyn BufRead>, String> {
    if path_str == "-" {
        return Ok(Box::new(BufReader::new(io::stdin())));
    }
    let f = File::open(path_str).map_err(|e| format!("{path_str}: {e}"))?;
    if path_str.to_ascii_lowercase().ends_with(".gz") {
        Ok(Box::new(BufReader::new(GzDecoder::new(f))))
    } else if path_str.to_ascii_lowercase().ends_with(".xz") {
        let mut out = Vec::new();
        lzma_rs::xz_decompress(&mut BufReader::new(f), &mut out).map_err(|e| e.to_string())?;
        Ok(Box::new(BufReader::new(Cursor::new(out))))
    } else {
        Ok(Box::new(BufReader::new(f)))
    }
}

// ============================================================================
// 1. tested-index (`LVTP2\0` binary CSR tested.bin)
// ============================================================================
fn cmd_tested_index(args: &[String]) -> Result<(), String> {
    if args.len() < 2 {
        return Err(
            "Usage: livia-atlas-cli tested-index <species_dir> <pairs.tsv|csv|.gz|- ...>".into(),
        );
    }
    let sp_dir = Path::new(&args[0]);
    let (keys, row_map) = load_species_keys(sp_dir)?;
    let n = keys.len();
    if n >= 65536 {
        return Err(format!("row numbers must fit in 16 bits (got N={n})"));
    }
    let sha = sha1_hex12(&keys);

    let mut best_by_pair: HashMap<(u16, u16), f64> = HashMap::new();
    let mut total_lost = 0usize;

    for file_arg in &args[1..] {
        let mut reader = open_line_reader(file_arg)?;
        let mut header_line = String::new();
        if reader.read_line(&mut header_line).map_err(|e| e.to_string())? == 0 {
            continue;
        }
        let trimmed_h = header_line.trim_end_matches(['\r', '\n']);
        let delim = if trimmed_h.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = trimmed_h.split(delim).collect();

        let ca = cols
            .iter()
            .position(|&c| matches!(c, "protein_a" | "chain_i" | "a" | "bait"))
            .unwrap_or(0);
        let cb = cols
            .iter()
            .position(|&c| matches!(c, "protein_b" | "chain_j" | "b" | "prey"))
            .unwrap_or(1);
        let cs = cols
            .iter()
            .position(|&c| c == "iLIS_best")
            .or_else(|| cols.iter().position(|&c| c == "iLIS"))
            .unwrap_or(2);

        let mut file_rows = 0usize;
        let mut file_lost = 0usize;
        let mut line = String::new();
        loop {
            line.clear();
            let nb = reader.read_line(&mut line).map_err(|e| e.to_string())?;
            if nb == 0 {
                break;
            }
            let l = line.trim_end_matches(['\r', '\n']);
            if l.is_empty() {
                continue;
            }
            file_rows += 1;
            let parts: Vec<&str> = l.split(delim).collect();
            let Some(&ka) = parts.get(ca) else {
                file_lost += 1;
                continue;
            };
            let Some(&kb) = parts.get(cb) else {
                file_lost += 1;
                continue;
            };
            let (Some(&ia), Some(&ib)) = (row_map.get(ka), row_map.get(kb)) else {
                file_lost += 1;
                continue;
            };
            let raw_s = parts
                .get(cs)
                .and_then(|s| s.parse::<f64>().ok())
                .filter(|v| !v.is_nan())
                .unwrap_or(0.0)
                .clamp(0.0, 1.0);
            let i = ia.min(ib) as u16;
            let j = ia.max(ib) as u16;
            best_by_pair
                .entry((i, j))
                .and_modify(|cur| {
                    if raw_s > *cur {
                        *cur = raw_s;
                    }
                })
                .or_insert(raw_s);
        }
        total_lost += file_lost;
        let bn = Path::new(file_arg)
            .file_name()
            .map(|s| s.to_string_lossy().into_owned())
            .unwrap_or_else(|| file_arg.clone());
        println!("{bn} {file_rows} pairs, {file_lost} with a key outside the index");
    }

    let mut sorted_pairs: Vec<((u16, u16), f64)> = best_by_pair.into_iter().collect();
    sorted_pairs.sort_unstable_by_key(|&(ij, _)| ij);

    let p_len = sorted_pairs.len();
    let mut offsets = vec![0u32; n + 1];
    let mut partners = Vec::with_capacity(p_len);
    let mut scores = Vec::with_capacity(p_len);
    let mut homodimers = 0usize;

    for &((i, j), s) in &sorted_pairs {
        offsets[(i as usize) + 1] += 1;
        partners.push(j);
        // Match NumPy `np.round(s * 255)` round-ties-to-even!
        scores.push((s * 255.0).round_ties_even() as u8);
        if i == j {
            homodimers += 1;
        }
    }
    for idx in 1..=n {
        offsets[idx] += offsets[idx - 1];
    }

    let out_bin = sp_dir.join("tested.bin");
    {
        let f = File::create(&out_bin).map_err(|e| e.to_string())?;
        let mut w = BufWriter::new(f);
        w.write_all(b"LVTP2\0").map_err(|e| e.to_string())?;
        w.write_all(&(n as u32).to_le_bytes()).map_err(|e| e.to_string())?;
        w.write_all(&(p_len as u32).to_le_bytes()).map_err(|e| e.to_string())?;
        w.write_all(sha.as_bytes()).map_err(|e| e.to_string())?;
        w.write_all(&[0u8; 6]).map_err(|e| e.to_string())?;
        for &o in &offsets {
            w.write_all(&o.to_le_bytes()).map_err(|e| e.to_string())?;
        }
        for &nb in &partners {
            w.write_all(&nb.to_le_bytes()).map_err(|e| e.to_string())?;
        }
        w.write_all(&scores).map_err(|e| e.to_string())?;
        w.flush().map_err(|e| e.to_string())?;
    }

    let man_path = sp_dir.join("manifest.json");
    if man_path.exists() {
        let mtxt = fs::read_to_string(&man_path).map_err(|e| e.to_string())?;
        let mut man: Value = serde_json::from_str(&mtxt).map_err(|e| e.to_string())?;
        if let Some(obj) = man.as_object_mut() {
            let files = obj
                .entry("files")
                .or_insert_with(|| Value::Object(Map::new()));
            if let Some(fo) = files.as_object_mut() {
                fo.insert("tested".to_string(), Value::String("tested.bin".to_string()));
            }
            let counts = obj
                .entry("counts")
                .or_insert_with(|| Value::Object(Map::new()));
            if let Some(co) = counts.as_object_mut() {
                co.insert("pairsTested".to_string(), json!(p_len));
            }
        }
        fs::write(&man_path, to_json_indent1(&man)).map_err(|e| e.to_string())?;
    }

    println!(
        "{}/tested.bin: N {n}, {p_len} pairs ({homodimers} homodimers), {total_lost} rows left out, index {sha}",
        sp_dir.display()
    );
    Ok(())
}

// ============================================================================
// 2. paralog-index (`paralogs.json`)
// ============================================================================
fn cmd_paralog_index(args: &[String]) -> Result<(), String> {
    if args.len() < 3 {
        return Err(
            "Usage: livia-atlas-cli paralog-index <species_dir> <cluster.tsv> <identity_note>"
                .into(),
        );
    }
    let sp_dir = Path::new(&args[0]);
    let tsv_path = &args[1];
    let note = &args[2];
    let (keys, row_map) = load_species_keys(sp_dir)?;
    let sha = sha1_hex12(&keys);

    let mut fam: HashMap<String, BTreeSet<usize>> = HashMap::new();
    let mut lost = 0usize;
    let reader = open_line_reader(tsv_path)?;
    for line_res in reader.lines() {
        let line = line_res.map_err(|e| e.to_string())?;
        let l = line.trim_end_matches(['\r', '\n']);
        if l.is_empty() {
            continue;
        }
        let Some((rep, mem)) = l.split_once('\t') else {
            lost += 1;
            continue;
        };
        if let (Some(&ir), Some(&im)) = (row_map.get(rep), row_map.get(mem)) {
            let set = fam.entry(rep.to_string()).or_default();
            set.insert(im);
            set.insert(ir);
        } else {
            lost += 1;
        }
    }

    let mut fams: Vec<Vec<usize>> = fam
        .into_values()
        .filter(|s| s.len() >= 2)
        .map(|s| s.into_iter().collect())
        .collect();
    fams.sort_by_key(|f| f[0]);

    let covered: usize = fams.iter().map(|f| f.len()).sum();
    let mut out_obj = Map::new();
    out_obj.insert("index".to_string(), Value::String(sha));
    out_obj.insert("identity".to_string(), Value::String(note.clone()));
    out_obj.insert("families".to_string(), json!(fams));

    let out_path = sp_dir.join("paralogs.json");
    let compact = serde_json::to_string(&Value::Object(out_obj)).map_err(|e| e.to_string())?;
    fs::write(&out_path, compact).map_err(|e| e.to_string())?;

    let man_path = sp_dir.join("manifest.json");
    if man_path.exists() {
        let mtxt = fs::read_to_string(&man_path).map_err(|e| e.to_string())?;
        let mut man: Value = serde_json::from_str(&mtxt).map_err(|e| e.to_string())?;
        if let Some(obj) = man.as_object_mut() {
            let files = obj
                .entry("files")
                .or_insert_with(|| Value::Object(Map::new()));
            if let Some(fo) = files.as_object_mut() {
                fo.insert(
                    "paralogs".to_string(),
                    Value::String("paralogs.json".to_string()),
                );
            }
        }
        fs::write(&man_path, to_json_indent1(&man)).map_err(|e| e.to_string())?;
    }

    println!(
        "{}/paralogs.json: {} families of 2+ covering {covered} rows; {lost} table rows outside the index",
        sp_dir.display(),
        fams.len()
    );
    Ok(())
}

// ============================================================================
// 3. extract-set (`tools/extract_set.py` equivalent)
// ============================================================================
fn parse_csv_line(line: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut chars = line.chars().peekable();
    let mut in_quotes = false;
    while let Some(c) = chars.next() {
        if in_quotes {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    cur.push('"');
                    chars.next();
                } else {
                    in_quotes = false;
                }
            } else {
                cur.push(c);
            }
        } else if c == '"' {
            in_quotes = true;
        } else if c == ',' {
            out.push(std::mem::take(&mut cur));
        } else {
            cur.push(c);
        }
    }
    out.push(cur);
    out
}

fn format_csv_fields(fields: &[String]) -> String {
    fields
        .iter()
        .map(|f| {
            if f.contains(',') || f.contains('"') || f.contains('\n') {
                format!("\"{}\"", f.replace('"', "\"\""))
            } else {
                f.clone()
            }
        })
        .collect::<Vec<_>>()
        .join(",")
}

fn fmt_thousands(mut n: u64) -> String {
    if n == 0 {
        return "0".to_string();
    }
    let mut groups = Vec::new();
    while n > 0 {
        groups.push(n % 1000);
        n /= 1000;
    }
    groups.reverse();
    let mut s = groups[0].to_string();
    for &g in &groups[1..] {
        s.push_str(&format!(",{g:03}"));
    }
    s
}

fn cmd_extract_set(args: &[String]) -> Result<(), String> {
    if args.is_empty() {
        return Err("Usage: livia-atlas-cli extract-set <dataset.zip> [set_id]".into());
    }
    let f = File::open(&args[0]).map_err(|e| format!("{}: {e}", args[0]))?;
    let mut z = zip::ZipArchive::new(BufReader::new(f)).map_err(|e| e.to_string())?;
    let sets_val: Value = {
        let mut sf = z.by_name("sets.json").map_err(|e| e.to_string())?;
        let mut s = String::new();
        sf.read_to_string(&mut s).map_err(|e| e.to_string())?;
        serde_json::from_str(&s).map_err(|e| e.to_string())?
    };
    let sets = sets_val
        .get("sets")
        .and_then(|v| v.as_array())
        .ok_or("sets.json missing 'sets'")?;

    if args.len() < 2 {
        for s in sets {
            let id = s.get("id").and_then(|v| v.as_str()).unwrap_or("");
            let tp = s.get("type").and_then(|v| v.as_str()).unwrap_or("");
            let cnt = s
                .get("counts")
                .and_then(|c| c.get("predictions"))
                .and_then(|v| v.as_u64())
                .unwrap_or(0);
            let title = s.get("title").and_then(|v| v.as_str()).unwrap_or("");
            println!("{id:26} {tp:9} {:>10} predictions  {title}", fmt_thousands(cnt));
        }
        return Ok(());
    }

    let target_id = &args[1];
    let Some(chosen) = sets
        .iter()
        .find(|s| s.get("id").and_then(|v| v.as_str()) == Some(target_id.as_str()))
    else {
        let all_ids: Vec<&str> = sets
            .iter()
            .filter_map(|s| s.get("id").and_then(|v| v.as_str()))
            .collect();
        return Err(format!(
            "No set {target_id:?}; the sets are: {}",
            all_ids.join(", ")
        ));
    };

    let runs: HashSet<String> = chosen
        .get("runs")
        .and_then(|v| v.as_array())
        .ok_or("set missing runs")?
        .iter()
        .map(|v| {
            v.as_str()
                .map(|s| s.to_string())
                .unwrap_or_else(|| v.to_string())
        })
        .collect();

    let prot_file = chosen
        .get("files")
        .and_then(|f| f.get("proteins"))
        .and_then(|v| v.as_str())
        .ok_or("set missing files.proteins")?;

    let genes: Vec<String> = {
        let mut pf = z.by_name(prot_file).map_err(|e| e.to_string())?;
        let mut txt = String::new();
        pf.read_to_string(&mut txt).map_err(|e| e.to_string())?;
        let pv: Value = serde_json::from_str(&txt).map_err(|e| e.to_string())?;
        pv.get("rows")
            .and_then(|r| r.as_array())
            .ok_or("proteins missing rows")?
            .iter()
            .filter_map(|r| r.get(0).and_then(|x| x.as_str()).map(|s| s.to_string()))
            .collect()
    };

    let stdout = io::stdout();
    let mut out = BufWriter::new(stdout.lock());
    let mut seen: HashSet<(String, String, String)> = HashSet::new();
    let mut wrote_header = false;
    let mut count = 0u64;

    for gene in &genes {
        let entry_name = format!("b/{gene}.zip");
        let mut inner_bytes = Vec::new();
        {
            let mut bf = z.by_name(&entry_name).map_err(|e| e.to_string())?;
            bf.read_to_end(&mut inner_bytes).map_err(|e| e.to_string())?;
        }
        let mut bz = zip::ZipArchive::new(Cursor::new(inner_bytes)).map_err(|e| e.to_string())?;
        let mut csv_name = None;
        for idx in 0..bz.len() {
            if let Ok(f) = bz.by_index(idx) {
                if f.name().ends_with(".csv") {
                    csv_name = Some(f.name().to_string());
                    break;
                }
            }
        }
        let Some(cname) = csv_name else { continue };
        let mut cf = bz.by_name(&cname).map_err(|e| e.to_string())?;
        let mut csv_txt = String::new();
        cf.read_to_string(&mut csv_txt).map_err(|e| e.to_string())?;

        let mut lines = csv_txt.lines();
        let Some(hline) = lines.next() else { continue };
        let h = parse_csv_line(hline.trim_end_matches('\r'));
        if !wrote_header {
            writeln!(out, "{}", format_csv_fields(&h)).map_err(|e| e.to_string())?;
            wrote_header = true;
        }
        let Some(bi) = h.iter().position(|c| c == "batch") else { continue };
        let Some(ni) = h.iter().position(|c| c == "name") else { continue };
        let Some(ri) = h.iter().position(|c| c == "rank") else { continue };

        for rline in lines {
            let trimmed = rline.trim_end_matches('\r');
            if trimmed.is_empty() {
                continue;
            }
            let row = parse_csv_line(trimmed);
            if row.len() <= bi.max(ni).max(ri) {
                continue;
            }
            if runs.contains(&row[bi]) {
                let key = (row[ni].clone(), row[ri].clone(), row[bi].clone());
                if seen.insert(key) {
                    writeln!(out, "{}", format_csv_fields(&row)).map_err(|e| e.to_string())?;
                    count += 1;
                }
            }
        }
    }
    out.flush().map_err(|e| e.to_string())?;
    let title = chosen.get("title").and_then(|v| v.as_str()).unwrap_or("");
    eprintln!(
        "{} predictions of {title} ({} proteins)",
        fmt_thousands(count),
        fmt_thousands(genes.len() as u64)
    );
    Ok(())
}

// ============================================================================
// 4. offsets-shards (`build/offsets_shards.py`)
// ============================================================================
fn shard_offsets_map(src: &Path, out_dir: &Path) -> Result<(usize, Vec<u64>), String> {
    let _ = fs::remove_dir_all(out_dir);
    fs::create_dir_all(out_dir).map_err(|e| e.to_string())?;
    let txt = fs::read_to_string(src).map_err(|e| e.to_string())?;
    let val: Value = serde_json::from_str(&txt).map_err(|e| e.to_string())?;
    let obj = val.as_object().ok_or("offsets.json must be object")?;

    let rx_ext = Regex::new(r"\.[^.]+$").unwrap();
    let rx_en = Regex::new(r"^([A-Z0-9]+)_[A-Z0-9]+$").unwrap();
    let rx_san = Regex::new(r"[^a-z0-9_-]").unwrap();

    let mut shard_order: Vec<String> = Vec::new();
    let mut shards: HashMap<String, Map<String, Value>> = HashMap::new();

    for (rel, at) in obj {
        let base = rel.rsplit('/').next().unwrap_or(rel);
        let mut name = rx_ext.replace(base, "").to_string();
        if let Some(c) = rx_en.captures(&name) {
            name = c[1].to_string();
        }
        let k_raw = if name.chars().count() >= 2 {
            let chars: Vec<char> = name.chars().collect();
            chars[chars.len() - 2..]
                .iter()
                .collect::<String>()
                .to_lowercase()
        } else {
            "_".to_string()
        };
        let k = rx_san.replace_all(&k_raw, "_").to_string();
        shards
            .entry(k.clone())
            .or_insert_with(|| {
                shard_order.push(k);
                Map::new()
            })
            .insert(rel.clone(), at.clone());
    }

    let mut sizes = Vec::with_capacity(shard_order.len());
    for k in shard_order {
        let m = shards.remove(&k).unwrap();
        let sp = out_dir.join(format!("{k}.json"));
        let s = serde_json::to_string(&Value::Object(m)).map_err(|e| e.to_string())?;
        fs::write(&sp, s.as_bytes()).map_err(|e| e.to_string())?;
        sizes.push(s.len() as u64);
    }
    Ok((obj.len(), sizes))
}

fn cmd_offsets_shards(args: &[String]) -> Result<(), String> {
    if args.is_empty() {
        return Err("Usage: livia-atlas-cli offsets-shards <site_repo>".into());
    }
    let site = Path::new(&args[0]);
    let reg_txt = fs::read_to_string(site.join("datasets.json")).map_err(|e| e.to_string())?;
    let reg: Value = serde_json::from_str(&reg_txt).map_err(|e| e.to_string())?;
    let datasets = reg
        .get("datasets")
        .and_then(|v| v.as_array())
        .ok_or("datasets.json missing datasets")?;
    let min_bytes: u64 = 100 * 1024;
    let mut tot = 0u64;

    for d in datasets {
        let id = d.get("id").and_then(|v| v.as_str()).unwrap_or("");
        let Some(z) = d.get("zip") else { continue };
        let Some(off) = z.get("offsets").and_then(|v| v.as_str()) else { continue };
        let src = site.join(off);
        if !src.exists() {
            println!("no offsets file for {id} {off}");
            continue;
        }
        let sd = src.parent().unwrap();
        let out = sd.join("offsets");
        let _ = fs::remove_dir_all(&out);
        let mp = sd.join("manifest.json");
        let mut man: Value =
            serde_json::from_str(&fs::read_to_string(&mp).map_err(|e| e.to_string())?)
                .map_err(|e| e.to_string())?;

        if fs::metadata(&src).map(|m| m.len()).unwrap_or(0) < min_bytes {
            let mut popped = false;
            if let Some(fo) = man.get_mut("files").and_then(|f| f.as_object_mut()) {
                popped = fo.shift_remove("offsetShards").is_some();
            }
            if popped {
                fs::write(&mp, to_json_indent1(&man)).map_err(|e| e.to_string())?;
            }
            continue;
        }

        let (n, sizes) = shard_offsets_map(&src, &out)?;
        tot += sizes.iter().sum::<u64>();
        if let Some(obj) = man.as_object_mut() {
            let files = obj
                .entry("files")
                .or_insert_with(|| Value::Object(Map::new()));
            if let Some(fo) = files.as_object_mut() {
                fo.insert("offsetShards".to_string(), Value::Bool(true));
            }
        }
        fs::write(&mp, to_json_indent1(&man)).map_err(|e| e.to_string())?;
        let max_kb = sizes.iter().copied().max().unwrap_or(0) / 1024;
        println!(
            "{id:<34} {n:>7} entries -> {:>4} shards, largest {max_kb} KB",
            sizes.len()
        );

        if let Some(fb) = z
            .get("fallback")
            .and_then(|f| f.get("offsets"))
            .and_then(|v| v.as_str())
        {
            let fsrc = site.join(fb);
            if !fsrc.exists() {
                println!("no fallback offsets file for {id} {fb}");
                continue;
            }
            let fout = PathBuf::from(fsrc.to_string_lossy().trim_end_matches(".json"));
            let (fn_cnt, fsizes) = shard_offsets_map(&fsrc, &fout)?;
            tot += fsizes.iter().sum::<u64>();
            let fmax_kb = fsizes.iter().copied().max().unwrap_or(0) / 1024;
            println!(
                "{:<34} {fn_cnt:>7} entries -> {:>4} shards, largest {fmax_kb} KB",
                format!("{id} (fallback)"),
                fsizes.len()
            );
        }
    }
    println!("total {} MB of shards", tot / 1_048_576);
    Ok(())
}

// ============================================================================
// 5. compact-index (`build/compact_index.py`)
// ============================================================================
fn encode_compact_proteins(p: &Value) -> Result<Value, String> {
    let cols_val = p
        .get("columns")
        .and_then(|v| v.as_array())
        .ok_or("missing columns")?;
    let col_idx = |name: &str| -> Result<usize, String> {
        cols_val
            .iter()
            .position(|c| c.as_str() == Some(name))
            .ok_or_else(|| format!("missing col {name}"))
    };
    let i_key = col_idx("key")?;
    let i_acc = col_idx("acc")?;
    let i_len = col_idx("len")?;
    let i_clen = col_idx("clen")?;
    let i_occ = col_idx("occ")?;

    let rows_val = p
        .get("rows")
        .and_then(|v| v.as_array())
        .ok_or("missing rows")?;
    let mut out_rows = Vec::with_capacity(rows_val.len());

    for r in rows_val {
        let r_arr = r.as_array().ok_or("row not array")?;
        let mut r2 = r_arr.clone();
        let key = r_arr[i_key].as_str().unwrap_or("");
        let occ_str = r_arr[i_occ].as_str().unwrap_or("");

        let mut toks: Vec<(&str, &str)> = Vec::new();
        let mut other_order: Vec<&str> = Vec::new();
        let mut other_counts: HashMap<&str, usize> = HashMap::new();

        for t in occ_str.split_whitespace() {
            let (d, n) = t.split_once(':').unwrap_or((t, ""));
            toks.push((d, n));
            if n != key {
                other_counts
                    .entry(n)
                    .and_modify(|cnt| *cnt += 1)
                    .or_insert_with(|| {
                        other_order.push(n);
                        1
                    });
            }
        }

        // First-seen tie-breaking matching Python `Counter(others).most_common(1)[0][0]`.
        let mut en = "";
        let mut best_cnt = 0usize;
        for cand in other_order {
            let cnt = other_counts[cand];
            if cnt > best_cnt {
                best_cnt = cnt;
                en = cand;
            }
        }

        let new_occ = toks
            .into_iter()
            .map(|(d, n)| {
                if n == key {
                    d.to_string()
                } else if n == en {
                    format!("{d}.")
                } else {
                    format!("{d}:{n}")
                }
            })
            .collect::<Vec<_>>()
            .join(" ");
        r2[i_occ] = Value::String(new_occ);
        if r2[i_acc].as_str() == Some(key) {
            r2[i_acc] = json!(0);
        }
        if r2[i_clen] == r2[i_len] {
            r2[i_clen] = json!(0);
        }
        r2.push(Value::String(en.to_string()));
        out_rows.push(Value::Array(r2));
    }

    let mut new_cols = cols_val.clone();
    new_cols.push(Value::String("en".to_string()));

    let mut q = Map::new();
    q.insert(
        "format".to_string(),
        Value::String("livia-atlas-index/1".to_string()),
    );
    q.insert("columns".to_string(), Value::Array(new_cols));
    q.insert("rows".to_string(), Value::Array(out_rows));
    Ok(Value::Object(q))
}

fn cmd_compact_index(args: &[String]) -> Result<(), String> {
    if args.is_empty() {
        return Err("Usage: livia-atlas-cli compact-index <site_repo> [--check]".into());
    }
    let site = Path::new(&args[0]);
    let check = args.iter().any(|a| a == "--check");
    let reg_path = site.join("datasets.json");
    let mut reg: Value =
        serde_json::from_str(&fs::read_to_string(&reg_path).map_err(|e| e.to_string())?)
            .map_err(|e| e.to_string())?;

    let mut bad = 0usize;
    let mut changed = false;

    if let Some(species) = reg.get_mut("species").and_then(|v| v.as_array_mut()) {
        for s in species.iter_mut() {
            let sid = s
                .get("id")
                .and_then(|v| v.as_str())
                .unwrap_or("")
                .to_string();
            let base = s
                .get("base")
                .and_then(|v| v.as_str())
                .unwrap_or("")
                .to_string();
            let pp = site.join(&base).join("proteins.json");
            let ip = site.join(&base).join("index.json");
            if !pp.exists() {
                continue;
            }
            let p_bytes = fs::read(&pp).map_err(|e| e.to_string())?;
            let md5 = md5_hex(&p_bytes);
            if check {
                if s.get("index").is_some() {
                    if !ip.exists() {
                        bad += 1;
                        println!("{sid} index.json missing");
                        continue;
                    }
                    let iv: Value = serde_json::from_str(
                        &fs::read_to_string(&ip).map_err(|e| e.to_string())?,
                    )
                    .map_err(|e| e.to_string())?;
                    let got = iv.get("source_md5").and_then(|v| v.as_str()).unwrap_or("");
                    if got != md5 {
                        bad += 1;
                        println!("{sid}: index.json was built from another proteins.json (md5 {got}, now {md5}); rerun compact_index.py");
                    }
                }
                continue;
            }

            let p_val: Value = serde_json::from_slice(&p_bytes).map_err(|e| e.to_string())?;
            let mut q = encode_compact_proteins(&p_val)?;
            if let Some(qo) = q.as_object_mut() {
                qo.insert("source_md5".to_string(), Value::String(md5));
            }
            let q_str = serde_json::to_string(&q).map_err(|e| e.to_string())?;
            fs::write(&ip, q_str.as_bytes()).map_err(|e| e.to_string())?;
            if s.get("index").and_then(|v| v.as_str()) != Some("index.json") {
                if let Some(so) = s.as_object_mut() {
                    so.insert("index".to_string(), Value::String("index.json".to_string()));
                    changed = true;
                }
            }
            if p_bytes.len() > 1_000_000 {
                println!(
                    "{sid:<34} {:>6} KB -> {:>6} KB",
                    p_bytes.len() / 1024,
                    q_str.len() / 1024
                );
            }
        }
    }

    if let Some(species) = reg.get("species").and_then(|v| v.as_array()) {
        for s in species {
            let sid = s.get("id").and_then(|v| v.as_str()).unwrap_or("");
            let base = s.get("base").and_then(|v| v.as_str()).unwrap_or("");
            let ep = site.join(base).join("edges.tsv");
            let gp = site.join(base).join("edges.grp");
            let mp = site.join(base).join("manifest.json");
            if !ep.exists() || !mp.exists() {
                continue;
            }
            let mut man: Value =
                serde_json::from_str(&fs::read_to_string(&mp).map_err(|e| e.to_string())?)
                    .map_err(|e| e.to_string())?;
            let e_bytes = fs::read(&ep).map_err(|e| e.to_string())?;
            let md5 = md5_hex(&e_bytes);

            if check {
                let has_grp = man
                    .get("files")
                    .and_then(|f| f.get("edgesGrp"))
                    .is_some();
                if has_grp {
                    let got = fs::read_to_string(&gp)
                        .ok()
                        .and_then(|t| {
                            t.lines()
                                .next()
                                .and_then(|l| l.split('\t').nth(1))
                                .map(|x| x.to_string())
                        })
                        .unwrap_or_default();
                    if got != md5 {
                        bad += 1;
                        println!(
                            "{sid}: edges.grp was built from another edges.tsv; rerun compact_index.py"
                        );
                    }
                }
                continue;
            }

            if e_bytes.len() < 1_000_000 {
                continue;
            }
            let e_txt = String::from_utf8(e_bytes).map_err(|e| e.to_string())?;
            let mut lines = e_txt.split('\n');
            let head = lines.next().unwrap_or("");
            if head != "i\tj\tiLIS_best\tipTM_best\tsrc" {
                println!("{sid} edges kept as they are: header {head}");
                continue;
            }

            let mut rows: Vec<(i64, i64, i64, Option<i64>, i64)> = Vec::new();
            let mut ok = true;
            for l in lines {
                if l.is_empty() {
                    continue;
                }
                let mut sp = l.split('\t');
                let (Some(a), Some(b), Some(c), Some(d), Some(e)) = (
                    sp.next(),
                    sp.next(),
                    sp.next(),
                    sp.next(),
                    sp.next(),
                ) else {
                    ok = false;
                    break;
                };
                let (Ok(ai), Ok(bi), Ok(cf), Ok(ei)) = (
                    a.parse::<i64>(),
                    b.parse::<i64>(),
                    c.parse::<f64>(),
                    e.parse::<i64>(),
                ) else {
                    ok = false;
                    break;
                };
                let ci = (cf * 1000.0).round_ties_even() as i64;
                if ((ci as f64) / 1000.0 - cf).abs() > 1e-9 {
                    ok = false;
                    break;
                }
                let di = if d.is_empty() {
                    None
                } else {
                    let Ok(df) = d.parse::<f64>() else {
                        ok = false;
                        break;
                    };
                    let dv = (df * 100.0).round_ties_even() as i64;
                    if ((dv as f64) / 100.0 - df).abs() > 1e-9 {
                        ok = false;
                        break;
                    }
                    Some(dv)
                };
                rows.push((ai, bi, ci, di, ei));
            }
            if !ok {
                println!("{sid} edges kept as they are: more decimals than 3 / 2");
                continue;
            }
            rows.sort_by_key(|&(a, b, _, _, _)| (a, b));

            let mut out = String::with_capacity(e_txt.len() / 2);
            out.push_str(&format!("livia-edges/1\t{md5}\n"));
            let mut pi = -1i64;
            let mut pj = 0i64;
            for &(a, b, c, d, e) in &rows {
                if a != pi {
                    out.push_str(&format!("#{}\n", a - pi));
                    pi = a;
                    pj = a;
                }
                if let Some(dv) = d {
                    out.push_str(&format!("{}\t{c}\t{dv}\t{e}\n", b - pj));
                } else {
                    out.push_str(&format!("{}\t{c}\t\t{e}\n", b - pj));
                }
                pj = b;
            }
            fs::write(&gp, out.as_bytes()).map_err(|e| e.to_string())?;
            let mut updated_man = false;
            if let Some(mo) = man.as_object_mut() {
                let fl = mo
                    .entry("files")
                    .or_insert_with(|| Value::Object(Map::new()));
                if let Some(fo) = fl.as_object_mut() {
                    if fo.get("edgesGrp").and_then(|v| v.as_str()) != Some("edges.grp") {
                        fo.insert("edgesGrp".to_string(), Value::String("edges.grp".to_string()));
                        updated_man = true;
                    }
                }
            }
            if updated_man {
                fs::write(&mp, to_json_indent1(&man)).map_err(|e| e.to_string())?;
            }
            println!(
                "{sid:<34} edges {:>6} KB -> {:>6} KB",
                e_txt.len() / 1024,
                out.len() / 1024
            );
        }
    }

    if check {
        println!("index check: {}", if bad > 0 { "FAIL" } else { "ok" });
        if bad > 0 {
            std::process::exit(1);
        }
        return Ok(());
    }
    if changed {
        fs::write(&reg_path, format!("{}\n", to_json_indent1(&reg)))
            .map_err(|e| e.to_string())?;
    }
    Ok(())
}

// ============================================================================
// 6. edges-shards (`build/edges_shards.py`)
// ============================================================================
type GrpRow = (i64, i64, String, String, String);

fn decode_edges_grp(txt: &str) -> Vec<GrpRow> {
    let mut out = Vec::new();
    let mut i = -1i64;
    let mut j = 0i64;
    for l in txt.split('\n').skip(1) {
        if l.is_empty() {
            continue;
        }
        if let Some(rest) = l.strip_prefix('#') {
            if let Ok(d) = rest.parse::<i64>() {
                i += d;
                j = i;
            }
            continue;
        }
        let mut sp = l.split('\t');
        let (Some(dj_s), Some(c), Some(d), Some(e)) = (sp.next(), sp.next(), sp.next(), sp.next())
        else {
            continue;
        };
        if let Ok(dj) = dj_s.parse::<i64>() {
            j += dj;
            out.push((i, j, c.to_string(), d.to_string(), e.to_string()));
        }
    }
    out
}

fn encode_edges_grp(rows: &[GrpRow], head: &str) -> String {
    let mut out = String::new();
    out.push_str(head);
    out.push('\n');
    let mut pi = -1i64;
    let mut pj = 0i64;
    for (a, b, c, d, e) in rows {
        if *a != pi {
            out.push_str(&format!("#{}\n", *a - pi));
            pi = *a;
            pj = *a;
        }
        out.push_str(&format!("{}\t{c}\t{d}\t{e}\n", *b - pj));
        pj = *b;
    }
    out
}

fn cmd_edges_shards(args: &[String]) -> Result<(), String> {
    if args.is_empty() {
        return Err("Usage: livia-atlas-cli edges-shards <site_repo> [--check]".into());
    }
    let site = Path::new(&args[0]);
    let check = args.iter().any(|a| a == "--check");
    let n_shards = 256usize;
    let reg: Value =
        serde_json::from_str(&fs::read_to_string(site.join("datasets.json")).map_err(|e| e.to_string())?)
            .map_err(|e| e.to_string())?;
    let species = reg
        .get("species")
        .and_then(|v| v.as_array())
        .ok_or("datasets.json missing species")?;

    let mut bad = 0usize;
    for s in species {
        let sid = s.get("id").and_then(|v| v.as_str()).unwrap_or("");
        let base = site.join(s.get("base").and_then(|v| v.as_str()).unwrap_or(""));
        let mp = base.join("manifest.json");
        let gp = base.join("edges.grp");
        if !mp.exists() {
            continue;
        }
        let mut man: Value =
            serde_json::from_str(&fs::read_to_string(&mp).map_err(|e| e.to_string())?)
                .map_err(|e| e.to_string())?;
        let has_grp = man.get("files").and_then(|f| f.get("edgesGrp")).is_some();
        if !has_grp || !gp.exists() {
            continue;
        }
        let txt = fs::read_to_string(&gp).map_err(|e| e.to_string())?;
        let head = txt.lines().next().unwrap_or("");

        if check {
            let has_shards = man.get("files").and_then(|f| f.get("edgesShards")).is_some();
            if !has_shards {
                continue;
            }
            for k in 0..n_shards {
                let p = base.join("edges").join(format!("{k}.grp"));
                let ok = fs::read_to_string(&p)
                    .ok()
                    .and_then(|t| t.lines().next().map(|l| l == head))
                    .unwrap_or(false);
                if !ok {
                    bad += 1;
                    println!("{sid}: edges/{k}.grp is missing or was built from another edges.grp; rerun edges_shards.py");
                    break;
                }
            }
            continue;
        }

        let rows = decode_edges_grp(&txt);
        let mut buckets: Vec<Vec<GrpRow>> = vec![Vec::new(); n_shards];
        for r in &rows {
            let k = (r.0.min(r.1).rem_euclid(n_shards as i64)) as usize;
            buckets[k].push(r.clone());
        }
        let mut shards = Vec::with_capacity(n_shards);
        let mut back = Vec::with_capacity(rows.len());
        for mut b in buckets {
            b.sort();
            let enc = encode_edges_grp(&b, head);
            back.extend(decode_edges_grp(&enc));
            shards.push(enc);
        }
        let mut orig_sorted = rows.clone();
        orig_sorted.sort();
        back.sort();
        if back != orig_sorted {
            return Err(format!("{sid}: round-trip verification failed!"));
        }

        let edges_dir = base.join("edges");
        fs::create_dir_all(&edges_dir).map_err(|e| e.to_string())?;
        let mut sz = 0usize;
        for (k, stxt) in shards.iter().enumerate() {
            sz += stxt.len();
            fs::write(edges_dir.join(format!("{k}.grp")), stxt.as_bytes())
                .map_err(|e| e.to_string())?;
        }
        let mut updated = false;
        if let Some(mo) = man.as_object_mut() {
            let fl = mo
                .entry("files")
                .or_insert_with(|| Value::Object(Map::new()));
            if let Some(fo) = fl.as_object_mut() {
                if fo.get("edgesShards").and_then(|v| v.as_u64()) != Some(n_shards as u64) {
                    fo.insert("edgesShards".to_string(), json!(n_shards));
                    updated = true;
                }
            }
        }
        if updated {
            fs::write(&mp, to_json_indent1(&man)).map_err(|e| e.to_string())?;
        }
        println!(
            "{sid:<20} {:>9} pairs -> {n_shards} shards, {:>6} KB ({} KB each)",
            fmt_thousands(rows.len() as u64),
            sz / 1024,
            sz / n_shards / 1024
        );
    }

    if check {
        println!("edge shards check: {}", if bad > 0 { "FAIL" } else { "ok" });
        if bad > 0 {
            std::process::exit(1);
        }
    }
    Ok(())
}

// ============================================================================
// 7. dedup-stream (100M+ model streaming warehouse pair deduplicator)
// ============================================================================
#[derive(Clone, Copy)]
struct AggPair {
    best_ilis: f32,
    best_iptm: f32,
    best_ipsae: f32,
    best_lipdockq: f32,
    models: u32,
}

fn ingest_table_stream<R: BufRead>(
    mut reader: R,
    intern: &mut HashMap<String, u32>,
    rev_intern: &mut Vec<String>,
    agg: &mut HashMap<u64, AggPair>,
) -> Result<u64, String> {
    let mut header = String::new();
    if reader.read_line(&mut header).map_err(|e| e.to_string())? == 0 {
        return Ok(0);
    }
    let h_trim = header.trim_end_matches(['\r', '\n']);
    let delim = if h_trim.contains('\t') { '\t' } else { ',' };
    let cols: Vec<&str> = h_trim.split(delim).collect();

    let col_pair_name = cols.iter().position(|&c| c == "name");
    let col_a = cols
        .iter()
        .position(|&c| matches!(c, "protein_a" | "gene_a" | "a"));
    let col_b = cols
        .iter()
        .position(|&c| matches!(c, "protein_b" | "gene_b" | "b"));
    let col_ilis = cols
        .iter()
        .position(|&c| matches!(c, "iLIS_best" | "iLIS"))
        .ok_or("table missing iLIS/iLIS_best column")?;
    let col_iptm = cols.iter().position(|&c| matches!(c, "ipTM_best" | "ipTM"));
    let col_ipsae = cols.iter().position(|&c| matches!(c, "ipSAE_best" | "ipSAE"));
    let col_lipdq = cols
        .iter()
        .position(|&c| matches!(c, "LIpDockQ_best" | "LIpDockQ"));

    let mut get_id = |k: &str| -> u32 {
        if let Some(&id) = intern.get(k) {
            id
        } else {
            let id = rev_intern.len() as u32;
            let owned = k.to_string();
            intern.insert(owned.clone(), id);
            rev_intern.push(owned);
            id
        }
    };

    let mut rows_read = 0u64;
    let mut line = String::new();
    loop {
        line.clear();
        if reader.read_line(&mut line).map_err(|e| e.to_string())? == 0 {
            break;
        }
        let l = line.trim_end_matches(['\r', '\n']);
        if l.is_empty() {
            continue;
        }
        rows_read += 1;
        // Fast split up to quote-free prefix columns (`name` .. `LIpDockQ` are all before quoted LIR_indices!).
        let parts: Vec<&str> = l.split(delim).collect();
        let (ka, kb) = if let (Some(ia), Some(ib)) = (col_a, col_b) {
            let Some(&sa) = parts.get(ia) else { continue };
            let Some(&sb) = parts.get(ib) else { continue };
            (sa, sb)
        } else if let Some(inm) = col_pair_name {
            let Some(&nm) = parts.get(inm) else { continue };
            if let Some((x, y)) = nm.split_once("___") {
                (x, y)
            } else if let Some((x, y)) = nm.split_once("__") {
                (x, y)
            } else if let Some((x, y)) = nm.split_once('_') {
                (x, y)
            } else {
                continue;
            }
        } else {
            continue;
        };

        let id_a = get_id(ka);
        let id_b = get_id(kb);
        let lo = id_a.min(id_b) as u64;
        let hi = id_a.max(id_b) as u64;
        let key64 = (lo << 32) | hi;

        let parse_f32 = |idx_opt: Option<usize>| -> f32 {
            idx_opt
                .and_then(|idx| parts.get(idx))
                .and_then(|s| s.parse::<f32>().ok())
                .filter(|v| !v.is_nan())
                .unwrap_or(f32::NEG_INFINITY)
        };
        let ilis = parse_f32(Some(col_ilis)).max(0.0);
        let iptm = parse_f32(col_iptm);
        let ipsae = parse_f32(col_ipsae);
        let lipdq = parse_f32(col_lipdq);

        agg.entry(key64)
            .and_modify(|e| {
                if ilis > e.best_ilis {
                    e.best_ilis = ilis;
                }
                if iptm > e.best_iptm {
                    e.best_iptm = iptm;
                }
                if ipsae > e.best_ipsae {
                    e.best_ipsae = ipsae;
                }
                if lipdq > e.best_lipdockq {
                    e.best_lipdockq = lipdq;
                }
                e.models += 1;
            })
            .or_insert(AggPair {
                best_ilis: ilis,
                best_iptm: iptm,
                best_ipsae: ipsae,
                best_lipdockq: lipdq,
                models: 1,
            });
    }
    Ok(rows_read)
}

fn cmd_dedup_stream(args: &[String]) -> Result<(), String> {
    let mut inputs = Vec::new();
    let mut out_path: Option<String> = None;
    let mut min_ilis = 0.0_f32;
    let mut i = 0usize;
    while i < args.len() {
        match args[i].as_str() {
            "-o" | "--output" => {
                i += 1;
                out_path = args.get(i).cloned();
            }
            "--min-ilis" => {
                i += 1;
                min_ilis = args.get(i).and_then(|s| s.parse().ok()).unwrap_or(0.0);
            }
            other => inputs.push(other.to_string()),
        }
        i += 1;
    }
    let Some(out_file) = out_path else {
        return Err(
            "Usage: livia-atlas-cli dedup-stream <inputs...> -o <out.tsv> [--min-ilis 0.0]".into(),
        );
    };
    if inputs.is_empty() {
        return Err("dedup-stream requires at least one input file".into());
    }

    let t0 = Instant::now();
    let mut intern: HashMap<String, u32> = HashMap::new();
    let mut rev_intern: Vec<String> = Vec::new();
    let mut agg: HashMap<u64, AggPair> = HashMap::new();
    let mut total_rows = 0u64;

    for inp in &inputs {
        let low = inp.to_ascii_lowercase();
        if low.ends_with(".tar.xz") || low.ends_with(".txz") {
            let f = File::open(inp).map_err(|e| e.to_string())?;
            let mut tar_bytes = Vec::new();
            lzma_rs::xz_decompress(&mut BufReader::new(f), &mut tar_bytes)
                .map_err(|e| e.to_string())?;
            let mut ar = tar::Archive::new(Cursor::new(tar_bytes));
            for item in ar.entries().map_err(|e| e.to_string())? {
                let entry = item.map_err(|e| e.to_string())?;
                if !entry.header().entry_type().is_file() {
                    continue;
                }
                let p = entry
                    .path()
                    .map_err(|e| e.to_string())?
                    .to_string_lossy()
                    .to_ascii_lowercase();
                if p.ends_with(".csv") || p.ends_with(".tsv") {
                    total_rows += ingest_table_stream(
                        BufReader::new(entry),
                        &mut intern,
                        &mut rev_intern,
                        &mut agg,
                    )?;
                }
            }
        } else {
            let rdr = open_line_reader(inp)?;
            total_rows += ingest_table_stream(rdr, &mut intern, &mut rev_intern, &mut agg)?;
        }
    }

    let mut pairs: Vec<(u32, u32, AggPair)> = agg
        .into_iter()
        .filter(|&(_, v)| v.best_ilis >= min_ilis)
        .map(|(k, v)| ((k >> 32) as u32, (k & 0xffff_ffff) as u32, v))
        .collect();
    pairs.sort_unstable_by(|a, b| {
        rev_intern[a.0 as usize]
            .cmp(&rev_intern[b.0 as usize])
            .then_with(|| rev_intern[a.1 as usize].cmp(&rev_intern[b.1 as usize]))
    });

    let f = File::create(&out_file).map_err(|e| e.to_string())?;
    let mut w = BufWriter::new(f);
    writeln!(
        w,
        "protein_a\tprotein_b\tiLIS_best\tipTM_best\tipSAE_best\tLIpDockQ_best\tmodels"
    )
    .map_err(|e| e.to_string())?;
    let fmt_opt = |x: f32| -> String {
        if x.is_finite() {
            format!("{x:.4}")
        } else {
            String::new()
        }
    };
    for (ia, ib, v) in &pairs {
        let mut ka = &rev_intern[*ia as usize];
        let mut kb = &rev_intern[*ib as usize];
        if ka > kb {
            std::mem::swap(&mut ka, &mut kb);
        }
        writeln!(
            w,
            "{ka}\t{kb}\t{:.4}\t{}\t{}\t{}\t{}",
            v.best_ilis,
            fmt_opt(v.best_iptm),
            fmt_opt(v.best_ipsae),
            fmt_opt(v.best_lipdockq),
            v.models
        )
        .map_err(|e| e.to_string())?;
    }
    w.flush().map_err(|e| e.to_string())?;

    println!(
        "dedup-stream: {total_rows} input rows -> {} unique pairs ({} proteins) in {:.2}s -> {out_file}",
        pairs.len(),
        rev_intern.len(),
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

fn print_usage() {
    eprintln!(
        "Usage: livia-atlas-cli <subcommand> [args...]\n\n\
         Subcommands:\n\
           tested-index <species_dir> <pairs.tsv|csv|.gz|- ...>\n\
           paralog-index <species_dir> <cluster.tsv> <identity_note>\n\
           extract-set <dataset.zip> [set_id]\n\
           offsets-shards <site_repo>\n\
           compact-index <site_repo> [--check]\n\
           edges-shards <site_repo> [--check]\n\
           dedup-stream <inputs...> -o <out.tsv> [--min-ilis 0.0]"
    );
}

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();
    let Some((cmd, rest)) = args.split_first() else {
        print_usage();
        std::process::exit(1);
    };
    let res = match cmd.as_str() {
        "tested-index" => cmd_tested_index(rest),
        "paralog-index" => cmd_paralog_index(rest),
        "extract-set" => cmd_extract_set(rest),
        "offsets-shards" => cmd_offsets_shards(rest),
        "compact-index" => cmd_compact_index(rest),
        "edges-shards" => cmd_edges_shards(rest),
        "dedup-stream" => cmd_dedup_stream(rest),
        "-h" | "--help" => {
            print_usage();
            Ok(())
        }
        other => Err(format!("Unknown subcommand: {other}")),
    };
    if let Err(e) = res {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
}
