"""
FlyPredictome → ChimeraX Agent
================================
Given a FlyPredictome famdb_details URL, this agent:
  1. Scrapes the rank table to extract cLIR/LIR residue indices and metadata
  2. Downloads the PDB structure file
  3. Generates a ChimeraX .cxc command script for visualization

Usage:
    python flypredictome_agent.py <URL> [--rank N] [--output-dir DIR]

Example:
    python flypredictome_agent.py \
        https://www.flyrnai.org/tools/fly_predictome/web/famdb_details/Egfr/spi/SET_69/ \
        --rank 1 --output-dir ./output
"""

import argparse
import json
import re
import sys
import urllib.request
import urllib.error
from html.parser import HTMLParser
from pathlib import Path


# ── HTML Table Parser ──────────────────────────────────────────────────────

class FlyPredictomeTableParser(HTMLParser):
    """Parse the rank table from a FlyPredictome famdb_details page."""

    # Column headers in expected order (0-indexed)
    EXPECTED_HEADERS = [
        "id", "Rank", "Protein 1", "Protein 2", "fbgn1", "fbgn2",
        "ncbi_gene_id1", "ncbi_gene_id2", "Symbol 1", "Symbol 2",
        "protein_1_size", "protein_2_size", "iLIS",
        "LIS (AB)", "LIS (BA)", "LIS",
        "LIA (AB)", "LIA (BA)", "LIA",
        "LIR (A)", "LIR (B)", "LIR",
        "LIpLDDT (A)", "LIpLDDT (B)", "LIpLDDT",
        "cLIS (AB)", "cLIS (BA)", "cLIS",
        "cLIA (AB)", "cLIA (BA)", "cLIA",
        "cLIR (A)", "cLIR (B)", "cLIR",
        "cLIpLDDT (A)", "cLIpLDDT (B)", "cLIpLDDT",
        "LIR Indice A", "LIR Indice B",
        "cLIR Indice A", "cLIR Indice B",
        "ipTM", "Confidence",
        "pae_file_name", "directory_name", "output_file_name",
        "enrichment_id1", "enrichment_id2",
        "pLDDT", "pTM", "Len A", "Len B",
    ]

    def __init__(self):
        super().__init__()
        self.in_table = False
        self.in_thead = False
        self.in_tbody = False
        self.in_th = False
        self.in_td = False
        self.in_tr = False
        self.headers = []
        self.rows = []
        self.current_row = []
        self.current_cell = ""
        self.table_id = None

    def handle_starttag(self, tag, attrs):
        attrs_dict = dict(attrs)
        if tag == "table" and attrs_dict.get("id") == "rank_table":
            self.in_table = True
            self.table_id = "rank_table"
        elif self.in_table:
            if tag == "thead":
                self.in_thead = True
            elif tag == "tbody":
                self.in_tbody = True
            elif tag == "tr":
                self.in_tr = True
                self.current_row = []
            elif tag == "th" and self.in_thead:
                self.in_th = True
                self.current_cell = ""
            elif tag == "td" and self.in_tbody:
                self.in_td = True
                self.current_cell = ""

    def handle_endtag(self, tag):
        if tag == "table" and self.in_table:
            self.in_table = False
        elif self.in_table:
            if tag == "thead":
                self.in_thead = False
            elif tag == "tbody":
                self.in_tbody = False
            elif tag == "tr":
                self.in_tr = False
                if self.in_tbody and self.current_row:
                    self.rows.append(self.current_row)
            elif tag == "th" and self.in_th:
                self.in_th = False
                self.headers.append(self.current_cell.strip())
            elif tag == "td" and self.in_td:
                self.in_td = False
                self.current_row.append(self.current_cell.strip())

    def handle_data(self, data):
        if self.in_th:
            self.current_cell += data
        elif self.in_td:
            self.current_cell += data

    def get_rows_as_dicts(self):
        """Return rows as list of dicts keyed by header name."""
        result = []
        for row in self.rows:
            d = {}
            for i, val in enumerate(row):
                if i < len(self.headers):
                    d[self.headers[i]] = val
            result.append(d)
        return result


# ── Residue Index Parsing ──────────────────────────────────────────────────

def parse_residue_ranges(raw: str) -> set:
    """Parse a JSON-encoded list of residue ranges into a set of integers.

    Input format: '["107-114","117","120-121"]'
    Returns: {107, 108, 109, 110, 111, 112, 113, 114, 117, 120, 121}
    """
    if not raw or raw in ("[]", "null", "None", ""):
        return set()

    # Try JSON parse first
    try:
        items = json.loads(raw)
    except (json.JSONDecodeError, TypeError):
        # Fallback: extract numbers/ranges from string
        items = re.findall(r'[\d]+-[\d]+|\d+', raw)

    positions = set()
    for item in items:
        item = str(item).strip()
        if "-" in item:
            parts = item.split("-")
            try:
                start, end = int(parts[0]), int(parts[1])
                positions.update(range(start, end + 1))
            except (ValueError, IndexError):
                continue
        else:
            try:
                positions.add(int(item))
            except ValueError:
                continue
    return positions


# ── Web Scraping ───────────────────────────────────────────────────────────

def fetch_page(url: str) -> str:
    """Fetch HTML content from a URL."""
    req = urllib.request.Request(url, headers={
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X 10_15_7) FlyPredictome-Agent/1.0"
    })
    with urllib.request.urlopen(req, timeout=30) as resp:
        return resp.read().decode("utf-8")


def extract_pdb_urls(html: str) -> list:
    """Extract PDB download URLs from page HTML/JavaScript."""
    # Pattern 1: loadStructureFromUrl or href containing .pdb
    patterns = [
        r'loadStructureFromUrl\(\s*["\']([^"\']+\.pdb)["\']',
        r'href=["\']([^"\']*\.pdb)["\']',
        r'["\'](/tools/fly_predictome/web/colabfold-output/[^"\']+\.pdb)["\']',
    ]
    urls = []
    seen = set()
    for pattern in patterns:
        for match in re.findall(pattern, html):
            if match not in seen:
                seen.add(match)
                urls.append(match)
    return urls


def scrape_flypredictome(url: str) -> tuple:
    """Scrape the rank table from a FlyPredictome famdb_details page.

    Returns (rows, pdb_urls) where rows is a list of dicts with parsed residue indices,
    and pdb_urls is a list of PDB download paths found on the page.
    """
    print(f"[Agent] Fetching: {url}")
    html = fetch_page(url)

    parser = FlyPredictomeTableParser()
    parser.feed(html)

    # Extract PDB URLs from page
    pdb_urls = extract_pdb_urls(html)
    if pdb_urls:
        print(f"[Agent] Found {len(pdb_urls)} PDB URL(s) on page")
        for u in pdb_urls:
            print(f"  {u}")

    if not parser.rows:
        print("[Agent] HTML parser found no rows — trying regex fallback...")
        return _regex_fallback(html), pdb_urls

    rows = parser.get_rows_as_dicts()
    print(f"[Agent] Parsed {len(rows)} rows from rank table")
    print(f"[Agent] Headers found: {parser.headers}")

    # Enrich each row with parsed residue sets
    for row in rows:
        row["cLIR_A_set"] = parse_residue_ranges(row.get("cLIR Indice A", ""))
        row["cLIR_B_set"] = parse_residue_ranges(row.get("cLIR Indice B", ""))
        row["LIR_A_set"] = parse_residue_ranges(row.get("LIR Indice A", ""))
        row["LIR_B_set"] = parse_residue_ranges(row.get("LIR Indice B", ""))

    return rows, pdb_urls


def _regex_fallback(html: str) -> list:
    """Fallback extraction using regex when HTML parser fails (e.g., JS-rendered tables)."""
    # Try to find table data in embedded JavaScript or inline JSON
    # Look for DataTable initialization data
    rows = []

    # Pattern: look for arrays that contain cLIR-like data
    # This is a best-effort fallback
    lir_pattern = re.findall(
        r'"LIR Indice [AB]"\s*:\s*"(\[.*?\])"', html, re.DOTALL
    )
    if lir_pattern:
        print(f"[Agent] Regex fallback found {len(lir_pattern)} LIR Indice matches")

    return rows


# ── PDB Download ───────────────────────────────────────────────────────────

def build_pdb_url(row: dict, base_url: str = "https://www.flyrnai.org") -> str:
    """Build the PDB download URL from table row data."""
    directory = row.get("directory_name", "")
    filename = row.get("output_file_name", "")

    if directory and filename:
        return f"{base_url}/tools/fly_predictome/web/colabfold-output/{directory}/{filename}"

    # Fallback: construct from gene names
    sym1 = row.get("Symbol 1", "")
    sym2 = row.get("Symbol 2", "")
    if sym1 and sym2:
        return f"{base_url}/tools/fly_predictome/web/colabfold-output/*/{sym1}___{sym2}_unrelaxed_rank_001_*.pdb"

    return ""


def download_pdb(row: dict, output_dir: Path) -> Path:
    """Download the PDB file for a given row."""
    url = build_pdb_url(row)
    if not url:
        print("[Agent] ERROR: Cannot determine PDB URL")
        return None

    filename = row.get("output_file_name", "")
    if not filename:
        sym1 = row.get("Symbol 1", "")
        sym2 = row.get("Symbol 2", "")
        filename = f"{sym1}___{sym2}_rank{row.get('Rank', '1')}.pdb"

    output_path = output_dir / filename
    if output_path.exists():
        print(f"[Agent] PDB already exists: {output_path}")
        return output_path

    print(f"[Agent] Downloading PDB: {url}")
    try:
        req = urllib.request.Request(url, headers={
            "User-Agent": "Mozilla/5.0 FlyPredictome-Agent/1.0"
        })
        with urllib.request.urlopen(req, timeout=60) as resp:
            data = resp.read()
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_bytes(data)
        print(f"[Agent] PDB saved: {output_path} ({len(data)} bytes)")
        return output_path
    except urllib.error.URLError as e:
        print(f"[Agent] ERROR downloading PDB: {e}")
        return None


# ── ChimeraX Script Generation ────────────────────────────────────────────

def _fill_gaps(positions: set, max_gap: int = 10) -> set:
    """Fill small gaps between segments for continuous cartoon display."""
    if not positions:
        return positions
    sorted_pos = sorted(positions)
    filled = set(sorted_pos)
    for i in range(len(sorted_pos) - 1):
        gap = sorted_pos[i + 1] - sorted_pos[i] - 1
        if 0 < gap <= max_gap:
            for j in range(sorted_pos[i] + 1, sorted_pos[i + 1]):
                filled.add(j)
    return filled


def _res_spec(positions: set, chain: str) -> str:
    """Build ChimeraX residue specifier: /A:384-486,490-500"""
    if not positions:
        return ""
    sorted_pos = sorted(positions)
    ranges = []
    start = end = sorted_pos[0]
    for p in sorted_pos[1:]:
        if p == end + 1:
            end = p
        else:
            ranges.append(f"{start}-{end}" if start != end else str(start))
            start = end = p
    ranges.append(f"{start}-{end}" if start != end else str(start))
    return f"/{chain}:" + ",".join(ranges)


def _ranges_str(positions: set) -> str:
    """Compact ranges string for display: '384-486, 490-575'."""
    if not positions:
        return "none"
    sorted_pos = sorted(positions)
    ranges = []
    start = end = sorted_pos[0]
    for p in sorted_pos[1:]:
        if p <= end + 15:
            end = p
        else:
            ranges.append(f"{start}-{end}")
            start = end = p
    ranges.append(f"{start}-{end}")
    return ", ".join(ranges)


def generate_chimerax_script(
    pdb_path: str,
    symbol_1: str,
    symbol_2: str,
    cLIR_A: set,
    cLIR_B: set,
    LIR_A: set,
    LIR_B: set,
    protein_len_A: int = None,
    protein_len_B: int = None,
    iLIS: float = None,
    ipTM: float = None,
) -> str:
    """Generate a ChimeraX .cxc script for visualizing a PPI prediction.

    Color scheme (blue/orange gradient):
      - Chain A LIR: light blue (#b3d4e8) → cLIR: dark blue (#2471A3)
      - Chain B LIR: light orange (#f5cba7) → cLIR: dark orange (#E67E22)
      - Non-LIR regions are hidden
    """
    iLIS_str = f"{iLIS:.3f}" if iLIS is not None else "N/A"
    ipTM_str = f"{ipTM:.2f}" if ipTM is not None else "N/A"
    lir_a_str = _ranges_str(LIR_A)
    lir_b_str = _ranges_str(LIR_B)

    lines = [
        f"# ChimeraX visualization: {symbol_1} — {symbol_2} (iLIS={iLIS_str}, ipTM={ipTM_str})",
        f"# Chain /A: {symbol_1} ({protein_len_A or '?'} residues, LIR: {lir_a_str})",
        f"# Chain /B: {symbol_2} ({protein_len_B or '?'} residues, LIR: {lir_b_str})",
        f"# Generated by FlyPredictome-ChimeraX Agent",
        f"",
        f"close",
        f"open {pdb_path}",
        f"",
        f"# ── Setup: white background, hide everything first ──",
        f"set bgColor white",
        f"graphics silhouettes true",
        f"lighting soft",
        f"color #1/A #b3d4e8",   # light blue for chain A (LIR base)
        f"color #1/B #f5cba7",   # light orange for chain B (LIR base)
        f"hide atoms",
        f"hide cartoons",
        f"",
    ]

    # ── Title labels (3 parts: protein A in blue, protein B in orange, metrics in black) ──
    part_a = f"{symbol_1} ({lir_a_str})"
    part_b = f"{symbol_2} ({lir_b_str})"
    part_c = f"iLIS: {iLIS_str}  ipTM: {ipTM_str}"
    xpos_b = 0.03 + len(part_a) * 0.009 + 0.015
    xpos_c = xpos_b + len(part_b) * 0.009 + 0.015
    lines.append(f"# ── Title ──")
    lines.append(f'2dlabels create title text "{part_a}" xpos 0.03 ypos 0.95 size 16 color #2471A3 bold true')
    lines.append(f'2dlabels create title2 text "{part_b}" xpos {xpos_b:.3f} ypos 0.95 size 16 color #E67E22 bold true')
    lines.append(f'2dlabels create title3 text "{part_c}" xpos {xpos_c:.3f} ypos 0.95 size 16 color black bold true')
    lines.append("")

    # ── Show ONLY LIR regions as cartoons ──
    LIR_A_filled = _fill_gaps(LIR_A, max_gap=20)
    LIR_B_filled = _fill_gaps(LIR_B, max_gap=20)
    lir_a_spec = _res_spec(LIR_A_filled, "A")
    lir_b_spec = _res_spec(LIR_B_filled, "B")

    lines.append(f"# ── Show ONLY LIR regions as cartoons ──")
    lines.append(f"# LIR_A: {symbol_1}, {len(LIR_A)} residues — light blue")
    if lir_a_spec:
        lines.append(f"show #1{lir_a_spec} cartoons")
    lines.append("")

    lines.append(f"# LIR_B: {symbol_2}, {len(LIR_B)} residues — light orange")
    if lir_b_spec:
        lines.append(f"show #1{lir_b_spec} cartoons")
    lines.append("")

    # ── Highlight cLIR residues — darker colors on cartoons ──
    lines.append(f"# ── Highlight cLIR residues (contact interface) — darker colors ──")
    if cLIR_A:
        clir_a_spec = _res_spec(cLIR_A, "A")
        lines.append(f"# cLIR_A: {symbol_1}, {len(cLIR_A)} residues — dark blue")
        lines.append(f"color #1{clir_a_spec} #2471A3")
        lines.append("")

    if cLIR_B:
        clir_b_spec = _res_spec(cLIR_B, "B")
        lines.append(f"# cLIR_B: {symbol_2}, {len(cLIR_B)} residues — dark orange")
        lines.append(f"color #1{clir_b_spec} #E67E22")
        lines.append("")

    # ── Final view ──
    lines.append("# ── Final view ──")
    lines.append("view")
    lines.append("lighting soft depthCue true")
    lines.append("")
    # ── Save figure ──
    png_name = f"{symbol_1}_{symbol_2}.png"
    lines.append(f"# ── Save figure ──")
    lines.append(f"save {png_name} transparentBackground true")
    lines.append("")

    # ── Color legend & metadata ──
    lines.append("# ── Color Legend ──")
    lines.append(f"# LIR  = Local Interaction Region (PAE <= 12)")
    lines.append(f"# cLIR = contact LIR (PAE <= 12 & C-beta distance <= 8 Å)")
    lines.append(f"#")
    lines.append(f"# Light blue cartoons (#b3d4e8) = {symbol_1} LIR (chain A)")
    lines.append(f"# Dark blue cartoons  (#2471A3) = {symbol_1} cLIR (chain A)")
    lines.append(f"# Light orange cartoons (#f5cba7) = {symbol_2} LIR (chain B)")
    lines.append(f"# Dark orange cartoons  (#E67E22) = {symbol_2} cLIR (chain B)")
    lines.append(f"# Non-LIR regions are HIDDEN")
    lines.append(f"#")
    if LIR_A:
        s = sorted(LIR_A)
        lines.append(f"# LIR_A: {symbol_1} {s[0]}-{s[-1]} ({len(LIR_A)} residues)")
    if LIR_B:
        s = sorted(LIR_B)
        lines.append(f"# LIR_B: {symbol_2} {s[0]}-{s[-1]} ({len(LIR_B)} residues)")
    if cLIR_A:
        lines.append(f"# cLIR_A: {symbol_1} {len(cLIR_A)} residues")
    if cLIR_B:
        lines.append(f"# cLIR_B: {symbol_2} {len(cLIR_B)} residues")

    return "\n".join(lines)


# ── Main Pipeline ──────────────────────────────────────────────────────────

def run(url: str, rank: int = 1, output_dir: str = "./output"):
    """Run the full pipeline: scrape → download PDB → generate ChimeraX script."""
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Step 1: Scrape the page
    rows, pdb_urls = scrape_flypredictome(url)
    if not rows:
        print("[Agent] ERROR: No data found in the rank table")
        sys.exit(1)

    # Step 2: Select the requested rank
    target_row = None
    for row in rows:
        try:
            if int(row.get("Rank", 0)) == rank:
                target_row = row
                break
        except (ValueError, TypeError):
            continue

    if target_row is None:
        print(f"[Agent] Rank {rank} not found, using first row")
        target_row = rows[0]

    sym1 = target_row.get("Symbol 1", "unknown1")
    sym2 = target_row.get("Symbol 2", "unknown2")
    print(f"\n[Agent] Selected: Rank {rank} — {sym1} / {sym2}")
    print(f"  iLIS: {target_row.get('iLIS', 'N/A')}")
    print(f"  ipTM: {target_row.get('ipTM', 'N/A')}")
    print(f"  cLIR_A: {len(target_row.get('cLIR_A_set', set()))} residues")
    print(f"  cLIR_B: {len(target_row.get('cLIR_B_set', set()))} residues")
    print(f"  LIR_A: {len(target_row.get('LIR_A_set', set()))} residues")
    print(f"  LIR_B: {len(target_row.get('LIR_B_set', set()))} residues")

    # Step 3: Download PDB (use URL extracted from page HTML)
    case_dir = output_path / f"{sym1}_{sym2}"
    case_dir.mkdir(parents=True, exist_ok=True)

    pdb_path = None
    if pdb_urls:
        # Use the first PDB URL found on the page (rank 1 structure)
        pdb_url = pdb_urls[0]
        if not pdb_url.startswith("http"):
            pdb_url = f"https://www.flyrnai.org{pdb_url}"
        pdb_filename = pdb_url.split("/")[-1]
        pdb_local = case_dir / pdb_filename
        if pdb_local.exists():
            print(f"[Agent] PDB already exists: {pdb_local}")
            pdb_path = pdb_local
        else:
            print(f"[Agent] Downloading PDB: {pdb_url}")
            try:
                req = urllib.request.Request(pdb_url, headers={
                    "User-Agent": "Mozilla/5.0 FlyPredictome-Agent/1.0"
                })
                with urllib.request.urlopen(req, timeout=60) as resp:
                    data = resp.read()
                pdb_local.write_bytes(data)
                print(f"[Agent] PDB saved: {pdb_local} ({len(data)} bytes)")
                pdb_path = pdb_local
            except urllib.error.URLError as e:
                print(f"[Agent] ERROR downloading PDB: {e}")

    if pdb_path is None:
        # Fallback to table-based URL construction
        pdb_path = download_pdb(target_row, case_dir)

    pdb_name = pdb_path.name if pdb_path else "STRUCTURE_NOT_FOUND.pdb"

    # Step 4: Generate ChimeraX script
    try:
        iLIS_val = float(target_row.get("iLIS", 0))
    except (ValueError, TypeError):
        iLIS_val = None
    try:
        ipTM_val = float(target_row.get("ipTM", 0))
    except (ValueError, TypeError):
        ipTM_val = None
    # protein_1_size may be "full" — use Len A/B columns instead
    try:
        plen_A = int(target_row.get("Len A", 0)) or None
    except (ValueError, TypeError):
        try:
            plen_A = int(target_row.get("protein_1_size", 0)) or None
        except (ValueError, TypeError):
            plen_A = None
    try:
        plen_B = int(target_row.get("Len B", 0)) or None
    except (ValueError, TypeError):
        try:
            plen_B = int(target_row.get("protein_2_size", 0)) or None
        except (ValueError, TypeError):
            plen_B = None

    script = generate_chimerax_script(
        pdb_path=pdb_name,
        symbol_1=sym1,
        symbol_2=sym2,
        cLIR_A=target_row.get("cLIR_A_set", set()),
        cLIR_B=target_row.get("cLIR_B_set", set()),
        LIR_A=target_row.get("LIR_A_set", set()),
        LIR_B=target_row.get("LIR_B_set", set()),
        protein_len_A=plen_A,
        protein_len_B=plen_B,
        iLIS=iLIS_val,
        ipTM=ipTM_val,
    )

    cxc_path = case_dir / f"{sym1}_{sym2}_interface.cxc"
    cxc_path.write_text(script)
    print(f"\n[Agent] ChimeraX script saved: {cxc_path}")
    print(f"[Agent] Open in ChimeraX: open {cxc_path}")

    return {
        "symbol_1": sym1,
        "symbol_2": sym2,
        "pdb_path": str(pdb_path) if pdb_path else None,
        "cxc_path": str(cxc_path),
        "cLIR_A": sorted(target_row.get("cLIR_A_set", set())),
        "cLIR_B": sorted(target_row.get("cLIR_B_set", set())),
        "LIR_A_count": len(target_row.get("LIR_A_set", set())),
        "LIR_B_count": len(target_row.get("LIR_B_set", set())),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="FlyPredictome → ChimeraX Agent")
    parser.add_argument("url", help="FlyPredictome famdb_details URL")
    parser.add_argument("--rank", type=int, default=1, help="Which rank to visualize (default: 1)")
    parser.add_argument("--output-dir", default="./output", help="Output directory (default: ./output)")
    args = parser.parse_args()

    result = run(args.url, rank=args.rank, output_dir=args.output_dir)
    print(f"\n[Agent] Done! Result: {json.dumps(result, indent=2)}")
