"""
FlyPredictome → ChimeraX Web App
Flask server that wraps the flypredictome_agent pipeline.
"""

from flask import Flask, render_template, request, jsonify, send_file
from pathlib import Path
import json
import zipfile
import io

import flypredictome_agent as agent

app = Flask(__name__, template_folder="templates", static_folder="static")
OUTPUT_DIR = Path("output")


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/generate", methods=["POST"])
def generate():
    """Run the agent pipeline and return results."""
    data = request.json
    url = data.get("url", "").strip()
    rank = int(data.get("rank", 1))
    colors = data.get("colors", {})

    if not url:
        return jsonify({"error": "URL is required"}), 400

    if "flyrnai.org" not in url or "famdb_details" not in url:
        return jsonify({"error": "Please enter a valid FlyPredictome famdb_details URL"}), 400

    try:
        # Step 1: Scrape
        rows, pdb_urls = agent.scrape_flypredictome(url)
        if not rows:
            return jsonify({"error": "No data found in the rank table"}), 404

        # Step 2: Select rank
        target_row = None
        for row in rows:
            try:
                if int(row.get("Rank", 0)) == rank:
                    target_row = row
                    break
            except (ValueError, TypeError):
                continue
        if target_row is None:
            target_row = rows[0]

        sym1 = target_row.get("Symbol 1", "unknown1")
        sym2 = target_row.get("Symbol 2", "unknown2")

        # Step 3: Download PDB
        case_dir = OUTPUT_DIR / f"{sym1}_{sym2}"
        case_dir.mkdir(parents=True, exist_ok=True)

        pdb_path = None
        if pdb_urls:
            pdb_url = pdb_urls[0]
            if not pdb_url.startswith("http"):
                pdb_url = f"https://www.flyrnai.org{pdb_url}"
            pdb_filename = pdb_url.split("/")[-1]
            pdb_local = case_dir / pdb_filename
            if pdb_local.exists():
                pdb_path = pdb_local
            else:
                import urllib.request
                req = urllib.request.Request(pdb_url, headers={
                    "User-Agent": "Mozilla/5.0 FlyPredictome-Agent/1.0"
                })
                with urllib.request.urlopen(req, timeout=60) as resp:
                    pdb_data = resp.read()
                pdb_local.write_bytes(pdb_data)
                pdb_path = pdb_local

        pdb_name = pdb_path.name if pdb_path else "STRUCTURE_NOT_FOUND.pdb"

        # Parse metrics
        try:
            iLIS_val = float(target_row.get("iLIS", 0))
        except (ValueError, TypeError):
            iLIS_val = None
        try:
            ipTM_val = float(target_row.get("ipTM", 0))
        except (ValueError, TypeError):
            ipTM_val = None
        try:
            plen_A = int(target_row.get("Len A", 0)) or None
        except (ValueError, TypeError):
            plen_A = None
        try:
            plen_B = int(target_row.get("Len B", 0)) or None
        except (ValueError, TypeError):
            plen_B = None

        # Step 4: Generate ChimeraX script with custom colors
        cxc_script = generate_cxc_with_colors(
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
            colors=colors,
        )

        cxc_path = case_dir / f"{sym1}_{sym2}_interface.cxc"
        cxc_path.write_text(cxc_script)

        # Build all-ranks summary
        all_ranks = []
        for row in rows:
            try:
                r = int(row.get("Rank", 0))
            except (ValueError, TypeError):
                r = 0
            all_ranks.append({
                "rank": r,
                "symbol_1": row.get("Symbol 1", ""),
                "symbol_2": row.get("Symbol 2", ""),
                "iLIS": row.get("iLIS", ""),
                "ipTM": row.get("ipTM", ""),
                "cLIR_A": len(row.get("cLIR_A_set", set())),
                "cLIR_B": len(row.get("cLIR_B_set", set())),
                "LIR_A": len(row.get("LIR_A_set", set())),
                "LIR_B": len(row.get("LIR_B_set", set())),
            })

        return jsonify({
            "success": True,
            "symbol_1": sym1,
            "symbol_2": sym2,
            "iLIS": iLIS_val,
            "ipTM": ipTM_val,
            "protein_len_A": plen_A,
            "protein_len_B": plen_B,
            "cLIR_A_count": len(target_row.get("cLIR_A_set", set())),
            "cLIR_B_count": len(target_row.get("cLIR_B_set", set())),
            "LIR_A_count": len(target_row.get("LIR_A_set", set())),
            "LIR_B_count": len(target_row.get("LIR_B_set", set())),
            "cxc_script": cxc_script,
            "cxc_filename": f"{sym1}_{sym2}_interface.cxc",
            "pdb_filename": pdb_name,
            "case_dir": str(case_dir),
            "all_ranks": all_ranks,
        })

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/download/<path:filename>")
def download_file(filename):
    """Download a generated file."""
    filepath = OUTPUT_DIR / filename
    if filepath.exists():
        return send_file(filepath, as_attachment=True)
    return jsonify({"error": "File not found"}), 404


@app.route("/api/download-zip/<sym1>/<sym2>")
def download_zip(sym1, sym2):
    """Download PDB + CXC as a zip."""
    case_dir = OUTPUT_DIR / f"{sym1}_{sym2}"
    if not case_dir.exists():
        return jsonify({"error": "Not found"}), 404

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w", zipfile.ZIP_DEFLATED) as zf:
        for f in case_dir.iterdir():
            if f.suffix in (".pdb", ".cxc"):
                zf.write(f, f.name)
    buf.seek(0)
    return send_file(buf, mimetype="application/zip",
                     as_attachment=True,
                     download_name=f"{sym1}_{sym2}_chimerax.zip")


def generate_cxc_with_colors(
    pdb_path, symbol_1, symbol_2,
    cLIR_A, cLIR_B, LIR_A, LIR_B,
    protein_len_A=None, protein_len_B=None,
    iLIS=None, ipTM=None,
    colors=None,
):
    """Generate CXC script with user-customizable colors."""
    if colors is None:
        colors = {}

    lir_a_color = colors.get("lir_a", "#b3d4e8")
    lir_b_color = colors.get("lir_b", "#f5cba7")
    clir_a_color = colors.get("clir_a", "#2471A3")
    clir_b_color = colors.get("clir_b", "#E67E22")

    iLIS_str = f"{iLIS:.3f}" if iLIS is not None else "N/A"
    ipTM_str = f"{ipTM:.2f}" if ipTM is not None else "N/A"
    lir_a_str = agent._ranges_str(LIR_A)
    lir_b_str = agent._ranges_str(LIR_B)

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
        f"color #1/A {lir_a_color}",
        f"color #1/B {lir_b_color}",
        f"hide atoms",
        f"hide cartoons",
        f"",
    ]

    # Title (3 parts: protein A in clir_a color, protein B in clir_b color, metrics in black)
    part_a = f"{symbol_1} ({lir_a_str})"
    part_b = f"{symbol_2} ({lir_b_str})"
    part_c = f"iLIS: {iLIS_str}  ipTM: {ipTM_str}"
    xpos_b = 0.03 + len(part_a) * 0.009 + 0.015
    xpos_c = xpos_b + len(part_b) * 0.009 + 0.015
    lines.append(f"# ── Title ──")
    lines.append(f'2dlabels create title text "{part_a}" xpos 0.03 ypos 0.95 size 16 color {clir_a_color} bold true')
    lines.append(f'2dlabels create title2 text "{part_b}" xpos {xpos_b:.3f} ypos 0.95 size 16 color {clir_b_color} bold true')
    lines.append(f'2dlabels create title3 text "{part_c}" xpos {xpos_c:.3f} ypos 0.95 size 16 color black bold true')
    lines.append("")

    # LIR regions
    LIR_A_filled = agent._fill_gaps(LIR_A, max_gap=20)
    LIR_B_filled = agent._fill_gaps(LIR_B, max_gap=20)
    lir_a_spec = agent._res_spec(LIR_A_filled, "A")
    lir_b_spec = agent._res_spec(LIR_B_filled, "B")

    lines.append(f"# ── Show ONLY LIR regions as cartoons ──")
    lines.append(f"# LIR_A: {symbol_1}, {len(LIR_A)} residues")
    if lir_a_spec:
        lines.append(f"show #1{lir_a_spec} cartoons")
    lines.append("")
    lines.append(f"# LIR_B: {symbol_2}, {len(LIR_B)} residues")
    if lir_b_spec:
        lines.append(f"show #1{lir_b_spec} cartoons")
    lines.append("")

    # cLIR highlights
    lines.append(f"# ── Highlight cLIR residues (contact interface) ──")
    if cLIR_A:
        clir_a_spec = agent._res_spec(cLIR_A, "A")
        lines.append(f"# cLIR_A: {symbol_1}, {len(cLIR_A)} residues")
        lines.append(f"color #1{clir_a_spec} {clir_a_color}")
        lines.append("")
    if cLIR_B:
        clir_b_spec = agent._res_spec(cLIR_B, "B")
        lines.append(f"# cLIR_B: {symbol_2}, {len(cLIR_B)} residues")
        lines.append(f"color #1{clir_b_spec} {clir_b_color}")
        lines.append("")

    # Final
    lines.append("# ── Final view ──")
    lines.append("view")
    lines.append("lighting soft depthCue true")
    lines.append("")
    # Legend
    lines.append("# ── Color Legend ──")
    lines.append(f"# LIR  = Local Interaction Residues (PAE <= 12)")
    lines.append(f"# cLIR = contact residues (PAE <= 12 & C-beta distance <= 8 Å)")
    lines.append(f"#")
    lines.append(f"# LIR cartoons  ({lir_a_color}) = {symbol_1} LIR (chain A)")
    lines.append(f"# cLIR cartoons ({clir_a_color}) = {symbol_1} cLIR (chain A)")
    lines.append(f"# LIR cartoons  ({lir_b_color}) = {symbol_2} LIR (chain B)")
    lines.append(f"# cLIR cartoons ({clir_b_color}) = {symbol_2} cLIR (chain B)")
    lines.append(f"# Non-LIR regions are HIDDEN")
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


if __name__ == "__main__":
    OUTPUT_DIR.mkdir(exist_ok=True)
    app.run(debug=True, port=5000)
