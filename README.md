# ChimeraX Visualization for LIS

Generate [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/) visualization scripts from protein-protein interaction predictions, with automatic [AFM-LIS](https://github.com/flyark/AFM-LIS) metric calculation.

<img src="images/upd1_dome_example.png" alt="Example: upd1-dome interaction visualized in ChimeraX" width="600">

> **Example:** upd1–dome interaction (iLIS: 0.529, ipTM: 0.57). Chain A (upd1) in blue, Chain B (dome) in orange. Light colors show LIR regions (confident interaction domain, PAE ≤ 12 Å). Dark colors highlight cLIR residues (direct physical contacts, PAE ≤ 12 Å & Cβ ≤ 8 Å). Non-interacting regions are hidden.

Two standalone web tools — **no installation required**. All analysis runs **locally in your browser** — no data is uploaded to any server, no Python needed:

**[Open Web Tools](https://flyark.github.io/chimerax_visualization_lis/)**

| Tool | Input |
|------|-------|
| **FlyPredictome → ChimeraX** | FlyPredictome URL |
| **AlphaFold3 → ChimeraX** | AF3 zip file |

## FlyPredictome → ChimeraX

Automatically fetch interaction data from [FlyPredictome](https://www.flyrnai.org/tools/fly_predictome) and generate ChimeraX scripts.

**Quick start:**
1. Paste a FlyPredictome result page URL
2. Choose rank and color scheme
3. Click **Generate**
4. Download `.cxc` script and `.pdb` structure file
5. Place both in the same folder, open the `.cxc` in ChimeraX

**Features:**
- Automatic data fetching via CORS proxy (parallel racing for speed)
- Pre-computed LIR/cLIR residue indices from FlyPredictome database
- Color presets: gradient (LIR light / cLIR dark), solid, ChimeraX defaults
- All 5 ranks shown with iLIS/ipTM scores

## AlphaFold3 → ChimeraX

Upload an AlphaFold3 prediction zip file to calculate AFM-LIS metrics and generate ChimeraX scripts. Supports multi-chain complexes (tested with 8-chain CCT complex, 4334 residues).

**Quick start:**
1. Download your prediction from [AlphaFold Server](https://alphafoldserver.com) as a zip
2. Drop the zip file on the page and click **Process**
3. Explore the generated maps and score matrices
4. Click a chain pair to generate a ChimeraX script
5. Download `.cxc` + `.cif`, place in the same folder, open the `.cxc` in ChimeraX

**Features:**
- Full AFM-LIS metric calculation locally in the browser (no data uploaded, no Python needed)
- **PAE Maps** — per-model Predicted Aligned Error (bwr colorscale)
- **Score Matrix** — combined iLIS (Oranges) / ipTM (Purples) per model
- **Residue Count Matrix** — combined LIR (Blues) / cLIR (Greens) per model
- **LIS Maps** — Local Interaction Score maps (Blues, dilated for large complexes)
- **cLIS Maps** — contact LIS maps (Greens, dilated for visibility)
- Multi-chain support with staggered axis labels
- Per-model and averaged metrics in the chain pair table
- Adjustable PAE cutoff and Cβ distance cutoff
- Color presets: gradient, solid, ChimeraX defaults, multi-chain (tab10 palette)
- **Download All Chains .cxc** — single script showing all chains with `color bychain`

## AFM-LIS Metrics

| Metric | Definition |
|--------|------------|
| **LIR** | Local Interaction Residues — residues in the confident interaction region (PAE ≤ 12 Å) |
| **cLIR** | contact-filtered LIR — residues in direct physical contact (PAE ≤ 12 Å & Cβ ≤ 8 Å) |
| **LIS** | Local Interaction Score — normalized PAE confidence across the interface (0–1) |
| **cLIS** | contact-filtered LIS — LIS restricted to direct physical contacts |
| **iLIS** | integrated LIS — geometric mean of LIS and cLIS: √(LIS × cLIS) |
| **ipTM** | interface predicted TM-score — AlphaFold-Multimer's global interface confidence |

See [AFM-LIS](https://github.com/flyark/AFM-LIS) for details.

## Color Scheme

Default blue/orange gradient:

| Color | Region |
|-------|--------|
| Light blue (`#b3d4e8`) | Chain A — LIR |
| Dark blue (`#2471A3`) | Chain A — cLIR |
| Light orange (`#f5cba7`) | Chain B — LIR |
| Dark orange (`#E67E22`) | Chain B — cLIR |

Non-LIR regions are hidden. Multiple preset schemes available including solid colors and ChimeraX defaults.

## Local Python Agents (Alternative)

For command-line usage or batch processing (requires Python + numpy + scipy):

```bash
# FlyPredictome
python flypredictome_agent.py "https://www.flyrnai.org/tools/fly_predictome/web/famdb_details/Egfr/spi/SET_69/" --rank 1

# AlphaFold3
python af3_agent.py /path/to/af3_prediction.zip --chain-pair A,B

# Flask web app (local server)
pip install flask
python app.py  # → http://127.0.0.1:5000
```

## References

**AFM-LIS**
- [Kim et al. 2024](https://www.biorxiv.org/content/10.1101/2024.02.19.580970) — Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer
- [Kim et al. 2025](https://www.biorxiv.org/content/10.1101/2025.10.10.681672) — A Structure-Guided Kinase–Transcription Factor Interactome Atlas

**AlphaFold**
- [Abramson et al. 2024](https://doi.org/10.1038/s41586-024-07487-w) — Accurate structure prediction of biomolecular interactions with AlphaFold 3. *Nature* 630, 493–500.

- [Evans et al. 2022](https://doi.org/10.1101/2021.10.04.463034) — Protein complex prediction with AlphaFold-Multimer. *bioRxiv*.

**ChimeraX**
- [Pettersen et al. 2021](https://doi.org/10.1002/pro.3943) — UCSF ChimeraX: Structure visualization for researchers, educators, and developers. *Protein Science* 30, 70–82.

## Acknowledgments

This tool was developed with assistance from Anthropic's Claude Code.

## License

MIT
