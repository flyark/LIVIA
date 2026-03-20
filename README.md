# *In silico* PPI Analysis & ChimeraX Script Generator

Analyze protein-protein interactions from structure predictions and generate [UCSF ChimeraX](https://www.cgl.ucsf.edu/chimerax/) visualization scripts, with automatic [AFM-LIS](https://github.com/flyark/AFM-LIS) metric calculation.

All analysis runs **locally in your browser** — no installation, no data uploaded, no Python needed.

**[Open Web Tools](https://flyark.github.io/in_silico_ppi_analysis/)**

| Tool | Input | Description |
|------|-------|-------------|
| **FlyPredictome Analysis** | FlyPredictome URL | Drosophila PPI predictions |
| **Universal Prediction Analysis** | AF2/AF3/ColabFold/Boltz/Chai-1/OpenFold files | Multi-platform structure prediction analysis |
| **AlphaFold DB Dimer Analysis** | UniProt ID or model ID | Pre-computed dimer predictions from AlphaFold Database |
| **AlphaFold DB Monomer Subdomain Analysis** | UniProt ID | Intramolecular domain interaction analysis |

## Universal Prediction Analysis

Upload prediction output from **any major platform** — the tool auto-detects the format:

- **AlphaFold3** — .cif + summary_confidences + full_data JSON
- **AlphaFold2** — ranked_*.pdb + PAE JSON
- **ColabFold** — *_unrelaxed_rank_*.pdb + *_scores_rank_*.json
- **Boltz-1/2** — .pdb/.cif + confidence_*.json + pae_*.npz
- **Chai-1** — pred.rank_*.cif + scores.rank_*.json + pae.rank_*.npy
- **OpenFold3** — result_sample_*_model.pdb + confidences JSON
- **Generic** — any .cif/.pdb with a PAE JSON

Accepts .zip, .gz, folders, or individual files. Handles .npz and .npy PAE formats.

## Features

- **PAE/LIS/cLIS maps** — bwr, matplotlib Blues, matplotlib Greens colormaps
- **Score matrix** — iLIS (Oranges) / ipTM (Purples) per chain pair
- **Residue count matrix** — LIR (Blues) / cLIR (Greens)
- **Sequence viewer** — amino acid letters with LIR/cLIR highlighting
- **Sortable tables** — click column headers to sort
- **ChimeraX presets** — gradient, solid, pLDDT coloring, color bychain, color bypolymer
- **Domain auto-detection** (monomer) — LIS/pLDDT-based segmentation with adjustable parameters
- **UniProt autocomplete** (monomer/dimer) — search suggestions as you type
- **CSV download** — full metrics with residue indices

## AFM-LIS Metrics

| Metric | Definition |
|--------|------------|
| **iLIS** | integrated LIS — √(LIS × cLIS) |
| **LIS** | Local Interaction Score — normalized PAE confidence (0–1) |
| **cLIS** | contact-filtered LIS — restricted to direct contacts |
| **LIR** | Local Interaction Residues (PAE ≤ 12 Å) |
| **cLIR** | contact-filtered LIR (PAE ≤ 12 Å & Cβ ≤ 8 Å) |
| **ipTM** | interface predicted TM-score |

Default cutoffs: PAE ≤ 12 Å, Cβ ≤ 8 Å (adjustable). PAE averaged in both directions per residue pair.

See [AFM-LIS](https://github.com/flyark/AFM-LIS) for details.

## References

- [Kim et al. 2024](https://www.biorxiv.org/content/10.1101/2024.02.19.580970) — Enhanced PPI Discovery via AlphaFold-Multimer
- [Kim et al. 2025](https://www.biorxiv.org/content/10.1101/2025.10.10.681672) — *In silico* Kinase-TF atlas for human and fly

## Acknowledgments

This tool was developed with assistance from Anthropic's Claude Code.

## License

MIT
