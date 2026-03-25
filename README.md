# LIVIA — **L**ocal **I**nteraction **VI**sualization and **A**nalysis

Browser-based tools for analyzing protein-protein interactions from structure predictions. All analysis runs locally in your browser — no data leaves your device and no installation is needed.

**https://flyark.github.io/LIVIA/**

## Features

- **Interaction residue detection** — LIR (Local Interaction Region, PAE ≤ 12Å) and cLIR (contact LIR, PAE ≤ 12Å & Cβ ≤ 8Å)
- **iLIS / iLIA / iLISA calculation** — interface confidence scoring ([Kim et al. 2024](https://doi.org/10.1101/2024.02.19.580970), [Kim et al. 2025](https://doi.org/10.1101/2025.10.10.681672))
- **3D structure viewer** — Mol* viewer with LIR/cLIR coloring
- **ChimeraX and PyMOL script generation** — multiple color modes (gradient, solid, high contrast, pLDDT, bychain)
- **Interactive visualizations** — sequence viewer, PAE/LIS/cLIS maps, interface contact maps, chord diagrams
- **Color presets** — gradient, solid, high contrast, and custom color options
- **Show complete structure / Gray non-LIR** — toggle full structure display with optional gray for non-interacting residues

## Pages

| Page | Description |
|------|-------------|
| **[Prediction Analysis](https://flyark.github.io/LIVIA/universal.html)** | General-purpose analysis for AlphaFold2, AlphaFold3, ColabFold, Boltz, Chai-1, and OpenFold |
| **[FlyPredictome](https://flyark.github.io/LIVIA/flypredictome.html)** | *Drosophila* protein-protein interaction analysis from [FlyPredictome](https://www.flyrnai.org/tools/fly_predictome) |
| **[Ortholog Predictome](https://flyark.github.io/LIVIA/ortholog_predictome.html)** | Non-fly predictions from [FlyPredictome](https://www.flyrnai.org/tools/fly_predictome) (including >240,000 human S/T kinase–TF predictions) |
| **[AlphaFold DB Dimer](https://flyark.github.io/LIVIA/dimer.html)** | Two-chain complex analysis from [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk) |
| **[Monomer Subdomain](https://flyark.github.io/LIVIA/monomer.html)** | Single-chain subdomain interaction analysis from [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk) |
| **[Guide](https://flyark.github.io/LIVIA/guide.html)** | Step-by-step instructions, metric definitions, and references |

## Supported Platforms

| Platform | Structure | PAE | Auto-detection |
|----------|-----------|-----|----------------|
| AlphaFold3 | `.cif` | `full_data_*.json` | Yes |
| AlphaFold2 | `.pdb` | PAE JSON | Yes |
| ColabFold | `.pdb` | `scores_*.json` | Yes |
| Boltz-1/2 | `.pdb` / `.cif` | `pae_*.npz` | Yes |
| Chai-1 | `.cif` | `scores.*.json` | Yes |
| OpenFold3 | `.pdb` | confidences JSON | Yes |

## Key Metrics

- **LIS** (Local Interaction Score) — average confidence of residue pairs with PAE ≤ 12Å ([Kim et al. 2024](https://doi.org/10.1101/2024.02.19.580970))
- **cLIS** (contact LIS) — same as LIS but restricted to pairs with Cβ–Cβ ≤ 8Å
- **iLIS** — `sqrt(LIS × cLIS)` — single integrated score ([Kim et al. 2025](https://doi.org/10.1101/2025.10.10.681672))
- **iLIA** — `sqrt(LIA × cLIA)` — integrated interaction area
- **iLISA** — `iLIS × iLIA` — combined score and area
- **ipSAE** — TM-score-like interface metric with adaptive normalization ([Dunbrack 2025](https://doi.org/10.1101/2025.01.17.633614))
- **LIR** — Local Interaction Region (residues with PAE ≤ 12Å to other chain)
- **cLIR** — contact LIR (LIR residues with Cβ–Cβ ≤ 8Å)

## Note

- Tested on Chrome and Safari (macOS/iOS).
- Each prediction platform may produce different confidence calibrations. The iLIS ≥ 0.223 threshold was established using ColabFold/AlphaFold-Multimer predictions. Other platforms may require adjusted thresholds.

## Related Resources

- **[AFM-LIS](https://github.com/flyark/AFM-LIS)** — Python framework for iLIS/LIS calculation
- **[FlyPredictome](https://www.flyrnai.org/tools/fly_predictome)** — Large-scale *Drosophila* PPI prediction database (>1.7 million predictions)
- **[AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk)** — Predicted protein structures

## References

- Kim, A.-R. et al. (2025). A Structure-Guided Kinase–Transcription Factor Interactome Atlas Reveals Docking Landscapes of the Kinome. *bioRxiv*. https://doi.org/10.1101/2025.10.10.681672
- Kim, A.-R. et al. (2024). Enhanced Protein-Protein Interaction Discovery via AlphaFold-Multimer. *bioRxiv*. https://doi.org/10.1101/2024.02.19.580970

## Citation

If you use LIVIA in your research, please cite:

```bibtex
@misc{livia,
  author = {Kim, Ah-Ram},
  title = {LIVIA: Local Interaction Visualization and Analysis},
  year = {2026},
  publisher = {GitHub},
  url = {https://github.com/flyark/LIVIA}
}
```

## License

MIT
