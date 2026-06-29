# Huang Teng — Academic Research Portfolio

> This repository contains the core code, reports, and documentation for my undergraduate research projects in computational biology, computational chemistry, and multi-omics analysis.
>
> **Contact**: alexanderh_t@163.com  
> **Affiliation**: School of Pharmacy, Xi'an Medical University (Expected graduation: July 2026)

## Projects Overview

| # | Project | Domain | Key Methods | Key Results |
|---|---------|--------|-------------|-------------|
| 01 | [6UJN Sequential Docking](./01_6UJN_docking) | Computational Chemistry / Drug Design | AutoDock Vina, Amber/OpenMM MD, MM-GBSA | Meropenem-first Vancomycin ΔG = −52.85 ± 4.20 kcal/mol; ligand-ligand COM 7.38 Å |
| 02 | [PirB Stroke Multi-omics](./02_PirB_stroke) | Single-cell & Spatial Transcriptomics | Scanpy, Seurat, Squidpy, scVI, CellChat | Pirb peaks at D3 in microglia (26.6%); enriched at ischemic penumbra (p = 3.3×10⁻¹⁷) |
| 03 | [Virtual Cell Modeling](./03_virtual_cell) | Multi-scale Modeling | HH neuron model, FBA (COBRApy), scVI | HH↔FBA interface R² = 0.9996; AD neuroinflammation-cerebrovascular cascade |
| 04 | [LC-MS Metabolomics Pipeline](./04_metabolomics) | Non-targeted Metabolomics | Python, OPLS-DA, PubChem API, KEGG ORA | 116 significant metabolite features; 12 pathways at FDR < 0.25 |
| 05 | [LLM Automation Workflow](./05_llm_automation) | Research Productivity | Python, LLM APIs, modular pipelines | Reusable modules for technical chains and boundary definition |

## Skills

- **Programming**: Python (pandas, numpy, scipy, scanpy, scvi-tools, cobrapy), R, Bash
- **Computational Chemistry**: AutoDock Vina, AmberTools, OpenMM, MMPBSA.py, PyMOL
- **Bioinformatics**: scRNA-seq/snRNA-seq analysis, Visium spatial transcriptomics, bulk RNA-seq, pathway enrichment
- **Modeling**: ODE-based electrophysiology, flux balance analysis (FBA), variational autoencoders (scVI)
- **AI-assisted R&D**: LLM API integration for code generation, automated report drafting, and pipeline modularization

## Repository Structure

```
.
├── 01_6UJN_docking/          # P-glycoprotein sequential docking & MD/MM-GBSA
├── 02_PirB_stroke/           # Brain ischemia PirB multi-omics atlas
├── 03_virtual_cell/          # HH + FBA + scVI multi-scale virtual cell
├── 04_metabolomics/          # LC-MS non-targeted metabolomics pipeline
├── 05_llm_automation/        # LLM-assisted research workflow modules
└── README.md
```

## Notes

- This is a curated portfolio. Raw sequencing/mass-spec data files (e.g., `.h5ad`, `.mzML`, FASTQ) are excluded due to size.
- Each project directory contains its own README with setup instructions and key outputs.
- For the complete AI Agent OS project (Kaelis), please see the private repository [Kaelis-main](https://github.com/Alex-conder/Kaelis-main) (access available upon request).

## License

Code in this repository is provided for academic review and demonstration purposes.
