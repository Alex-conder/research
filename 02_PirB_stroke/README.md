# PirB/LILRB2 Post-Ischemic Neuroinflammation: Multi-omics Atlas

## Background

Paired immunoglobulin-like receptor B (PirB, mouse ortholog of human LILRB2) is an immune inhibitory receptor implicated in neuroinflammation and axon regeneration failure. This project integrates seven public GEO datasets to map the cell-type-specific, temporal, and spatial dynamics of PirB-expressing cells after experimental stroke.

## Methods

- **Data**: 7 GEO datasets (GSE174574, GSE171169, GSE225948, GSE233815, GSE233812, GSE233813, GSE233814)
- **Technologies**: scRNA-seq, snRNA-seq, bulk RNA-seq, Visium spatial transcriptomics
- **Tools**: Scanpy, Seurat, Harmony, ComBat, Squidpy, scVI-tools, gseapy, CellChat, LIANA
- **Languages**: Python, R

## Key Results

- **D3 microglial peak**: 26.6% Pirb+ microglia in GSE233812 scRNA-seq; 14.47% in cross-dataset integration.
- **A1-like astrocytes**: Pirb+ fraction 7.04% in MCAO vs 0.39% in sham (GSE174574).
- **Peripheral myeloid cells**: PB monocytes 51.72% and neutrophils 52.64% Pirb+ at acute D02 (GSE225948).
- **Spatial localization**: Pirb+ spots enriched at ischemic penumbra boundary (distance ratio Pirb+/Pirb− = 0.60, p = 3.3 × 10⁻¹⁷).

## Repository Contents

```
.
├── README.md
├── 03_analysis_code/         # Python/R analysis scripts
├── 04_reports/               # Reports, figures, and supplementary tables
└── README.md
```

## How to Run

1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
2. Start with `03_analysis_code/01_qc_pirb_GSE174574.py` for QC and annotation.
3. See `04_reports/` for generated reports including a Nature Communications-style manuscript.

## Requirements

See [`requirements.txt`](./requirements.txt).
