# Plant LC-MS Non-targeted Metabolomics Pipeline

## Background

A command-line Python pipeline for plant/cell LC-MS non-targeted metabolomics starting from an SPSS `.sav` abundance matrix. The pipeline performs QC, normalization, univariate and multivariate differential analysis (OPLS-DA), metabolite annotation, and pathway enrichment.

## Methods

- **Input**: SPSS `.sav` abundance matrix with m/z and RT metadata
- **Annotation**: Local standard library matching + PubChem API
- **Statistics**: t-test/Mann-Whitney, OPLS-DA, VIP, FDR correction
- **Pathway**: KEGG ORA, Mummichog m/z-to-pathway
- **Tools**: Python, pandas, scipy, scikit-learn, pyreadstat

## Key Results

- **H₂O₂-induced Neuro-2a cells**: 116 significant features in POS mode; log2FC up to 14.27.
- **Annotation**: 116/116 POS features annotated via PubChem; 20 unique KEGG IDs at 5 ppm.
- **Pathways**: 39 Mummichog pathways, 12 at FDR < 0.25.

## Repository Contents

```
.
├── README.md
├── run_analysis.py           # Main pipeline entry point
├── *.py                      # Analysis modules
├── library/                  # Compound libraries
├── output/                   # Example outputs
└── requirements.txt
```

## How to Run

```bash
python run_analysis.py \
  --input ./data/raw_data.sav \
  --lib ./library/plant_library.csv \
  --output ./output
```

## Requirements

See [`requirements.txt`](./requirements.txt).
