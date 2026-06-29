# Multi-scale Virtual Cell Modeling

## Background

This project builds a multi-scale virtual cell framework that couples three modeling paradigms:
- **Mechanism-driven**: Hodgkin-Huxley (HH) neuronal electrophysiology
- **Constraint-driven**: Flux Balance Analysis (FBA) of genome-scale metabolism
- **Data-driven**: scVI single-cell variational inference and virtual perturbation

The framework is applied to Alzheimer's disease (AD) and extended to cerebral ischemia via the PirB/LILRB2 project.

## Methods

- **HH Model**: 2D excitability phase diagram, 5-state Markov sodium channel, two-compartment cable model
- **FBA**: COBRApy, E. coli core/iML1515, three-cell neuron-microglia-astrocyte model
- **scVI**: scvi-tools, PBMC 3k, AD snRNA-seq (GSE138852)
- **Integration**: Linear and non-linear interfaces between HH, FBA, and scVI modules

## Key Results

| Interface | R² | Description |
|---|---|---|
| HH → FBA | 0.9996 | Firing rate → ATP demand |
| FBA → HH | 0.6587 | ATP capacity → K-ATP channel → excitability |
| scVI → HH | 1.0000 | Gene perturbation latent shift → Kv block |
| scVI → FBA | 1.0000 | Gene perturbation latent shift → ETC inhibition |

**AD Cascade**: Cerebrovascular hypoperfusion → Aβ clearance reduction → neuroinflammation → neuronal hyperexcitability → mitochondrial dysfunction.

## Repository Contents

```
.
├── README.md
├── notebooks/                # HH, FBA, scVI, and integration notebooks
├── src/                      # Reusable functions
├── data/                     # Data download instructions
├── outputs/                  # Expected output directory
├── run_all.py                # One-click run all notebooks
└── requirements.txt
```

## How to Run

```bash
python -m venv venv
source venv/bin/activate  # Windows: venv\Scripts\activate
pip install -r requirements.txt
python run_all.py
```

## Requirements

See [`requirements.txt`](./requirements.txt).
