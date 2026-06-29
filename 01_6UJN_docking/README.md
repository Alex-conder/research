# 6UJN (P-glycoprotein/ABCB1) Sequential Docking & MD/MM-GBSA

## Background

P-glycoprotein (P-gp/ABCB1/MDR1) is a key efflux transporter mediating multidrug resistance. This project investigates the sequential binding of two clinical antibiotics — Meropenem (carbapenem) and Vancomycin (glycopeptide) — to the 6UJN (P-gp C952A mutant) structure using molecular docking, molecular dynamics relaxation, and MM-GBSA binding free energy calculations.

## Methods

- **Docking**: AutoDock Vina 1.2.5
- **MD Relaxation**: AmberTools + OpenMM 8.5.2 (ff14SB/GAFF2/TIP3P, 100 ps NPT)
- **Free Energy**: MMPBSA.py (igb=2, saltcon=0.15 M)
- **Analysis**: Python (pandas, matplotlib, seaborn)

## Key Results

| Order | 2nd Ligand | Docking E (kcal/mol) | Ligand-Ligand COM (Å) | H-bonds | MM-GBSA ΔG (kcal/mol) |
|---|---|---|---|---|---|
| Meropenem → Vancomycin | Vancomycin | −8.502 | **7.38** | **30** | **−52.85 ± 4.20** |
| Vancomycin → Meropenem | Meropenem | −6.738 | **51.05** | 2 | **−23.03 ± 2.48** |

**Conclusion**: Meropenem pre-binding creates a ligand-dependent pre-organization effect that enables Vancomycin to co-bind in a crowded ternary complex.

## Repository Contents

```
.
├── README.md
├── scripts/                  # Docking analysis and PDB preparation scripts
├── input/                    # Receptor/ligand PDBQT files and config
├── output/                   # Docking logs, CSV results, and final poses
├── 04_reports/               # Validation reports (MD/MM-GBSA)
└── 6UJN_Meropenem_Vancomycin.csv
```

## How to Run

1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
2. Run sequential docking:
   ```bash
   # Example: analyze docking results
   python scripts/analyze_results.py
   ```
3. View reports in `04_reports/`.

## Requirements

See [`requirements.txt`](./requirements.txt).
