# Vancomycin-first -> Meropenem Molecular Docking Supplement Report

## 1. Objective

Dock Vancomycin onto 6UJN first, then dock Meropenem onto the resulting holo receptor. This report documents the Vancomycin-first sequential docking results and visualization.

## 2. Input Files and Parameters

| Item | Path / Value |
|---|---|
| Receptor (holo, visualization) | `D:\6UJN_Meropenem_Vancomycin\md\input_van_first\complex.pdb` |
| Docked ligand poses | `D:\Tools\Docking\docking_project\data\receptors\6UJN_Vancomycin\Meropenem_out.pdbqt` |
| Docking result CSV | `D:\6UJN_Meropenem_Vancomycin\input\6UJN_Vancomycin_Meropenem.csv` |
| Docking center (A) | (56.0328, 17.1666, 25.8216) |
| Box size (A) | 60.68 x 60.68 x 60.68 |
| Exhaustiveness | 8 |
| num_modes | 10 |
| seed | 42 |

## 3. Docking Results

| Mode | Binding Affinity (kcal/mol) | RMSD/lb | RMSD/ub | Distance to VAN centroid (A) | vs apo Meropenem RMSD (A) |
|---|---|---|---|---|---|
| 1 | -6.738 | 0.000 | 0.000 | 51.05 | 52.809 |
| 2 | -6.713 | 30.680 | 33.330 | 51.71 | 40.233 |
| 3 | -6.685 | 36.790 | 40.410 | 53.70 | 42.406 |
| 4 | -6.630 | 3.079 | 4.411 | 18.39 | 4.720 |
| 5 | -6.459 | 4.192 | 6.486 | 16.71 | 7.188 |
| 6 | -6.459 | 36.480 | 40.060 | 44.02 | 46.083 |
| 7 | -6.418 | 29.430 | 32.080 | 51.16 | 40.228 |
| 8 | -6.415 | 3.711 | 8.083 | 26.83 | 16.456 |
| 9 | -6.376 | 50.870 | 55.590 | 53.24 | 41.737 |

**Best affinity**: -6.738 kcal/mol (Mode 1).

## 4. Key Observations

1. **Meropenem stays far from Vancomycin**: the best pose (Mode 1) centroid is **51.05 A** from the Vancomycin centroid.
2. **Dispersed poses**: the 10 poses span **16.7-53.7 A** from Vancomycin, indicating Meropenem does not cluster near the Vancomycin site.
3. **Large difference from apo pose**: heavy-atom RMSD vs apo Meropenem is **52.81 A**, showing a distinct binding mode.
4. **Energetic implication**: the best affinity (-6.738 kcal/mol) is close to the apo Meropenem reference (-6.788 kcal/mol), so Vancomycin pre-binding does not promote Meropenem binding to its own site.

## 5. Visualizations

### 5.1 Complex Overview (best pose)

![Vancomycin-first complex overview](figures/van_first_complex_overview.png)

*Figure 1. Protein CA (gray), Vancomycin (blue), and best Meropenem pose (red).*

### 5.2 All 10 Meropenem Docking Poses

![Vancomycin-first docking poses](figures/van_first_docking_poses.png)

*Figure 2. Ten Meropenem poses colored by affinity; Vancomycin in blue.*

### 5.3 Meropenem-Vancomycin Distance vs Affinity

![Distance vs affinity](figures/van_first_distance_vs_affinity.png)

*Figure 3. Distance from each Meropenem pose centroid to Vancomycin centroid versus docking affinity.*

### 5.4 XY Plane Projection

![XY projection](figures/van_first_xy_projection.png)

*Figure 4. XY projection showing Meropenem poses distributed on the protein surface, away from Vancomycin.*

## 6. Conclusion

> In the Vancomycin-first sequential docking, Meropenem does not occupy the pocket near Vancomycin; instead it binds at sites remote from Vancomycin. The best affinity (-6.738 kcal/mol) is almost identical to the apo reference, indicating no significant promoting effect of Vancomycin pre-binding on Meropenem. This is consistent with the subsequent symmetric MD relaxation (COM ~46.76 A, MM-GBSA dG = -23.03 kcal/mol).

---

*Report generated: 2026-06-26*
