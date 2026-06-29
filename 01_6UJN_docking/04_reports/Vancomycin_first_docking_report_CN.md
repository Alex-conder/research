# 万古霉素先结合（Vancomycin-first）→ 美罗培南分子对接补充报告

## 1. 实验目的

在 6UJN 受体上先对接 Vancomycin，再在其形成的 holo 受体上对接 Meropenem，以获取 Vancomycin-first 序贯对接的详细结果与可视化。

## 2. 输入文件与参数

| 项目 | 路径 / 数值 |
|---|---|
| 受体（holo，可视化用） | `D:\6UJN_Meropenem_Vancomycin\md\input_van_first\complex.pdb` |
| 对接配体 poses | `D:\Tools\Docking\docking_project\data\receptors\6UJN_Vancomycin\Meropenem_out.pdbqt` |
| 对接结果 CSV | `D:\6UJN_Meropenem_Vancomycin\input\6UJN_Vancomycin_Meropenem.csv` |
| 对接中心 (Å) | (56.0328, 17.1666, 25.8216) |
| 盒子大小 (Å) | 60.68 × 60.68 × 60.68 |
| Exhaustiveness | 8 |
| num_modes | 10 |
| seed | 42 |

## 3. 对接结果

| Mode | Binding Affinity (kcal/mol) | RMSD/lb | RMSD/ub | 到 VAN 质心距离 (Å) | vs apo Meropenem RMSD (Å) |
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

**最佳亲和力**: -6.738 kcal/mol（Mode 1）。

## 4. 关键观察

1. **Meropenem 远离 Vancomycin**：最佳 pose（Mode 1）的质心距 Vancomycin 质心约 **51.05 Å**。
2. **多模式分散**：10 个 pose 距 Vancomycin 的范围约 **16.7–53.7 Å**，说明 Meropenem 未聚集于 Vancomycin 位点。
3. **与 apo Meropenem 差异大**：最佳 pose 与 apo Meropenem（6UJN 直接对接）的重原子 RMSD 为 **52.81 Å**，结合模式显著不同。
4. **能量意义**：最佳亲和力 −6.738 kcal/mol 与 apo Meropenem 参考能 −6.788 kcal/mol 接近（ΔE ≈ +0.050 kcal/mol），表明 Vancomycin 预结合几乎未提升 Meropenem 在其位点的结合。

## 5. 可视化

### 5.1 复合物总览（最佳 pose）

![Vancomycin-first complex overview](figures/van_first_complex_overview.png)

*图 1. 蛋白 CA（灰）、Vancomycin（蓝）与最佳 Meropenem pose（红）。*

### 5.2 全部 10 个 Meropenem 对接 pose

![Vancomycin-first docking poses](figures/van_first_docking_poses.png)

*图 2. 10 个 Meropenem pose 按亲和力着色；Vancomycin（蓝）。*

### 5.3 Meropenem–Vancomycin 距离 vs 亲和力

![Distance vs affinity](figures/van_first_distance_vs_affinity.png)

*图 3. 各 pose 的 Meropenem 质心到 Vancomycin 质心距离与其对接亲和力。*

### 5.4 XY 平面投影

![XY projection](figures/van_first_xy_projection.png)

*图 4. XY 平面投影显示 Meropenem pose 分布于蛋白表面，远离 Vancomycin。*

## 6. 结论

> Vancomycin-first 序贯对接中，Meropenem 未占据 Vancomycin 附近的口袋，而是结合在远离 Vancomycin 的位点。最佳亲和力 −6.738 kcal/mol 与 apo 状态几乎相同，说明 Vancomycin 预结合对 Meropenem 结合无显著促进作用。这与后续对称 MD 松弛（COM ~46.76 Å，MM-GBSA ΔG = −23.03 kcal/mol）一致。

---

*报告生成时间：2026-06-26*
