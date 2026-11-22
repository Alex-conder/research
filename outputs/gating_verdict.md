# gating_verdict.md — 空位 3：可靠性门控（Day5，2026-10-13）

## 问题

GEARS 对某个留出扰动会不会答错，能否用**模型自身的廉价特征**提前门控？（CoT benchmark ECE 置信度探针思路的迁移检验）

## 方法

- 数据：seed 0，217 个留出扰动（predictions_seed0.csv + pred_expr_seed0.npz）。
- 标签：**error = 1 − pearson_delta**（top-50 DEG 口径）；**y = error > 中位数（=0.4141）即判「答错」**——定义性选择，对应「预测质量落在全体留出扰动后一半」。
- 特征（仅模型侧）：预测谱 8 项（delta_l2 / delta_l1 / delta_max / delta_top50_mean / n_de_025 / frac_de_025 / neg_frac / delta_entropy，delta = 预测 − ctrl 训练均值）+ log_n_cells。
- 排除（泄漏防护）：**strength 由留出组真实表达算出（step1_run.py:96-105），属真值派生，已排除**；pearson_delta_all 同理；in_pert_graph 为常数列自动剔除。
- 探针：L2 logistic 回归（C=1.0），fold 内标准化，分层 5-fold CV（shuffle, seed=0），纯 CPU。

## 结果

- **5-fold CV AUROC = 0.6556 ± 0.0480**（逐 fold: 0.657, 0.682, 0.578, 0.656, 0.706；pooled = 0.6547）
- 单特征消融（CV AUROC 均值）：
    - delta_max: 0.6274 ± 0.0542
    - delta_top50_mean: 0.6186 ± 0.0518
    - n_de_025: 0.6176 ± 0.0506
    - frac_de_025: 0.6176 ± 0.0506
    - delta_l2: 0.6139 ± 0.0526
    - delta_l1: 0.6118 ± 0.0504
    - neg_frac: 0.6020 ± 0.0516
    - delta_entropy: 0.5743 ± 0.0527
    - log_n_cells: 0.4721 ± 0.0467
- 最好单特征：**delta_max** = 0.6274

## 判据（定义性阈值，写死）

**CV AUROC 均值 > 0.6 → 「门控信号存在」**（0.5 = 随机；与 CoT benchmark 中 ECE 探针达标线同一量级，此为项目内约定阈值，非外部标准）。

## 结论

**AUROC = 0.6556 > 0.6 → 门控信号存在。** 模型自身的廉价特征（以 delta_max 为主）能显著优于随机地提前识别 GEARS 会答错的留出扰动——**支持空位 3 成立**：可靠性门控可作为虚拟细胞预测的部署前过滤层。

## 对照与后续

- GEARS 原生 uncertainty 列全 NaN（单种子无 ensemble 方差）。**种子 1/2 训练完成后，用 3 种子预测方差补全 uncertainty，重跑本脚本做对照**——若原生不确定性有效，应优于或互补本探针。
- 判阴性/阳性的稳健性：可换 pearson_delta_all 口径重跑标签，看中位数二分结论是否翻转。

## 复现

- 脚本：step6_gating_probe.py（特征 + 探针核心），step6_run.py（运行器/产物）
- 产物：outputs/gating_probe.json, gating_probe_cv.png, 本文件