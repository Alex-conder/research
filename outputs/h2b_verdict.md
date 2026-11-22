# h2b_verdict.md — H2b 补全检验：ensemble 不确定性 × 违例分数（Day5，2026-10-13）

> **路线关闭后的补全检验**：守恒约束路线已按预注册判阴性归档（verdict.md），本检验只为档案完整性，**不改变已归档结论**。

## 方法

- ensemble 不确定性：3 种子（seed 0/1/2，SPLIT_SEED=0 同考场）预测谱在共同 ctrl 训练均值基线上取 delta，逐基因跨种子方差（ddof=1）再对基因取均值。
- 违例分数：outputs/violation_score.csv（Day3 绝对表达口径，主判据文件）。
- 判据（preregistration.md 原文）：r < 0.3 通过；≥ 0.3 信号冗余。0.3 为判定阈值非对标值。

## 结果

- n = 217，**r = -0.6168，p = 3.862e-24** → **通过（违例提供独立信号）**

## 结论

H2b 通过。但因守恒路线已关闭，此结果仅入档，不触发 H1 或任何后续动作。

## 复现

- 脚本：step7_ensemble_gating.py, step7_run.py；产物：outputs/h2b_ensemble.json, 本文件