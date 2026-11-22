# audit_verdict.md — 门控信号机制审计（Step9，2026-10-14）

## 四个问题的答案（replogle_rpe1_essential，n=306，seed 0）

1. **幅度家族消融**：全特征 CV AUROC = 0.7965；只留幅度家族 =
   0.8021；**非幅度残余 = 0.7970**。
2. **置换 null**（1000 次）：null 均值 0.4922 ±
   0.0461，观察到 0.7894 的 p = 0.0000。
3. **Bootstrap 95% CI**（pooled OOF，1000 次重采样）：
   [0.7398, 0.8370]。
4. **覆盖-风险**：过滤最可疑 20% 查询 → 捕获
   37.9% 的总误差，保留预测均误差
   0.2028（全体 0.2614）。
5. **强度分解**（真值诊断）：rho(error, strength) =
   -0.330；rho(error, delta_max) =
   -0.504；strength 单独 AUROC =
   0.6323（不可部署，仅诊断）。

## 判读（与预注册判据无关，不改写 step6 归档结论）

见 audit.json 的 interpretation_rule。产物：audit.json、
audit_coverage_risk.png、本文件。
