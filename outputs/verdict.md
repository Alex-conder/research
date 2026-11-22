# verdict.md — 第一阶段 H2a 检验结论（Day3，2026-10-13）

## 结果

- n = 217 个留出测试扰动（seed 0，单种子；3 种子为预注册锁定设置，本轮时间盒偏离）
- **H2a：r = 0.0226，p = 0.7402**（判据：r > 0.3 且 p < 0.05）→ **未成立**
- H2a-null：置换 200 次，观测 r 的 empirical p = 0.3333（null 95 分位 = 0.0951）→ 未通过置换检验
- H2a-partial（控制 strength）：r = -0.2181，t(214) = -3.270，p = 0.001255
- **H2b：本轮跳过** —— uncertainty 全 NaN（种子 1/2 未训），待 ensemble 补齐后检验。

## 三条件分流

1. H2a 不成立（r < 0.1 或 p ≥ 0.05）→ 按预注册**归档阴性**，做归因分析（映射距离 / 守恒律不适用性），方案终止。
2. 归档格式与阳性结论相同。

## 归因线索（阴性时参考）

- strength（扰动效应量）与 error 的相关：r = -0.4994
- strength 与 violation 的相关：r = -0.3933
- 代理通量为表达均值粗近似，一阶段噪声预期内（预注册边界声明）。
## 协议偏离记录（按偏离条款并列报告）

1. OOD 划分：scPerturBench 官方留出名单未取得，用 seed=0 自行留出 20%（D2-0）。
2. 数据：GEARS 过滤版 replogle_k562_essential 子集，非完整 Replogle K562。
3. pearson_delta 差异基因数：top-50 为自定（GEARS 内部 top-20，scPerturBench 未确证）；全基因口径 pearson_delta_all 已在 CSV 并列。
4. uncertainty 全 NaN：3 种子 ensemble 未跑（时间盒偏离）。
5. 高变子集解释：表达方差基于 217 个**预测**谱（非原始数据）；保留反应 167/3943。

## 复现

- 脚本：step1_run.py / step2_run.py / step3_run.py / step4_run.py
- 中间产物：outputs/predictions.csv, violation_score.csv, pred_expr_seed0.npz, proxy_setup.json
- 日志：setup_log.md（Day3 小节）