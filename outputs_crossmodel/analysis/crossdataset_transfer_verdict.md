# step9b 判据落档：同一 K562 探针跨数据集迁移到 RPE1

- 迁移数字：K562 就地探针（deployment 版，217 扰动 in-sample 拟合，未重训、scaler 未 refit）原样打分 RPE1 全部 306 扰动，**AUROC = 0.758，bootstrap 95% CI [0.705, 0.809]（B=1000，扰动级，seed 20261022）；单侧置换 p = 0.0010**。
- 对照：K562 就地 0.656（归档 cv_auroc_mean）、RPE1 就地 0.797（归档 0.7965）、delta_max 系数符号迁移 0.768。
- 预注册判据落档 **transfer_ok**（≥0.6 且 CI 下界 > 0.5）：**检测跨数据集可迁移，签名方向不可迁移**——同一探针跨数据集检测性能几乎不损（0.758 vs RPE1 就地 0.797），而误差签名方向余弦仅 0.01–0.044。
- Abstract 措辞建议：可将 "error detection transfers across datasets while signature direction does not" 作为核心发现的补充句，与 "error signature is an emergent model×dataset property" 并列，形成"检测通用、签名特异"的完整叙事。
- 产物：`outputs_crossmodel/analysis/crossdataset_transfer.{json,png}`（300 dpi 三柱图）。
