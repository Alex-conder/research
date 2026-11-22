# Conformal 门控接入预测流水线总结

## 1. 做法
- 用模型侧 9 特征训练 Ridge 回归器预测 `err = 1 - pearson_delta`。
- Split Conformal：train/calibration/test = 60/20/20，在 calibration 上取残差 90% 分位数作为 margin `q`。
- 对 test 每个扰动，若 `err_pred + q < tau` 则保留 GEARS 预测；否则退回 `ctrl_mean`。
- 在真实表达谱上重新计算保留前后的 `pearson_delta`。

## 2. 结果（K562 essential，217 个留出扰动，10 次随机划分）
- **GEARS 全集 mean pearson_delta ≈ 0.9854**（预测质量已经很高）。
- **Conformal 区间覆盖 ≈ 0.886 @ α=0.1**，平均宽度 ≈ 0.59，说明区间统计有效但较宽。
- **err_pred 对高误差扰动的 AUROC ≈ 0.665**，与 gating probe 的 0.654 同量级——有信号但不够强。
- **门控后平均质量提升极小**：最佳 tau≈0.775 时，gated mean = 0.9860，比全集高 0.0006。
- **accepted set 的实际质量反而略低于 rejected**：accepted mean 0.981 vs rejected GEARS mean 0.987。说明当前误差预测模型会把一些其实预测得好的扰动判为“低置信”，把部分难样本判为“高置信”。

## 3. 解释
- 这不是 conformal 方法本身失败：覆盖率和宽度都正常。
- 失败的是**误差预测模型**：R2_train ≈ 0.11，AUROC ≈ 0.66，模型只能解释 10% 的误差方差。
- 在 GEARS 已经很强的数据上（pearson≈0.98），ctrl_mean 退守策略的边际收益也非常有限。

## 4. 对论文的写法
- 把 conformal gating 作为**不确定性量化的技术路线**，而不是**精度提升的黑箱**。
- 强调：门控的有效性上限由误差预测模型的质量决定；当前特征只能给出“宽区间、弱排序”，未来需要更强的模型侧特征（如集成方差、注意力熵、扰动图结构）或更大的校准集。
- 这一阴性结果反而支持论文的核心论点：门控不是“加一个阈值就能发”，需要诚实地报告 coverage/width/AUROC 三件套。

产物：
- `outputs/conformal_gating_pipeline.csv`
- `outputs/conformal_gating_pipeline_summary.md`
- `figures/conformal_gating_pipeline.png`
