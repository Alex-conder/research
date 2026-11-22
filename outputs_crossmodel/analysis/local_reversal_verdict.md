# step9c 审稿人补充检验 verdict

## 检验 A：16/217 局部优胜置换检验（主口径 = step7 归档标签）
- wins_obs = 16/217（7.4%），零分布均值 = 33.73（范围 21–47），上侧 p = 1.000，下侧 p = 1.00e-04，双侧 p = 2.00e-04（B=10000）。
- 判据落档：**no_excess_wins_below_independence**。
- observed 16 wins lies FAR BELOW the independence-null mean (~34): the two models' per-perturbation scores are positively associated, so random overlap would predict MORE wins than observed. The 16 wins are neither an excess of structure nor an artifact of overlap noise.
- 灵敏度（step7b 统一 4813 基因空间重算）：wins = 14/217，零分布均值 31.61，下侧 p = 1.00e-04，方向与主口径一致。
- 定量锚：oracle−全 GEARS headroom = 0.0152±0.0069，仅占全 GEARS（0.561）的 2.7%——即使结构显著也不影响主结论。

## 检验 B：迁移 vs 就地配对 bootstrap（RPE1 306）
- AUROC：迁移 0.758 vs 就地（in-sample）0.802（归档 CV 参照 0.797）。
- 差值均值 = -0.0441，95% CI [-0.0776, -0.0099]，P(<0) = 9.95e-01。
- 判据落档：**small_degradation**——CI 全 < 0 → detection transfers with a small, significant degradation
