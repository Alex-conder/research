"""Step 3: 守恒违例分数。

一阶段锁定：violation(g) = ||S . v_hat(g)||_2   （不做归一化，先看分布）

技术决策留痕（2026-10-12 评审意见）：
  - 之前版本除以 sqrt(n_reactions)，已移除——S 矩阵系数与 v_hat 单位均为近似，
    伪归一化会掩盖分布信息；先看 raw 分布，再决定归一化方式（预注册已锁定此条）。
  - v_hat 的三种可选近似（一阶段锁定用 (a)）：
    (a) v_hat_r = 酶表达代理值直接当通量代理（量纲不严谨但单调，零成本）★锁定
    (b) v_hat_r = 代理值 × 反应方向先验（可逆反应取正）
    (c) v_hat = argmin ||v - proxy|| s.t. S v = b（最小范数解，半软约束）
产出 outputs/violation_score.csv: perturbation, violation_score
"""
import numpy as np
import pandas as pd

import cobra
model = cobra.io.read_sbml_model("data/Human-GEM.xml")
S = cobra.util.create_stoichiometric_matrix(model)      # (m_metabolites, n_reactions)
S = np.asarray(S, dtype=float)

proxy = pd.read_csv("outputs/enzyme_proxy.csv").set_index("perturbation")
v = proxy[S_cols_order].values                          # 列须与 S 的反应列一致

violation = np.linalg.norm(S @ v.T, axis=0)             # raw，不归一化（锁定）
pd.DataFrame({"perturbation": proxy.index, "violation_score": violation}).to_csv(
    "outputs/violation_score.csv", index=False)
print("Step 3 done -> outputs/violation_score.csv")
print(pd.Series(violation).describe())                  # 先看分布
