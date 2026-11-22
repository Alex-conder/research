"""Step 4: H2a / H2a-null / H2a-partial / H2b 判据检验（2026-10-12 增补版）。

H2a:       corr(violation_score, error)                r > 0.3, p < 0.05
H2a-null:  置换检验——打乱基因-反应映射 ≥100 次，真实 r 须 > null 95 分位
H2a-partial: 偏相关——控制扰动强度（敲低效率/全基因组效应量）后仍显著
H2b:       corr(violation_score, uncertainty) < 0.3

产出 outputs/scatter_h2a.png, null_distribution.png, corr_matrix_h2b.txt, verdict.md
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr, perm_test_stub  # placeholder

pred = pd.read_csv("outputs/predictions.csv")
# 列契约（Day2 须满足）：perturbation, pearson_delta, uncertainty, strength
#   strength = 扰动强度协变量（敲低效率或扰动的全基因组效应量），
#   用于 H2a-partial 偏相关（预注册 2026-10-12 增补）
strength = pred[["perturbation", "strength"]]
df = pred.merge(viol, on="perturbation")

error = 1 - df["pearson_delta"]          # 误差口径（锁定）

def pearson(a, b):
    m = np.isfinite(a) & np.isfinite(b)
    return pearsonr(a[m], b[m])

# ---- H2a ----
r_h2a, p_h2a = pearson(df["violation_score"].values, error.values)

# ---- H2a-null：置换基因-反应映射 ----
# 实现方式：对 enzyme_proxy 的列（反应）随机重排 n_perm=200 次，重算违例分数与 r
# 注意：严格做法应固定表达、打乱 GPR 映射；此处等价于打乱反应列标签
n_perm = 200
rng = np.random.default_rng(42)
proxy = pd.read_csv("outputs/enzyme_proxy.csv").set_index("perturbation")
import cobra
model = cobra.io.read_sbml_model("data/Human-GEM.xml")
S = np.asarray(cobra.util.create_stoichiometric_matrix(model), dtype=float)
null_r = []
for i in range(n_perm):
    idx = rng.permutation(S.shape[1])
    v_null = proxy[S_cols_order].values[:, idx]
    vs = np.linalg.norm(S @ v_null.T, axis=0)
    r_null, _ = pearson(vs, error.values)
    null_r.append(r_null)
null_95 = np.percentile(null_r, 95)

plt.figure(figsize=(5, 3.5))
plt.hist(null_r, bins=30, alpha=0.7, label="permutation null")
plt.axvline(r_h2a, color="red", label=f"observed r={r_h2a:.3f}")
plt.axvline(null_95, color="gray", ls="--", label=f"null 95th={null_95:.3f}")
plt.legend(); plt.tight_layout(); plt.savefig("outputs/null_distribution.png", dpi=150)

# ---- H2a-partial：控制扰动强度 ----
# 偏相关 r(viol, err · strength) = (r_ve - r_vs*r_es) / sqrt((1-r_vs^2)(1-r_es^2))
r_vs, _ = pearson(df["violation_score"].values, df["strength"].values)
r_es, _ = pearson(error.values, df["strength"].values)
r_partial = (r_h2a - r_vs * r_es) / np.sqrt(max(1e-12, (1 - r_vs**2) * (1 - r_es**2)))
# 偏相关 p 值近似：n-3 自由度下的 t 检验
n = len(df)
t_partial = r_partial * np.sqrt((n - 3) / max(1e-12, 1 - r_partial**2))
from scipy.stats import t as t_dist
p_partial = 2 * (1 - t_dist.cdf(abs(t_partial), df=n - 3))

# ---- H2b ----
r_h2b, p_h2b = pearson(df["violation_score"].values, df["uncertainty"].values)

# ---- H2a 图 ----
plt.figure(figsize=(5, 4))
plt.scatter(df["violation_score"], error, s=18)
plt.xlabel("conservation violation score"); plt.ylabel("prediction error (1 - Pearson delta)")
plt.title(f"H2a: r = {r_h2a:.3f}, p = {p_h2a:.4f}")
plt.tight_layout(); plt.savefig("outputs/scatter_h2a.png", dpi=150)

corr = df[["violation_score", "uncertainty", "pearson_delta", "strength"]].corr()
with open("outputs/corr_matrix_h2b.txt", "w") as f:
    f.write(str(corr))

h2a_pass = (r_h2a > 0.3) and (p_h2a < 0.05)
null_pass = r_h2a > null_95
partial_pass = p_partial < 0.05
h2b_pass = r_h2b < 0.3

verdict = f"""# 判据结果（自动生成）
- H2a:       r = {r_h2a:.3f}, p = {p_h2a:.4f} -> {'PASS' if h2a_pass else 'FAIL'}
- H2a-null:  observed r vs null 95th ({null_95:.3f}) -> {'PASS' if null_pass else 'FAIL'}
- H2a-partial: r_partial = {r_partial:.3f}, p = {p_partial:.4f} -> {'PASS' if partial_pass else 'FAIL'}
- H2b:       r = {r_h2b:.3f} (阈值 < 0.3) -> {'PASS' if h2b_pass else 'FAIL'}
n = {n}；功效声明见 preregistration.md（n<30 时 p 值仅供辅助）
处置（按 preregistration.md）: {
  '进入二阶段（H1 投影校正）' if (h2a_pass and null_pass and partial_pass and h2b_pass) else
  '转仪表盘路线' if (h2a_pass and null_pass and partial_pass) else
  '归档阴性/伪相关/强度代理，归因分析，方案终止'}
"""
with open("outputs/verdict.md", "w") as f:
    f.write(verdict)
print(verdict)
