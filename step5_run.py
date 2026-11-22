"""Day4 E1（探索性，已预登记于 preregistration.md）：delta 表达代理的违例分数。

背景：Day3 绝对表达代理 H2a 阴性、违例 CV 仅 5.9%；E1 改用
delta = 预测表达 − ctrl 均值 作为酶表达代理输入，其余口径与 Day3 完全同
（同一 GPR AND-only、同一高变 67 基因子集、同一 167 反应、S_sub、
‖S·v̂‖₂ 不归一化），保证可比。

纪律：探索性——不复活 H2a 主判据，p 值标注"探索性"；报告效应量+CI。

输入：
    outputs/pred_expr_seed0.npz                       (217×5000 预测表达)
    data/gears_demo/replogle_k562_essential/perturb_processed.h5ad  (ctrl 均值)
    data/human_gem/parsed_model_*.{npz,pkl}           (step2)
    outputs/predictions.csv                            (Day2 error/strength，同 step4)
输出：
    outputs/violation_delta.csv
    outputs/scatter_e1_delta.png
    outputs/verdict_delta.md
"""
import json
import pickle

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats

OUT = "outputs"
N_PERM = 200
rng = np.random.default_rng(0)   # 与 step4 同种子，置换方法一致

# ---- 1. 预测谱 + ctrl 均值 → delta ----
z = np.load(f"{OUT}/pred_expr_seed0.npz")
pred_expr = z["pred_expr"]
perts = z["pred_perturbations"].tolist()
pred_genes = z["pred_genes"].tolist()
gidx = {g: i for i, g in enumerate(pred_genes)}

adata = ad.read_h5ad(
    "data/gears_demo/replogle_k562_essential/perturb_processed.h5ad")
var_genes = adata.var["gene_name"].tolist()
ctrl = adata[adata.obs.condition == "ctrl"]
ctrl_mean = np.asarray(ctrl.X.mean(axis=0)).ravel()     # 按 adata.var 顺序
del adata, ctrl
# 对齐到 pred_genes 顺序（同序则直接可用）
if var_genes == pred_genes:
    ctrl_aligned = ctrl_mean
else:
    m = {g: v for g, v in zip(var_genes, ctrl_mean)}
    missing = [g for g in pred_genes if g not in m]
    assert not missing, f"{len(missing)} 个面板基因不在 h5ad var 中"
    ctrl_aligned = np.array([m[g] for g in pred_genes])

delta = pred_expr - ctrl_aligned.astype(np.float32)     # (217, 5000)
print(f"[E1] delta: {delta.shape}, ctrl_mean[:3]={ctrl_aligned[:3]}")

# ---- 2. 模型与映射（同 step3）----
S = sparse.load_npz("data/human_gem/parsed_model_S.npz").tocsr()
meta = np.load("data/human_gem/parsed_model_meta.npz", allow_pickle=True)
rxn_ids = meta["rxn_ids"].tolist()
with open("data/human_gem/parsed_model_rxn_genes.pkl", "rb") as f:
    rxn_genes = pickle.load(f)
model_genes = set(meta["gene_symbols"].tolist())
mapped_genes = sorted(model_genes & set(pred_genes))

# 高变子集：与 Day3 同一算法（绝对预测表达跨 217 方差 top 25%）→ 同一 67 基因
var = pred_expr[:, [gidx[g] for g in mapped_genes]].var(axis=0)
n_top = max(1, int(0.25 * len(mapped_genes)))
hv_genes = set(np.array(mapped_genes)[np.argsort(-var)[:n_top]])
keep = [j for j, gl in enumerate(rxn_genes)
        if all(g in gidx for g in gl) and all(g in hv_genes for g in gl)]
S_k = S[:, keep]
print(f"[E1] 高变子集 {n_top}/{len(mapped_genes)}，保留反应 {len(keep)}/{len(rxn_ids)}"
      f"（与 Day3 同口径）")

# ---- 3. delta 代理 v̂ 与违例分数 ----
hv_list = sorted(hv_genes)
col = {g: i for i, g in enumerate(hv_list)}
rows_r, cols_r, vals_r = [], [], []
for i, j in enumerate(keep):
    genes_r = [g for g in rxn_genes[j] if g in col]
    for g in genes_r:
        rows_r.append(i)
        cols_r.append(col[g])
        vals_r.append(1.0 / len(genes_r))
R = sparse.csr_matrix((vals_r, (rows_r, cols_r)),
                      shape=(len(keep), len(hv_list)))
E = delta[:, [gidx[g] for g in hv_list]]               # (217, n_hv) delta 谱
V = (R @ E.T).T
Sv = np.asarray(S_k @ V.T)
violation = np.sqrt((Sv ** 2).sum(axis=0)).ravel()

df_v = pd.DataFrame({"perturbation": perts, "violation_delta": violation})
df_v.to_csv(f"{OUT}/violation_delta.csv", index=False)
cv = df_v.violation_delta.std() / df_v.violation_delta.mean()   # ddof=1，与 describe 一致
print(f"[E1] violation_delta describe():")
print(df_v.violation_delta.describe())
print(f"[E1] CV = {cv*100:.2f}%（Day3 绝对表达代理 CV = 5.86%）")

# ---- 4. 四项读出 ----
pred = pd.read_csv(f"{OUT}/predictions.csv")
df = pred.merge(df_v, on="perturbation", how="inner").dropna(
    subset=["pearson_delta", "violation_delta", "strength"])
df["error"] = 1.0 - df["pearson_delta"]
n = len(df)
assert n >= 20, f"样本太少 n={n}"
v, e, s = df.violation_delta.values, df.error.values, df.strength.values

def fisher_ci(r, n):
    zr = np.arctanh(np.clip(r, -0.999999, 0.999999))
    lo, hi = zr - 1.96 / np.sqrt(n - 3), zr + 1.96 / np.sqrt(n - 3)
    return np.tanh(lo), np.tanh(hi)

# a. 分布/CV（上方已打印）
# b. r（效应量 + 95% CI）
r_obs, p_obs = stats.pearsonr(v, e)
r_lo, r_hi = fisher_ci(r_obs, n)
print(f"[E1] r = {r_obs:.4f} (95% CI [{r_lo:.4f}, {r_hi:.4f}]), "
      f"p = {p_obs:.4g} (探索性)")

# c. 置换 null（200 次，同 Day3）
r_null = np.array([stats.pearsonr(v, e[rng.permutation(n)])[0]
                   for _ in range(N_PERM)])
null_p95 = np.percentile(r_null, 95)
emp_p = (1 + np.sum(r_null >= r_obs)) / (1 + N_PERM)
print(f"[E1] null 95% = {null_p95:.4f}, empirical p = {emp_p:.4f} (探索性)")

# d. 偏相关（控制 strength，残差法，t df=n-3）
def resid(y, x):
    b = np.polyfit(x, y, 1)
    return y - np.polyval(b, x)

r_par = stats.pearsonr(resid(v, s), resid(e, s))[0]
t_par = r_par * np.sqrt((n - 3) / (1 - r_par ** 2))
p_par = 2 * stats.t.sf(abs(t_par), n - 3)
par_lo, par_hi = fisher_ci(r_par, n)
print(f"[E1] partial|r(strength) = {r_par:.4f} (95% CI [{par_lo:.4f}, {par_hi:.4f}]),"
      f" t({n-3}) = {t_par:.3f}, p = {p_par:.4g} (探索性)")

# ---- 5. 图 ----
fig, ax = plt.subplots(figsize=(6, 4.5))
ax.scatter(v, e, alpha=0.6)
b = np.polyfit(v, e, 1)
xs = np.linspace(v.min(), v.max(), 50)
ax.plot(xs, np.polyval(b, xs), "r--", lw=1)
ax.set_xlabel("conservation violation ||S*vhat_delta||_2 (unnormalized)")
ax.set_ylabel("error = 1 - pearson_delta (top50 DEG)")
ax.set_title(f"E1 exploratory: r={r_obs:.3f}, 95%CI [{r_lo:.3f},{r_hi:.3f}], n={n}")
fig.tight_layout()
fig.savefig(f"{OUT}/scatter_e1_delta.png", dpi=150)
plt.close(fig)

# ---- 6. verdict_delta.md（含与 Day3 对照表）----
lines = []
lines.append("# verdict_delta.md — E1 探索性检验：delta 表达代理（Day4，2026-10-13）\n")
lines.append("> 本报告为**探索性分析**（E1，出数据前登记于 preregistration.md"
             "「探索性分析登记」小节）。按多重比较锁定条款：**不复活 H2a 主判据、"
             "不进主判据推断**；下述 p 值一律为探索性 p，结论只看效应量与置信区间。\n")
lines.append("## 方法对照（与 Day3 唯一差异）\n")
lines.append("| 环节 | Day3（绝对表达代理） | E1（delta 代理） |")
lines.append("|---|---|---|")
lines.append("| 酶表达输入 | 预测表达绝对值 | 预测表达 − ctrl 训练均值 |")
lines.append("| GPR 清洗 | AND-only，3,943 反应 | 同左（不变） |")
lines.append("| 高变子集 | 67 基因（top 25%） | **同一 67 基因**（不变） |")
lines.append("| 保留反应 | 167，S_sub 4,308×167 | 同左（不变） |")
lines.append("| 违例分数 | ‖S·v̂‖₂ 不归一化 | 同左（不变） |")
lines.append("| error 来源 | outputs/predictions.csv | 同左（不变） |\n")
lines.append("## 四项读出\n")
lines.append(f"n = {n}（合并后，与 Day3 相同）\n")
lines.append(f"1. **违例分数分布**：mean {df_v.violation_delta.mean():.1f}，"
             f"std {df_v.violation_delta.std():.1f}，"
             f"min {df_v.violation_delta.min():.1f}，"
             f"中位 {df_v.violation_delta.median():.1f}，"
             f"max {df_v.violation_delta.max():.1f}；**CV = {cv*100:.2f}%**"
             f"（Day3 = 5.86%）。")
lines.append(f"2. **r（violation_delta × error）= {r_obs:.4f}**，"
             f"95% CI [{r_lo:.4f}, {r_hi:.4f}]，探索性 p = {p_obs:.4g}。")
lines.append(f"3. **置换 null（200 次）**：null 95 分位 {null_p95:.4f}，"
             f"empirical p = {emp_p:.4f}（探索性）；观测 r "
             f"{'高于' if r_obs > null_p95 else '未超'} null 95 分位。")
lines.append(f"4. **偏相关（控制 strength）= {r_par:.4f}**，"
             f"95% CI [{par_lo:.4f}, {par_hi:.4f}]，t({n-3}) = {t_par:.3f}，"
             f"探索性 p = {p_par:.4g}。\n")
lines.append("## 与 Day3 对照表\n")
lines.append("| 读出 | Day3 绝对表达 | E1 delta |")
lines.append("|---|---|---|")
lines.append(f"| CV | 5.86% | {cv*100:.2f}% |")
lines.append(f"| r | 0.0226 | {r_obs:.4f} |")
lines.append(f"| null 95% 分位 | 0.0951 | {null_p95:.4f} |")
lines.append(f"| empirical p | 0.3333 | {emp_p:.4f} |")
lines.append("| partial r | -0.2181 | " + f"{r_par:.4f} |\n")
lines.append("## 解读与决策边界\n")
cv_recovered = cv > 2 * 0.0586   # 恢复判据：明显高于 Day3（内部参考线，非判据）
lines.append(f"- CV 恢复（明显高于 5.86%）？**{'是' if cv_recovered else '否'}**。")
lines.append(f"- 若 CV 恢复但偏相关为零/反号：信号仍可能是强度代理，"
             f"不构成升级证据。")
lines.append("- **E1 不触发任何自动后续动作**；是否升级二阶段由用户决策。\n")
lines.append("## 复现\n")
lines.append("- 脚本：step5_run.py（delta 代理，同管线）")
lines.append("- 产物：outputs/violation_delta.csv, scatter_e1_delta.png, "
             "verdict_delta.md")
open(f"{OUT}/verdict_delta.md", "w", encoding="utf-8").write("\n".join(lines))
print(f"[E1] -> {OUT}/violation_delta.csv, scatter_e1_delta.png, verdict_delta.md")

json.dump({
    "formula": "violation_delta_i = ||S_sub @ vhat_delta_i||_2; "
               "vhat_delta_r = mean delta expression of AND genes; "
               "delta = pred_expr - ctrl_train_mean",
    "n_perts": len(perts), "n_mapped_genes": len(mapped_genes),
    "n_hv_genes": n_top, "n_rxns_kept": len(keep),
    "cv_percent": float(cv * 100),
    "r": float(r_obs), "r_ci": [float(r_lo), float(r_hi)], "p_exploratory": float(p_obs),
    "null_p95": float(null_p95), "emp_p_exploratory": float(emp_p),
    "partial_r": float(r_par), "partial_ci": [float(par_lo), float(par_hi)],
    "partial_p_exploratory": float(p_par),
}, open(f"{OUT}/e1_delta_setup.json", "w"), indent=2, ensure_ascii=False)
print(f"[E1] -> {OUT}/e1_delta_setup.json")
