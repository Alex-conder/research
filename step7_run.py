"""Day5 Step7 运行器：3 种子 ensemble 不确定性 → H2b 补全检验 + 门控探针对照。

两件事（详见 step7_ensemble_gating.py docstring 的公式与判据）：
1. H2b（预注册判据）：corr(violation_score.csv 的违例, ensemble 不确定性)，
   判据 r < 0.3 通过 / ≥ 0.3 冗余。守恒路线已关闭，本检验只为档案完整性，
   结论标注"不改变已归档结论"。输出 outputs/h2b_ensemble.json + h2b_verdict.md。
2. 门控探针对照：step6 特征集 + log_ens_var，重跑 5-fold CV，
   与 step6 的 0.6556 对照；报集成特征 AUROC、单特征 AUROC、
   与 delta_max 的相关性。输出 outputs/gating_probe_ensemble.json。

文件缺失处理：pred_expr_seed{1,2}.npz / predictions_seed{1,2}.csv 任一缺失
（种子 1/2 仍在训练）→ 打印明确提示并以退出码 2 退出，不产出任何文件。

用法：
    source venv/Scripts/activate
    python step7_run.py
"""
import json
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd

from step6_gating_probe import SPECTRUM_FEATURES, compute_features, run_probe
from step7_ensemble_gating import SEEDS, corr_features, ensemble_uncertainty, run_h2b

OUT = "outputs"
GATING_CRITERION = 0.6          # step6 定义性阈值，沿用
H2B_CRITERION = 0.3             # preregistration.md 原文判据

# ---- 0. 输入文件齐备性检查（种子 1/2 训练中 → 退出码 2）----
need = [f"{OUT}/pred_expr_seed{s}.npz" for s in SEEDS] + \
       [f"{OUT}/predictions_seed{s}.csv" for s in SEEDS]
missing = [f for f in need if not os.path.exists(f)]
if missing:
    print("[step7] 以下文件尚不存在（种子 1/2 可能仍在训练）：")
    for f in missing:
        print(f"    - {f}")
    print("[step7] 待训练完成后重跑：python step7_run.py")
    sys.exit(2)

# ---- 1. 读 3 种子预测谱（校验扰动顺序一致）+ ctrl 训练均值 ----
mats, perts_ref = [], None
for s in SEEDS:
    z = np.load(f"{OUT}/pred_expr_seed{s}.npz")
    perts = z["pred_perturbations"].tolist()
    if perts_ref is None:
        perts_ref = perts
    assert perts == perts_ref, f"seed{s} 扰动顺序与 seed0 不一致，无法 ensemble"
    mats.append(z["pred_expr"])
    print(f"[step7] seed{s}: pred_expr {z['pred_expr'].shape}")

# backed 模式只物化 ctrl 子集：全量 X ~GB 级，训练并行期内存紧张会 OOM（同 step8 修法）
adata = ad.read_h5ad(
    "data/gears_demo/replogle_k562_essential/perturb_processed.h5ad",
    backed="r")
var_genes = adata.var["gene_name"].tolist()
ctrl_mean = np.asarray(
    adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
adata.file.close()
del adata
pred_genes = np.load(f"{OUT}/pred_expr_seed0.npz")["pred_genes"].tolist()
if var_genes != pred_genes:                      # 防御：按基因名重排（同 step5/6）
    m = dict(zip(var_genes, ctrl_mean))
    ctrl_mean = np.array([m[g] for g in pred_genes])

ens_var = ensemble_uncertainty(mats, ctrl_mean)
log_ens_var = np.log(ens_var)
print(f"[step7] ens_var: min={ens_var.min():.4g}, med={np.median(ens_var):.4g}, "
      f"max={ens_var.max():.4g}")

# ---- 2. H2b：违例分数 × ensemble 不确定性 ----
viol = pd.read_csv(f"{OUT}/violation_score.csv").set_index("perturbation")
viol = viol.loc[perts_ref]                       # 对齐扰动顺序
h2b = run_h2b(viol["violation"].values, ens_var)
h2b_pass = h2b["r"] < H2B_CRITERION
print(f"[step7] H2b: r={h2b['r']:.4f}, p={h2b['p']:.4g} "
      f"（判据 r < {H2B_CRITERION}）→ {'通过' if h2b_pass else '冗余'}")

json.dump({
    "task": "H2b post-hoc completeness check (conservation route closed)",
    "note": "路线关闭后的补全检验，不改变已归档结论",
    "formula": "ens_var_i = mean_g var_s(pred_expr_{s,i,g} - ctrl_mean), "
               "s in {0,1,2}, ddof=1",
    "violation_source": "outputs/violation_score.csv (Day3 绝对表达口径)",
    "criterion": "r < 0.3 pass / >= 0.3 redundant (preregistration.md)",
    **h2b, "r_pass": h2b_pass,
}, open(f"{OUT}/h2b_ensemble.json", "w", encoding="utf-8"),
    indent=2, ensure_ascii=False)

lines = []
lines.append("# h2b_verdict.md — H2b 补全检验：ensemble 不确定性 × 违例分数"
             "（Day5，2026-10-13）\n")
lines.append("> **路线关闭后的补全检验**：守恒约束路线已按预注册判阴性归档"
             "（verdict.md），本检验只为档案完整性，**不改变已归档结论**。\n")
lines.append("## 方法\n")
lines.append("- ensemble 不确定性：3 种子（seed 0/1/2，SPLIT_SEED=0 同考场）"
             "预测谱在共同 ctrl 训练均值基线上取 delta，逐基因跨种子方差"
             "（ddof=1）再对基因取均值。")
lines.append("- 违例分数：outputs/violation_score.csv（Day3 绝对表达口径，"
             "主判据文件）。")
lines.append(f"- 判据（preregistration.md 原文）：r < {H2B_CRITERION} 通过；"
             f"≥ {H2B_CRITERION} 信号冗余。0.3 为判定阈值非对标值。\n")
lines.append("## 结果\n")
lines.append(f"- n = {h2b['n']}，**r = {h2b['r']:.4f}，p = {h2b['p']:.4g}** → "
             f"**{'通过（违例提供独立信号）' if h2b_pass else '冗余（违例与不确定性同源）'}**\n")
lines.append("## 结论\n")
lines.append(f"H2b {'通过' if h2b_pass else '不通过'}。但因守恒路线已关闭，"
             f"此结果仅入档，不触发 H1 或任何后续动作。\n")
lines.append("## 复现\n")
lines.append("- 脚本：step7_ensemble_gating.py, step7_run.py；"
             "产物：outputs/h2b_ensemble.json, 本文件")
open(f"{OUT}/h2b_verdict.md", "w", encoding="utf-8").write("\n".join(lines))

# ---- 3. 门控探针集成对照（step6 特征 + log_ens_var）----
pred0 = pd.read_csv(f"{OUT}/predictions_seed0.csv").set_index("perturbation")
pred0 = pred0.loc[perts_ref]
feat = compute_features(mats[0], perts_ref, pred_genes, ctrl_mean)
feat["log_n_cells"] = np.log1p(pred0["n_cells"].values)
feat["log_ens_var"] = log_ens_var
feature_names = SPECTRUM_FEATURES + ["log_n_cells", "log_ens_var"]
X = feat[feature_names].values.astype(np.float64)
assert not np.isnan(X).any()

error = 1.0 - pred0["pearson_delta"].values
y = (error > np.median(error)).astype(int)

res = run_probe(X, y, feature_names, seed=0)
auroc6 = GATING_CRITERION  # 占位，下方尝试读 step6 实际值
try:
    prev = json.load(open(f"{OUT}/gating_probe.json", encoding="utf-8"))
    auroc6 = prev["cv_auroc_mean"]
except (OSError, KeyError):
    print("[step7] 警告：读不到 step6 的 gating_probe.json，对照基线用 0.6 线")
delta = res["auroc_mean"] - auroc6
corr_ev_dm = corr_features(log_ens_var, feat["delta_max"].values)
print(f"[step7] 集成探针 AUROC = {res['auroc_mean']:.4f} ± "
      f"{res['auroc_std']:.4f}（step6 = {auroc6:.4f}，Δ = {delta:+.4f}）")
print(f"[step7] log_ens_var 单特征 AUROC = "
      f"{res['single_feature']['log_ens_var']['auroc_mean']:.4f}")
print(f"[step7] corr(log_ens_var, delta_max) = {corr_ev_dm['r']:.4f} "
      f"(p={corr_ev_dm['p']:.3g})")

json.dump({
    "task": "gating probe with ensemble-uncertainty feature "
            "(vs step6 cheap-features-only)",
    "criterion": "CV AUROC mean > 0.6 (definitional, same as step6)",
    "step6_baseline_auroc": auroc6,
    "n_perts": int(len(y)),
    "features": feature_names,
    "cv_auroc_mean": res["auroc_mean"],
    "cv_auroc_std": res["auroc_std"],
    "cv_auroc_per_fold": res["auroc_per_fold"],
    "pooled_cv_auroc": res["pooled_auroc"],
    "delta_vs_step6": delta,
    "ens_var_single_feature_auroc": res["single_feature"]["log_ens_var"],
    "delta_max_single_feature_auroc": res["single_feature"]["delta_max"],
    "corr_log_ens_var_vs_delta_max": corr_ev_dm,
    "verdict_pass": bool(res["auroc_mean"] > GATING_CRITERION),
    "single_feature_auroc": res["single_feature"],
}, open(f"{OUT}/gating_probe_ensemble.json", "w", encoding="utf-8"),
    indent=2, ensure_ascii=False)
print(f"[step7] -> {OUT}/h2b_ensemble.json, h2b_verdict.md, "
      f"gating_probe_ensemble.json")
