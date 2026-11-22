"""Day5 Step6 运行器：门控探针（空位 3：可靠性门控）—— 特征 → 探针 → 判据 → 产物。

目的：检验"GEARS 对某个留出扰动会不会答错"能否被模型自身的廉价特征
预测（CoT benchmark ECE 置信度探针思路迁移到虚拟细胞语境）。

标签（定义性选择）：error = 1 - pearson_delta（top-50 DEG 口径），
二分类 y = error 是否高于全体留出扰动的中位数（"答错" = 后一半）。

特征（模型侧，泄漏纪律见 step6_gating_probe.py docstring）：
- 预测谱派生：delta_l2 / delta_l1 / delta_max / delta_top50_mean /
  n_de_025 / frac_de_025 / neg_frac / delta_entropy
  （delta = 预测表达 − ctrl 训练均值）
- 协变量：log_n_cells；in_pert_graph 全为 True（常数）自动剔除留痕
- **strength 排除**：step1_run.py 第 96-105 行证实 strength 由留出组
  真实表达 sub.X 计算，属真值派生，进特征即泄漏。

探针：sklearn LogisticRegression(L2, C=1.0)，Pipeline 内 StandardScaler，
分层 5-fold CV（shuffle, seed=0），纯 CPU。

判据（定义性阈值，写死）：CV AUROC 均值 > 0.6 → "门控信号存在"
（0.5 = 随机；CoT benchmark 的 ECE 探针达标线同一量级，此为项目内约定）。

对照说明：GEARS 原生 uncertainty 列全 NaN（单种子无 ensemble 方差）。
种子 1/2 训练完成后可用 3 种子 ensemble 预测方差补全 uncertainty，
届时重跑本脚本做对照（预期原生不确定性若有效应优于/互补本探针）。

用法：
    source venv/Scripts/activate
    python step6_run.py                       # 默认 K562，输出 outputs/
    DATA=replogle_rpe1_essential OUT_DIR=outputs_rpe1 python step6_run.py
输出（默认 outputs/，跨数据集时 $OUT_DIR/）：
    gating_probe.json       全部数字（CV AUROC、逐 fold、单特征表）
    gating_probe_cv.png     5-fold pooled ROC 曲线
    gating_verdict.md       判据结论
"""
import json
import os

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from step6_gating_probe import (COVARIATES, SPECTRUM_FEATURES,
                                compute_features, run_probe)

OUT = os.environ.get("OUT_DIR", "outputs")   # 跨数据集：OUT_DIR=outputs_rpe1（同 step1_run.py 模式）
DATA = os.environ.get("DATA", "replogle_k562_essential")
SEED = 0
os.makedirs(OUT, exist_ok=True)

# ---- 1. 读预测谱 + ctrl 训练均值 + 标签 ----
z = np.load(f"{OUT}/pred_expr_seed{SEED}.npz")
pred_expr = z["pred_expr"]
perts = z["pred_perturbations"].tolist()
pred_genes = z["pred_genes"].tolist()

# backed 模式只物化 ctrl 子集：全量 X ~GB 级，训练并行期内存紧张会 OOM（同 step8 修法）
adata = ad.read_h5ad(
    f"data/gears_demo/{DATA}/perturb_processed.h5ad", backed="r")
var_genes = adata.var["gene_name"].tolist()
ctrl_mean = np.asarray(
    adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
adata.file.close()
del adata
if var_genes == pred_genes:                      # 同序则直接可用（与 step5 同）
    ctrl_aligned = ctrl_mean
else:                                            # 防御：按基因名重排
    m = dict(zip(var_genes, ctrl_mean))
    ctrl_aligned = np.array([m[g] for g in pred_genes])

pred_df = pd.read_csv(f"{OUT}/predictions_seed{SEED}.csv").set_index(
    "perturbation")
# 顺序对齐到 npz（双保险；两文件当前同序）
pred_df = pred_df.loc[perts]

# ---- 2. 特征矩阵 ----
feat = compute_features(pred_expr, perts, pred_genes, ctrl_aligned)
feat["log_n_cells"] = np.log1p(pred_df["n_cells"].values)

excluded = {}
if pred_df["in_pert_graph"].nunique() <= 1:      # 常数无方差，剔除留痕
    excluded["in_pert_graph"] = (
        f"常数列（全为 {bool(pred_df['in_pert_graph'].iloc[0])}），无方差剔除")
else:
    feat["in_pert_graph"] = pred_df["in_pert_graph"].astype(int).values
excluded["strength"] = (
    "真值派生（step1_run.py:96-105 由留出组真实表达算出），防泄漏排除")
excluded["pearson_delta_all"] = "真值派生（标签来源），防泄漏排除"

feature_names = SPECTRUM_FEATURES + \
    [c for c in COVARIATES if c in feat.columns]
X = feat[feature_names].values.astype(np.float64)
assert not np.isnan(X).any(), "特征矩阵含 NaN"
print(f"[step6] n={X.shape[0]}, 特征={len(feature_names)}: {feature_names}")

# ---- 3. 标签：error 高于中位数 = "答错"（定义性选择）----
error = 1.0 - pred_df["pearson_delta"].values
med = float(np.median(error))
y = (error > med).astype(int)
print(f"[step6] error 中位数={med:.4f}，正类（答错）={y.sum()}/{len(y)}")

# ---- 4. 探针 ----
res = run_probe(X, y, feature_names, seed=SEED)
print(f"[step6] CV AUROC = {res['auroc_mean']:.4f} ± {res['auroc_std']:.4f}"
      f"（pooled={res['pooled_auroc']:.4f}）")
best_single = max(res["single_feature"],
                  key=lambda k: res["single_feature"][k]["auroc_mean"])
print(f"[step6] 最好单特征: {best_single} = "
      f"{res['single_feature'][best_single]['auroc_mean']:.4f}")
for name, s in sorted(res["single_feature"].items(),
                      key=lambda kv: -kv[1]["auroc_mean"]):
    print(f"    {name:18s} {s['auroc_mean']:.4f} ± {s['auroc_std']:.4f}")

# ---- 5. 产物：JSON + ROC 图 ----
json.dump({
    "task": "gating probe (slot 3: reliability gating)",
    "seed": SEED,
    "n_perts": int(len(y)),
    "label_definition": "y = (1 - pearson_delta) > median; "
                        "error 高于中位数即判'答错'（定义性选择）",
    "features": feature_names,
    "excluded_features": excluded,
    "model": "LogisticRegression(L2, C=1.0) + StandardScaler, "
             "stratified 5-fold CV (shuffle, seed=0)",
    "criterion": "CV AUROC mean > 0.6 -> gating signal present "
                 "(0.5 = chance; definitional threshold)",
    "error_median": med,
    "n_positive": int(y.sum()),
    "cv_auroc_mean": res["auroc_mean"],
    "cv_auroc_std": res["auroc_std"],
    "cv_auroc_per_fold": res["auroc_per_fold"],
    "pooled_cv_auroc": res["pooled_auroc"],
    "best_single_feature": best_single,
    "single_feature_auroc": res["single_feature"],
    "verdict_pass": bool(res["auroc_mean"] > 0.6),
}, open(f"{OUT}/gating_probe.json", "w", encoding="utf-8"),
    indent=2, ensure_ascii=False)

fig, ax = plt.subplots(figsize=(6, 4.5))
ax.plot(res["roc_fpr"], res["roc_tpr"],
        label=f"pooled ROC (AUROC={res['pooled_auroc']:.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=1, label="chance (0.5)")
ax.set_xlabel("false positive rate")
ax.set_ylabel("true positive rate")
ax.set_title(f"Gating probe 5-fold CV: AUROC = {res['auroc_mean']:.3f} "
             f"± {res['auroc_std']:.3f}, n={len(y)}")
ax.legend()
fig.tight_layout()
fig.savefig(f"{OUT}/gating_probe_cv.png", dpi=300)
plt.close(fig)

# ---- 6. verdict ----
passed = res["auroc_mean"] > 0.6
lines = []
lines.append("# gating_verdict.md — 空位 3：可靠性门控（Day5，2026-10-13）\n")
lines.append("## 问题\n")
lines.append("GEARS 对某个留出扰动会不会答错，能否用**模型自身的廉价特征**"
             "提前门控？（CoT benchmark ECE 置信度探针思路的迁移检验）\n")
lines.append("## 方法\n")
lines.append(f"- 数据：seed {SEED}，{len(y)} 个留出扰动（predictions_seed{SEED}.csv "
             f"+ pred_expr_seed{SEED}.npz）。")
lines.append("- 标签：**error = 1 − pearson_delta**（top-50 DEG 口径）；"
             f"**y = error > 中位数（={med:.4f}）即判「答错」**——定义性选择，"
             "对应「预测质量落在全体留出扰动后一半」。")
lines.append("- 特征（仅模型侧）：预测谱 8 项（delta_l2 / delta_l1 / "
             "delta_max / delta_top50_mean / n_de_025 / frac_de_025 / "
             "neg_frac / delta_entropy，delta = 预测 − ctrl 训练均值）"
             "+ log_n_cells。")
lines.append("- 排除（泄漏防护）：**strength 由留出组真实表达算出"
             "（step1_run.py:96-105），属真值派生，已排除**；"
             "pearson_delta_all 同理；in_pert_graph 为常数列自动剔除。")
lines.append("- 探针：L2 logistic 回归（C=1.0），fold 内标准化，"
             "分层 5-fold CV（shuffle, seed=0），纯 CPU。\n")
lines.append("## 结果\n")
lines.append(f"- **5-fold CV AUROC = {res['auroc_mean']:.4f} ± "
             f"{res['auroc_std']:.4f}**（逐 fold: "
             + ", ".join(f"{a:.3f}" for a in res["auroc_per_fold"])
             + f"；pooled = {res['pooled_auroc']:.4f}）")
lines.append("- 单特征消融（CV AUROC 均值）：")
for name, s in sorted(res["single_feature"].items(),
                      key=lambda kv: -kv[1]["auroc_mean"]):
    lines.append(f"    - {name}: {s['auroc_mean']:.4f} ± {s['auroc_std']:.4f}")
lines.append(f"- 最好单特征：**{best_single}** = "
             f"{res['single_feature'][best_single]['auroc_mean']:.4f}\n")
lines.append("## 判据（定义性阈值，写死）\n")
lines.append("**CV AUROC 均值 > 0.6 → 「门控信号存在」**"
             "（0.5 = 随机；与 CoT benchmark 中 ECE 探针达标线同一量级，"
             "此为项目内约定阈值，非外部标准）。\n")
lines.append(f"## 结论\n")
if passed:
    lines.append(f"**AUROC = {res['auroc_mean']:.4f} > 0.6 → 门控信号存在。** "
                 f"模型自身的廉价特征（以 {best_single} 为主）能显著优于随机地"
                 f"提前识别 GEARS 会答错的留出扰动——**支持空位 3 成立**："
                 f"可靠性门控可作为虚拟细胞预测的部署前过滤层。")
else:
    lines.append(f"**AUROC = {res['auroc_mean']:.4f} ≤ 0.6 → 门控信号不存在（未达标）。** "
                 f"模型侧廉价特征无法以可用精度提前识别 GEARS 的答错扰动——"
                 f"**本轮不支持空位 3**；不升级为部署门控，按阴性归档。")
lines.append("\n## 对照与后续\n")
lines.append("- GEARS 原生 uncertainty 列全 NaN（单种子无 ensemble 方差）。"
             "**种子 1/2 训练完成后，用 3 种子预测方差补全 uncertainty，"
             "重跑本脚本做对照**——若原生不确定性有效，应优于或互补本探针。")
lines.append("- 判阴性/阳性的稳健性：可换 pearson_delta_all 口径重跑标签，"
             "看中位数二分结论是否翻转。\n")
lines.append("## 复现\n")
lines.append("- 脚本：step6_gating_probe.py（特征 + 探针核心），"
             "step6_run.py（运行器/产物）")
lines.append("- 产物：outputs/gating_probe.json, gating_probe_cv.png, "
             "本文件")
open(f"{OUT}/gating_verdict.md", "w", encoding="utf-8").write("\n".join(lines))
print(f"[step6] -> {OUT}/gating_probe.json, gating_probe_cv.png, "
      f"gating_verdict.md  (pass={passed})")
