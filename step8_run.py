"""Day5 Step8 运行器：门控探针种子稳定性检查（种子 2 落地收尾）。

对 seed ∈ {0,1,2} 分别用 outputs/predictions_seed{s}.csv +
outputs/pred_expr_seed{s}.npz 重跑 step6 完整门控流程
（特征构造 → L2 logistic 探针 → 分层 5-fold CV AUROC），检验
step6 信号（seed0 AUROC 0.6556）对模型初始化是否稳健。

不变量检查（[修复 D5] 预期）：三种子留出扰动集合必须完全一致
（predictions_seed{s}.csv 第一列全等），不等即报错退出。

稳定性判据（三条，层级：存在性 > 方向 > 量级；量级不过不否定存在性，
详见 step8_seed_stability.py docstring）：
- 方向稳健：三种子 AUROC 均 > 0.6 且 best 单特征一致；
- 量级稳健：三种子 AUROC 极差 < 0.05（阈值依据：≈ step6 fold std 0.048；
  n=217 小样本 fold 级 std 约 0.05-0.06，极差落在该区间优先归因小样本）；
- 信号存在：至少两粒 > 0.6 且三粒均值 > 0.6。
组合结论写入 outputs/seed_stability.json 的 stability.verdict。

文件缺失：任一 predictions_seed{s}.csv / pred_expr_seed{s}.npz 不存在
（种子仍在训练）→ 打印缺失清单、退出码 2，不产出任何文件。

用法：
    source venv/Scripts/activate
    python step8_run.py
输出：outputs/seed_stability.json
"""
import json
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd

from step8_seed_stability import gating_for_seed, stability_summary

OUT = "outputs"
SEEDS = (0, 1, 2)

# ---- 0. 输入齐备性检查 ----
need = [f"{OUT}/predictions_seed{s}.csv" for s in SEEDS] + \
       [f"{OUT}/pred_expr_seed{s}.npz" for s in SEEDS]
missing = [f for f in need if not os.path.exists(f)]
if missing:
    print("[step8] 以下文件尚不存在（种子可能仍在训练）：")
    for f in missing:
        print(f"    - {f}")
    print("[step8] 待训练完成后重跑：python step8_run.py")
    sys.exit(2)

# ---- 1. 读三种子预测 + ctrl 训练均值（共同基线，与 step6 同口径）----
per_seed = {}
pert_sets = {}
for s in SEEDS:
    z = np.load(f"{OUT}/pred_expr_seed{s}.npz")
    perts = z["pred_perturbations"].tolist()
    pred_genes = z["pred_genes"].tolist()
    pdf = pd.read_csv(f"{OUT}/predictions_seed{s}.csv")
    perts_csv = pdf["perturbation"].tolist()
    pert_sets[s] = perts_csv
    assert set(perts_csv) == set(perts), \
        f"seed{s}: CSV 与 npz 扰动集合不一致"
    per_seed[s] = (z["pred_expr"], perts, pred_genes,
                   pdf.set_index("perturbation"))
    print(f"[step8] seed{s}: {len(perts)} 个留出扰动, "
          f"pred_expr {z['pred_expr'].shape}")

# ---- 2. 不变量：三种子留出扰动集合全等（[修复 D5] 预期）----
ref = pert_sets[0]
for s in SEEDS[1:]:
    if pert_sets[s] != ref:                      # 顺序级全等（同考场同划分）
        raise SystemExit(
            f"[step8] 不变量违反：seed{s} 扰动列表与 seed0 不完全一致"
            f"（[修复 D5] 后 SPLIT_SEED=0 应完全一致），拒绝继续。")
print(f"[step8] 不变量通过：3 种子扰动集合完全一致（n={len(ref)}）")

# backed 模式只物化 ctrl 子集：全量 X ~GB 级，训练并行期内存紧张会 OOM
adata = ad.read_h5ad(
    "data/gears_demo/replogle_k562_essential/perturb_processed.h5ad",
    backed="r")
var_genes = adata.var["gene_name"].tolist()
ctrl_mean = np.asarray(
    adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
adata.file.close()
del adata
if var_genes != pred_genes:                      # 防御：按基因名重排（同 step5-7）
    m = dict(zip(var_genes, ctrl_mean))
    ctrl_mean = np.array([m[g] for g in pred_genes])

# ---- 3. 逐种子重跑 step6 流程 ----
results = {}
for s in SEEDS:
    pred_expr, perts, _, pdf = per_seed[s]
    results[s] = gating_for_seed(pred_expr, perts, pred_genes, ctrl_mean, pdf)
    r = results[s]
    print(f"[step8] seed{s}: AUROC = {r['cv_auroc_mean']:.4f} ± "
          f"{r['cv_auroc_std']:.4f}（pooled={r['pooled_cv_auroc']:.4f}, "
          f"best={r['best_single_feature']}）")

# ---- 4. 稳定性判据 ----
summ = stability_summary(results)
cr = summ["criteria"]
print(f"[step8] 判据：存在 {'通过' if cr['existence_at_least_2_above_0.6_and_mean_above_0.6'] else '未过'}"
      f"｜方向 {'通过' if cr['direction_all_above_0.6_and_same_best_feature'] else '未过'}"
      f"（best: {summ['best_single_feature_by_seed']}）"
      f"｜量级 {'通过' if cr['magnitude_range_below_0.05'] else '未过'}"
      f"（极差 {summ['range']:.4f}）")
print(f"[step8] 组合结论：{summ['verdict']}")

json.dump({
    "task": "gating probe seed stability (step6 signal replication)",
    "n_perts": len(ref),
    "perturbation_sets_identical": True,
    "per_seed": {str(s): results[s] for s in SEEDS},
    "stability": summ,          # 含 criteria（三判据逐项 true/false + 阈值 +
                                # 量级归因注释）与 verdict 组合结论
}, open(f"{OUT}/seed_stability.json", "w", encoding="utf-8"),
    indent=2, ensure_ascii=False)
print(f"[step8] -> {OUT}/seed_stability.json")
