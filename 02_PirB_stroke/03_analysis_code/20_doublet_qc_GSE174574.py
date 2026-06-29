"""
GSE174574 Doublet 检测与 Pirb+ 星形胶质细胞纯度评估
- 使用 Scrublet 预测 doublet
- 重新计算 Pirb+ 比例，排除潜在 doublet
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import scrublet as scr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

IN_H5AD = "../04_reports/figures/GSE174574/GSE174574_annotated.h5ad"
OUT_DIR = "../04_reports/figures/GSE174574/doublet_qc"
os.makedirs(OUT_DIR, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)
print(f"[INFO] Loaded: {adata.n_obs} cells")

# 按样本分别做 Scrublet（推荐）
adata.obs["doublet_score"] = 0.0
adata.obs["predicted_doublet"] = False

for sample in adata.obs["sample"].unique():
    mask = adata.obs["sample"] == sample
    sub = adata[mask].copy()
    # 使用原始 counts 或 normalize_total 后的数据；Scrublet 期望 counts
    counts = sub.X.copy()
    scrub = scr.Scrublet(counts, expected_doublet_rate=0.06)
    doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    adata.obs.loc[mask, "doublet_score"] = doublet_scores
    adata.obs.loc[mask, "predicted_doublet"] = predicted_doublets
    print(f"  {sample}: {predicted_doublets.sum()} / {len(predicted_doublets)} predicted doublets")

# 保存
adata.obs.to_csv(os.path.join(OUT_DIR, "cell_metadata_with_doublet.csv"))
adata.write_h5ad(os.path.join(OUT_DIR, "GSE174574_with_doublet.h5ad"))

# 总体统计
total_doublets = adata.obs["predicted_doublet"].sum()
print(f"\n[INFO] Total predicted doublets: {total_doublets} / {adata.n_obs} ({total_doublets/adata.n_obs*100:.2f}%)")

# 按细胞类型统计 doublet 率
doublet_by_ct = adata.obs.groupby("cell_type").agg(
    n_cells=("cell_type", "size"),
    n_doublets=("predicted_doublet", "sum"),
    doublet_rate=("predicted_doublet", "mean")
).reset_index()
doublet_by_ct["doublet_rate_pct"] = doublet_by_ct["doublet_rate"] * 100
doublet_by_ct.to_csv(os.path.join(OUT_DIR, "doublet_rate_by_celltype.csv"), index=False)
print("\nDoublet rate by cell type:")
print(doublet_by_ct)

# Pirb+ 星形胶质细胞中 doublet 比例
astro = adata[adata.obs["cell_type"] == "Astrocyte"]
pirb_pos_astro = astro[astro.obs["Pirb_detected"] == 1]
pirb_neg_astro = astro[astro.obs["Pirb_detected"] == 0]

print(f"\n[Astrocyte Pirb+] n={len(pirb_pos_astro)}, doublets={pirb_pos_astro.obs['predicted_doublet'].sum()} ({pirb_pos_astro.obs['predicted_doublet'].mean()*100:.2f}%)")
print(f"[Astrocyte Pirb-] n={len(pirb_neg_astro)}, doublets={pirb_neg_astro.obs['predicted_doublet'].sum()} ({pirb_neg_astro.obs['predicted_doublet'].mean()*100:.2f}%)")

# 排除 doublet 后重新计算 Pirb+ 比例
adata_nodbl = adata[~adata.obs["predicted_doublet"]].copy()
pirb_stats_nodbl = adata_nodbl.obs.groupby(["cell_type", "condition"]).agg(
    n_cells=("cell_type", "size"),
    pirb_positive=("Pirb_detected", "sum"),
    pirb_positive_pct=("Pirb_detected", "mean"),
).reset_index()
pirb_stats_nodbl["pirb_positive_pct"] *= 100
pirb_stats_nodbl.to_csv(os.path.join(OUT_DIR, "pirb_by_celltype_condition_no_doublet.csv"), index=False)
print("\nPirb+ % after removing predicted doublets:")
print(pirb_stats_nodbl)

# 与原始对比
original = pd.read_csv("../04_reports/figures/GSE174574/pirb_by_celltype_condition.csv")
compare = original.merge(pirb_stats_nodbl, on=["cell_type", "condition"], suffixes=("_orig", "_nodbl"))
compare["pct_diff"] = compare["pirb_positive_pct_nodbl"] - compare["pirb_positive_pct_orig"]
compare.to_csv(os.path.join(OUT_DIR, "pirb_comparison_doublet_removal.csv"), index=False)
print("\nComparison (original vs no-doublet):")
print(compare[["cell_type", "condition", "pirb_positive_pct_orig", "pirb_positive_pct_nodbl", "pct_diff"]])

# 可视化
fig, axes = plt.subplots(1, 2, figsize=(14, 5))
sns.barplot(data=doublet_by_ct, x="cell_type", y="doublet_rate_pct", ax=axes[0])
axes[0].set_title("Predicted doublet rate by cell type")
axes[0].set_ylabel("Doublet rate (%)")
axes[0].tick_params(axis='x', rotation=45)

# Pirb+ astro doublet score 分布
sns.histplot(pirb_pos_astro.obs["doublet_score"], bins=50, kde=True, color="red", label="Pirb+", ax=axes[1])
sns.histplot(pirb_neg_astro.obs["doublet_score"], bins=50, kde=True, color="blue", label="Pirb-", ax=axes[1])
axes[1].set_title("Doublet score distribution in astrocytes")
axes[1].set_xlabel("Scrublet doublet score")
axes[1].legend()
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "doublet_qc_overview.png"), dpi=150)
plt.close()

print(f"[DONE] Results saved to {OUT_DIR}")
