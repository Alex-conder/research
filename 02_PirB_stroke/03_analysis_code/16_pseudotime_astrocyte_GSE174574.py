"""
GSE174574 星形胶质细胞伪时间轨迹分析（scanpy DPT）
- 以 Sham 星形胶质细胞为根，观察 Pirb 表达沿伪时间变化
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

IN_H5AD = "D:/Pirb_stroke_project/04_reports/figures/GSE174574/GSE174574_annotated.h5ad"
OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE174574/pseudotime"
os.makedirs(OUT_DIR, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)
astro = adata[adata.obs["cell_type"] == "Astrocyte"].copy()
print(f"[INFO] Astrocytes: {astro.n_obs}")

# 重新高变基因、PCA、邻域、扩散图（仅星形胶质细胞）
sc.pp.highly_variable_genes(astro, min_mean=0.0125, max_mean=3, min_disp=0.5)
astro_hvg = astro[:, astro.var.highly_variable].copy()
sc.pp.scale(astro_hvg, max_value=10)
sc.tl.pca(astro_hvg, svd_solver="arpack")
sc.pp.neighbors(astro_hvg, n_neighbors=15, n_pcs=30)
sc.tl.diffmap(astro_hvg)

# 以 Sham 星形胶质细胞为根
sham_idx = np.where(astro_hvg.obs["condition"] == "Sham")[0]
astro_hvg.uns["iroot"] = sham_idx[0]
sc.tl.dpt(astro_hvg)

# 将坐标和伪时间写回 astro
astro.obs["dpt_pseudotime"] = astro_hvg.obs["dpt_pseudotime"].values
for key in ["X_pca", "X_diffmap"]:
    astro.obsm[key] = astro_hvg.obsm[key]

# Pirb 表达
astro.obs["Pirb_counts"] = astro[:, "Pirb"].X.toarray().flatten()

# 保存
astro.obs.to_csv(os.path.join(OUT_DIR, "astrocyte_pseudotime_metadata.csv"))
astro.write_h5ad(os.path.join(OUT_DIR, "astrocyte_pseudotime.h5ad"))

# 可视化
fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.diffmap(astro_hvg, color="dpt_pseudotime", ax=axes[0], show=False, title="Diffusion pseudotime")
sc.pl.diffmap(astro_hvg, color="condition", ax=axes[1], show=False, title="Condition")
sc.pl.diffmap(astro_hvg, color="Pirb_counts", ax=axes[2], show=False, title="Pirb expression")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "diffmap_overview.png"), dpi=150)
plt.close()

# 伪时间 vs Pirb 表达散点图 + 平滑曲线
plot_df = astro.obs[["dpt_pseudotime", "Pirb_counts", "condition"]].copy()
plt.figure(figsize=(10, 6))
sns.scatterplot(data=plot_df, x="dpt_pseudotime", y="Pirb_counts", hue="condition", alpha=0.3, s=10)
# 局部平滑
plot_df_sorted = plot_df.sort_values("dpt_pseudotime")
window = max(50, int(len(plot_df_sorted) * 0.05))
plot_df_sorted["smooth"] = plot_df_sorted["Pirb_counts"].rolling(window=window, center=True, min_periods=1).mean()
sns.lineplot(data=plot_df_sorted, x="dpt_pseudotime", y="smooth", color="black", linewidth=2, label="Smoothed")
plt.xlabel("Diffusion pseudotime")
plt.ylabel("Pirb expression (log1p norm)")
plt.title("Pirb expression along astrocyte pseudotime")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pirb_vs_pseudotime.png"), dpi=150)
plt.close()

# 分 bins 统计 Pirb 阳性率
plot_df["pt_bin"] = pd.qcut(plot_df["dpt_pseudotime"], q=10, duplicates="drop")
bin_stats = plot_df.groupby("pt_bin").agg(
    mean_pirb=("Pirb_counts", "mean"),
    pct_pirb_pos=("Pirb_counts", lambda x: (x > 0).mean() * 100),
    mid_pt=("dpt_pseudotime", "mean")
).reset_index()
bin_stats.to_csv(os.path.join(OUT_DIR, "pirb_by_pseudotime_bin.csv"), index=False)

fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sns.lineplot(data=bin_stats, x="mid_pt", y="mean_pirb", ax=axes[0], marker="o")
axes[0].set_title("Mean Pirb expression by pseudotime bin")
axes[0].set_xlabel("Pseudotime")
sns.lineplot(data=bin_stats, x="mid_pt", y="pct_pirb_pos", ax=axes[1], marker="o", color="red")
axes[1].set_title("Pirb+ % by pseudotime bin")
axes[1].set_xlabel("Pseudotime")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pirb_trend_pseudotime.png"), dpi=150)
plt.close()

print(f"[DONE] Results saved to {OUT_DIR}")
