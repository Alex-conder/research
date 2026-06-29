"""
GSE227651 单细胞数据质控与 Pirb 时间序列表达分析
- 4 个样本：1d, 3d, 7d, sham
- 10x 标准格式，样本位于子目录
"""
import os, glob
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = "D:/Pirb_stroke_project/01_raw_data/GSE227651"
OUT_DIR  = "D:/Pirb_stroke_project/04_reports/figures/GSE227651"
os.makedirs(OUT_DIR, exist_ok=True)

sample_dirs = sorted([d for d in glob.glob(os.path.join(DATA_DIR, "*")) if os.path.isdir(d)])
print(f"[INFO] Found {len(sample_dirs)} samples")

adatas = []
for d in sample_dirs:
    sname = os.path.basename(d)
    # 标准 10x 子目录
    ad = sc.read_10x_mtx(d, var_names="gene_symbols", cache=False)
    ad.obs["sample"] = sname
    if "sham" in sname.lower():
        ad.obs["condition"] = "Sham"
        ad.obs["day"] = "sham"
    else:
        ad.obs["condition"] = "MCAO"
        m = sname.split("_")[-1]
        ad.obs["day"] = m
    adatas.append(ad)
    print(f"  {sname}: {ad.n_obs} cells × {ad.n_vars} genes")

adata = sc.concat(adatas, label="batch", keys=[a.obs["sample"][0] for a in adatas], index_unique="-")
adata.obs_names_make_unique()
print(f"[INFO] Merged: {adata.n_obs} cells × {adata.n_vars} genes")

# 质控
adata.var["mt"] = adata.var_names.str.startswith("mt-")
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo"], percent_top=None, log1p=False, inplace=True)

n_before = adata.n_obs
adata = adata[adata.obs.n_genes_by_counts >= 500]
adata = adata[adata.obs.total_counts >= 1000]
adata = adata[adata.obs.pct_counts_mt < 20]
n_after = adata.n_obs
print(f"[INFO] QC: {n_before} -> {n_after} cells")

qc = adata.obs.groupby("sample").agg(
    condition=("condition", "first"),
    day=("day", "first"),
    n_cells=("sample", "size"),
    median_n_genes=("n_genes_by_counts", "median"),
    median_total_counts=("total_counts", "median"),
    median_pct_mt=("pct_counts_mt", "median")
).reset_index()
qc.to_csv(os.path.join(OUT_DIR, "qc_summary.csv"), index=False)
print(qc)

# 标准化
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

# Pirb
pirb = "Pirb"
adata.obs["Pirb_counts"] = adata[:, pirb].X.toarray().flatten()
adata.obs["Pirb_detected"] = adata.obs["Pirb_counts"] > 0

pirb_stats = adata.obs.groupby(["day", "condition"]).agg(
    n_cells=("day", "size"),
    pirb_positive=("Pirb_detected", "sum"),
    pirb_positive_pct=("Pirb_detected", "mean"),
    pirb_mean_counts=("Pirb_counts", "mean"),
    pirb_median_counts=("Pirb_counts", "median")
).reset_index()
pirb_stats["pirb_positive_pct"] *= 100
pirb_stats.to_csv(os.path.join(OUT_DIR, "pirb_time_summary.csv"), index=False)
print("\nPirb expression over time:")
print(pirb_stats)

# 小提琴图
fig, ax = plt.subplots(figsize=(8, 5))
order = ["sham", "1d", "3d", "7d"]
sc.pl.violin(adata, keys=["Pirb_counts"], groupby="day", rotation=45, show=False, ax=ax)
ax.set_title("Pirb expression across timepoints")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pirb_violin_time.png"), dpi=150)
plt.close()

# 降维聚类
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata_hvg = adata[:, adata.var.highly_variable].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, svd_solver="arpack")
sc.pp.neighbors(adata_hvg, n_neighbors=10, n_pcs=30)
sc.tl.umap(adata_hvg)
sc.tl.leiden(adata_hvg, resolution=0.5)

adata.obs["leiden"] = adata_hvg.obs["leiden"].values
for key in ["X_pca", "X_umap"]:
    adata.obsm[key] = adata_hvg.obsm[key]
adata.obs.to_csv(os.path.join(OUT_DIR, "cell_metadata.csv"))
adata.write_h5ad(os.path.join(OUT_DIR, "GSE227651_qc.h5ad"))

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.umap(adata_hvg, color="leiden", ax=axes[0], show=False, title="Leiden clusters")
sc.pl.umap(adata_hvg, color="day", ax=axes[1], show=False, title="Timepoint")
sc.pl.umap(adata_hvg, color="Pirb_counts", ax=axes[2], show=False, title="Pirb expression")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "umap_overview.png"), dpi=150)
plt.close()

print(f"[DONE] Results saved to {OUT_DIR}")
