"""
GSE171169 CD45high 免疫细胞数据质控与 Pirb 表达分析
- 4 个样本：5d_N1, 5d_N2, 14d_N1, 14d_N2
- 文件平铺，legacy genes.tsv.gz 命名
"""
import os, glob
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = "D:/Pirb_stroke_project/01_raw_data/GSE171169"
OUT_DIR  = "D:/Pirb_stroke_project/04_reports/figures/GSE171169"
os.makedirs(OUT_DIR, exist_ok=True)

barcode_files = sorted(glob.glob(os.path.join(DATA_DIR, "*_barcodes.tsv.gz")))
prefixes = sorted(set(os.path.basename(f).replace("_barcodes.tsv.gz", "") for f in barcode_files))
print(f"[INFO] Found {len(prefixes)} samples")

adatas = []
for prefix in prefixes:
    bc_file = os.path.join(DATA_DIR, f"{prefix}_barcodes.tsv.gz")
    gene_file = os.path.join(DATA_DIR, f"{prefix}_genes.tsv.gz")
    mtx_file = os.path.join(DATA_DIR, f"{prefix}_matrix.mtx.gz")
    barcodes = pd.read_csv(bc_file, header=None, sep="\t")[0].values
    genes = pd.read_csv(gene_file, header=None, sep="\t")
    gene_ids = genes[0].values
    gene_symbols = genes[1].values if genes.shape[1] > 1 else gene_ids
    X = mmread(mtx_file).T.tocsr()
    ad = sc.AnnData(X=X, obs=pd.DataFrame(index=barcodes), var=pd.DataFrame(index=gene_symbols))
    ad.var_names_make_unique()
    ad.var["gene_ids"] = gene_ids
    ad.obs["sample"] = prefix
    if "5d" in prefix:
        ad.obs["day"] = "5d"
    elif "14d" in prefix:
        ad.obs["day"] = "14d"
    else:
        ad.obs["day"] = "unknown"
    ad.obs["condition"] = "MCAO"
    adatas.append(ad)
    print(f"  {prefix}: {ad.n_obs} cells × {ad.n_vars} genes")

adata = sc.concat(adatas, label="batch", keys=[a.obs["sample"][0] for a in adatas], index_unique="-")
adata.obs_names_make_unique()
print(f"[INFO] Merged: {adata.n_obs} cells × {adata.n_vars} genes")

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
    day=("day", "first"),
    n_cells=("sample", "size"),
    median_n_genes=("n_genes_by_counts", "median"),
    median_total_counts=("total_counts", "median"),
    median_pct_mt=("pct_counts_mt", "median")
).reset_index()
qc.to_csv(os.path.join(OUT_DIR, "qc_summary.csv"), index=False)
print(qc)

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

pirb = "Pirb"
adata.obs["Pirb_counts"] = adata[:, pirb].X.toarray().flatten()
adata.obs["Pirb_detected"] = adata.obs["Pirb_counts"] > 0

pirb_stats = adata.obs.groupby("day").agg(
    n_cells=("day", "size"),
    pirb_positive=("Pirb_detected", "sum"),
    pirb_positive_pct=("Pirb_detected", "mean"),
    pirb_mean_counts=("Pirb_counts", "mean"),
).reset_index()
pirb_stats["pirb_positive_pct"] *= 100
pirb_stats.to_csv(os.path.join(OUT_DIR, "pirb_summary.csv"), index=False)
print("\nPirb expression by day:")
print(pirb_stats)

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
adata.write_h5ad(os.path.join(OUT_DIR, "GSE171169_qc.h5ad"))

fig, axes = plt.subplots(1, 3, figsize=(18, 5))
sc.pl.umap(adata_hvg, color="leiden", ax=axes[0], show=False, title="Leiden clusters")
sc.pl.umap(adata_hvg, color="day", ax=axes[1], show=False, title="Day")
sc.pl.umap(adata_hvg, color="Pirb_counts", ax=axes[2], show=False, title="Pirb expression")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "umap_overview.png"), dpi=150)
plt.close()

print(f"[DONE] Results saved to {OUT_DIR}")
