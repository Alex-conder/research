"""
GSE227651 部分样本分析（目前 1d + 3d 完整；7d/sham 下载中）
"""
import os, glob, gzip
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = "../01_raw_data/GSE227651"
OUT_DIR  = "../04_reports/figures/GSE227651"
os.makedirs(OUT_DIR, exist_ok=True)

sample_dirs = sorted([d for d in glob.glob(os.path.join(DATA_DIR, "*")) if os.path.isdir(d)])
print(f"[INFO] Found {len(sample_dirs)} sample dirs")

adatas = []
for d in sample_dirs:
    sname = os.path.basename(d)
    mtx = os.path.join(d, "matrix.mtx.gz")
    # 检查 gzip 完整性
    try:
        with gzip.open(mtx, 'rb') as f:
            f.read(1)
    except Exception as e:
        print(f"[SKIP] {sname}: matrix corrupted ({e})")
        continue
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

if len(adatas) == 0:
    print("[WARN] No valid samples")
    exit(0)

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
pirb_stats.to_csv(os.path.join(OUT_DIR, "pirb_partial_time_summary.csv"), index=False)
print("\nPirb expression over time (partial):")
print(pirb_stats)

# 可视化
fig, ax = plt.subplots(figsize=(8, 5))
order = sorted(adata.obs["day"].unique(), key=lambda x: (x != "sham", x))
sc.pl.violin(adata, keys=["Pirb_counts"], groupby="day", rotation=45, show=False, ax=ax)
ax.set_title("Pirb expression across available timepoints")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pirb_violin_partial_time.png"), dpi=150)
plt.close()

adata.obs.to_csv(os.path.join(OUT_DIR, "cell_metadata_partial.csv"))
adata.write_h5ad(os.path.join(OUT_DIR, "GSE227651_partial_qc.h5ad"))
print(f"[DONE] Results saved to {OUT_DIR}")
