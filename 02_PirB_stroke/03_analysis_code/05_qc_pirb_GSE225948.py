"""
GSE225948 Brain + Blood 单细胞数据 Pirb 表达分析
- 29 个样本：CSV counts + metadata（已标准化）
- 读取并合并，按 treatment / tissue / parent cell type 分析 Pirb
"""
import os, glob
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

DATA_DIR = "../01_raw_data/GSE225948"
OUT_DIR  = "../04_reports/figures/GSE225948"
os.makedirs(OUT_DIR, exist_ok=True)

count_files = sorted(glob.glob(os.path.join(DATA_DIR, "*_counts.csv.gz")))
print(f"[INFO] Found {len(count_files)} count files")

adatas = []
for cf in count_files:
    sname = os.path.basename(cf).replace("_counts.csv.gz", "")
    meta_file = cf.replace("_counts.csv.gz", "_metadata.csv.gz")
    
    # 读取 counts：行为基因，列为细胞
    df = pd.read_csv(cf, index_col=0)
    # 转置为细胞×基因
    X = df.T.values
    # 处理重复基因名
    var = pd.DataFrame(index=df.index)
    var.index = var.index.astype(str)
    obs = pd.DataFrame(index=df.columns)
    ad = sc.AnnData(X=X, obs=obs, var=var)
    ad.var_names_make_unique()
    
    # 读取 metadata
    meta = pd.read_csv(meta_file, index_col=0)
    meta.index = meta.index.astype(str)
    # 确保索引一致
    ad = ad[ad.obs.index.isin(meta.index)].copy()
    ad.obs = ad.obs.join(meta)
    
    ad.obs["sample"] = sname
    ad.obs["tissue_raw"] = sname.split("_")[1]  # Brain / Blood
    ad.obs_names_make_unique()
    adatas.append(ad)
    print(f"  {sname}: {ad.n_obs} cells × {ad.n_vars} genes")

adata = sc.concat(adatas, label="batch", keys=[a.obs["sample"][0] for a in adatas], index_unique="-")
adata.obs_names_make_unique()
print(f"[INFO] Merged: {adata.n_obs} cells × {adata.n_vars} genes")

# 数据已经是标准化后的，直接用于 Pirb 分析
pirb = "Pirb"
if pirb not in adata.var_names:
    print(f"[WARN] {pirb} not found")
else:
    adata.obs["Pirb_counts"] = adata[:, pirb].X.toarray().flatten()
    adata.obs["Pirb_detected"] = adata.obs["Pirb_counts"] > 0

# metadata 列名去掉引号空格
adata.obs.columns = [c.strip('"') for c in adata.obs.columns]

# 按 treatment 和 tissue 统计
for grp in ["treatment", "tissue_raw", "parent", "sub.celltype"]:
    if grp not in adata.obs.columns:
        continue
    stats = adata.obs.groupby(grp).agg(
        n_cells=(grp, "size"),
        pirb_positive=("Pirb_detected", "sum"),
        pirb_positive_pct=("Pirb_detected", "mean"),
        pirb_mean_counts=("Pirb_counts", "mean"),
    ).reset_index()
    stats["pirb_positive_pct"] *= 100
    stats.to_csv(os.path.join(OUT_DIR, f"pirb_by_{grp}.csv"), index=False)
    print(f"\nPirb by {grp}:")
    print(stats)

# 按 treatment + tissue
stats2 = adata.obs.groupby(["treatment", "tissue_raw"]).agg(
    n_cells=("treatment", "size"),
    pirb_positive=("Pirb_detected", "sum"),
    pirb_positive_pct=("Pirb_detected", "mean"),
    pirb_mean_counts=("Pirb_counts", "mean"),
).reset_index()
stats2["pirb_positive_pct"] *= 100
stats2.to_csv(os.path.join(OUT_DIR, "pirb_by_treatment_tissue.csv"), index=False)
print("\nPirb by treatment + tissue:")
print(stats2)

# 简单可视化
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
sc.pl.violin(adata, keys=["Pirb_counts"], groupby="treatment", ax=axes[0], show=False)
axes[0].set_title("Pirb by treatment")
sc.pl.violin(adata, keys=["Pirb_counts"], groupby="tissue_raw", ax=axes[1], show=False)
axes[1].set_title("Pirb by tissue")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pirb_violin.png"), dpi=150)
plt.close()

# 保存
adata.write_h5ad(os.path.join(OUT_DIR, "GSE225948_merged.h5ad"))
adata.obs.to_csv(os.path.join(OUT_DIR, "cell_metadata.csv"))

print(f"[DONE] Results saved to {OUT_DIR}")
