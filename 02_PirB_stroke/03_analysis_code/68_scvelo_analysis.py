"""
基于 velocyto 生成的 loom 文件运行 scVelo RNA 速率分析。
需要前置步骤：在 WSL/Linux 中运行 01_raw_data/GSE233812_velocity/scripts/run_velocity_pipeline.sh
生成合并后的 loom 文件。
"""
import os
import sys
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# 尝试导入 scvelo
try:
    import scvelo as scv
except ImportError:
    print("[ERROR] scvelo not installed. Run: pip install scvelo")
    sys.exit(1)

# 路径
PROJECT_ROOT = ".."
VELO_DIR = os.path.join(PROJECT_ROOT, "01_raw_data/GSE233812_velocity/velocity")
LOOM_PATH = os.path.join(VELO_DIR, "GSE233812_merged.loom")
H5AD_IN = os.path.join(PROJECT_ROOT, "04_reports/figures/GSE233812_processed.h5ad")
OUT_DIR = os.path.join(PROJECT_ROOT, "04_reports/figures/rna_velocity")
H5AD_OUT = os.path.join(OUT_DIR, "GSE233812_velocity.h5ad")
os.makedirs(OUT_DIR, exist_ok=True)

if not os.path.exists(LOOM_PATH):
    print(f"[ERROR] Loom file not found: {LOOM_PATH}")
    print("Please run the WSL/Linux pipeline first:")
    print("  bash 01_raw_data/GSE233812_velocity/scripts/run_velocity_pipeline.sh")
    sys.exit(1)

print(f"[INFO] Loading loom: {LOOM_PATH}")
ldata = scv.read(LOOM_PATH, cache=True)
print(f"[INFO] Loom shape: {ldata.n_obs} cells × {ldata.n_vars} genes")

# 加载已有 h5ad
print(f"[INFO] Loading processed h5ad: {H5AD_IN}")
adata = sc.read_h5ad(H5AD_IN)
print(f"[INFO] H5AD shape: {adata.n_obs} cells × {adata.n_vars} genes")

# velocyto 的 cell barcode 通常是 ACCAGCAAGGTCTTCC-1-0 格式，需要匹配
# 常见做法：取 loom 的 obs_names 与 adata 的 obs_names 交集
print("[INFO] Matching cell barcodes...")
print(f"  Loom barcodes (first 3): {list(ldata.obs_names[:3])}")
print(f"  H5AD barcodes (first 3): {list(adata.obs_names[:3])}")

# 尝试清理 barcode 名称
# loom barcode 可能为 "sample:ACGT...-1" 或 "ACGT...-1-0"
# adata barcode 可能为 "ACGT...-1"
loom_bc_clean = pd.Series(ldata.obs_names).str.replace(r'-\d+$', '', regex=True)
adata_bc_clean = pd.Series(adata.obs_names).str.replace(r'-\d+$', '', regex=True)

# 如果 adata 有 sample/library 前缀，需要相应处理
common = set(loom_bc_clean) & set(adata_bc_clean)
print(f"[INFO] Common cells after cleaning: {len(common)}")

if len(common) == 0:
    print("[WARN] No common cells found with simple cleaning. Trying prefix matching...")
    # 备用方案：根据样本名手动匹配
    sys.exit(1)

# 创建 loom 的 clean index
ldata.obs['clean_bc'] = loom_bc_clean.values
adata.obs['clean_bc'] = adata_bc_clean.values

# 合并 spliced/unspliced 到 adata
# scvelo 期望 layers: spliced, unspliced
ldata.obs_names_make_unique()
adata.obs_names_make_unique()

# 使用 clean barcode 作为临时 key
ldata.obs_names = ldata.obs['clean_bc'].values
adata.obs_names = adata.obs['clean_bc'].values

# 取交集
common_list = list(common)
ldata_sub = ldata[common_list].copy()
adata_sub = adata[common_list].copy()

# 合并 layers
adata_sub.layers['spliced'] = ldata_sub.layers['spliced']
adata_sub.layers['unspliced'] = ldata_sub.layers['unspliced']
adata_sub.layers['ambiguous'] = ldata_sub.layers['ambiguous']

print(f"[INFO] Merged shape: {adata_sub.n_obs} cells × {adata_sub.n_vars} genes")
print(f"[INFO] Spliced layer non-zero: {(adata_sub.layers['spliced'] > 0).sum()}")
print(f"[INFO] Unspliced layer non-zero: {(adata_sub.layers['unspliced'] > 0).sum()}")

# scVelo 标准流程
scv.settings.figdir = OUT_DIR
scv.settings.set_figure_params('scvelo')

print("[INFO] scVelo preprocessing...")
scv.pp.filter_and_normalize(adata_sub, min_shared_counts=20, n_top_genes=2000)
scv.pp.moments(adata_sub, n_pcs=30, n_neighbors=30)

print("[INFO] Recovering dynamics...")
scv.tl.recover_dynamics(adata_sub, n_jobs=4)

print("[INFO] Computing velocity...")
scv.tl.velocity(adata_sub, mode='dynamical')
scv.tl.velocity_graph(adata_sub)

# 按时间点和 Pirb 表达可视化
print("[INFO] Plotting velocity streams...")
fig_kwargs = dict(basis='umap', color='time_point', show=False, dpi=300)
scv.pl.velocity_embedding_stream(adata_sub, save='_time_point.png', **fig_kwargs)

fig_kwargs = dict(basis='umap', color='Pirb', show=False, dpi=300)
scv.pl.velocity_embedding_stream(adata_sub, save='_Pirb.png', **fig_kwargs)

# 计算 velocity pseudotime
scv.tl.velocity_pseudotime(adata_sub)
scv.pl.scatter(adata_sub, color='velocity_pseudotime', save='_pseudotime.png', show=False, dpi=300)

# 识别速度基因
scv.tl.rank_velocity_genes(adata_sub, groupby='time_point', n_genes=10)
velocity_genes = adata_sub.uns['rank_velocity_genes']['names']
velocity_genes_df = pd.DataFrame(velocity_genes)
velocity_genes_df.to_csv(os.path.join(OUT_DIR, 'velocity_genes_by_timepoint.csv'), index=False)

# 计算 Pirb 相关的速度基因
if 'Pirb' in adata_sub.var_names:
    scv.tl.rank_velocity_genes(adata_sub, groupby='Pirb_positive_str', n_genes=20)
    pirb_velocity_genes = pd.DataFrame(adata_sub.uns['rank_velocity_genes']['names'])
    pirb_velocity_genes.to_csv(os.path.join(OUT_DIR, 'velocity_genes_by_Pirb_status.csv'), index=False)

# 保存结果
adata_sub.write(H5AD_OUT)
print(f"[SAVE] {H5AD_OUT}")

# 更新 Pirb 阳性细胞的来源/去向
pirb_pos = adata_sub.obs['Pirb_positive_str'] == 'Pirb+'
if pirb_pos.sum() > 0:
    print(f"[INFO] Pirb+ cells: {pirb_pos.sum()}")
    # 计算每个时间点的 Pirb+ 细胞平均速度
    velocity_df = pd.DataFrame({
        'time_point': adata_sub.obs['time_point'].values,
        'Pirb_positive': adata_sub.obs['Pirb_positive_str'].values,
        'velocity_pseudotime': adata_sub.obs['velocity_pseudotime'].values
    })
    velocity_df.to_csv(os.path.join(OUT_DIR, 'velocity_summary.csv'), index=False)

print("[DONE] scVelo analysis complete.")
