"""
尝试使用 scVelo 推断 GSE233812 小胶质细胞 RNA 速率。
注意：GSE233812 处理后 h5ad 仅含 counts 层，缺乏 spliced/unspliced counts，
      因此 scVelo 无法直接运行。本脚本记录该限制并输出说明文件。
"""
import os
import sys
import scanpy as sc
import numpy as np

H5AD_PATH = "../04_reports/figures/GSE233812_processed.h5ad"
OUT_DIR = "../04_reports/figures/rna_velocity"
os.makedirs(OUT_DIR, exist_ok=True)

print("[INFO] RNA velocity analysis for GSE233812 microglia")

# 检查 scvelo 是否可用
try:
    import scvelo as scv
    has_scvelo = True
    print("[INFO] scvelo version:", scv.__version__)
except ImportError as e:
    has_scvelo = False
    print("[WARN] scvelo import failed:", e)

# 读取数据
if not os.path.exists(H5AD_PATH):
    print(f"[ERROR] H5AD not found: {H5AD_PATH}")
    sys.exit(1)

adata = sc.read_h5ad(H5AD_PATH)
print(f"[INFO] Loaded {adata.n_obs} cells × {adata.n_vars} genes")
print("[INFO] Layers:", list(adata.layers.keys()))
print("[INFO] Pirb positive cells:", int((adata[:, 'Pirb'].X.toarray() > 0).sum()))

# 判断是否有 spliced/unspliced
has_spliced = 'spliced' in adata.layers and 'unspliced' in adata.layers
print("[INFO] Has spliced/unspliced layers:", has_spliced)

if not has_scvelo or not has_spliced:
    msg = (
        "RNA velocity analysis could not be completed.\n"
        "Reasons:\n"
        "1. scvelo is not installed in the current environment, OR\n"
        "2. The processed h5ad lacks 'spliced' and 'unspliced' count layers.\n\n"
        "Current layers: " + str(list(adata.layers.keys())) + "\n\n"
        "To run RNA velocity, please:\n"
        "- Install scvelo: pip install scvelo\n"
        "- Generate spliced/unspliced counts from raw fastq using velocyto or STARsolo,\n"
        "  then merge into the h5ad object.\n"
    )
    out_txt = os.path.join(OUT_DIR, "rna_velocity_status.txt")
    with open(out_txt, 'w', encoding='utf-8') as f:
        f.write(msg)
    print(f"[SAVE] {out_txt}")
    print("[DONE] RNA velocity analysis skipped due to missing data/dependency.")
    sys.exit(0)

# 如果未来数据可用，则继续标准 scVelo 流程
print("[INFO] Running scvelo standard pipeline...")
scv.pp.filter_and_normalize(adata)
scv.pp.moments(adata, n_pcs=30, n_neighbors=30)
scv.tl.velocity(adata)
scv.tl.velocity_graph(adata)
scv.tl.umap(adata)

# 按时间着色并绘制速率
adata.obs['time_point'] = adata.obs.get('time_point', 'unknown')
fig, ax = plt.subplots(figsize=(8, 6))
scv.pl.velocity_embedding_stream(adata, basis='umap', color='time_point', ax=ax, show=False)
fig.savefig(os.path.join(OUT_DIR, "velocity_stream_time_point.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# 按 Pirb 表达着色
fig, ax = plt.subplots(figsize=(8, 6))
scv.pl.velocity_embedding_stream(adata, basis='umap', color='Pirb', ax=ax, show=False)
fig.savefig(os.path.join(OUT_DIR, "velocity_stream_Pirb.png"), dpi=300, bbox_inches='tight')
plt.close(fig)

# 保存 velocity 结果
adata.write(os.path.join(OUT_DIR, "microglia_velocity.h5ad"))
print("[DONE] RNA velocity analysis completed.")
