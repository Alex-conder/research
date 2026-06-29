"""
GSE174574 细胞类型注释与 Pirb 跨细胞类型表达分析
- 读取 01_qc_pirb_GSE174574.py 生成的 h5ad
- 基于已知 marker 对 Leiden 簇进行注释
- 计算各细胞类型 Pirb 阳性率 / 平均表达
- 比较 MCAO vs Sham
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

IN_H5AD = "../04_reports/figures/GSE174574/GSE174574_qc.h5ad"
OUT_DIR = "../04_reports/figures/GSE174574"
os.makedirs(OUT_DIR, exist_ok=True)

# 读取
adata = sc.read_h5ad(IN_H5AD)
print(f"[INFO] Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# 已知 marker 字典（小鼠脑）
markers = {
    "Astrocyte": ["Aqp4", "Gja1", "Aldh1l1"],
    "Microglia": ["P2ry12", "Tmem119", "Cx3cr1"],
    "Neuron": ["Rtn1", "Snap25", "Syn1", "Rbfox3"],
    "Oligodendrocyte": ["Mog", "Mbp", "Plp1"],
    "OPC": ["Pdgfra", "Cspg4"],
    "Endothelial": ["Cldn5", "Ly6c1", "Pecam1"],
    "Pericyte": ["Pdgfrb", "Acta2", "Myh11"],
    "Immune": ["Ptprc", "Lyz2", "Cd68"],
    "Ependymal": ["Foxj1", "Cfap126"],
}

# 为每个簇计算 marker 得分（平均表达）
# 取 log1p 标准化后的表达
cluster_means = {}
for ctype, genes in markers.items():
    genes_present = [g for g in genes if g in adata.var_names]
    if not genes_present:
        print(f"[WARN] No marker genes found for {ctype}")
        continue
    # 计算每个细胞该细胞类型的得分
    score = adata[:, genes_present].X.mean(axis=1).A1 if hasattr(adata.X, 'A1') else adata[:, genes_present].X.mean(axis=1)
    adata.obs[f"score_{ctype}"] = score

# 对每个簇，根据得分最高的细胞类型注释
scores = adata.obs[[c for c in adata.obs.columns if c.startswith("score_")]]
cluster_score = scores.groupby(adata.obs["leiden"]).mean()
cluster_ct = cluster_score.idxmax(axis=1).str.replace("score_", "")
cluster_ct = cluster_ct.to_dict()
print("\nCluster annotation:")
for cl, ct in sorted(cluster_ct.items(), key=lambda x: int(x[0])):
    print(f"  Cluster {cl}: {ct}")

adata.obs["cell_type"] = adata.obs["leiden"].map(cluster_ct)

# 保存注释结果
adata.obs.to_csv(os.path.join(OUT_DIR, "cell_metadata_annotated.csv"))

# 各细胞类型数量
counts = adata.obs.groupby(["cell_type", "condition"]).size().unstack(fill_value=0)
counts.to_csv(os.path.join(OUT_DIR, "celltype_counts.csv"))
print("\nCell type counts:")
print(counts)

# Pirb 表达分析
pirb = "Pirb"
adata.obs["Pirb_detected"] = (adata[:, pirb].X.toarray().flatten() > 0).astype(int)
adata.obs["Pirb_counts"] = adata[:, pirb].X.toarray().flatten()

pirb_by_ct = adata.obs.groupby(["cell_type", "condition"]).agg(
    n_cells=("cell_type", "size"),
    pirb_positive=("Pirb_detected", "sum"),
    pirb_positive_pct=("Pirb_detected", "mean"),
    pirb_mean_expr=("Pirb_counts", "mean"),
).reset_index()
pirb_by_ct["pirb_positive_pct"] *= 100
pirb_by_ct.to_csv(os.path.join(OUT_DIR, "pirb_by_celltype_condition.csv"), index=False)
print("\nPirb expression by cell type and condition:")
print(pirb_by_ct)

# 可视化
fig, axes = plt.subplots(2, 2, figsize=(14, 12))
sc.pl.umap(adata, color="cell_type", ax=axes[0, 0], show=False, title="Cell type annotation")
sc.pl.umap(adata, color="leiden", ax=axes[0, 1], show=False, title="Leiden clusters")
sc.pl.umap(adata, color="Pirb_counts", ax=axes[1, 0], show=False, title="Pirb expression")
sc.pl.umap(adata, color="condition", ax=axes[1, 1], show=False, title="Condition")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "umap_celltype.png"), dpi=150)
plt.close()

# 点图：marker 基因
marker_genes = [g for genes in markers.values() for g in genes if g in adata.var_names]
sc.pl.dotplot(adata, marker_genes, groupby="cell_type", save=False, show=False)
plt.savefig(os.path.join(OUT_DIR, "marker_dotplot.png"), dpi=150, bbox_inches='tight')
plt.close()

# 小提琴图：Pirb 按细胞类型
sc.pl.violin(adata, keys=["Pirb_counts"], groupby="cell_type", rotation=45, show=False)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pirb_violin_celltype.png"), dpi=150)
plt.close()

# 柱状图：Pirb 阳性率按细胞类型和条件
pivot_pct = pirb_by_ct.pivot(index="cell_type", columns="condition", values="pirb_positive_pct")
pivot_pct.plot(kind="bar", figsize=(10, 6))
plt.ylabel("Pirb+ cells (%)")
plt.title("Pirb detection rate by cell type and condition")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "pirb_positive_rate_celltype.png"), dpi=150)
plt.close()

# 保存带注释的 h5ad
adata.write_h5ad(os.path.join(OUT_DIR, "GSE174574_annotated.h5ad"))
print(f"\n[DONE] Results saved to {OUT_DIR}")
