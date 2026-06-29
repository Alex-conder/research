"""
GSE174574 配体-受体互作分析（liana + CellChatDB）
- 预测缺血后脑内各细胞类型间的信号通讯
- 重点关注指向 Astrocyte 的 Pirb 相关配体
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import liana as li
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

IN_H5AD = "D:/Pirb_stroke_project/04_reports/figures/GSE174574/GSE174574_annotated.h5ad"
OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE174574/cellchat"
os.makedirs(OUT_DIR, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)
print(f"[INFO] Loaded: {adata.n_obs} cells × {adata.n_vars} genes")

# 仅保留 condition == MCAO（缺血后的互作更相关）
adata_mcao = adata[adata.obs["condition"] == "MCAO"].copy()
print(f"[INFO] MCAO cells: {adata_mcao.n_obs}")

# 运行 liana rank_aggregate
li.mt.rank_aggregate(adata_mcao, groupby="cell_type", resource_name='mouseconsensus', expr_prop=0.1, min_cells=10, use_raw=False, verbose=True)

# 结果在 adata.uns['liana_res']
lr_res = adata_mcao.uns["liana_res"]
lr_res.to_csv(os.path.join(OUT_DIR, "liana_rank_aggregate_all.csv"), index=False)
print(f"[INFO] Total LR pairs: {len(lr_res)}")

# 筛选指向 Astrocyte 的互作
lr_to_astro = lr_res[lr_res["target"] == "Astrocyte"].sort_values("magnitude_rank")
lr_to_astro.to_csv(os.path.join(OUT_DIR, "liana_to_astrocyte.csv"), index=False)
print("\nTop 20 interactions to Astrocyte:")
print(lr_to_astro.head(20)[["source", "target", "ligand", "receptor", "magnitude_rank", "specificity_rank"]])

# 筛选 Pirb 相关配体-受体对
pirb_ligands = ["Rtn4", "Mag", "Mog", "Omgp", "H2-D1", "H2-K1", "H2-M3"]  # Nogo-A, MAG, MOG, OMgp, MHC I
pirb_related = lr_res[
    (lr_res["receptor"] == "Pirb") |
    (lr_res["ligand"].isin(pirb_ligands)) |
    (lr_res["ligand"].str.contains("H2-", na=False))
]
pirb_related.to_csv(os.path.join(OUT_DIR, "liana_pirb_related.csv"), index=False)
print(f"\nPirb-related LR pairs: {len(pirb_related)}")
print(pirb_related.head(20)[["source", "target", "ligand", "receptor", "magnitude_rank"]])

# 热图：source × target 的互作强度
agg = lr_res.groupby(["source", "target"])["magnitude_rank"].mean().unstack(fill_value=0)
plt.figure(figsize=(10, 8))
import seaborn as sns
sns.heatmap(agg, cmap="viridis_r", annot=False)
plt.title("Mean magnitude rank of LR interactions (MCAO)")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "interaction_heatmap.png"), dpi=150)
plt.close()

# 点图：top 30 指向 Astrocyte 的配体-受体对
top30 = lr_to_astro.head(30)
plt.figure(figsize=(10, 8))
sns.scatterplot(data=top30, x="ligand", y="source", size="magnitude_rank", hue="specificity_rank", palette="coolwarm_r")
plt.xticks(rotation=90, ha="right", fontsize=8)
plt.title("Top 30 ligands to Astrocyte (MCAO)")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "top_ligands_to_astrocyte.png"), dpi=150)
plt.close()

print(f"[DONE] Results saved to {OUT_DIR}")
