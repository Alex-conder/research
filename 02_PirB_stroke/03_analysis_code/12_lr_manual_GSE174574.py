"""
GSE174574 简化版配体-受体互作分析（基于 mouseconsensus）
- 使用稀疏矩阵计算细胞类型平均表达，避免内存爆炸
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

IN_H5AD = "../04_reports/figures/GSE174574/GSE174574_annotated.h5ad"
LR_FILE = "../04_reports/gene_sets/mouseconsensus_lr.csv"
OUT_DIR = "../04_reports/figures/GSE174574/lr_manual"
os.makedirs(OUT_DIR, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)
lr = pd.read_csv(LR_FILE)
print(f"[INFO] Loaded adata: {adata.n_obs} cells, LR pairs: {len(lr)}")

# 只保留 LR 数据库中出现的基因
lr_genes = set(lr["ligand"]).union(set(lr["receptor"]))
lr_genes = [g for g in lr_genes if g in adata.var_names]
adata_lr = adata[:, lr_genes].copy()
print(f"[INFO] Subset to {len(lr_genes)} LR-related genes")

results = []
for condition in ["MCAO", "Sham"]:
    ad_cond = adata_lr[adata_lr.obs["condition"] == condition]
    cell_types = ad_cond.obs["cell_type"].unique()
    
    # 计算每种细胞类型的平均表达（稀疏）
    ct_mean = {}
    for ct in cell_types:
        mask = ad_cond.obs["cell_type"] == ct
        sub = ad_cond[mask]
        mean_expr = np.array(sub.X.mean(axis=0)).flatten()
        ct_mean[ct] = pd.Series(mean_expr, index=sub.var_names)
    
    for source in cell_types:
        for target in cell_types:
            if source == target:
                continue
            for _, row in lr.iterrows():
                ligand = row["ligand"]
                receptor = row["receptor"]
                if ligand not in ct_mean[source].index or receptor not in ct_mean[target].index:
                    continue
                score = ct_mean[source][ligand] * ct_mean[target][receptor]
                if score > 0:
                    results.append({
                        "condition": condition,
                        "source": source,
                        "target": target,
                        "ligand": ligand,
                        "receptor": receptor,
                        "score": score
                    })

res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(OUT_DIR, "lr_scores_all.csv"), index=False)
print(f"[INFO] Total interactions: {len(res_df)}")

# Top 互作对
top_lr = res_df.sort_values("score", ascending=False).head(100)
top_lr.to_csv(os.path.join(OUT_DIR, "top100_lr.csv"), index=False)
print("\nTop 20 LR interactions:")
print(top_lr.head(20)[["condition", "source", "target", "ligand", "receptor", "score"]])

# 指向 Astrocyte 的 top 互作
to_astro = res_df[res_df["target"] == "Astrocyte"].sort_values("score", ascending=False)
to_astro.to_csv(os.path.join(OUT_DIR, "lr_to_astrocyte.csv"), index=False)
print("\nTop 20 interactions to Astrocyte:")
print(to_astro.head(20)[["condition", "source", "ligand", "receptor", "score"]])

# Pirb 相关配体-受体对
pirb_ligands = ["Rtn4", "Mag", "Mog", "Omgp"]
pirb_related = res_df[
    (res_df["receptor"] == "Pirb") |
    (res_df["ligand"].isin(pirb_ligands)) |
    (res_df["ligand"].str.startswith("H2-", na=False))
].sort_values("score", ascending=False)
pirb_related.to_csv(os.path.join(OUT_DIR, "lr_pirb_related.csv"), index=False)
print(f"\nPirb-related LR pairs: {len(pirb_related)}")
print(pirb_related.head(20)[["condition", "source", "target", "ligand", "receptor", "score"]])

# 可视化：MCAO 条件下 source × target 互作强度热图
mcao = res_df[res_df["condition"] == "MCAO"]
agg = mcao.groupby(["source", "target"])["score"].sum().unstack(fill_value=0)
plt.figure(figsize=(10, 8))
sns.heatmap(np.log1p(agg), cmap="viridis", annot=False, cbar_kws={"label": "log1p(sum score)"})
plt.title("LR interaction strength (MCAO)")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "interaction_heatmap_mcao.png"), dpi=150)
plt.close()

# 可视化：top 30 指向 Astrocyte 的互作（点图）
top30 = to_astro.head(30)
plt.figure(figsize=(10, 8))
sns.scatterplot(data=top30, x="ligand", y="source", size="score", hue="score", palette="coolwarm", sizes=(20, 200))
plt.xticks(rotation=90, ha="right", fontsize=8)
plt.title("Top 30 ligands to Astrocyte (MCAO + Sham)")
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "top_ligands_to_astrocyte.png"), dpi=150)
plt.close()

print(f"[DONE] Results saved to {OUT_DIR}")
