"""
GSE174574 Pirb 特异性配体-受体互作分析
- 手动构建 Pirb 已知配体数据库
- 计算各细胞类型到目标细胞类型的 Pirb 信号强度
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
OUT_DIR = "../04_reports/figures/GSE174574/lr_pirb_specific"
os.makedirs(OUT_DIR, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)
print(f"[INFO] Loaded adata: {adata.n_obs} cells × {adata.n_vars} genes")

# Pirb 已知配体（小鼠 gene symbol）
pirb_lr = pd.DataFrame([
    {"ligand": "Rtn4", "receptor": "Pirb", "pathway": "Nogo-A/RTN4 → PirB"},
    {"ligand": "Mag", "receptor": "Pirb", "pathway": "MAG → PirB"},
    {"ligand": "Omgp", "receptor": "Pirb", "pathway": "OMgp → PirB"},
    {"ligand": "Mog", "receptor": "Pirb", "pathway": "MOG → PirB"},
    {"ligand": "H2-D1", "receptor": "Pirb", "pathway": "MHC I → PirB"},
    {"ligand": "H2-K1", "receptor": "Pirb", "pathway": "MHC I → PirB"},
    {"ligand": "H2-M3", "receptor": "Pirb", "pathway": "MHC I → PirB"},
])

# 检查基因是否存在
available_ligands = [g for g in pirb_lr["ligand"].unique() if g in adata.var_names]
available_receptors = [g for g in pirb_lr["receptor"].unique() if g in adata.var_names]
print(f"[INFO] Available ligands: {available_ligands}")
print(f"[INFO] Available receptors: {available_receptors}")

results = []
for condition in ["MCAO", "Sham"]:
    ad_cond = adata[adata.obs["condition"] == condition]
    cell_types = ad_cond.obs["cell_type"].unique()
    
    ct_mean = {}
    for ct in cell_types:
        mask = ad_cond.obs["cell_type"] == ct
        sub = ad_cond[mask]
        mean_expr = pd.Series(np.array(sub.X.mean(axis=0)).flatten(), index=sub.var_names)
        ct_mean[ct] = mean_expr
    
    for target in cell_types:
        if "Pirb" not in ct_mean[target].index:
            continue
        pirb_expr = ct_mean[target]["Pirb"]
        for source in cell_types:
            for _, row in pirb_lr.iterrows():
                ligand = row["ligand"]
                if ligand not in ct_mean[source].index:
                    continue
                lig_expr = ct_mean[source][ligand]
                score = lig_expr * pirb_expr
                if score > 0:
                    results.append({
                        "condition": condition,
                        "source": source,
                        "target": target,
                        "ligand": ligand,
                        "pathway": row["pathway"],
                        "ligand_expr": lig_expr,
                        "pirb_expr": pirb_expr,
                        "score": score
                    })

res_df = pd.DataFrame(results)
res_df.to_csv(os.path.join(OUT_DIR, "pirb_lr_scores.csv"), index=False)
print(f"[INFO] Total Pirb-related interactions: {len(res_df)}")

# 汇总：每个 source → target 的 Pirb 信号总强度
summary = res_df.groupby(["condition", "source", "target"]).agg(
    total_score=("score", "sum"),
    n_pairs=("score", "size"),
    top_pathway=("pathway", lambda x: x.value_counts().index[0])
).reset_index().sort_values("total_score", ascending=False)
summary.to_csv(os.path.join(OUT_DIR, "pirb_lr_summary.csv"), index=False)
print("\nTop 20 Pirb signaling source → target:")
print(summary.head(20))

# 按 pathway 汇总
pathway_summary = res_df.groupby(["condition", "pathway", "target"]).agg(
    total_score=("score", "sum"),
    max_source=("source", lambda x: x.value_counts().index[0])
).reset_index().sort_values("total_score", ascending=False)
pathway_summary.to_csv(os.path.join(OUT_DIR, "pirb_pathway_summary.csv"), index=False)
print("\nTop 20 Pirb pathway contributions by target:")
print(pathway_summary.head(20))

# 热图：source × target 的 total Pirb 信号（MCAO）
mcao = summary[summary["condition"] == "MCAO"]
if len(mcao) > 0:
    pivot = mcao.pivot(index="source", columns="target", values="total_score").fillna(0)
    plt.figure(figsize=(10, 8))
    sns.heatmap(pivot, cmap="Reds", annot=True, fmt=".2f", cbar_kws={"label": "Pirb signal score"})
    plt.title("Pirb ligand-receptor signal strength (MCAO)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "pirb_signal_heatmap_mcao.png"), dpi=150)
    plt.close()

# 条形图：指向 Astrocyte 的各 pathway 贡献
astro = res_df[(res_df["target"] == "Astrocyte") & (res_df["condition"] == "MCAO")]
if len(astro) > 0:
    pathway_astro = astro.groupby("pathway")["score"].sum().sort_values(ascending=True)
    plt.figure(figsize=(8, 5))
    pathway_astro.plot(kind="barh", color="darkred")
    plt.xlabel("Pirb signal score")
    plt.title("Pirb ligand contributions to Astrocyte (MCAO)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "pirb_pathways_to_astrocyte.png"), dpi=150)
    plt.close()

print(f"[DONE] Results saved to {OUT_DIR}")
