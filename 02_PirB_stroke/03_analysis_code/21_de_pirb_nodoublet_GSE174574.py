"""
排除 predicted doublet 后，重新分析 Pirb+ 星形胶质细胞
- 差异表达
- 手动通路富集
- 与原始结果比较
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import hypergeom
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

IN_H5AD = "../04_reports/figures/GSE174574/doublet_qc/GSE174574_with_doublet.h5ad"
OUT_DIR = "../04_reports/figures/GSE174574/doublet_qc"
os.makedirs(OUT_DIR, exist_ok=True)

adata = sc.read_h5ad(IN_H5AD)
adata_nodbl = adata[~adata.obs["predicted_doublet"]].copy()
print(f"[INFO] After removing doublets: {adata_nodbl.n_obs} cells")

# 星形胶质细胞
astro = adata_nodbl[adata_nodbl.obs["cell_type"] == "Astrocyte"].copy()
print(f"[INFO] Astrocytes after doublet removal: {astro.n_obs}")

# Pirb 分组
astro.obs["Pirb_group"] = pd.Categorical((astro[:, "Pirb"].X.toarray().flatten() > 0).astype(int).astype(str))
n_pos = (astro.obs["Pirb_group"] == "1").sum()
n_neg = (astro.obs["Pirb_group"] == "0").sum()
print(f"[DE] Astrocyte Pirb+ = {n_pos}, Pirb- = {n_neg}")

sc.tl.rank_genes_groups(astro, groupby="Pirb_group", method="wilcoxon", pts=True)
result = astro.uns["rank_genes_groups"]
df = pd.DataFrame({
    "gene": result["names"]['1'],
    "logfoldchange": result["logfoldchanges"]['1'],
    "pvals": result["pvals"]['1'],
    "pvals_adj": result["pvals_adj"]['1'],
    "scores": result["scores"]['1'],
})
df.to_csv(os.path.join(OUT_DIR, "DE_PirbPos_vs_Neg_Astrocyte_no_doublet.csv"), index=False)

# 火山图
df["-log10_padj"] = -np.log10(df["pvals_adj"].replace(0, 1e-300))
df["significant"] = (df["pvals_adj"] < 0.05) & (abs(df["logfoldchange"]) > 0.5)
plt.figure(figsize=(8, 6))
plt.scatter(df["logfoldchange"], df["-log10_padj"], c="gray", s=5, alpha=0.5)
sig = df[df["significant"]]
plt.scatter(sig["logfoldchange"], sig["-log10_padj"], c="red", s=10, alpha=0.7)
top = pd.concat([sig.nlargest(8, "logfoldchange"), sig.nsmallest(8, "logfoldchange")])
for _, row in top.iterrows():
    plt.text(row["logfoldchange"], row["-log10_padj"], row["gene"], fontsize=7)
plt.xlabel("Log2 fold change (Pirb+ vs Pirb-)")
plt.ylabel("-log10 adjusted p-value")
plt.title("Astrocyte Pirb+ vs Pirb- (doublets removed)")
plt.axhline(-np.log10(0.05), color="blue", linestyle="--", linewidth=0.8)
plt.axvline(0, color="black", linestyle="-", linewidth=0.8)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "volcano_Astrocyte_no_doublet.png"), dpi=150)
plt.close()

# 手动通路富集
all_genes = df["gene"].tolist()
up_genes = df[(df["logfoldchange"] > 0) & (df["pvals_adj"] < 0.05)]["gene"].tolist()
down_genes = df[(df["logfoldchange"] < 0) & (df["pvals_adj"] < 0.05)]["gene"].tolist()
N = len(all_genes)

pathways = {
    "Complement_activation": ["C1qa", "C1qb", "C1qc", "C3", "C4a", "C4b", "C5", "C5ar1", "C3ar1"],
    "Lysosome": ["Ctsb", "Ctsd", "Ctsz", "Ctss", "Ctsc", "Ctsl", "Lamp1", "Lamp2", "Hexa", "Hexb"],
    "Antigen_presentation_MHC": ["H2-D1", "H2-K1", "H2-M3", "Cd74", "H2-Ab1", "H2-Aa", "B2m"],
    "Inflammasome": ["Nlrp3", "Pycard", "Casp1", "Casp4", "Il1b", "Il18"],
    "Myelin_inhibition_Pirb_ligands": ["Rtn4", "Mag", "Omgp", "Mog"],
    "NFkB_signaling": ["Nfkb1", "Nfkb2", "Rela", "Relb", "Rel", "Nfkbia", "Ikbkb", "Tnfaip3"],
    "JAK_STAT_signaling": ["Stat1", "Stat3", "Stat6", "Jak1", "Jak2", "Socs1", "Socs3"],
    "IRF_signaling": ["Irf1", "Irf3", "Irf5", "Irf7", "Irf8", "Irf9"],
    "RhoA_ROCK": ["Rhoa", "Rock1", "Rock2", "Limk1", "Limk2", "Cfl1", "Arhgef1"],
    "A1_astrocyte_reactive": ["C3", "Gfap", "Vim", "S100a10", "Serping1", "H2-D1", "H2-T23"],
    "Microglia_induced_A1": ["Il1a", "Tnf", "C1qa", "C1qb", "C1qc"],
}

results = []
for name, genes in pathways.items():
    genes_in_background = [g for g in genes if g in all_genes]
    K = len(genes_in_background)
    if K == 0:
        continue
    for direction, deg_set in [("up", up_genes), ("down", down_genes)]:
        k = len(set(deg_set) & set(genes_in_background))
        n = len(deg_set)
        if k > 0:
            pval = hypergeom.sf(k - 1, N, K, n)
            results.append({
                "pathway": name,
                "direction": direction,
                "n_deg": n,
                "n_pathway_genes_in_background": K,
                "overlap": k,
                "overlap_genes": ", ".join(sorted(set(deg_set) & set(genes_in_background))),
                "p_value": pval,
                "log10_p": -np.log10(pval),
            })

res_df = pd.DataFrame(results)
if len(res_df) > 0:
    res_df = res_df.sort_values("p_value")
    res_df.to_csv(os.path.join(OUT_DIR, "manual_pathway_enrichment_no_doublet.csv"), index=False)
    print("\nPathway enrichment (no doublet):")
    print(res_df)

print(f"[DONE]")
