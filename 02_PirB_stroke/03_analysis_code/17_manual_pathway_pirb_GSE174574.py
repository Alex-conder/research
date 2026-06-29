"""
GSE174574 Pirb+ 星形胶质细胞差异基因手动通路富集
- 对预定义的免疫/炎症/髓鞘抑制通路基因集做超几何检验
"""
import os
import numpy as np
import pandas as pd
from scipy.stats import hypergeom

OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE174574/de_pirb"
DEG_FILE = os.path.join(OUT_DIR, "DE_PirbPos_vs_Neg_Astrocyte.csv")

# 读取 DEG
deg = pd.read_csv(DEG_FILE)
all_genes = deg["gene"].tolist()
up_genes = deg[(deg["logfoldchange"] > 0) & (deg["pvals_adj"] < 0.05)]["gene"].tolist()
down_genes = deg[(deg["logfoldchange"] < 0) & (deg["pvals_adj"] < 0.05)]["gene"].tolist()
N = len(all_genes)

# 预定义通路基因集（小鼠）
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
    res_df.to_csv(os.path.join(OUT_DIR, "manual_pathway_enrichment.csv"), index=False)
    print(res_df)
else:
    print("[WARN] No significant overlaps")

print(f"[DONE]")
