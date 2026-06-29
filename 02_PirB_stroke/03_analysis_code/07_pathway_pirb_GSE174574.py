"""
GSE174574 Pirb+ vs Pirb- 星形胶质细胞差异基因通路富集
- 使用 gseapy 本地基因集进行 GO / KEGG 富集
"""
import os
import numpy as np
import pandas as pd
import gseapy as gp

OUT_DIR = "../04_reports/figures/GSE174574/de_pirb"
DEG_FILE = os.path.join(OUT_DIR, "DE_PirbPos_vs_Neg_Astrocyte.csv")

deg = pd.read_csv(DEG_FILE)
deg_up = deg[(deg["logfoldchange"] > 0) & (deg["pvals_adj"] < 0.05)].sort_values("logfoldchange", ascending=False)
deg_down = deg[(deg["logfoldchange"] < 0) & (deg["pvals_adj"] < 0.05)].sort_values("logfoldchange")

print(f"[INFO] Up-regulated genes: {len(deg_up)}")
print(f"[INFO] Down-regulated genes: {len(deg_down)}")

# 准备 ranked list（排除 Pirb 自身）
deg["rank_score"] = deg["logfoldchange"] * (-np.log10(deg["pvals_adj"].replace(0, 1e-300)))
deg_nopirb = deg[deg["gene"] != "Pirb"]
rnk = deg_nopirb.set_index("gene")["rank_score"].sort_values(ascending=False)

# GSEA 富集（本地基因集）
for db in ["GO_Biological_Process_2021", "KEGG_2019_Mouse"]:
    outdir = os.path.join(OUT_DIR, f"gsea_{db}")
    try:
        res = gp.prerank(rnk=rnk, gene_sets=db, outdir=outdir, permutation_num=100, seed=42)
        print(f"\n[ENRICHMENT] {db} top terms:")
        print(res.res2d.head(10)[["Term", "NES", "FDR q-val"]])
    except Exception as e:
        print(f"[WARN] {db} enrichment failed: {e}")

# ORA 富集：上调基因
top_up = deg_up["gene"].tolist()
try:
    ora = gp.enrichr(gene_list=top_up, gene_sets="GO_Biological_Process_2021", organism="mouse", outdir=os.path.join(OUT_DIR, "ora_up_GO_BP"))
    print("\n[ORA] GO BP top terms for up-regulated genes:")
    print(ora.res2d.head(10)[["Term", "Adjusted P-value"]])
except Exception as e:
    print(f"[WARN] ORA failed: {e}")

print(f"[DONE] Results saved to {OUT_DIR}")
