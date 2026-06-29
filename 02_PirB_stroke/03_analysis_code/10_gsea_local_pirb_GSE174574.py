"""
GSE174574 Pirb+ vs Pirb- 星形胶质细胞 GSEA（本地基因集）
- 使用所有基因的 ranking（logFC * -log10(padj)）
- 本地 GO_BP 和 KEGG
"""
import os
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

OUT_DIR = "D:/Pirb_stroke_project/04_reports/figures/GSE174574/de_pirb"
DEG_FILE = os.path.join(OUT_DIR, "DE_PirbPos_vs_Neg_Astrocyte.csv")
GS_DIR = "D:/Pirb_stroke_project/04_reports/gene_sets"

deg = pd.read_csv(DEG_FILE)
# 用 logFC 作为 ranking；添加微小随机扰动打破平局
deg["rank_score"] = deg["logfoldchange"] + np.random.normal(0, 1e-6, size=len(deg))
deg = deg[deg["gene"] != "Pirb"]  # 排除 Pirb 自身
rnk = deg.set_index("gene")["rank_score"].sort_values(ascending=False)
print(f"[INFO] Ranked genes: {len(rnk)}")

for db in ["GO_Biological_Process_2021", "KEGG_2019_Mouse"]:
    gmt = os.path.join(GS_DIR, f"{db}.gmt")
    if not os.path.exists(gmt):
        continue
    outdir = os.path.join(OUT_DIR, f"gsea_{db}")
    try:
        res = gp.prerank(rnk=rnk, gene_sets=gmt, outdir=outdir, permutation_num=100, seed=42, format="png")
        df = res.res2d
        df.to_csv(os.path.join(OUT_DIR, f"gsea_{db}.csv"), index=False)
        print(f"\n[GSEA {db}] Top 10 terms:")
        print(df.head(10)[["Term", "NES", "FDR q-val"]])
        
        # 条形图 top 10 正/负 NES
        sig = df[df["FDR q-val"] < 0.25]
        if len(sig) > 0:
            top = sig.head(20).sort_values("NES", ascending=True)
            plt.figure(figsize=(8, 8))
            colors = ["blue" if x < 0 else "red" for x in top["NES"]]
            plt.barh(range(len(top)), top["NES"], color=colors)
            plt.yticks(range(len(top)), top["Term"], fontsize=7)
            plt.xlabel("Normalized Enrichment Score (NES)")
            plt.title(f"GSEA {db} (FDR < 0.25)")
            plt.axvline(0, color="black", linewidth=0.8)
            plt.tight_layout()
            plt.savefig(os.path.join(OUT_DIR, f"gsea_barplot_{db}.png"), dpi=150)
            plt.close()
    except Exception as e:
        print(f"[WARN] {db} GSEA failed: {e}")

print(f"[DONE]")
