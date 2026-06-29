"""
GSE174574 Pirb+ vs Pirb- 星形胶质细胞差异基因本地通路富集
- 使用本地 GMT 文件（GO_BP, KEGG）
- 输出富集结果 CSV 和条形图
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
deg_up = deg[(deg["logfoldchange"] > 0) & (deg["pvals_adj"] < 0.05)]["gene"].tolist()
deg_down = deg[(deg["logfoldchange"] < 0) & (deg["pvals_adj"] < 0.05)]["gene"].tolist()

print(f"[INFO] Up genes: {len(deg_up)}, Down genes: {len(deg_down)}")

# 本地 ORA
for db in ["GO_Biological_Process_2021", "KEGG_2019_Mouse"]:
    gmt = os.path.join(GS_DIR, f"{db}.gmt")
    if not os.path.exists(gmt):
        print(f"[SKIP] {gmt} not found")
        continue
    
    for direction, genes in [("up", deg_up), ("down", deg_down)]:
        if len(genes) < 5:
            continue
        try:
            res = gp.enrichr(gene_list=genes, gene_sets=gmt, outdir=os.path.join(OUT_DIR, f"ora_{direction}_{db}"), no_plot=True)
            df = res.res2d
            df.to_csv(os.path.join(OUT_DIR, f"ora_{direction}_{db}.csv"), index=False)
            print(f"\n[{direction.upper()} {db}] Top 10 terms:")
            print(df.head(10)[["Term", "Overlap", "Adjusted P-value"]])
            
            # 条形图（top 10）
            top = df.head(10).sort_values("Adjusted P-value", ascending=True)
            top["-log10_padj"] = -np.log10(top["Adjusted P-value"].replace(0, 1e-300))
            plt.figure(figsize=(8, 6))
            plt.barh(range(len(top)), top["-log10_padj"])
            plt.yticks(range(len(top)), top["Term"], fontsize=8)
            plt.xlabel("-log10 adjusted p-value")
            plt.title(f"{direction.upper()}-regulated: {db}")
            plt.gca().invert_yaxis()
            plt.tight_layout()
            plt.savefig(os.path.join(OUT_DIR, f"barplot_{direction}_{db}.png"), dpi=150)
            plt.close()
        except Exception as e:
            print(f"[WARN] {db} {direction} failed: {e}")

print(f"[DONE] Results saved to {OUT_DIR}")
