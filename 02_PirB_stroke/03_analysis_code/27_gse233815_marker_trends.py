"""
GSE233815 bulk 中关键 marker 基因的时间趋势
- Pirb, C3, Gfap, C1qa, Spp1, Tnf, Il1a 等
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

OUT_DIR = "../04_reports/figures/GSE233815"
cpm_symbol = pd.read_csv(os.path.join(OUT_DIR, "cpm_matrix_symbol.csv"), index_col=0) if os.path.exists(os.path.join(OUT_DIR, "cpm_matrix_symbol.csv")) else None

# 如果没有预保存，重新生成
if cpm_symbol is None:
    cpm = pd.read_csv(os.path.join(OUT_DIR, "cpm_matrix.csv"), index_col=0)
    feature_file = "../01_raw_data/GSE174574/GSM5319987_sham1_genes.tsv.gz"
    gene_map = pd.read_csv(feature_file, sep="\t", header=None, names=["gene_id", "gene_symbol"])
    gene_map["gene_id"] = gene_map["gene_id"].str.split(".").str[0]
    gene_map = gene_map.drop_duplicates("gene_id").set_index("gene_id")["gene_symbol"]
    cpm_symbol = cpm.copy()
    cpm_symbol.index = cpm_symbol.index.map(gene_map)
    cpm_symbol = cpm_symbol[~cpm_symbol.index.isna()]
    cpm_symbol = cpm_symbol.groupby(cpm_symbol.index).sum()
    cpm_symbol.to_csv(os.path.join(OUT_DIR, "cpm_matrix_symbol.csv"))

meta = pd.read_csv(os.path.join(OUT_DIR, "sample_metadata.csv"))

markers = {
    "Pirb": "Pirb",
    "Astrocyte": ["Aqp4", "Gja1", "Aldh1l1", "Gfap"],
    "A1_reactive": ["C3", "Gfap", "Vim", "Serping1"],
    "Microglia_A1_inducer": ["C1qa", "C1qb", "Tnf", "Il1a"],
    "Inflammation": ["Spp1", "Ccl2", "Ccl12", "Lyz2"],
    "Endothelial": ["Cldn5", "Abcb1a", "Cxcl12"],
}

# 计算每个 marker 在每个时间点的平均表达
time_order = ["3h", "12h", "24h", "3D", "7D"]
records = []
for category, genes in markers.items():
    if isinstance(genes, str):
        genes = [genes]
    for gene in genes:
        if gene not in cpm_symbol.index:
            continue
        expr = cpm_symbol.loc[gene].to_frame(name="CPM").merge(meta.set_index("sample"), left_index=True, right_index=True)
        summary = expr.groupby(["timepoint", "group"])["CPM"].mean().reset_index()
        for _, row in summary.iterrows():
            records.append({
                "category": category,
                "gene": gene,
                "timepoint": row["timepoint"],
                "group": row["group"],
                "CPM": row["CPM"]
            })

df = pd.DataFrame(records)
df.to_csv(os.path.join(OUT_DIR, "marker_time_trends.csv"), index=False)

# 画图
n_categories = len(markers)
fig, axes = plt.subplots(n_categories, 1, figsize=(10, 3 * n_categories), sharex=True)
if n_categories == 1:
    axes = [axes]

for ax, (category, genes) in zip(axes, markers.items()):
    if isinstance(genes, str):
        genes = [genes]
    sub = df[(df["category"] == category) & (df["gene"].isin(genes))]
    if len(sub) == 0:
        continue
    for gene in genes:
        gsub = sub[sub["gene"] == gene]
        if len(gsub) == 0:
            continue
        gsub_m = gsub[gsub["group"] == "MCAO"].set_index("timepoint").reindex(time_order)
        gsub_s = gsub[gsub["group"] == "SH"].set_index("timepoint").reindex(time_order)
        ax.plot(time_order, gsub_m["CPM"], marker="o", label=f"{gene} (MCAO)", linewidth=2)
        ax.plot(time_order, gsub_s["CPM"], marker="s", linestyle="--", label=f"{gene} (SH)", linewidth=1.5)
    ax.set_title(category)
    ax.set_ylabel("CPM")
    ax.legend(fontsize=7, ncol=2)
    ax.grid(True, alpha=0.3)

axes[-1].set_xlabel("Timepoint")
plt.suptitle("Marker gene expression dynamics in GSE233815 bulk", y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(OUT_DIR, "marker_time_trends.png"), dpi=150)
plt.close()

print(f"[DONE] Results saved to {OUT_DIR}")
