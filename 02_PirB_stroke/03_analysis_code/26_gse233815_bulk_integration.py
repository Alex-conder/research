"""
GSE233815 bulk 时间序列与 GSE174574 单细胞结果整合
- 找出 MCAO vs SH 的 bulk 差异基因
- 与 Pirb+ 星形胶质细胞 DEG 取交集
- 热图展示时间动态
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

OUT_DIR = "../04_reports/figures/GSE233815"
os.makedirs(OUT_DIR, exist_ok=True)

# 读取 GSE233815 CPM
cpm = pd.read_csv(os.path.join(OUT_DIR, "cpm_matrix.csv"), index_col=0)
meta = pd.read_csv(os.path.join(OUT_DIR, "sample_metadata.csv"))

# 基因 symbol 映射
feature_file = "../01_raw_data/GSE174574/GSM5319987_sham1_genes.tsv.gz"
gene_map = pd.read_csv(feature_file, sep="\t", header=None, names=["gene_id", "gene_symbol"])
gene_map["gene_id"] = gene_map["gene_id"].str.split(".").str[0]
gene_map = gene_map.drop_duplicates("gene_id").set_index("gene_id")["gene_symbol"]

cpm_symbol = cpm.copy()
cpm_symbol.index = cpm_symbol.index.map(gene_map)
cpm_symbol = cpm_symbol[~cpm_symbol.index.isna()]
cpm_symbol = cpm_symbol.groupby(cpm_symbol.index).sum()

# 读取 GSE174574 Pirb+ astro DEG
pirb_de = pd.read_csv("../04_reports/figures/GSE174574/doublet_qc/DE_PirbPos_vs_Neg_Astrocyte_no_doublet.csv")
pirb_up = pirb_de[(pirb_de["logfoldchange"] > 0) & (pirb_de["pvals_adj"] < 0.05)]["gene"].tolist()
pirb_down = pirb_de[(pirb_de["logfoldchange"] < 0) & (pirb_de["pvals_adj"] < 0.05)]["gene"].tolist()

# 对每个时间点做 t-test，找出 bulk DEG
from scipy import stats
from statsmodels.stats.multitest import multipletests

time_order = ["3h", "12h", "24h", "3D", "7D"]
bulk_de_results = []
for tp in time_order:
    samples_m = meta[(meta["timepoint"] == tp) & (meta["group"] == "MCAO")]["sample"].tolist()
    samples_s = meta[(meta["timepoint"] == tp) & (meta["group"] == "SH")]["sample"].tolist()
    if len(samples_m) < 2 or len(samples_s) < 2:
        continue
    for gene in cpm_symbol.index:
        m_vals = cpm_symbol.loc[gene, samples_m].values
        s_vals = cpm_symbol.loc[gene, samples_s].values
        t, p = stats.ttest_ind(m_vals, s_vals)
        log2fc = np.log2((m_vals.mean() + 1) / (s_vals.mean() + 1))
        bulk_de_results.append({
            "timepoint": tp,
            "gene": gene,
            "log2FC": log2fc,
            "p_value": p,
            "mcao_mean": m_vals.mean(),
            "sh_mean": s_vals.mean()
        })

bulk_de = pd.DataFrame(bulk_de_results)
bulk_de["p_adj"] = multipletests(bulk_de["p_value"], method="fdr_bh")[1]
bulk_de.to_csv(os.path.join(OUT_DIR, "bulk_de_all_genes.csv"), index=False)

# 找出每个时间点显著上调的基因（与 Pirb+ astro 上调基因交集）
overlap_results = []
for tp in time_order:
    tp_de = bulk_de[(bulk_de["timepoint"] == tp) & (bulk_de["p_adj"] < 0.1) & (bulk_de["log2FC"] > 0)]
    overlap = set(tp_de["gene"]) & set(pirb_up)
    overlap_results.append({
        "timepoint": tp,
        "n_bulk_up": len(tp_de),
        "n_pirb_up": len(pirb_up),
        "n_overlap": len(overlap),
        "overlap_genes": ", ".join(sorted(overlap))
    })
    print(f"\n{tp}: {len(overlap)} / {len(pirb_up)} Pirb+ astro genes upregulated in bulk MCAO vs SH")
    print(f"  Genes: {sorted(overlap)[:20]}")

pd.DataFrame(overlap_results).to_csv(os.path.join(OUT_DIR, "pirb_astro_bulk_overlap.csv"), index=False)

# 热图：Pirb+ astro 上调基因在 bulk 时间序列中的表达
overlap_genes_all = set()
for r in overlap_results:
    if r["overlap_genes"]:
        overlap_genes_all.update(r["overlap_genes"].split(", "))
overlap_genes_all = sorted(overlap_genes_all)
print(f"\nTotal overlapping genes: {len(overlap_genes_all)}")

if len(overlap_genes_all) > 0:
    # 计算每个时间点的 MCAO/SH 平均表达
    heatmap_data = []
    for tp in time_order:
        samples_m = meta[(meta["timepoint"] == tp) & (meta["group"] == "MCAO")]["sample"].tolist()
        samples_s = meta[(meta["timepoint"] == tp) & (meta["group"] == "SH")]["sample"].tolist()
        for gene in overlap_genes_all:
            m_mean = cpm_symbol.loc[gene, samples_m].mean() if gene in cpm_symbol.index else 0
            s_mean = cpm_symbol.loc[gene, samples_s].mean() if gene in cpm_symbol.index else 0
            heatmap_data.append({
                "gene": gene,
                "timepoint": tp,
                "MCAO": m_mean,
                "SH": s_mean,
                "MCAO_vs_SH": m_mean - s_mean
            })
    
    heatmap_df = pd.DataFrame(heatmap_data)
    pivot = heatmap_df.pivot(index="gene", columns="timepoint", values="MCAO_vs_SH")
    pivot = pivot[time_order]
    
    plt.figure(figsize=(10, max(6, len(pivot) * 0.3)))
    sns.heatmap(pivot, cmap="RdBu_r", center=0, linewidths=0.5, cbar_kws={"label": "MCAO - SH CPM"})
    plt.title("Pirb+ astrocyte genes in GSE233815 bulk time series")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "pirb_genes_bulk_heatmap.png"), dpi=150)
    plt.close()

#  also plot Pirb itself
if "Pirb" in cpm_symbol.index:
    pirb_expr = cpm_symbol.loc["Pirb"].to_frame(name="Pirb_CPM")
    pirb_expr = pirb_expr.merge(meta.set_index("sample"), left_index=True, right_index=True)
    pirb_summary = pirb_expr.groupby(["timepoint", "group"])["Pirb_CPM"].mean().unstack()
    pirb_summary = pirb_summary.reindex(time_order)
    pirb_summary.plot(kind="bar", figsize=(8, 5), color={"MCAO": "red", "SH": "blue"})
    plt.title("Pirb CPM in GSE233815 bulk (MCAO vs SH)")
    plt.ylabel("Pirb CPM")
    plt.xlabel("Timepoint")
    plt.xticks(rotation=0)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "pirb_barplot_bulk.png"), dpi=150)
    plt.close()

print(f"[DONE] Results saved to {OUT_DIR}")
