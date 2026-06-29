"""
GSE233815 bulk RNA-seq 时间序列分析
- 读取 57 个 count.txt.gz 文件
- 基因 ID 转换为 gene symbol
- CPM 标准化
- PCA、热图、Pirb 时间趋势
"""
import os, glob, re
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.io import mmread
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

DATA_DIR = "../01_raw_data/GSE233815"
OUT_DIR = "../04_reports/figures/GSE233815"
os.makedirs(OUT_DIR, exist_ok=True)

# 读取所有 count 文件
files = sorted(glob.glob(os.path.join(DATA_DIR, "*.count.txt.gz")))
print(f"[INFO] Found {len(files)} count files")

counts_dict = {}
sample_meta = []
for f in files:
    sname = os.path.basename(f).replace("_union_name.count.txt.gz", "")
    # 解析样本信息
    parts = sname.split("_")
    # 例如 GSM7437165_MCAO_12h_A 或 GSM7437165_MCAO_12h_SH_A
    gsm = parts[0]
    condition = parts[1]  # MCAO or SH
    timepoint = parts[2]  # 3h, 12h, 24h, 3D, 7D
    replicate = parts[-1]  # A, B, C...
    # 判断是否为 SH 样本
    is_sh = "SH" in sname
    group = "SH" if is_sh else "MCAO"
    
    df = pd.read_csv(f, sep="\t", header=None, names=["gene_id", "count"])
    df = df.set_index("gene_id")
    counts_dict[sname] = df["count"]
    sample_meta.append({
        "sample": sname,
        "gsm": gsm,
        "condition": condition,
        "timepoint": timepoint,
        "replicate": replicate,
        "group": group,
        "full_time": f"{timepoint}_{group}"
    })

meta = pd.DataFrame(sample_meta)
counts = pd.DataFrame(counts_dict)
counts = counts.fillna(0)
print(f"[INFO] Count matrix: {counts.shape}")

# 简单基因 ID 转换：去除版本号 ENSMUSG00000000001.4 -> ENSMUSG00000000001
counts.index = counts.index.str.split(".").str[0]
# 去重：取平均
counts = counts.groupby(counts.index).mean()
print(f"[INFO] After dedup: {counts.shape}")

# CPM 标准化
cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
logcpm = np.log2(cpm + 1)

# 保存
meta.to_csv(os.path.join(OUT_DIR, "sample_metadata.csv"), index=False)
counts.to_csv(os.path.join(OUT_DIR, "raw_counts_matrix.csv"))
cpm.to_csv(os.path.join(OUT_DIR, "cpm_matrix.csv"))
logcpm.to_csv(os.path.join(OUT_DIR, "logcpm_matrix.csv"))

# Pirb 表达趋势
# 首先找到 Pirb 对应的 ENSG ID
# 由于我们没有完整的注释，可以尝试在 GSE174574 的 features 中找对应关系
# 先读取 GSE174574 的 features 文件获取 gene_id -> symbol 映射
feature_file = "../01_raw_data/GSE174574/GSM5319987_sham1_genes.tsv.gz"
gene_map = pd.read_csv(feature_file, sep="\t", header=None, names=["gene_id", "gene_symbol"])
gene_map["gene_id"] = gene_map["gene_id"].str.split(".").str[0]
gene_map = gene_map.drop_duplicates("gene_id").set_index("gene_id")["gene_symbol"]

# 将 counts 索引转为 gene symbol
symbol_counts = counts.copy()
symbol_counts.index = symbol_counts.index.map(gene_map)
symbol_counts = symbol_counts[~symbol_counts.index.isna()]
symbol_counts = symbol_counts.groupby(symbol_counts.index).sum()
cpm_symbol = symbol_counts.div(symbol_counts.sum(axis=0), axis=1) * 1e6

if "Pirb" in cpm_symbol.index:
    pirb_expr = cpm_symbol.loc["Pirb"].to_frame(name="Pirb_CPM")
    pirb_expr = pirb_expr.merge(meta.set_index("sample"), left_index=True, right_index=True)
    pirb_expr.to_csv(os.path.join(OUT_DIR, "pirb_expression.csv"))
    
    # 时间顺序
    time_order = ["3h", "12h", "24h", "3D", "7D"]
    pirb_expr["timepoint"] = pd.Categorical(pirb_expr["timepoint"], categories=time_order, ordered=True)
    
    # 画图
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=pirb_expr, x="timepoint", y="Pirb_CPM", hue="group", palette={"MCAO": "red", "SH": "blue"})
    sns.stripplot(data=pirb_expr, x="timepoint", y="Pirb_CPM", hue="group", dodge=True, color="black", size=4, legend=False)
    plt.title("Pirb expression across timepoints (GSE233815 bulk RNA-seq)")
    plt.ylabel("Pirb CPM")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "pirb_time_trend.png"), dpi=150)
    plt.close()
    
    print("\nPirb expression by timepoint and group:")
    print(pirb_expr.groupby(["timepoint", "group"])["Pirb_CPM"].agg(["mean", "std", "count"]))
else:
    print("[WARN] Pirb not found in symbol-mapped matrix")

# PCA
try:
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    
    # 取高变基因
    gene_var = logcpm.var(axis=1)
    top_genes = gene_var.nlargest(2000).index
    pca_data = logcpm.loc[top_genes].T
    
    pca = PCA(n_components=2)
    pca_res = pca.fit_transform(StandardScaler().fit_transform(pca_data))
    pca_df = pd.DataFrame(pca_res, columns=["PC1", "PC2"], index=pca_data.index)
    pca_df = pca_df.merge(meta.set_index("sample"), left_index=True, right_index=True)
    pca_df.to_csv(os.path.join(OUT_DIR, "pca_results.csv"))
    
    plt.figure(figsize=(10, 7))
    sns.scatterplot(data=pca_df, x="PC1", y="PC2", hue="timepoint", style="group", s=100)
    plt.title(f"PCA of GSE233815 bulk samples (top 2000 variable genes)\nPC1={pca.explained_variance_ratio_[0]:.1%}, PC2={pca.explained_variance_ratio_[1]:.1%}")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "pca.png"), dpi=150)
    plt.close()
except Exception as e:
    print(f"[WARN] PCA failed: {e}")

# 差异表达：MCAO vs SH 在每个时间点（简单 t-test）
from scipy import stats
de_results = []
for tp in time_order:
    if tp not in meta["timepoint"].values:
        continue
    samples_m = meta[(meta["timepoint"] == tp) & (meta["group"] == "MCAO")]["sample"].tolist()
    samples_s = meta[(meta["timepoint"] == tp) & (meta["group"] == "SH")]["sample"].tolist()
    if len(samples_m) < 2 or len(samples_s) < 2:
        continue
    for gene in cpm_symbol.index:
        m_vals = cpm_symbol.loc[gene, samples_m].values
        s_vals = cpm_symbol.loc[gene, samples_s].values
        t, p = stats.ttest_ind(m_vals, s_vals)
        log2fc = np.log2((m_vals.mean() + 1) / (s_vals.mean() + 1))
        de_results.append({
            "timepoint": tp,
            "gene": gene,
            "log2FC": log2fc,
            "p_value": p,
            "mcao_mean": m_vals.mean(),
            "sh_mean": s_vals.mean()
        })

de_df = pd.DataFrame(de_results)
if len(de_df) > 0:
    # BH correction per timepoint
    from statsmodels.stats.multitest import multipletests
    de_df["p_adj"] = de_df.groupby("timepoint")["p_value"].transform(lambda x: multipletests(x, method="fdr_bh")[1])
    de_df.to_csv(os.path.join(OUT_DIR, "differential_expression_mcaovssh.csv"), index=False)
    
    # Pirb 的 DE 结果
    pirb_de = de_df[de_df["gene"] == "Pirb"]
    print("\nPirb differential expression (MCAO vs SH):")
    print(pirb_de[["timepoint", "log2FC", "p_value", "p_adj", "mcao_mean", "sh_mean"]])

print(f"[DONE] Results saved to {OUT_DIR}")
