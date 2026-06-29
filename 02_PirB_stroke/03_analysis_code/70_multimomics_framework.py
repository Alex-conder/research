"""
70_multimomics_framework.py
多组学整合框架搭建：蛋白质组数据检索 + MOFA 整合 + Pseudo-bulk 矩阵生成
"""
import os, json, requests, gzip, numpy as np, pandas as pd
import scanpy as sc
from pathlib import Path
import warnings; warnings.filterwarnings('ignore')

BASE = Path("..")
FIG = BASE / "04_reports/figures/multiomics"
FIG.mkdir(parents=True, exist_ok=True)

# ============================================================
# 1. 搜索 ProteomeXchange / PRIDE 中脑缺血蛋白质组数据集
# ============================================================
print("=" * 60)
print("1. 搜索 ProteomeXchange / PRIDE 蛋白质组数据集")
print("=" * 60)

# 使用 PRIDE API 搜索 ischemia/stroke + brain proteomics
pride_queries = [
    "cerebral ischemia AND mouse AND brain",
    "stroke AND mouse brain proteomics",
    "MCAO AND proteomics",
    "ischemia reperusion brain proteome",
]

all_pride_results = []
for query in pride_queries:
    url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects?keyword={query.replace(' ', '%20')}&pageSize=20&sortDirection=DESC&sortFields=submissionDate"
    try:
        r = requests.get(url, timeout=30, headers={"Accept": "application/json"})
        if r.status_code == 200:
            data = r.json()
            projects = data.get("_embedded", {}).get("projects", [])
            for p in projects:
                all_pride_results.append({
                    "accession": p.get("accession", ""),
                    "title": p.get("title", ""),
                    "submissionDate": p.get("submissionDate", ""),
                    "publicationDate": p.get("publicationDate", ""),
                    "query": query,
                })
            print(f"  Query '{query}': {len(projects)} results")
        else:
            print(f"  Query '{query}': HTTP {r.status_code}")
    except Exception as e:
        print(f"  Query '{query}': Error {e}")

# 去重
if all_pride_results:
    pride_df = pd.DataFrame(all_pride_results).drop_duplicates(subset="accession")
    pride_df = pride_df.sort_values("publicationDate", ascending=False)
    pride_df.to_csv(FIG / "pride_ischemia_proteomics_search.csv", index=False)
    print(f"\n  总共找到 {len(pride_df)} 个不重复 PRIDE 项目")
    print("  Top 10:")
    for _, row in pride_df.head(10).iterrows():
        print(f"    {row['accession']}: {row['title'][:80]}")
else:
    pride_df = pd.DataFrame()
    print("  未找到 PRIDE 结果")

# ============================================================
# 2. 搜索 GEO 中脑缺血蛋白质组 / DIA / TMT 数据
# ============================================================
print("\n" + "=" * 60)
print("2. 搜索 GEO 中蛋白质组数据")
print("=" * 60)

geo_proteomics = []
geo_search_terms = [
    "ischemia[Title] AND proteomics[Title] AND mouse",
    "stroke[Title] AND mass spectrometry[Title] AND brain",
    "MCAO[Title] AND TMT[Title]",
    "ischemia[Title] AND DIA[Title] AND brain",
]
for term in geo_search_terms:
    url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&term={term.replace(' ', '+')}&retmode=json&retmax=20"
    try:
        r = requests.get(url, timeout=30)
        if r.status_code == 200:
            data = r.json()
            ids = data.get("esearchresult", {}).get("idlist", [])
            print(f"  '{term[:50]}...': {len(ids)} results")
            if ids:
                # 获取详情
                id_str = ",".join(ids[:10])
                detail_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&id={id_str}&retmode=json"
                dr = requests.get(detail_url, timeout=30)
                if dr.status_code == 200:
                    for gid, info in dr.json().get("result", {}).items():
                        if gid == "uids":
                            continue
                        geo_proteomics.append({
                            "uid": gid,
                            "title": info.get("title", ""),
                            "summary": str(info.get("summary", ""))[:200],
                            "gse": info.get("accession", ""),
                        })
    except Exception as e:
        print(f"  Error: {e}")

if geo_proteomics:
    geo_df = pd.DataFrame(geo_proteomics).drop_duplicates(subset="uid")
    geo_df.to_csv(FIG / "geo_proteomics_search.csv", index=False)
    print(f"\n  GEO 蛋白质组: {len(geo_df)} 个")
    for _, row in geo_df.iterrows():
        print(f"    {row['gse']}: {row['title'][:80]}")
else:
    geo_df = pd.DataFrame()
    print("  GEO 无蛋白质组结果")

# ============================================================
# 3. 生成 Pseudo-bulk 矩阵
# ============================================================
print("\n" + "=" * 60)
print("3. 生成 Pseudo-bulk 矩阵 (GSE174574)")
print("=" * 60)

h5ad_path = BASE / "04_reports/figures/GSE174574/GSE174574_annotated.h5ad"
if h5ad_path.exists():
    adata = sc.read_h5ad(str(h5ad_path))
    print(f"  Loaded: {adata.shape}")
    
    # 按 cell_type × condition 聚合
    if 'counts' in adata.layers:
        adata.X = adata.layers['counts']
    
    # 计算 pseudo-bulk: 每个 cell_type × sample 的平均表达
    obs_cols = ['cell_type', 'sample'] if 'sample' in adata.obs.columns else ['cell_type', 'condition']
    
    # 创建伪样本标识
    if 'sample' in adata.obs.columns:
        adata.obs['pseudobulk_key'] = adata.obs['cell_type'].astype(str) + '_' + adata.obs['sample'].astype(str)
    else:
        adata.obs['pseudobulk_key'] = adata.obs['cell_type'].astype(str) + '_' + adata.obs['condition'].astype(str)
    
    pseudobulk = {}
    for key in adata.obs['pseudobulk_key'].unique():
        mask = adata.obs['pseudobulk_key'] == key
        sub = adata[mask]
        pseudobulk[key] = np.asarray(sub.X.sum(axis=0)).flatten()
    
    pb_df = pd.DataFrame(pseudobulk, index=adata.var_names)
    pb_df.index.name = 'gene'
    pb_path = FIG / "pseudobulk_GSE174574_celltype_sample.csv"
    pb_df.to_csv(str(pb_path))
    print(f"  Pseudo-bulk 矩阵: {pb_df.shape} (genes × pseudobulk_keys)")
    print(f"  Saved: {pb_path}")
    
    # 也生成按 cell_type 的聚合
    pb_ct = {}
    for ct in adata.obs['cell_type'].unique():
        mask = adata.obs['cell_type'] == ct
        sub = adata[mask]
        pb_ct[ct] = np.asarray(sub.X.sum(axis=0)).flatten()
    
    pb_ct_df = pd.DataFrame(pb_ct, index=adata.var_names)
    pb_ct_df.index.name = 'gene'
    pb_ct_path = FIG / "pseudobulk_GSE174574_celltype.csv"
    pb_ct_df.to_csv(str(pb_ct_path))
    print(f"  Cell-type 级 pseudo-bulk: {pb_ct_df.shape}")
    
    # 按 cell_type × condition
    pb_ct_cond = {}
    for ct in adata.obs['cell_type'].unique():
        for cond in adata.obs['condition'].unique():
            mask = (adata.obs['cell_type'] == ct) & (adata.obs['condition'] == cond)
            n = mask.sum()
            if n >= 10:
                sub = adata[mask]
                pb_ct_cond[f"{ct}_{cond}"] = np.asarray(sub.X.sum(axis=0)).flatten()
    
    pb_ctc_df = pd.DataFrame(pb_ct_cond, index=adata.var_names)
    pb_ctc_df.index.name = 'gene'
    pb_ctc_path = FIG / "pseudobulk_GSE174574_celltype_condition.csv"
    pb_ctc_df.to_csv(str(pb_ctc_path))
    print(f"  Cell-type × Condition pseudo-bulk: {pb_ctc_df.shape}")
else:
    print("  GSE174574 h5ad not found!")

# 也对 GSE233812 生成 pseudo-bulk
h5ad_812 = BASE / "04_reports/figures/GSE233812_processed.h5ad"
if h5ad_812.exists():
    ad812 = sc.read_h5ad(str(h5ad_812))
    print(f"\n  GSE233812: {ad812.shape}")
    if 'counts' in ad812.layers:
        ad812.X = ad812.layers['counts']
    
    ct_col = 'cell_type' if 'cell_type' in ad812.obs.columns else 'leiden'
    tp_col = 'time_point' if 'time_point' in ad812.obs.columns else 'condition'
    
    pb812 = {}
    for ct in ad812.obs[ct_col].unique():
        for tp in ad812.obs[tp_col].unique():
            mask = (ad812.obs[ct_col] == ct) & (ad812.obs[tp_col] == tp)
            n = mask.sum()
            if n >= 10:
                sub = ad812[mask]
                pb812[f"{ct}_{tp}"] = np.asarray(sub.X.sum(axis=0)).flatten()
    
    pb812_df = pd.DataFrame(pb812, index=ad812.var_names)
    pb812_df.index.name = 'gene'
    pb812_path = FIG / "pseudobulk_GSE233812_celltype_timepoint.csv"
    pb812_df.to_csv(str(pb812_path))
    print(f"  GSE233812 pseudo-bulk: {pb812_df.shape}")

# ============================================================
# 4. MOFA 整合分析框架
# ============================================================
print("\n" + "=" * 60)
print("4. MOFA 整合分析框架")
print("=" * 60)

try:
    import mofapy2
    print(f"  mofapy2 已安装: {mofapy2.__version__}")
    mofa_available = True
except ImportError:
    print("  mofapy2 未安装, 尝试安装...")
    import subprocess, sys
    subprocess.run([sys.executable, "-m", "pip", "install", "mofapy2"], capture_output=True)
    try:
        import mofapy2
        print(f"  mofapy2 安装成功: {mofapy2.__version__}")
        mofa_available = True
    except:
        print("  mofapy2 安装失败, 使用手动 PCA 替代")
        mofa_available = False

# 准备 MOFA 输入: 多视图数据
# View 1: Gene expression (pseudo-bulk)
# View 2: Pathway scores (from GSVA/ssGSEA)
# View 3: TF activities (from SCENIC/DoRothEA)

# 选取 Pirb 相关通路基因集做 pathway score 视图
pirb_pathways = {
    'lysosome': ['Ctsd','Ctsb','Ctsl','Lamp1','Lamp2','Cd68','Atp6v0d1','Atp6v1b2','Gba','Hexa','Hexb','Man2b1','Ctss','Ctsz'],
    'complement': ['C1qa','C1qb','C1qc','C3','C3ar1','C4b','Cfb','Cfd','Cd55','Cd59a','Serping1','Cfh'],
    'nfkb': ['Nfkb1','Nfkb2','Rel','Rela','Relb','Nfkbia','Nfkbiz','Tnfaip3','Tnip1','Ikbkg','Chuk','Ikbkb'],
    'a1_astro': ['C3','Serpina3n','H2-T23','H2-D1','Ggta1','Gbp2','Fkbp5','Ifitm3','Ugt8a','Amigo2','Tgm1','Psmb8'],
    'interferon': ['Irf1','Irf3','Irf5','Irf7','Irf8','Stat1','Stat2','Ifit1','Ifit3','Isg15','Mx1','Oas1a'],
}

if mofa_available and (FIG / "pseudobulk_GSE174574_celltype_condition.csv").exists():
    from mofapy2.run.entry_point import entry_point
    
    pb = pd.read_csv(FIG / "pseudobulk_GSE174574_celltype_condition.csv", index_col=0)
    
    # 计算 pathway score matrix
    genes_in_data = set(pb.index)
    pw_scores = {}
    for pw_name, pw_genes in pirb_pathways.items():
        valid_genes = [g for g in pw_genes if g in genes_in_data]
        if valid_genes:
            pw_scores[pw_name] = pb.loc[valid_genes].mean()
    
    pw_df = pd.DataFrame(pw_scores).T
    print(f"  Pathway score 视图: {pw_df.shape}")
    
    # 标准化
    from sklearn.preprocessing import StandardScaler
    view1 = StandardScaler().fit_transform(pb.T)  # samples x genes (top HVG)
    view2 = StandardScaler().fit_transform(pw_df.T)  # samples x pathways
    
    # 取 top 500 HVG for view1
    var_genes = pb.var(axis=1).nlargest(500).index
    view1_sub = pb.loc[var_genes].T
    view1_norm = StandardScaler().fit_transform(view1_sub)
    
    # MOFA 输入准备
    print("  Running MOFA...")
    try:
        ent = entry_point()
        ent.set_data_options(scale_views=True)
        ent.set_data_matrix(
            [[view1_norm, view2.T.values]],  # list of groups, each group has list of views
            views_names=["GeneExpression", "PathwayScores"],
            groups_names=["GSE174574"],
            samples_names=[list(pb.columns)],
            features_names=[list(var_genes), list(pw_df.index)]
        )
        ent.set_model_options(factors=5)
        ent.set_train_options(iter=500, convergence_mode="medium", seed=42)
        ent.build()
        ent.run()
        
        # 提取因子
        factors = ent.model.getExpectations("Z")["GSE174574"]
        factor_df = pd.DataFrame(factors, index=pb.columns, columns=[f"Factor{i+1}" for i in range(5)])
        factor_df.to_csv(FIG / "mofa_factors_GSE174574.csv")
        print(f"  MOFA 因子: {factor_df.shape}")
        print(factor_df.head())
        
    except Exception as e:
        print(f"  MOFA 运行失败: {e}")
        print("  使用 PCA 替代...")
        # PCA 替代
        from sklearn.decomposition import PCA
        pca = PCA(n_components=5)
        pca_result = pca.fit_transform(view1_norm)
        factor_df = pd.DataFrame(pca_result, index=pb.columns, columns=[f"PC{i+1}" for i in range(5)])
        factor_df.to_csv(FIG / "pca_factors_GSE174574.csv")
        print(f"  PCA 因子: {factor_df.shape}")
        print(f"  方差解释: {pca.explained_variance_ratio_}")
else:
    print("  MOFA 不可用，生成 PCA 替代")
    if (FIG / "pseudobulk_GSE174574_celltype_condition.csv").exists():
        from sklearn.decomposition import PCA
        from sklearn.preprocessing import StandardScaler
        pb = pd.read_csv(FIG / "pseudobulk_GSE174574_celltype_condition.csv", index_col=0)
        var_genes = pb.var(axis=1).nlargest(500).index
        view1_norm = StandardScaler().fit_transform(pb.loc[var_genes].T)
        pca = PCA(n_components=5)
        pca_result = pca.fit_transform(view1_norm)
        factor_df = pd.DataFrame(pca_result, index=pb.columns, columns=[f"PC{i+1}" for i in range(5)])
        factor_df.to_csv(FIG / "pca_factors_GSE174574.csv")
        print(f"  PCA 因子: {factor_df.shape}")
        print(f"  方差解释: {pca.explained_variance_ratio_}")

# ============================================================
# 5. 生成 MOFA 整合分析结果可视化
# ============================================================
print("\n" + "=" * 60)
print("5. 可视化")
print("=" * 60)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import font_manager as fm

# 设置中文字体
plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False

# 5a. MOFA/PCA 因子图
factor_csv = FIG / "mofa_factors_GSE174574.csv"
if not factor_csv.exists():
    factor_csv = FIG / "pca_factors_GSE174574.csv"

if factor_csv.exists():
    fdf = pd.read_csv(str(factor_csv), index_col=0)
    
    # 解析 cell_type 和 condition
    fdf['cell_type'] = [idx.rsplit('_', 1)[0] for idx in fdf.index]
    fdf['condition'] = [idx.rsplit('_', 1)[1] if '_' in idx else 'unknown' for idx in fdf.index]
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    
    # PC1 vs PC2 by cell_type
    for ct in fdf['cell_type'].unique():
        mask = fdf['cell_type'] == ct
        axes[0].scatter(fdf.loc[mask, fdf.columns[0]], fdf.loc[mask, fdf.columns[1]], label=ct, s=80, alpha=0.7)
    axes[0].set_xlabel(fdf.columns[0])
    axes[0].set_ylabel(fdf.columns[1])
    axes[0].set_title('MOFA/PCA: Cell Type')
    axes[0].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=8)
    
    # PC1 vs PC2 by condition
    for cond in fdf['condition'].unique():
        mask = fdf['condition'] == cond
        axes[1].scatter(fdf.loc[mask, fdf.columns[0]], fdf.loc[mask, fdf.columns[1]], label=cond, s=80, alpha=0.7)
    axes[1].set_xlabel(fdf.columns[0])
    axes[1].set_ylabel(fdf.columns[1])
    axes[1].set_title('MOFA/PCA: Condition')
    axes[1].legend()
    
    plt.tight_layout()
    fig.savefig(str(FIG / "mofa_pca_factors.png"), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: mofa_pca_factors.png")

# 5b. Pirb 通路评分热图
if (FIG / "pseudobulk_GSE174574_celltype_condition.csv").exists():
    pb = pd.read_csv(FIG / "pseudobulk_GSE174574_celltype_condition.csv", index_col=0)
    
    # 计算通路评分
    pw_scores_all = {}
    for pw_name, pw_genes in pirb_pathways.items():
        valid_genes = [g for g in pw_genes if g in pb.index]
        if valid_genes:
            pw_scores_all[pw_name] = pb.loc[valid_genes].mean()
    
    pw_heatmap = pd.DataFrame(pw_scores_all)
    
    # 标准化
    from sklearn.preprocessing import StandardScaler
    pw_norm = pd.DataFrame(
        StandardScaler().fit_transform(pw_heatmap),
        index=pw_heatmap.index,
        columns=pw_heatmap.columns
    )
    
    fig, ax = plt.subplots(figsize=(10, 12))
    sns.heatmap(pw_norm, cmap='RdBu_r', center=0, ax=ax, 
                xticklabels=True, yticklabels=True, 
                linewidths=0.5, linecolor='white')
    ax.set_title('Pirb-Related Pathway Scores\n(Pseudo-bulk, GSE174574)')
    ax.set_xlabel('Pathway')
    ax.set_ylabel('Cell Type × Condition')
    plt.tight_layout()
    fig.savefig(str(FIG / "pirb_pathway_heatmap.png"), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: pirb_pathway_heatmap.png")

# 5c. Proteomics 搜索结果汇总图
if not pride_df.empty:
    fig, ax = plt.subplots(figsize=(12, max(4, len(pride_df.head(15)) * 0.4)))
    top_pride = pride_df.head(15).copy()
    top_pride['short_title'] = top_pride['title'].str[:60]
    ax.barh(range(len(top_pride)), [1]*len(top_pride), color='steelblue', alpha=0.7)
    ax.set_yticks(range(len(top_pride)))
    ax.set_yticklabels([f"{r['accession']}: {r['short_title']}" for _, r in top_pride.iterrows()], fontsize=8)
    ax.set_xlabel('PRIDE Projects Found')
    ax.set_title('Brain Ischemia Proteomics Datasets (PRIDE)')
    ax.invert_yaxis()
    plt.tight_layout()
    fig.savefig(str(FIG / "pride_proteomics_search.png"), dpi=150, bbox_inches='tight')
    plt.close()
    print("  Saved: pride_proteomics_search.png")

# ============================================================
# 6. 汇总
# ============================================================
print("\n" + "=" * 60)
print("6. 完成汇总")
print("=" * 60)

outputs = list(FIG.glob("*"))
print(f"  输出目录: {FIG}")
print(f"  文件数: {len(outputs)}")
for f in sorted(outputs):
    size = f.stat().st_size
    unit = "KB" if size < 1024*1024 else "MB"
    val = size/1024 if unit == "KB" else size/(1024*1024)
    print(f"    {f.name}: {val:.1f} {unit}")

print("\n  多组学框架搭建完成!")