"""Day3 Step3: 酶表达代理 → 守恒违例分数（预注册：不归一化）。

输入：
    outputs/pred_expr_seed0.npz  (GEARS 预测表达谱, 217 测试扰动 × 5000 基因符号)
    data/human_gem/parsed_model_S.npz / _meta.npz / _rxn_genes.pkl  (step2 产出)
输出：
    outputs/violation_score.csv  (perturbation, violation, ...)
    outputs/proxy_setup.json     （口径留痕：映射率、子集大小、公式）

口径（预注册锁定）：
- 酶表达代理：AND 反应 v̂_r = 其基因表达（GEARS 预测谱）的均值
- 高变代谢基因子集：映射成功的代谢基因中，预测表达跨 217 扰动的方差 top 25%
- 保留反应：基因全部落在该子集内的 AND 反应
- 违例分数 ‖S·v̂‖₂，不归一化
"""
import json
import pickle

import numpy as np
import pandas as pd
from scipy import sparse

OUT = "outputs"

# ---- 1. 加载预测表达谱 ----
z = np.load(f"{OUT}/pred_expr_seed0.npz")
pred_expr = z["pred_expr"]                    # (217, 5000) float32
perts = z["pred_perturbations"].tolist()      # 217 个扰动基因符号
pred_genes = z["pred_genes"].tolist()         # 5000 个基因符号
gidx = {g: i for i, g in enumerate(pred_genes)}
print(f"[step3] 预测谱: {pred_expr.shape}, 扰动数={len(perts)}")

# ---- 2. 加载模型 ----
S = sparse.load_npz("data/human_gem/parsed_model_S.npz").tocsr()
meta = np.load("data/human_gem/parsed_model_meta.npz", allow_pickle=True)
rxn_ids = meta["rxn_ids"].tolist()
with open("data/human_gem/parsed_model_rxn_genes.pkl", "rb") as f:
    rxn_genes = pickle.load(f)
print(f"[step3] S: {S.shape}, AND 反应={len(rxn_ids)}")

# ---- 3. 基因映射（Human-GEM symbol → GEARS 预测谱 5000 基因）----
model_genes = set(meta["gene_symbols"].tolist())
mapped_genes = sorted(model_genes & set(pred_genes))
map_rate_panel = len(mapped_genes) / len(pred_genes)      # 面板基因是代谢基因的比例
map_rate_model = len(mapped_genes) / len(model_genes)     # 模型基因被面板覆盖的比例
print(f"[step3] 映射: Human-GEM 基因 {len(model_genes)}，命中预测谱 {len(mapped_genes)} "
      f"(模型覆盖 {map_rate_model*100:.1f}%)")

# ---- 4. 高变代谢基因子集（top 25%，预注册锁定）----
var = pred_expr[:, [gidx[g] for g in mapped_genes]].var(axis=0)
n_top = max(1, int(0.25 * len(mapped_genes)))
top_idx = np.argsort(-var)[:n_top]
hv_genes = set(np.array(mapped_genes)[top_idx])
print(f"[step3] 高变子集: {n_top}/{len(mapped_genes)} 基因 (top 25%)")

# ---- 5. 保留反应：基因 ⊆ 高变子集 ----
keep = [j for j, gl in enumerate(rxn_genes)
        if all(g in gidx for g in gl) and all(g in hv_genes for g in gl)]
S_k = S[:, keep]
rxn_kept = [rxn_ids[j] for j in keep]
print(f"[step3] 保留反应: {len(keep)}/{len(rxn_ids)} "
      f"(S_sub = {S_k.shape[0]} mets × {S_k.shape[1]} rxns)")

# ---- 6. 酶表达代理 v̂ 与违例分数 ----
hv_list = sorted(hv_genes)
col = {g: i for i, g in enumerate(hv_list)}           # 高变基因 → R/E 的局部列号
# 反应 × 基因 均值矩阵 R：(n_rxns_kept, n_hv_genes)，每行对其基因列取均值
rows_r, cols_r, vals_r = [], [], []
for i, j in enumerate(keep):
    genes_r = [g for g in rxn_genes[j] if g in col]
    for g in genes_r:
        rows_r.append(i)
        cols_r.append(col[g])
        vals_r.append(1.0 / len(genes_r))
R = sparse.csr_matrix((vals_r, (rows_r, cols_r)),
                      shape=(len(keep), len(hv_list)))
E = pred_expr[:, [gidx[g] for g in hv_list]]         # (217, n_hv)
V = (R @ E.T).T                                       # (217, n_rxns) 代理通量
Sv = S_k @ V.T                                        # (n_mets, 217)
Sv = np.asarray(Sv)
violation = np.sqrt((Sv ** 2).sum(axis=0)).ravel()

df = pd.DataFrame({"perturbation": perts, "violation": violation})
df.to_csv(f"{OUT}/violation_score.csv", index=False)
print("[step3] violation describe():")
print(df.violation.describe())

json.dump({
    "formula": "violation_i = ||S_sub @ vhat_i||_2, unnormalized; "
               "vhat_r = mean expression of AND genes (GEARS pred, seed0)",
    "human_gem": "2.0.1 (2026-03-26)",
    "n_perts": len(perts), "n_panel_genes": len(pred_genes),
    "n_model_genes": len(model_genes), "n_mapped_genes": len(mapped_genes),
    "model_gene_coverage": map_rate_model,
    "hv_gene_frac": 0.25, "n_hv_genes": n_top,
    "n_rxns_total": len(rxn_ids), "n_rxns_kept": len(keep),
    "n_mets": S_k.shape[0],
}, open(f"{OUT}/proxy_setup.json", "w"), indent=2, ensure_ascii=False)
print(f"[step3] -> {OUT}/violation_score.csv, {OUT}/proxy_setup.json")
