"""Step 2: 表达谱 -> 酶表达代理（Human-GEM GPR 粗映射）。

产出 outputs/enzyme_proxy.csv：每行一个测试扰动，列为反应 ID 的酶代理值。

映射规则（锁定，见 preregistration.md）：
  每条代谢反应 r 的 GPR 取其所涉基因的表达均值。
  一阶段【只用 AND 关系反应】（多基因同时表达才激活），跳过 OR 关系——
  降低映射噪声（2026-10-12 评审意见，已留痕预注册）。
  只用高变代谢基因子集（表达变异 top 25%）。
  技术隐患备忘：Human-GEM 的 GPR 非标准格式，需手动清洗；若清洗卡壳超过半天，
  先取糖酵解等单条通路的子网络跑通链路，再扩展（README Day3 计划）。
"""
import numpy as np
import pandas as pd

# ---- 1. 加载 Human-GEM 与 GPR ----
import cobra
model = cobra.io.read_sbml_model("data/Human-GEM.xml")   # 或 .mat via scipy.io
gpr = {r.id: [g.id for g in r.genes] for r in model.reactions if r.gene_reaction_rule}

# ---- 2. 高变代谢基因子集（锁定：表达变异 top 25%）----
# all_expr = 加载全部细胞表达（train+test）
# var = all_expr.var(axis=0); high_var_genes = var.nlargest(int(0.25 * len(var))).index

# ---- 3. 对每个测试扰动的预测表达计算酶代理 ----
pred = pd.read_csv("outputs/predictions.csv")
rows = []
for _, row in pred.iterrows():
    expr = np.load(f"outputs/pred_{row.perturbation}_seed0.npy")  # 用 ensemble 均值更稳，可改
    gene2expr = dict(zip(gene_names, expr))                        # gene_names 需从 adata 导出
    proxy = {rid: np.mean([gene2expr[g] for g in genes if g in gene2expr])
             for rid, genes in gpr.items()
             if any(g in high_var_genes for g in genes)}
    rows.append({"perturbation": row.perturbation, **proxy})

pd.DataFrame(rows).to_csv("outputs/enzyme_proxy.csv", index=False)
print("Step 2 done -> outputs/enzyme_proxy.csv")
