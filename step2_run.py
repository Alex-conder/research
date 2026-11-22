"""Day3 Step2: Human-GEM 解析 + GPR 清洗（仅 AND 反应）+ Ensembl→symbol 映射。

预注册锁定：只用 AND 关系反应（含单基因规则，跳过任何含 or 的规则）。
Human-GEM yml 的 genes 段自带 id(Ensembl)→name(symbol) 映射，无需外部 API。

产出 data/human_gem/parsed_model.npz：
    reactions:     AND 反应 id 列表（全部基因成功映射到 symbol）
    S:             (n_mets, n_reactions) stoichiometric 矩阵（float32，scipy sparse 存 csr→npz）
    met_ids:       代谢物 id
    rxn_genes:     object 数组，每个反应的基因 symbol 列表
    gene_symbols:  全部映射成功 symbol（去重）
打印：映射率、AND 反应数等，供 setup_log.md 留痕。
"""
import re
import sys

import numpy as np
import pandas as pd
import yaml
from scipy import sparse

YML = "data/human_gem/Human-GEM_2.0.1.yml"
OUT = "data/human_gem/parsed_model.npz"


class OMapLoader(yaml.SafeLoader):
    pass


def _omap(loader, node):
    # !!omap 在 PyYAML 里是「单键 mapping 的序列」；收成 dict
    result = {}
    for sub in node.value:
        k_node, v_node = sub.value[0]
        result[loader.construct_object(k_node)] = loader.construct_object(v_node)
    return result


OMapLoader.add_constructor("tag:yaml.org,2002:omap", _omap)

print("[step2] 解析 YAML（~8MB，需数十秒）...")
with open(YML, encoding="utf-8") as f:
    model = yaml.load(f, Loader=OMapLoader)
print(f"[step2] Human-GEM version: {model['metaData'].get('version')}, "
      f"date: {model['metaData'].get('date')}")

genes = [dict(g) for g in model["genes"]]
reactions = [dict(r) for r in model["reactions"]]
mets = [dict(m) for m in model["metabolites"]]
print(f"[step2] genes={len(genes)}, reactions={len(reactions)}, mets={len(mets)}")

# Ensembl -> symbol（模型自带）
ens2sym = {g["id"]: g.get("name") for g in genes}
n_noname = sum(1 for v in ens2sym.values() if not v)
print(f"[step2] 模型内无 symbol 的 Ensembl 基因: {n_noname}")

ENSG = re.compile(r"ENSG\d+")

rxn_ids, rxn_gene_syms, rxn_mets, rxn_coefs = [], [], [], []
n_no_gene = n_or = n_unmap = 0
for r in reactions:
    rule = r.get("gene_reaction_rule")
    if not rule or rule.strip() in ("", "nan", "None"):
        n_no_gene += 1
        continue
    if " or " in rule.lower():
        n_or += 1
        continue
    ens = ENSG.findall(rule)
    if not ens:
        n_no_gene += 1
        continue
    syms = []
    ok = True
    for e in ens:
        s = ens2sym.get(e)
        if not s:
            ok = False
            break
        syms.append(s)
    if not ok:
        n_unmap += 1
        continue
    rxn_ids.append(r["id"])
    rxn_gene_syms.append(sorted(set(syms)))
    # 注意：yml 中 metabolites 是 omap → dict {metId: coeff}
    met_d = r["metabolites"]
    rxn_mets.append(list(met_d.keys()))
    rxn_coefs.append([float(c) for c in met_d.values()])

print(f"[step2] 反应清洗: 总={len(reactions)}, 跳过含OR={n_or}, "
      f"无GPR基因={n_no_gene}, 基因未映射={n_unmap}, 保留AND反应={len(rxn_ids)}")

# 构建 S 矩阵
met_ids = sorted({m for ml in rxn_mets for m in ml})
met_index = {m: i for i, m in enumerate(met_ids)}
rows, cols, data = [], [], []
for j, (ml, cl) in enumerate(zip(rxn_mets, rxn_coefs)):
    for m, c in zip(ml, cl):
        rows.append(met_index[m])
        cols.append(j)
        data.append(c)
S = sparse.csr_matrix(
    (np.array(data, dtype=np.float32), (np.array(rows), np.array(cols))),
    shape=(len(met_ids), len(rxn_ids)))
print(f"[step2] S 矩阵: {S.shape[0]} mets × {S.shape[1]} rxns, "
      f"nnz={S.nnz}")

all_syms = sorted({s for gl in rxn_gene_syms for s in gl})
sparse.save_npz(OUT.replace(".npz", "_S.npz"), S)
np.savez(
    OUT.replace(".npz", "_meta.npz"),
    rxn_ids=np.array(rxn_ids),
    met_ids=np.array(met_ids),
    gene_symbols=np.array(all_syms),
)
# 反应→基因列表用 pickle（ragged）
import pickle
with open(OUT.replace(".npz", "_rxn_genes.pkl"), "wb") as f:
    pickle.dump(rxn_gene_syms, f)

# 与 Replogle 基因符号的映射率
pred_genes_available = False
try:
    z = np.load("outputs/pred_expr_seed0.npz")
    expr_syms = set(z["pred_genes"].tolist())
    pred_genes_available = True
except FileNotFoundError:
    expr_syms = set()
if pred_genes_available:
    mapped = expr_syms & set(all_syms)
    print(f"[step2] 映射率: GEARS预测谱基因 {len(expr_syms)} 个，"
          f"命中 Human-GEM AND 反应基因 {len(mapped)} 个 "
          f"({len(mapped)/len(expr_syms)*100:.1f}%)")

print(f"[step2] -> {OUT}")
