# -*- coding: utf-8 -*-
"""本机对齐校验：服务器 scGenePT 产物 vs 本地 step6 口径（GEARS 产物）。

支持 K562 与 RPE1 两个数据集；校验项：
1. 扰动名单与顺序：scgenept npz.pred_perturbations == 本地 pred_expr_seed0.npz
2. 基因空间：scgenept full_gene_space == 本地 npz pred_genes
3. 真值侧对齐：scgenept_labels.csv 的 n_cells == 本地 predictions_seed0.csv 的 n_cells
4. GEARS vs scGenePT pearson_delta 并列（全集 / target_in_input_space 子集）
"""
import json

import numpy as np
import pandas as pd

BASE = "../../"          # virtual_cell_constraint_mvp/
SC = "."                 # outputs_crossmodel/scgenept/

DATASETS = {
    "k562": {"prefix": "scgenept_k562", "local_dir": BASE + "outputs/"},
    "rpe1": {"prefix": "scgenept_rpe1", "local_dir": BASE + "outputs_rpe1/"},
}

report = {}
for name, cfg in DATASETS.items():
    print(f"===== {name} =====")
    z = np.load(SC + cfg["prefix"] + "_predictions.npz")
    lab = pd.read_csv(SC + cfg["prefix"] + "_labels.csv")
    local_npz = np.load(cfg["local_dir"] + "pred_expr_seed0.npz")
    local_csv = pd.read_csv(cfg["local_dir"] + "predictions_seed0.csv")

    # 1) 扰动名单/顺序
    sc_perts = [str(g) for g in z["pred_perturbations"]]
    lc_perts = [str(g) for g in local_npz["pred_perturbations"]]
    print("1) 扰动名单一致:", sc_perts == lc_perts, f"({len(sc_perts)} 个)")

    # 2) 基因空间
    sc_full = [str(g) for g in z["full_gene_space"]]
    lc_genes = [str(g) for g in local_npz["pred_genes"]]
    print("2) 全基因空间一致:", sc_full == lc_genes, f"({len(sc_full)} 个)")

    # 3) n_cells 对齐
    m = lab.set_index("perturbation").loc[lc_perts]
    lc = local_csv.set_index("perturbation").loc[lc_perts]
    nc_match = (m["n_cells"].values == lc["n_cells"].values).all()
    print("3) n_cells 全部一致:", bool(nc_match),
          f"(scgenept 合计 {int(m['n_cells'].sum())} vs 本地 {int(lc['n_cells'].sum())})")

    # 4) GEARS vs scGenePT
    lab2 = lab.dropna(subset=["pearson_delta"]).set_index("perturbation")
    cmp = pd.DataFrame({
        "gears": local_csv.set_index("perturbation")["pearson_delta"],
        "scgenept": lab2["pearson_delta"],
        "scgenept_all": lab2["pearson_delta_all"],
        "target_in_input_space": lab2["target_in_input_space"],
    }).dropna()
    print(f"4) 可比扰动 n={len(cmp)}")
    print(f"   GEARS      pearson_delta: mean={cmp.gears.mean():.4f} "
          f"median={cmp.gears.median():.4f}")
    print(f"   scGenePT   pearson_delta: mean={cmp.scgenept.mean():.4f} "
          f"median={cmp.scgenept.median():.4f}")
    sub = cmp[cmp.target_in_input_space]
    print(f"   目标在输入空间子集 n={len(sub)}:")
    print(f"   GEARS      pearson_delta: mean={sub.gears.mean():.4f} "
          f"median={sub.gears.median():.4f}")
    print(f"   scGenePT   pearson_delta: mean={sub.scgenept.mean():.4f} "
          f"median={sub.scgenept.median():.4f}")

    # 5) 预测谱健全性
    pe = z["pred_expr"]
    pf = z["pred_expr_full"]
    print("5) pred_expr(matched) 形状:", pe.shape, "NaN 数:", int(np.isnan(pe).sum()))
    print("   pred_expr_full 形状:", pf.shape, "NaN 比例:", f"{np.isnan(pf).mean():.4f}")

    report[name] = {
        "perts_match": sc_perts == lc_perts,
        "n_perts": len(sc_perts),
        "gene_space_match": sc_full == lc_genes,
        "n_cells_match": bool(nc_match),
        "n_compared": int(len(cmp)),
        "pearson_delta_mean": {"gears": float(cmp.gears.mean()),
                               "scgenept": float(cmp.scgenept.mean())},
        "subset_flagged_n": int(len(sub)),
        "pearson_delta_mean_flagged": {
            "gears": float(sub.gears.mean()),
            "scgenept": float(sub.scgenept.mean())},
        "pred_expr_nan": int(np.isnan(pe).sum()),
        "pred_expr_full_nan_ratio": float(np.isnan(pf).mean()),
    }
    print()

with open(SC + "alignment_check.json", "w", encoding="utf-8") as f:
    json.dump(report, f, indent=2, ensure_ascii=False)
print("-> alignment_check.json")
