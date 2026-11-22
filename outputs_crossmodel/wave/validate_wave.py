# -*- coding: utf-8 -*-
"""Wave 产物校验：扰动数、嵌入维度、labels 探针列、pearson 全 NaN、值有限。

只读校验，不改数据。用法: python validate_wave.py <wave_dir>
实际 npz 结构（以 scgpt_k562 为参照）:
  pred_perturbations (n,) / pred_emb_cells (n,32,dim) / pred_emb_mean (n,dim) /
  pred_delta_emb (n,dim) / true_emb_cells (n,<=32,dim) / true_emb_mean (n,dim) /
  true_n_cells (n,) / ctrl_emb (dim,) / ctrl_emb_cells (32,dim) /
  target_in_input_space (n,) / target_masked|target_deleted (n,) / output_type
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

EXPECTED_N = {"k562": 217, "rpe1": 306}
EMB_DIM = {"scgpt": 512, "geneformer": 768, "uce": 1280}
BASE_COLS = ["perturbation", "pearson_delta", "pearson_delta_all", "n_cells",
             "target_in_input_space"]
ENGINES = ["scgpt", "geneformer", "uce"]
ENGINE_FLAG = {"scgpt": "target_masked", "geneformer": "target_deleted", "uce": "target_deleted"}
PROBE_COLS = ["pred_emb_norm", "true_emb_norm", "delta_emb_norm", "err_emb_dist"]
REQUIRED_NPZ = ["pred_perturbations", "pred_emb_cells", "pred_emb_mean", "pred_delta_emb",
                "true_emb_cells", "true_emb_mean", "true_n_cells", "ctrl_emb",
                "ctrl_emb_cells", "target_in_input_space", "output_type"]


def is_nan(x):
    s = str(x).strip()
    if not s:
        return True  # 空串 == NaN 口径（表达谱列对 embedding 产物设计为空）
    try:
        return math.isnan(float(s))
    except (TypeError, ValueError):
        return False


def main(wave_dir):
    wave = Path(wave_dir)
    ok_all = True
    for eng in ["scgpt", "geneformer", "uce"]:
        for ds, n_exp in EXPECTED_N.items():
            tag = f"{eng}/{ds}"
            npz_p = wave / eng / f"{eng}_wave_{ds}_predictions.npz"
            lab_p = wave / eng / f"{eng}_wave_{ds}_labels.csv"
            meta_p = wave / eng / f"{eng}_wave_{ds}_meta.json"
            checks = []

            def chk(name, cond, detail=""):
                checks.append((name, bool(cond), str(detail)))

            with open(lab_p, newline="", encoding="utf-8") as f:
                rows = list(csv.DictReader(f))
            cols = list(rows[0].keys()) if rows else []
            chk("labels_rows", len(rows) == n_exp, f"{len(rows)} vs {n_exp}")
            need = BASE_COLS + [ENGINE_FLAG[eng]] + PROBE_COLS
            chk("labels_cols", all(c in cols for c in need),
                "missing=" + str([c for c in need if c not in cols]))
            chk("pearson_all_nan",
                all(is_nan(r["pearson_delta"]) and is_nan(r["pearson_delta_all"]) for r in rows))
            chk("probe_cols_finite",
                all(all(not is_nan(r[c]) for c in PROBE_COLS) for r in rows))

            z = np.load(npz_p, allow_pickle=False)
            keys = set(z.files)
            dim = EMB_DIM[eng]
            miss = [k for k in REQUIRED_NPZ + [ENGINE_FLAG[eng]] if k not in keys]
            chk("npz_keys", not miss, "missing=" + str(miss))
            if "pred_perturbations" in keys:
                chk("npz_n", len(z["pred_perturbations"]) == n_exp,
                    len(z["pred_perturbations"]))
            for k in ["pred_emb_cells", "pred_emb_mean", "pred_delta_emb",
                      "true_emb_mean", "ctrl_emb", "ctrl_emb_cells"]:
                if k in keys:
                    a = z[k]
                    got_dim = a.shape[-1]
                    chk(f"dim[{k}]", got_dim == dim, f"{got_dim} vs {dim}")
                    chk(f"finite[{k}]", bool(np.isfinite(a).all()), a.shape)
            if "true_emb_cells" in keys:
                a = z["true_emb_cells"]
                chk("dim[true_emb_cells]", a.shape[-1] == dim, a.shape[-1])
                nan_cell = ~np.isfinite(a).any(axis=2)      # 整细胞 NaN（设计：条件细胞数<32 的填充）
                nan_part = np.isfinite(a).any(axis=2) & ~np.isfinite(a).all(axis=2)
                chk("true_cells_no_partial_nan", not bool(nan_part.any()),
                    f"partial-nan cells={int(nan_part.sum())}")
                if "true_n_cells" in keys:
                    tn = z["true_n_cells"]
                    chk("true_padding_rule", bool((nan_cell.sum(axis=1) == 32 - tn).all()),
                        f"rows_with_padding={int((nan_cell.sum(axis=1) > 0).sum())}/{len(tn)}")
                    valid = np.isfinite(a).all(axis=2)
                    chk("finite[true_emb_cells_valid]",
                        bool(np.isfinite(a[valid]).all()),
                        f"valid_cells={int(valid.sum())}")
            if "pred_emb_cells" in keys and "true_n_cells" in keys:
                chk("ctrl_pool_32", z["pred_emb_cells"].shape[1] == 32,
                    z["pred_emb_cells"].shape)
                chk("true_cells_pos", bool((z["true_n_cells"] > 0).all()),
                    f"min={z['true_n_cells'].min()}")
            if "output_type" in keys:
                ot = str(z["output_type"][0]).lower()
                chk("output_type_emb", "embedding" in ot, ot)
            meta = json.loads(meta_p.read_text(encoding="utf-8"))
            chk("meta_emb", "embedding" in meta["model"]["output_type"].lower(),
                meta["model"]["output_type"])
            chk("meta_ds", ds in meta["data"]["h5ad"].lower(), meta["data"]["h5ad"])

            bad = [(n, d) for n, ok, d in checks if not ok]
            ok_all = ok_all and not bad
            short = "; ".join(n for n, ok, _ in checks if ok)
            line = f"PASS {tag} [{len(checks)-len(bad)}/{len(checks)}] {short}"
            if bad:
                line += " || FAIL: " + "; ".join(f"{n}({d})" for n, d in bad)
            print(line)
    print("OVERALL:", "PASS" if ok_all else "FAIL")
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "."))
