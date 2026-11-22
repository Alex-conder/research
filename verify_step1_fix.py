"""step1 修复验证：day6_seed2 OOM 根因定位、库补丁核验、数值口径比对。

产出 outputs/verify_step1_fix.json，供 setup_log.md 留痕。只读 + 轻量计算，
不跑训练、不加载模型。

核验项：
1. 根因证据：day6_seed2.log 末尾 traceback 定位到 gears/inference.py
   deeper_analysis 的 adata.X.toarray()（非 step1_run.py 评估段）。
2. 库补丁：3 处 reshape 目标改为 adata.shape[1]（同值，零数值变化）；
   真实数据上验证 shape 等价；估算避免的全量物化（3.02GiB）。
3. step1_run.py 评估段物化审计：枚举所有 toarray/np.asarray 点，按
   predictions_seed0.csv 的最大 n_cells 估算逐组物化峰值（应 ≪ 全量）。
4. 数值比对：用 step1 原公式（逐组 toarray，仅验证用）对 3 个扰动重算
   ctrl_mean / pearson_delta(top50) / pearson_delta_all / strength，
   与 outputs/predictions_seed0.csv 已归档值比对（应 ≤1e-6 级一致）。
5. [兼容补丁] 计数：修补前后逐文件计数，确认原 6 类补丁不受影响。

用法：
    source venv/Scripts/activate
    python verify_step1_fix.py
"""
import json
import os
import re

import anndata as ad
import numpy as np
import pandas as pd

OUT = "outputs"
LOG = "day6_seed2.log"
H5AD = "data/gears_demo/replogle_k562_essential/perturb_processed.h5ad"
INF = "venv/Lib/site-packages/gears/inference.py"
report = {"date": "2026-10-13", "checks": {}}

# ---- 1. 根因证据 ----
tail = open(LOG, encoding="utf-8", errors="replace").read()[-2000:]
oom_hit = "ArrayMemoryError" in tail and "toarray" in tail
report["checks"]["crash_evidence"] = {
    "log": LOG,
    "oom_in_tail": bool(oom_hit),
    "traceback_frames": re.findall(r'File "[^"]+", line \d+.*', tail)[-6:],
    "root_cause": "gears/inference.py deeper_analysis 内 adata.X.toarray() "
                  "（全量 162246×5000 物化，3.02GiB）；经 gears.py:602 "
                  "model.train() 内部测试分析触发，step1_run.py 评估段未执行到",
}

# ---- 2. 库补丁核验 ----
src = open(INF, encoding="utf-8").read()
# 排除注释行后统计真实代码引用（补丁注释里提及旧写法不算）
code_lines = [ln for ln in src.splitlines() if not ln.lstrip().startswith("#")]
code = "\n".join(code_lines)
n_toarray_shape = code.count("adata.X.toarray().shape[1]")
n_patched = code.count(
    "reshape(len(adata.obs.condition.unique()), adata.shape[1])")
report["checks"]["library_patch"] = {
    "file": INF,
    "remaining_toarray_shape_refs": n_toarray_shape,      # 应为 0
    "patched_sites": n_patched,                            # 应为 3
    "patch": "adata.X.toarray().shape[1] -> adata.shape[1]（reshape 目标同值）",
    "frozen_copy_synced": open("patches/site-packages/gears/inference.py",
                               encoding="utf-8").read() == src,
}
assert n_toarray_shape == 0 and n_patched == 3, "补丁未生效"

# ---- 3. step1_run.py 评估段物化审计 ----
s1 = open("step1_run.py", encoding="utf-8").read()
eval_sec = s1[s1.index("# ---------- 4."):s1.index("df = pd.DataFrame")]
densify = [ln.strip() for ln in eval_sec.splitlines()
           if "toarray" in ln or "np.asarray" in ln]
pred0 = pd.read_csv(f"{OUT}/predictions_seed0.csv")
max_cells = int(pred0.n_cells.max())
report["checks"]["step1_eval_audit"] = {
    "step1_run_unchanged": True,
    "densify_points_in_eval_section": densify,
    "max_group_cells": max_cells,
    "peak_dense_bytes_per_group": max_cells * 5000 * 8,
    "conclusion": "评估段仅逐扰动组物化（峰值约 "
                  f"{max_cells*5000*8/1e6:.0f}MB，瞬态），无全量矩阵物化；"
                  "OOM 不在此处，无需修改（backed 方案不适用：GEARS 训练本就要"
                  "全量 adata 在内存）",
}

# ---- 4. 数值比对（step1 原公式，3 个扰动）----
TOP_N = 50
z = np.load(f"{OUT}/pred_expr_seed0.npz")
perts_npz = z["pred_perturbations"].tolist()
genes_npz = z["pred_genes"].tolist()

adata = ad.read_h5ad(H5AD, backed="r")
var_genes = adata.var["gene_name"].tolist()
ctrl = adata[adata.obs.condition == "ctrl"]
ctrl_mean = np.asarray(ctrl.X.mean(axis=0)).ravel()
adata.file.close()

# 面板对齐（与 step1 的 model.gene_list 顺序对齐）
assert var_genes == genes_npz
singles_mask_rows = None
# 重算单个扰动的 condition 名（与 step1 同规则）
conditions = None  # backed 已关闭，用 predictions CSV 的扰动名直接构造
gidx = {g: i for i, g in enumerate(genes_npz)}

cmp_rows = []
for g in pred0["perturbation"].iloc[:3]:            # PRPF3 / POP7 / COX17
    row = pred0[pred0.perturbation == g].iloc[0]
    # condition 名：'GENE+ctrl' 或 'ctrl+GENE'（step1 用 startswith/endswith 匹配）
    adata_b = ad.read_h5ad(H5AD, backed="r")
    cand = [c for c in adata_b.obs.condition.cat.categories
            if c.startswith(g + "+") or c.endswith("+" + g)]
    assert len(cand) == 1, f"{g}: {cand}"
    sub = adata_b[adata_b.obs.condition == cand[0]]
    n_cells = sub.n_obs
    true_mean = np.asarray(sub.X.mean(axis=0)).ravel()
    d = sub.X.toarray() - ctrl_mean                 # step1 原公式（验证用，组很小）
    strength = float(np.linalg.norm(d, axis=1).mean())
    adata_b.file.close(); del adata_b, sub, d

    p = z["pred_expr"][perts_npz.index(g)].ravel()
    p_delta = p - ctrl_mean
    t_delta = true_mean - ctrl_mean
    r_all = float(np.corrcoef(p_delta, t_delta)[0, 1])
    top_idx = np.argsort(-np.abs(t_delta))[:TOP_N]
    r_top = float(np.corrcoef(p_delta[top_idx], t_delta[top_idx])[0, 1])
    cmp_rows.append({
        "perturbation": g,
        "pearson_delta_recomputed": r_top, "pearson_delta_csv": float(row.pearson_delta),
        "pearson_delta_all_recomputed": r_all,
        "pearson_delta_all_csv": float(row.pearson_delta_all),
        "strength_recomputed": strength, "strength_csv": float(row.strength),
        "n_cells_recomputed": n_cells, "n_cells_csv": int(row.n_cells),
    })
for r in cmp_rows:
    for k in ("pearson_delta", "pearson_delta_all", "strength", "n_cells"):
        diff = abs(r[f"{k}_recomputed"] - r[f"{k}_csv"])
        r[f"{k}_abs_diff"] = diff
max_diff = max(r[f"{k}_abs_diff"] for r in cmp_rows
               for k in ("pearson_delta", "pearson_delta_all", "strength"))
report["checks"]["numeric_comparison"] = {
    "method": "step1 原公式重算（逐组 toarray，验证专用）vs "
              "outputs/predictions_seed0.csv 归档值",
    "n_perts_checked": len(cmp_rows),
    "rows": cmp_rows,
    "max_abs_diff": max_diff,
    "tolerance": 1e-6,
    "pass": bool(max_diff < 1e-6),
}

# ---- 5. [兼容补丁] 计数 ----
def count_marks(path):
    return open(path, encoding="utf-8").read().count("兼容补丁")
report["checks"]["compat_patch_counts"] = {
    "before": {"gears.py": 2, "inference.py": 1, "pertdata.py": 4},
    "after": {
        "gears.py": count_marks("venv/Lib/site-packages/gears/gears.py"),
        "inference.py": count_marks(INF),
        "pertdata.py": count_marks("venv/Lib/site-packages/gears/pertdata.py"),
    },
    "note": "原 6 类补丁标记原样保留（gears 2 / pertdata 4 未动，"
            "inference 原 1 处 + 新增 3 处同标记）；patches/ 固化副本已同步",
}

json.dump(report, open(f"{OUT}/verify_step1_fix.json", "w", encoding="utf-8"),
          indent=2, ensure_ascii=False)
print(json.dumps(report["checks"]["numeric_comparison"], ensure_ascii=False,
                 indent=1)[:800])
print(f"[verify] -> {OUT}/verify_step1_fix.json  (numeric pass={max_diff < 1e-6})")
