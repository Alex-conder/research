# step22_gb2_residual.py — G-B2：错误家族结构的定量互补性分析（残差 AUROC）
# 预注册：prereg_gb2_residual.md（2026-10-31，运行前锁定）。
#
# 设计（锁定，详见预注册）：
#   1) 共享难度成分提取：每组合误差秩 r_B（组合内 rank(err)/n）对"除源引擎 A
#      外全部引擎"的 OOF 探针分数 OLS（含截距）→ 残差 resid_B^(−A)
#      （留源口径：A 若入回归元则正交性使 corr(probe_A, resid)=0 按构造成立，
#      检验 vacuous；预注册 §2 已记录）。
#   2) 互补性检验：AUROC(probe_A, 1[resid_B^(−A) > 0])，扰动级 bootstrap
#      B=1000 的 95% CI；CI 完全在 0.5 之上 = 正向显著，完全在下 = 负向显著
#      （"互反"的定量形态）。逐家族对矩阵：家族内（表达谱/embedding）+
#      跨家族（双向）。
#   3) 判定（锁定）：主分析（K562/RPE1 全覆盖 5 引擎，跨家族有序对 24 对）中
#      存在 ≥1 对 AUROC ≥ 0.6 且 CI 完全在 0.5 之上 → 跨家族互补定量成立；
#      全无 → 家族结构即全部故事（如实收敛，R2 结论不变）。
#   描述性（不判定）：共享难度占比 = 留自回归 R²；scf_k562 45 扰动子集诊断。
#
# 防泄漏：探针分数全部 5-fold OOF（复用 G-B 规格）；残差 OLS in-sample
# （p ≤ 5 + 截距，n ≥ 45，空间可忽略，预注册 §3 记录）；bootstrap 只重采样
# 扰动索引。
#
# 产物：outputs_rebuild/gb2_residual.{json,md,png}
#       outputs_rebuild/smoke/gb2_residual.smoke.{json,png}（--smoke）
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step22_gb2_residual.py --smoke
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step22_gb2_residual.py
import argparse
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import roc_auc_score

for _s in ("stdout", "stderr"):  # GBK 控制台纪律：不崩、不丢字符
    getattr(sys, _s).reconfigure(encoding="utf-8", errors="replace")
os.environ.setdefault("OMP_NUM_THREADS", "1")

REPO = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(REPO, "outputs_rebuild")
SMOKE_DIR = os.path.join(OUT_DIR, "smoke")
WAVE_DIR = os.path.join(REPO, "outputs_crossmodel", "wave")

BOOT_MASTER_SEED = 20261031
FULL_B = 1000
SMOKE_B = 100

CRITERION = ("primary (K562+RPE1 full-coverage cross-family ordered pairs): "
             "exists >=1 pair with residual AUROC >= 0.6 and 95% CI entirely "
             "above 0.5 -> cross-family complementarity quantitatively "
             "established; else family structure is the whole story "
             "(honest convergence, R2 conclusion unchanged)")

from step10_signature_resolution import COMBOS
from step10_two_layer_gate import oof_probe_scores

FULL_KEYS = {
    "K562": ["gears_k562", "scgenept_k562", "scgpt_k562",
             "geneformer_k562", "uce_k562"],
    "RPE1": ["gears_rpe1", "scgenept_rpe1", "scgpt_rpe1",
             "geneformer_rpe1", "uce_rpe1"],
}
SCF_KEY = "scf_k562"
ENGINE = {k: COMBOS[k]["engine"] for k in COMBOS}
CLS = {k: COMBOS[k]["cls"] for k in COMBOS}


# ---------------------------------------------------------------- 加载（含连续误差）
def load_combo_full(key):
    """镜像 step17.load_combo_xy_perts，额外返回连续误差 err。
    返回 (X, y, err, feature_names, perts)，行序与 perts 一一对应。"""
    c = COMBOS[key]
    if c["kind"] == "gears":
        from step6d_naive_baselines import load_dataset
        feat, feature_names, y, pred_df, _ = load_dataset(c["data"], c["out"])
        X = feat[feature_names].values.astype(np.float64)
        err = 1.0 - pred_df["pearson_delta"].values.astype(np.float64)
        perts = [str(g) for g in pred_df.index]
    elif c["kind"] == "crossmodel":
        from step6f_crossmodel import load_crossmodel
        X, y, feature_names, perts, _extra = load_crossmodel(c["npz"],
                                                             c["labels"])
        lab = pd.read_csv(c["labels"]).set_index("perturbation")
        err = 1.0 - lab.loc[[str(p) for p in perts]]["pearson_delta"].values
        err = err.astype(np.float64)
        perts = [str(p) for p in perts]
    elif c["kind"] == "wave":
        lab = pd.read_csv(os.path.join(
            WAVE_DIR, c["wave"],
            f"{c['wave']}_wave_{'k562' if c['dataset'] == 'K562' else 'rpe1'}"
            f"_labels.csv"))
        err = lab["err_emb_dist"].values.astype(np.float64)
        assert np.isfinite(err).all(), f"{key}: err_emb_dist 含非有限值"
        y = (err > np.median(err)).astype(int)
        in_space = lab["target_in_input_space"].astype(bool).astype(np.float64)
        X = np.column_stack([
            lab["pred_emb_norm"].values.astype(np.float64),
            lab["delta_emb_norm"].values.astype(np.float64),
            np.log1p(lab["n_cells"].values.astype(np.float64)),
            in_space,
        ])
        feature_names = ["pred_emb_norm", "delta_emb_norm", "log_n_cells",
                         "target_in_input_space"]
        perts = [str(g) for g in lab["perturbation"]]
    else:
        raise FileNotFoundError(f"{key}: 产物缺失（注册为 missing）")
    assert len(perts) == len(y) == len(X) == len(err), f"{key}: 行不对齐"
    assert np.isfinite(err).all(), f"{key}: 误差含非有限值"
    return X, y, err, list(feature_names), perts


def common_perts(perts_by_key):
    """全部组合扰动名交集（排序），矩阵口径统一。"""
    sets = [set(v) for v in perts_by_key.values()]
    return sorted(set.intersection(*sets))


def align(vec, perts, common):
    m = {p: v for p, v in zip(perts, vec)}
    return np.array([m[p] for p in common], dtype=np.float64)


def ols_residual(y_vec, X_reg):
    """OLS（含截距）残差与 R²。X_reg 不含截距列。"""
    X1 = np.column_stack([np.ones(len(y_vec)), X_reg])
    beta, _, _, _ = np.linalg.lstsq(X1, y_vec, rcond=None)
    fit = X1 @ beta
    resid = y_vec - fit
    ss_res = float((resid ** 2).sum())
    ss_tot = float(((y_vec - y_vec.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return resid, r2


def boot_auroc_ci(scores, z, b_boot, rng):
    """扰动级 bootstrap（向量固定，只重采样索引），95% percentile CI。"""
    n = len(z)
    vals = np.empty(b_boot)
    for b in range(b_boot):
        idx = rng.integers(0, n, n)
        while len(np.unique(z[idx])) < 2:
            idx = rng.integers(0, n, n)
        vals[b] = roc_auc_score(z[idx], scores[idx])
    lo, hi = np.percentile(vals, [2.5, 97.5])
    return {"ci_lo": float(lo), "ci_hi": float(hi),
            "above_0_5": bool(lo > 0.5), "below_0_5": bool(hi < 0.5)}


def family_pair(src_cls, tgt_cls):
    if src_cls == tgt_cls:
        return f"within_{src_cls}"
    return "cross_expression_to_embedding" \
        if src_cls == "expression" else "cross_embedding_to_expression"


# ---------------------------------------------------------------- 数据集块
def run_dataset_block(ds_tag, keys, loaded, err_rank, probes, b_boot,
                      rng_master):
    """keys: 该数据集全部组合；loaded/err_rank/probes 已按 common 对齐。"""
    common = common_perts({k: loaded[k][4] for k in keys})
    n = len(common)
    eng_keys = list(keys)
    P = {k: align(probes[k], loaded[k][4], common) for k in eng_keys}
    R = {k: align(err_rank[k], loaded[k][4], common) for k in eng_keys}

    # 描述性：共享难度占比（留自回归：除 B 自身外全部引擎）
    shared = {}
    for k in eng_keys:
        others = [P[o] for o in eng_keys if o != k]
        _, r2 = ols_residual(R[k], np.column_stack(others))
        shared[k] = {"r2_leave_self_out": r2,
                     "specific_fraction": 1.0 - r2}

    pairs = {}
    for b_key in eng_keys:
        resid_cache = {}
        for a_key in eng_keys:
            if a_key == b_key:
                continue
            if a_key not in resid_cache:
                regs = [P[o] for o in eng_keys if o != a_key]
                resid_cache[a_key], _ = ols_residual(
                    R[b_key], np.column_stack(regs))
            resid = resid_cache[a_key]
            z = (resid > 0).astype(int)
            pk = f"{a_key}__{b_key}"
            if len(np.unique(z)) < 2:
                pairs[pk] = {"source": a_key, "target": b_key,
                             "source_engine": ENGINE[a_key],
                             "target_engine": ENGINE[b_key],
                             "source_cls": CLS[a_key], "target_cls": CLS[b_key],
                             "family_pair": family_pair(CLS[a_key], CLS[b_key]),
                             "n": n, "auroc": None,
                             "note": "resid binary single-class"}
                continue
            auroc = float(roc_auc_score(z, P[a_key]))
            rng = np.random.default_rng(rng_master.integers(0, 2**31))
            ci = boot_auroc_ci(P[a_key], z, b_boot, rng)
            pairs[pk] = {"source": a_key, "target": b_key,
                         "source_engine": ENGINE[a_key],
                         "target_engine": ENGINE[b_key],
                         "source_cls": CLS[a_key], "target_cls": CLS[b_key],
                         "family_pair": family_pair(CLS[a_key], CLS[b_key]),
                         "n": n, "auroc": auroc,
                         "deviation_from_0_5": auroc - 0.5,
                         "b_boot": b_boot, **ci}
    print(f"[GB2][{ds_tag}] n_common={n}, {len(pairs)} ordered pairs, "
          f"shared-difficulty R2 (leave-self): "
          + ", ".join(f"{k}={shared[k]['r2_leave_self_out']:.3f}"
                      for k in eng_keys))
    return {"engines": [ENGINE[k] for k in eng_keys], "keys": eng_keys,
            "n_common": n, "shared_difficulty": shared, "pairs": pairs}


def family_block_summary(block):
    agg = {}
    for pk, r in block["pairs"].items():
        if r["auroc"] is None:
            continue
        fp = r["family_pair"]
        d = agg.setdefault(fp, {"pairs": 0, "mean_auroc": 0.0,
                                "mean_abs_dev": 0.0, "n_above": 0,
                                "n_below": 0, "n_pass_locked": 0})
        d["pairs"] += 1
        d["mean_auroc"] += r["auroc"]
        d["mean_abs_dev"] += abs(r["deviation_from_0_5"])
        d["n_above"] += int(r["above_0_5"])
        d["n_below"] += int(r["below_0_5"])
        d["n_pass_locked"] += int(r["auroc"] >= 0.6 and r["above_0_5"])
    for fp, d in agg.items():
        d["mean_auroc"] /= d["pairs"]
        d["mean_abs_dev"] /= d["pairs"]
    return agg


# ---------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true",
                    help="迷你规格（B=100，仅 K562）")
    args = ap.parse_args()

    smoke = args.smoke
    b_boot = SMOKE_B if smoke else FULL_B
    ds_list = ["K562"] if smoke else ["K562", "RPE1"]

    t_all = time.perf_counter()
    need = set()
    for ds in ds_list:
        need.update(FULL_KEYS[ds])
    if not smoke:
        need.add(SCF_KEY)
    loaded = {k: load_combo_full(k) for k in sorted(need)}
    for k in sorted(loaded):
        print(f"[LOAD] {k}: n={len(loaded[k][1])} cls={CLS[k]}")

    probes = {k: oof_probe_scores(v[0], v[1]) for k, v in loaded.items()}
    err_rank = {k: rankdata(v[2]) / len(v[2]) for k, v in loaded.items()}
    print(f"[BASE] OOF probe scores + error ranks: {len(probes)} combos")

    boot_master = np.random.default_rng(BOOT_MASTER_SEED)
    res = {"task": "G-B2 quantitative complementarity analysis "
                   "(step22, prereg_gb2_residual.md)",
           "prereg": "prereg_gb2_residual.md",
           "smoke": smoke,
           "spec": {
               "error_rank": "within-combo rank percentile rank(err)/n; "
                             "err = 1 - pearson_delta (expression), "
                             "err_emb_dist (embedding)",
               "residual": "leave-source-out OLS (intercept): err_rank_B ~ "
                           "probe scores of all dataset engines except A; "
                           "target B's own probe included",
               "test": "AUROC(probe_A, 1[resid_B(-A) > 0]), perturbation "
                       "bootstrap, 95% percentile CI",
               "b_boot": b_boot,
               "full_run_spec": {"b_boot": FULL_B}},
           "boot_master_seed": BOOT_MASTER_SEED,
           "criterion": CRITERION,
           "datasets": {}}

    for ds in ds_list:
        keys = FULL_KEYS[ds]
        block = run_dataset_block(ds, keys, loaded, err_rank, probes,
                                  b_boot, boot_master)
        block["family_blocks"] = family_block_summary(block)
        res["datasets"][ds] = block

    if not smoke:
        # scf_k562 诊断：45 扰动子集，6 引擎
        keys6 = FULL_KEYS["K562"] + [SCF_KEY]
        block = run_dataset_block("scf_k562_diag", keys6, loaded, err_rank,
                                  probes, b_boot, boot_master)
        block["family_blocks"] = family_block_summary(block)
        block["diagnostic_only"] = True
        res["datasets"]["scf_k562_diag"] = block

    # ---- 判定（锁定）----
    primary = []
    for ds in ds_list:
        for pk, r in res["datasets"][ds]["pairs"].items():
            if r["family_pair"].startswith("cross"):
                primary.append((ds, pk, r))
    n_pass = sum(1 for _, _, r in primary
                 if r["auroc"] is not None and r["auroc"] >= 0.6
                 and r["above_0_5"])
    n_sig_above = sum(1 for _, _, r in primary if r.get("above_0_5"))
    n_sig_below = sum(1 for _, _, r in primary if r.get("below_0_5"))
    n_valid = sum(1 for _, _, r in primary if r["auroc"] is not None)
    if smoke:
        verdict = {"tier": "smoke_only_not_adjudicated",
                   "note": "smoke numbers do not enter the locked criteria",
                   "n_primary_pairs": len(primary), "n_pass_smoke": n_pass}
    else:
        verdict = {"rule": CRITERION,
                   "n_primary_pairs": len(primary),
                   "n_valid_pairs": n_valid,
                   "n_pass_auroc_ge_0.6_and_ci_above_0.5": n_pass,
                   "n_ci_above_0_5": n_sig_above,
                   "n_ci_below_0_5": n_sig_below,
                   "threshold": "auroc >= 0.6 AND ci_lo > 0.5 (locked)"}
        verdict["tier"] = (
            "cross_family_complementarity_residual" if n_pass >= 1
            else "family_structure_is_the_whole_story")
    res["verdict"] = verdict
    res["timing_total_s"] = round(time.perf_counter() - t_all, 2)

    out_dir = SMOKE_DIR if smoke else OUT_DIR
    os.makedirs(out_dir, exist_ok=True)
    stem = "gb2_residual.smoke" if smoke else "gb2_residual"
    out_json = os.path.join(out_dir, stem + ".json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("[SAVED]", out_json)
    make_figure(res, os.path.join(out_dir, stem + ".png"))
    if not smoke:
        write_verdict_md(res, os.path.join(OUT_DIR, "gb2_residual.md"))
    print(f"[DONE] total {res['timing_total_s']}s | tier {verdict.get('tier')}")


def auroc_mat(block):
    """有序对 AUROC 方阵：行 = 源，列 = 目标，对角 None。"""
    keys = block["keys"]
    mat = np.full((len(keys), len(keys)), np.nan)
    for pk, r in block["pairs"].items():
        i, j = keys.index(r["source"]), keys.index(r["target"])
        if r["auroc"] is not None:
            mat[i, j] = r["auroc"]
    return mat


def make_figure(res, path):
    ds_blocks = [(d, res["datasets"][d]) for d in ("K562", "RPE1")
                 if d in res["datasets"]]
    if not res.get("smoke") and "scf_k562_diag" in res["datasets"]:
        ds_blocks.append(("scf_k562 (diag)", res["datasets"]["scf_k562_diag"]))
    n_p = len(ds_blocks)
    if not n_p:
        return
    fig, axes = plt.subplots(1, n_p, figsize=(6.2 * n_p, 5.2))
    if n_p == 1:
        axes = [axes]
    for ax, (tag, block) in zip(axes, ds_blocks):
        mat = auroc_mat(block)
        keys = block["keys"]
        im = ax.imshow(mat, vmin=0.3, vmax=0.7, cmap="RdYlGn_r")
        ax.set_xticks(range(len(keys)))
        ax.set_yticks(range(len(keys)))
        ax.set_xticklabels([k.replace("_k562", "").replace("_rpe1", "")
                            for k in keys], fontsize=7, rotation=45,
                           ha="right")
        ax.set_yticklabels([k.replace("_k562", "").replace("_rpe1", "")
                            for k in keys], fontsize=7)
        for i in range(len(keys)):
            for j in range(len(keys)):
                if i == j:
                    ax.text(j, i, "-", ha="center", va="center", fontsize=8)
                elif np.isfinite(mat[i, j]):
                    r = block["pairs"][f"{keys[i]}__{keys[j]}"]
                    mark = "*" if (r.get("above_0_5") or r.get("below_0_5")) \
                        else ""
                    ax.text(j, i, f"{mat[i, j]:.2f}{mark}", ha="center",
                            va="center", fontsize=6.5)
        ax.set_title(f"residual AUROC (source rows -> target cols), {tag} "
                     f"(n={block['n_common']})", fontsize=9)
        ax.set_xlabel("target engine", fontsize=8)
        ax.set_ylabel("source engine", fontsize=8)
    fig.colorbar(im, ax=axes[-1], fraction=0.03, pad=0.02,
                 label="AUROC (0.5 = no residual info)")
    fig.suptitle("G-B2 complementarity: predicting residual errors "
                 f"(B={res['spec']['b_boot']}, smoke={res['smoke']})",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 0.99, 0.95])
    fig.savefig(path, dpi=300)
    print("[SAVED]", path)


def fmt(v, nd=3):
    return "n/a" if v is None else f"{v:.{nd}f}"


def write_verdict_md(res, path):
    v = res["verdict"]
    lines = ["# G-B2 定量互补性分析 verdict（step22，prereg_gb2_residual.md）",
             ""]
    if v.get("tier") == "smoke_only_not_adjudicated":
        lines.append("smoke 运行，不进锁定判定。")
    else:
        lines.append(
            f"判定：主分析跨家族有序对 {v['n_primary_pairs']} 对中，"
            f"AUROC ≥ 0.6 且 95% CI 完全在 0.5 之上的对数 = "
            f"{v['n_pass_auroc_ge_0.6_and_ci_above_0.5']}（判据 ≥1）-> "
            f"**{v['tier']}**。")
        lines.append("")
        lines.append(f"- CI 完全在 0.5 之上（正向显著）：{v['n_ci_above_0_5']} 对；"
                     f"完全在下（负向显著，'互反'定量形态）："
                     f"{v['n_ci_below_0_5']} 对。")
        lines.append("- 负向显著的对按预注册 §4 记为描述性'互反可预测'，"
                     "不进主判定。")
        lines.append("")
        for ds in ("K562", "RPE1", "scf_k562_diag"):
            if ds not in res["datasets"]:
                continue
            block = res["datasets"][ds]
            tag = "（诊断，不进判定）" if ds == "scf_k562_diag" else ""
            lines.append(f"## 残差 AUROC 矩阵（{ds}{tag}，n={block['n_common']}）")
            lines.append("")
            lines.append("行 = 源引擎探针，列 = 目标引擎残差；"
                         "* = 95% CI 不含 0.5。")
            lines.append("")
            keys = block["keys"]
            short = [k.replace("_k562", "").replace("_rpe1", "") for k in keys]
            lines.append("| src\\\\tgt | " + " | ".join(short) + " |")
            lines.append("|" + "---|" * (len(keys) + 1))
            for a in keys:
                row = [a.replace("_k562", "").replace("_rpe1", "")]
                for b in keys:
                    if a == b:
                        row.append("-")
                        continue
                    r = block["pairs"][f"{a}__{b}"]
                    if r["auroc"] is None:
                        row.append("n/a")
                    else:
                        mark = "*" if (r.get("above_0_5")
                                       or r.get("below_0_5")) else ""
                        row.append(f"{r['auroc']:.3f}{mark}")
                lines.append("| " + " | ".join(row) + " |")
            lines.append("")
            fb = block["family_blocks"]
            lines.append("| family pair | pairs | mean AUROC | mean |dev| | "
                         "sig(+) | sig(-) | pass(locked) |")
            lines.append("|---|---|---|---|---|---|---|")
            for fp, d in sorted(fb.items()):
                lines.append(
                    f"| {fp} | {d['pairs']} | {d['mean_auroc']:.3f} | "
                    f"{d['mean_abs_dev']:.3f} | {d['n_above']} | "
                    f"{d['n_below']} | {d['n_pass_locked']} |")
            lines.append("")
            sd = block["shared_difficulty"]
            lines.append("共享难度占比（留自回归 R²）：" + "；".join(
                f"{k}={d['r2_leave_self_out']:.3f}" for k, d in sd.items()))
            lines.append("")
        if v["tier"] == "cross_family_complementarity_residual":
            lines.append("## 升级表述（锁定条款触发）")
            lines.append("")
            lines.append("存在跨家族组合残差 AUROC ≥ 0.6 且 CI 不含 0.5："
                         "G-B 节升级为'集成受家族结构限制但互补信息存在"
                         "（残差证据）'。")
        else:
            lines.append("## 收敛表述（锁定阴性条款触发）")
            lines.append("")
            lines.append("主分析跨家族对无一达 AUROC ≥ 0.6 且 CI 在 0.5 之上："
                         "家族结构即全部故事，如实收敛；R2 结论"
                         "（产品 = 符号层）不变，与 error_signatures 深潜"
                         "'签名噪声'互证。上表负向显著格（若有）为 kappa 结构"
                         "'跨家族互反'的定量回声，按预注册为描述性记录。")
        lines.append("")
        lines.append("注：主检验共 24 对（每数据集 12 = 表达谱 2 × embedding 3 × "
                     "双向 2），多重性未进判据（任务书锁定 '>=1 对'），"
                     "显著对数与家族块均值已在上方如实列出。")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    print("[SAVED]", path)


if __name__ == "__main__":
    main()
