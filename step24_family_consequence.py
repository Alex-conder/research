# step24_family_consequence.py -- Family-Consequence: 三条有后果的设计原则
# 预注册：prereg_family_consequence.md（2026-10-30，运行前锁定）。
#
#   P1 省钱原则：家族内全引擎集成门 vs 家族代表单引擎门的效用差（G-C 既定
#      成本口径，回放于固定目标引擎 err_T）；主判定单元 = (D, T, R) 上 4 个
#      预算档 Δ 的均值，配对 bootstrap 95% CI 全部含 0 -> 冗余推论成立。
#   P2 异质原则：表达谱家族 vs embedding 家族门分数的秩分歧度 d 作为第三
#      信号，AUROC(d, 家族平均误差秩最差一半) 点估计 >= 0.6（K562 与 RPE1）。
#   P3 选择原则：on-file 就地 AUROC 档案（ensemble_gate.json）+ 数据集属性，
#      LODO 回放家族感知规则 R1 vs 家族无感基线 R0（元学习雏形）。
#
# 防泄漏：全部 on-file 回放，分数固定（OOF 探针确定性重算，同 G-B2），
# bootstrap 只重采样扰动索引；P3 不 bootstrap。
#
# 产物：outputs_rebuild/family_consequence.{json,png} + family_consequence_verdict.md
#       smoke -> outputs_rebuild/smoke/family_consequence.smoke.*
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step24_family_consequence.py --smoke
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step24_family_consequence.py
import argparse
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import rankdata, spearmanr
from sklearn.metrics import roc_auc_score

for _s in ("stdout", "stderr"):  # GBK 控制台纪律：不崩、不丢字符
    getattr(sys, _s).reconfigure(encoding="utf-8", errors="replace")
os.environ.setdefault("OMP_NUM_THREADS", "1")

REPO = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(REPO, "outputs_rebuild")
SMOKE_DIR = os.path.join(OUT_DIR, "smoke")

BOOT_MASTER_SEED = 20261101
FULL_B = 1000
SMOKE_B = 100
BUDGETS = (5, 10, 20, 50)

from step22_gb2_residual import (  # noqa: E402
    CLS, ENGINE, SCF_KEY, align, common_perts, load_combo_full)
from step10_two_layer_gate import oof_probe_scores  # noqa: E402

FULL_KEYS = {
    "K562": ["gears_k562", "scgenept_k562", "scgpt_k562",
             "geneformer_k562", "uce_k562"],
    "RPE1": ["gears_rpe1", "scgenept_rpe1", "scgpt_rpe1",
             "geneformer_rpe1", "uce_rpe1"],
}
FAMILY_KEYS = {  # 主分析家族划分（scf 不入主分析，同 G-B2 先例）
    "K562": {"expression": ["gears_k562", "scgenept_k562"],
             "embedding": ["scgpt_k562", "geneformer_k562", "uce_k562"]},
    "RPE1": {"expression": ["gears_rpe1", "scgenept_rpe1"],
             "embedding": ["scgpt_rpe1", "geneformer_rpe1", "uce_rpe1"]},
}
SCF_FAMILY_KEYS = {"expression": ["gears_k562", "scgenept_k562", "scf_k562"],
                   "embedding": ["scgpt_k562", "geneformer_k562", "uce_k562"]}

P1_CRIT = ("primary cells (D in {K562,RPE1} x embedding family x T x R!=T; "
           "delta = mean over B in {5,10,20,50} of topk(err_T,s_R,B) - "
           "topk(err_T,s_{F\\T},B); paired bootstrap 95% CI of delta "
           "contains 0 in ALL 12 cells -> within-family redundancy has a "
           "cost consequence (deploy any single representative); any cell "
           "CI entirely above 0 -> P1 fails")
P2_CRIT = ("AUROC(d, z) point estimate >= 0.6 in BOTH K562 and RPE1 "
           "(d = |rank_pct(mean probe of expression family) - "
           "rank_pct(mean probe of embedding family)|; z = 1[family-mean "
           "error rank >= median]) -> family disagreement is a predictive "
           "third signal")
P3_CRIT = ("LODO over {K562,RPE1}: family-aware rule R1 top-1 hit in both "
           "folds AND R1 regret < family-blind R0 regret in both folds -> "
           "selection principle supported (prototype-grade, low power)")


# ---------------------------------------------------------------- 工具
def topk_sum(values, scores, B):
    """按 scores 降序取 top-B（k = min(B, n)）的 values 和。"""
    n = len(values)
    k = min(B, n)
    if k >= n:
        return float(values.sum())
    idx = np.argpartition(scores, -k)[-k:]
    return float(values[idx].sum())


def boot_p1_delta(err, s_full, s_r, b_boot, rng):
    """配对 bootstrap：4 档预算 Δ 的均值与 95% CI（重采样扰动索引）。"""
    n = len(err)
    bars = np.empty(b_boot)
    for b in range(b_boot):
        idx = rng.integers(0, n, n)
        e, sf, sr = err[idx], s_full[idx], s_r[idx]
        d = [topk_sum(e, sr, B) - topk_sum(e, sf, B) for B in BUDGETS]
        bars[b] = float(np.mean(d))
    lo, hi = np.percentile(bars, [2.5, 97.5])
    return {"mean_delta": float(bars.mean()), "ci_lo": float(lo),
            "ci_hi": float(hi), "ci_contains_zero": bool(lo <= 0.0 <= hi)}


def boot_auroc_ci(scores, z, b_boot, rng):
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


def short(k):
    return k.replace("_k562", "").replace("_rpe1", "")


# ---------------------------------------------------------------- P1
def run_p1_block(ds_tag, fam_keys, P, E, b_boot, rng_master, diagnostic):
    """fam_keys: {family: [keys]}；P/E 已按公共扰动对齐。

    返回 cells 列表。结构性（F∖{T} 单台）不跑 bootstrap，Δ≡0 记录。"""
    cells = []
    for fam, keys in fam_keys.items():
        for t in keys:
            others = [k for k in keys if k != t]
            s_full = np.mean([P[k] for k in others], axis=0)
            for r in others:
                s_r = P[r]
                k_eff = {B: min(B, len(E[t])) for B in BUDGETS}
                util_full, util_r = {}, {}
                for B in BUDGETS:
                    rand = k_eff[B] / len(E[t]) * float(E[t].sum())
                    util_full[B] = rand - topk_sum(E[t], s_full, B)
                    util_r[B] = rand - topk_sum(E[t], s_r, B)
                if len(others) == 1:
                    cell = {"dataset": ds_tag, "family": fam, "target": t,
                            "rep": r, "structural": True,
                            "note": "family has 2 members: fam_full == single_R",
                            "delta_by_budget": {str(B): 0.0 for B in BUDGETS},
                            "mean_delta": 0.0, "ci_lo": 0.0, "ci_hi": 0.0,
                            "ci_contains_zero": True,
                            "util_fam_full_vs_random": util_full,
                            "util_single_vs_random": util_r}
                else:
                    rng = np.random.default_rng(rng_master.integers(0, 2**31))
                    summ = boot_p1_delta(E[t], s_full, s_r, b_boot, rng)
                    cell = {"dataset": ds_tag, "family": fam, "target": t,
                            "rep": r, "structural": False, **summ,
                            "delta_by_budget": {
                                str(B): topk_sum(E[t], s_r, B)
                                - topk_sum(E[t], s_full, B) for B in BUDGETS},
                            "util_fam_full_vs_random": util_full,
                            "util_single_vs_random": util_r}
                cell["diagnostic"] = diagnostic
                cells.append(cell)
    return cells


def p1_verdict(primary_cells):
    n = len(primary_cells)
    n_pass = sum(1 for c in primary_cells if c["ci_contains_zero"])
    fails = [f"{c['dataset']}|{short(c['target'])}<-{short(c['rep'])}"
             for c in primary_cells if not c["ci_contains_zero"]]
    return {"rule": P1_CRIT, "n_primary_cells": n,
            "n_ci_contains_zero": n_pass,
            "all_contain_zero": n_pass == n, "failing_cells": fails,
            "tier": ("within_family_redundancy_has_cost_consequence"
                     if n_pass == n else "significant_utility_loss_detected")}


# ---------------------------------------------------------------- P2
def run_p2_block(ds_tag, fam_keys, all_keys, P, R, b_boot, rng_master,
                 diagnostic):
    s_expr = np.mean([P[k] for k in fam_keys["expression"]], axis=0)
    s_emb = np.mean([P[k] for k in fam_keys["embedding"]], axis=0)
    n = len(s_expr)
    d = np.abs(rankdata(s_expr) / n - rankdata(s_emb) / n)
    r_bar = np.mean([R[k] for k in all_keys], axis=0)
    z = (r_bar >= np.median(r_bar)).astype(int)
    auroc = float(roc_auc_score(z, d))
    rng = np.random.default_rng(rng_master.integers(0, 2**31))
    ci = boot_auroc_ci(d, z, b_boot, rng)
    rho = float(spearmanr(d, r_bar).statistic)
    # 描述性变体：各家族均值探针 vs 全体均值探针的 |秩差|
    s_all = np.mean([P[k] for k in all_keys], axis=0)
    r_all = rankdata(s_all) / n
    var_expr = np.abs(rankdata(s_expr) / n - r_all)
    var_emb = np.abs(rankdata(s_emb) / n - r_all)
    return {"dataset": ds_tag, "n": n, "auroc": auroc, **ci,
            "spearman_d_vs_rbar": rho,
            "variant_expr_vs_all_auroc": float(roc_auc_score(z, var_expr)),
            "variant_emb_vs_all_auroc": float(roc_auc_score(z, var_emb)),
            "diagnostic": diagnostic}


def p2_verdict(blocks):
    prim = {b["dataset"]: b for b in blocks if not b["diagnostic"]}
    ok = all(b["auroc"] >= 0.6 for b in prim.values())
    return {"rule": P2_CRIT,
            "auroc": {d: round(b["auroc"], 4) for d, b in prim.items()},
            "all_ge_0_6": ok,
            "tier": ("family_disagreement_is_predictive_third_signal" if ok
                     else "family_disagreement_not_predictive")}


# ---------------------------------------------------------------- P3
def run_p3():
    """on-file 档案 + LODO 回放。返回 (rows, folds, scf_note)。"""
    eg = json.load(open(os.path.join(OUT_DIR, "ensemble_gate.json"),
                        encoding="utf-8"))
    gb2 = json.load(open(os.path.join(OUT_DIR, "gb2_residual.json"),
                         encoding="utf-8"))
    r2_by_key = {}
    for blk in gb2["datasets"].values():
        for k, v in blk.get("shared_difficulty", {}).items():
            r2_by_key[k] = v["r2_leave_self_out"]
    rows = []
    for key, c in eg["combos"].items():
        a = c.get("auroc")
        if not isinstance(a, dict) or a.get("own_in_situ") is None:
            continue  # scf_rpe1 stub 等无记录行不进档案
        ds = c["dataset"]
        # scf 行 dataset 字段为 'K562'（n=45 子集）：按 key 归入诊断单元，
        # 不进主 LODO 折（引擎池/样本量不对齐）
        unit = "scf_k562_diag" if key == "scf_k562" else ds
        rows.append({"key": key, "engine": c["engine"], "dataset": ds,
                     "unit": unit, "family": c["class"], "n": c["n"],
                     "label_balance": c["label_balance"],
                     "shared_r2": r2_by_key.get(key),
                     "own_in_situ_auroc": a["own_in_situ"]})
    units = {"K562": [r for r in rows if r["unit"] == "K562"],
             "RPE1": [r for r in rows if r["unit"] == "RPE1"],
             "scf_k562_diag": [r for r in rows if r["unit"] == "scf_k562_diag"]}

    def predict(train_rows, rule):
        if rule == "R1":
            fam_mean = {}
            for fam in ("expression", "embedding"):
                fs = [r["own_in_situ_auroc"] for r in train_rows
                      if r["family"] == fam]
                if fs:
                    fam_mean[fam] = float(np.mean(fs))
            best_fam = max(fam_mean, key=fam_mean.get)
            pool = [r for r in train_rows if r["family"] == best_fam]
        else:
            pool = list(train_rows)
        eng_mean = {}
        for r in pool:
            eng_mean.setdefault(r["engine"], []).append(r["own_in_situ_auroc"])
        return max(eng_mean, key=lambda e: float(np.mean(eng_mean[e])))

    folds = []
    for test_u, train_u in (("K562", "RPE1"), ("RPE1", "K562")):
        test_rows, train_rows = units[test_u], units[train_u]
        best = max(test_rows, key=lambda r: r["own_in_situ_auroc"])
        fold = {"train_unit": train_u, "test_unit": test_u,
                "test_best_engine": best["engine"],
                "test_best_auroc": best["own_in_situ_auroc"]}
        for rule in ("R1", "R0"):
            pred_eng = predict(train_rows, rule)
            pred_row = next(r for r in test_rows if r["engine"] == pred_eng)
            fold[rule] = {
                "predicted_engine": pred_eng,
                "top1_hit": pred_eng == best["engine"],
                "regret": best["own_in_situ_auroc"]
                - pred_row["own_in_situ_auroc"]}
        fold["random_reference_regret"] = best["own_in_situ_auroc"] - float(
            np.mean([r["own_in_situ_auroc"] for r in test_rows]))
        folds.append(fold)
    # 描述性：scf 行是否会被各训练单元选中
    scf_note = {"scf_auroc": units["scf_k562_diag"][0]["own_in_situ_auroc"]
                if units["scf_k562_diag"] else None}
    for train_u in ("K562", "RPE1"):
        for rule in ("R1", "R0"):
            scf_note[f"{rule}_train_{train_u}_picks"] = predict(
                units[train_u], rule)
    ok = all(f["R1"]["top1_hit"] and
             f["R1"]["regret"] < f["R0"]["regret"] for f in folds)
    verdict = {"rule": P3_CRIT,
               "prototype_low_power": True,
               "r1_top1_hits": sum(int(f["R1"]["top1_hit"]) for f in folds),
               "n_folds": len(folds),
               "r1_regret_lt_r0_both_folds": all(
                   f["R1"]["regret"] < f["R0"]["regret"] for f in folds),
               "tier": ("selection_principle_supported" if ok
                        else "selection_principle_not_supported")}
    return rows, folds, scf_note, verdict


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
    run_scf = not smoke
    if run_scf:
        need.update(FULL_KEYS["K562"])
        need.add(SCF_KEY)
    loaded = {k: load_combo_full(k) for k in sorted(need)}
    probes = {k: oof_probe_scores(v[0], v[1]) for k, v in loaded.items()}
    err = {k: v[2].astype(np.float64) for k, v in loaded.items()}
    err_rank = {k: rankdata(v[2]) / len(v[2]) for k, v in loaded.items()}
    print(f"[BASE] loaded {len(loaded)} combos, OOF probes ready")

    rng_master = np.random.default_rng(BOOT_MASTER_SEED)
    res = {"task": "Family-Consequence: three consequential design "
                   "principles from the architecture-family error structure "
                   "(step24, prereg_family_consequence.md)",
           "prereg": "prereg_family_consequence.md",
           "smoke": smoke,
           "spec": {"p1": "fam_full (mean probe of F\\{T}) vs single_R on "
                          "fixed target err_T; utility = G-C cost savings vs "
                          "random at B in {5,10,20,50}; delta per (T,R) = "
                          "mean over B; paired perturbation bootstrap "
                          f"B={b_boot}, master seed {BOOT_MASTER_SEED}",
                    "p2": "d = |rank_pct(mean expression-family probe) - "
                          "rank_pct(mean embedding-family probe)|; "
                          "z = 1[family-mean error rank >= median]; "
                          "AUROC + bootstrap CI",
                    "p3": "on-file own_in_situ AUROC registry "
                          "(ensemble_gate.json) + G-B2 shared-difficulty R2; "
                          "LODO over {K562,RPE1}, rules R1 (family-aware) "
                          "vs R0 (family-blind)",
                    "full_run_spec": {"b_boot": FULL_B}},
           "boot_master_seed": BOOT_MASTER_SEED,
           "budgets": list(BUDGETS)}

    # ---- P1 ----
    p1_cells = []
    for ds in ds_list:
        keys = FULL_KEYS[ds]
        common = common_perts({k: loaded[k][4] for k in keys})
        P = {k: align(probes[k], loaded[k][4], common) for k in keys}
        E = {k: align(err[k], loaded[k][4], common) for k in keys}
        p1_cells += run_p1_block(ds, FAMILY_KEYS[ds], P, E, b_boot,
                                 rng_master, diagnostic=False)
        print(f"[P1][{ds}] {sum(1 for c in p1_cells if c['dataset'] == ds)} "
              f"cells done")
    if run_scf:
        keys6 = FULL_KEYS["K562"] + [SCF_KEY]
        common = common_perts({k: loaded[k][4] for k in keys6})
        P = {k: align(probes[k], loaded[k][4], common) for k in keys6}
        E = {k: align(err[k], loaded[k][4], common) for k in keys6}
        p1_cells += run_p1_block("scf_k562_diag", SCF_FAMILY_KEYS, P, E,
                                 b_boot, rng_master, diagnostic=True)
        print(f"[P1][scf_k562_diag] cells done")
    p1_prim = [c for c in p1_cells
               if not c["diagnostic"] and c["family"] == "embedding"]
    res["p1"] = {"criterion": P1_CRIT, "cells": p1_cells,
                 "verdict": p1_verdict(p1_prim)}
    print(f"[P1] verdict: {res['p1']['verdict']['tier']}")

    # ---- P2 ----
    p2_blocks = []
    for ds in ds_list:
        keys = FULL_KEYS[ds]
        common = common_perts({k: loaded[k][4] for k in keys})
        P = {k: align(probes[k], loaded[k][4], common) for k in keys}
        R = {k: align(err_rank[k], loaded[k][4], common) for k in keys}
        p2_blocks.append(run_p2_block(ds, FAMILY_KEYS[ds], keys, P, R,
                                      b_boot, rng_master, diagnostic=False))
        print(f"[P2][{ds}] AUROC={p2_blocks[-1]['auroc']:.3f}")
    if run_scf:
        keys6 = FULL_KEYS["K562"] + [SCF_KEY]
        common = common_perts({k: loaded[k][4] for k in keys6})
        P = {k: align(probes[k], loaded[k][4], common) for k in keys6}
        R = {k: align(err_rank[k], loaded[k][4], common) for k in keys6}
        p2_blocks.append(run_p2_block("scf_k562_diag", SCF_FAMILY_KEYS, keys6,
                                      P, R, b_boot, rng_master,
                                      diagnostic=True))
    res["p2"] = {"criterion": P2_CRIT, "blocks": p2_blocks,
                 "verdict": p2_verdict(p2_blocks)}
    print(f"[P2] verdict: {res['p2']['verdict']['tier']}")

    # ---- P3 ----
    rows, folds, scf_note, p3_verdict = run_p3()
    res["p3"] = {"criterion": P3_CRIT, "registry_rows": rows, "folds": folds,
                 "scf_descriptive": scf_note, "verdict": p3_verdict}
    print(f"[P3] verdict: {p3_verdict['tier']}")

    res["timing_total_s"] = round(time.perf_counter() - t_all, 2)

    out_dir = SMOKE_DIR if smoke else OUT_DIR
    os.makedirs(out_dir, exist_ok=True)
    stem = "family_consequence.smoke" if smoke else "family_consequence"
    out_json = os.path.join(out_dir, stem + ".json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("[SAVED]", out_json)
    make_figure(res, os.path.join(out_dir, stem + ".png"))
    if not smoke:
        write_verdict_md(res, os.path.join(OUT_DIR,
                                           "family_consequence_verdict.md"))
    print(f"[DONE] total {res['timing_total_s']}s")


# ---------------------------------------------------------------- 图
def make_figure(res, path):
    fig, axes = plt.subplots(1, 3, figsize=(16.5, 5.4))

    # P1: 每 (D,F) 的 Δ 均值 + CI（主判定格聚合）
    ax = axes[0]
    cells = [c for c in res["p1"]["cells"] if not c["structural"]]
    groups = {}
    for c in cells:
        g = f"{c['dataset']}|{c['family']}"
        groups.setdefault(g, []).append(c["mean_delta"])
    order = [g for g in groups if "scf" not in g] + \
            [g for g in groups if "scf" in g]
    means = [float(np.mean(groups[g])) for g in order]
    # 组内 CI：用格 CI 的最保守包络（min lo, max hi）作展示带
    los = [min(c["ci_lo"] for c in cells if f"{c['dataset']}|{c['family']}"
                == g) for g in order]
    his = [max(c["ci_hi"] for c in cells if f"{c['dataset']}|{c['family']}"
                == g) for g in order]
    yerr = [np.array(means) - np.array(los),
            np.array(his) - np.array(means)]
    colors = ["#4878cf" if "embedding" in g else "#d65f5f" for g in order]
    ax.bar(range(len(order)), means, color=colors, alpha=0.85,
           yerr=yerr, capsize=4,
           error_kw=dict(lw=1.2, ecolor="#333333"))
    ax.axhline(0.0, color="k", lw=1)
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([g.replace("_k562_diag", "(diag)")
                        .replace("|", "\n") for g in order],
                       fontsize=8)
    n_struct = sum(1 for c in res["p1"]["cells"] if c["structural"])
    ax.set_title(f"P1 savings principle: utility loss of single "
                 f"representative vs full-family gate\n(mean delta over "
                 f"(T,R) cells + cell-CI envelope; {n_struct} structural "
                 f"cells delta=0 not shown)", fontsize=9)
    ax.set_ylabel("delta utility (cost_single - cost_full), >0 = loss")
    v = res["p1"]["verdict"]
    ax.text(0.02, 0.98, f"verdict: {v['tier']}\n"
            f"{v['n_ci_contains_zero']}/{v['n_primary_cells']} cells CI~0",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(fc="white", ec="#888888", alpha=0.9))

    # P2: 分歧度 AUROC
    ax = axes[1]
    blocks = res["p2"]["blocks"]
    tags = [b["dataset"] for b in blocks]
    aucs = [b["auroc"] for b in blocks]
    cols = ["#4878cf" if not b["diagnostic"] else "#999999" for b in blocks]
    ax.bar(range(len(tags)), aucs, color=cols, alpha=0.85)
    ax.axhline(0.5, color="k", ls="--", lw=1, label="chance")
    ax.axhline(0.6, color="#d65f5f", ls="-", lw=1.4,
               label="locked threshold 0.6")
    for i, b in enumerate(blocks):
        ax.text(i, b["auroc"] + 0.008, f"{b['auroc']:.3f}\n"
                f"[{b['ci_lo']:.3f},{b['ci_hi']:.3f}]",
                ha="center", fontsize=7.5)
    ax.set_xticks(range(len(tags)))
    ax.set_xticklabels([t.replace("_diag", "\n(diag)") for t in tags],
                       fontsize=8)
    ax.set_ylim(0.4, 0.75)
    ax.set_title("P2 heterogeneity principle: family-disagreement "
                 "signal\nAUROC(d, family-mean-error worst half)",
                 fontsize=9)
    ax.legend(fontsize=8, loc="lower right")
    v = res["p2"]["verdict"]
    ax.text(0.02, 0.02, f"verdict: {v['tier']}", transform=ax.transAxes,
            fontsize=8, bbox=dict(fc="white", ec="#888888", alpha=0.9))

    # P3: LODO regret
    ax = axes[2]
    folds = res["p3"]["folds"]
    x = np.arange(len(folds))
    w = 0.25
    for j, (key, lab, col) in enumerate(
            [("R1", "R1 family-aware", "#4878cf"),
             ("R0", "R0 family-blind", "#d65f5f"),
             ("random_reference_regret", "random reference", "#999999")]):
        vals = [f[key]["regret"] if isinstance(f[key], dict) else f[key]
                for f in folds]
        ax.bar(x + (j - 1) * w, vals, width=w, label=lab, color=col,
               alpha=0.85)
    for i, f in enumerate(folds):
        ax.text(i, -0.012,
                f"train {f['train_unit']} -> test {f['test_unit']}\n"
                f"best={f['test_best_engine']}\n"
                f"R1 {'hit' if f['R1']['top1_hit'] else 'miss'}",
                ha="center", fontsize=7)
    ax.set_xticks(x)
    ax.set_xticklabels([f"fold {i+1}" for i in x], fontsize=8)
    ax.set_ylim(bottom=-0.03)
    ax.axhline(0.0, color="k", lw=1)
    ax.set_title("P3 selection principle: LODO regret\n(best - predicted "
                 "in-situ AUROC; lower is better)", fontsize=9)
    ax.set_ylabel("regret (AUROC points)")
    ax.legend(fontsize=8)
    v = res["p3"]["verdict"]
    ax.text(0.02, 0.98, f"verdict: {v['tier']}\n(prototype, low power)",
            transform=ax.transAxes, va="top", fontsize=8,
            bbox=dict(fc="white", ec="#888888", alpha=0.9))

    fig.suptitle("Family-consequence analysis (on-file replay, zero "
                 f"training; bootstrap B={res['spec']['full_run_spec']['b_boot']}, "
                 f"smoke={res['smoke']})", fontsize=11)
    fig.tight_layout(rect=[0, 0, 0.99, 0.94])
    fig.savefig(path, dpi=300)
    print("[SAVED]", path)


# ---------------------------------------------------------------- verdict
def fmt(v, nd=3):
    return "n/a" if v is None else f"{v:.{nd}f}"


def write_verdict_md(res, path):
    p1, p2, p3 = res["p1"], res["p2"], res["p3"]
    L = ["# Family-Consequence verdict（step24，prereg_family_consequence.md）",
         ""]
    L.append("三条判定互相独立；全部 on-file 回放，零基础模型训练。")
    L.append("")

    # ---- P1 ----
    v1 = p1["verdict"]
    L.append("## P1 省钱原则（家族内全引擎集成 vs 单代表门的效用差）")
    L.append("")
    L.append(f"判定：主判定单元 {v1['n_primary_cells']} 个（K562/RPE1 × "
             f"embedding 家族 × (T,R)），Δ 的 95% CI 含 0 的 "
             f"{v1['n_ci_contains_zero']} 个 -> **{v1['tier']}**。")
    L.append("")
    L.append("| dataset | family | target | rep | mean Δ | 95% CI | 含0 |")
    L.append("|---|---|---|---|---|---|---|")
    for c in p1["cells"]:
        mark = "结构性" if c["structural"] else (
            "是" if c["ci_contains_zero"] else "**否**")
        tag = "（诊断）" if c["diagnostic"] else ""
        L.append(f"| {c['dataset']}{tag} | {c['family']} | "
                 f"{short(c['target'])} | {short(c['rep'])} | "
                 f"{fmt(c['mean_delta'])} | "
                 f"[{fmt(c['ci_lo'])},{fmt(c['ci_hi'])}] | {mark} |")
    L.append("")
    uf = [np.mean(list(c["util_fam_full_vs_random"].values()))
          for c in p1["cells"] if not c["structural"]]
    ur = [np.mean(list(c["util_single_vs_random"].values()))
          for c in p1["cells"] if not c["structural"]]
    L.append(f"- 两臂绝对效用（对 random 的成本节省，4 档预算均值，非结构性格）："
             f"fam_full mean={fmt(float(np.mean(uf)))}，single mean="
             f"{fmt(float(np.mean(ur)))}。P1 判据只依赖 Δ 的 CI；"
             f"绝对效用水平用于解读'省钱'含义（两臂皆近 0 时，结论读作"
             f"'加家族成员不买来门价值'）。")
    L.append("")

    # ---- P2 ----
    v2 = p2["verdict"]
    L.append("## P2 异质原则（家族分歧度作为第三信号）")
    L.append("")
    L.append(f"判定：AUROC(d, 家族平均误差最差一半) 点估计 — "
             + "，".join(f"{d}={fmt(a)}" for d, a in v2["auroc"].items())
             + f"；判据双数据集 ≥ 0.6 -> **{v2['tier']}**。")
    L.append("")
    L.append("| dataset | n | AUROC | 95% CI | sig>0.5 | sig<0.5 | "
             "Spearman(d, r̄) | 变体 expr/all | 变体 emb/all |")
    L.append("|---|---|---|---|---|---|---|---|---|")
    for b in p2["blocks"]:
        tag = "（诊断）" if b["diagnostic"] else ""
        L.append(f"| {b['dataset']}{tag} | {b['n']} | {fmt(b['auroc'])} | "
                 f"[{fmt(b['ci_lo'])},{fmt(b['ci_hi'])}] | "
                 f"{'是' if b['above_0_5'] else '否'} | "
                 f"{'是' if b['below_0_5'] else '否'} | "
                 f"{fmt(b['spearman_d_vs_rbar'])} | "
                 f"{fmt(b['variant_expr_vs_all_auroc'])} | "
                 f"{fmt(b['variant_emb_vs_all_auroc'])} |")
    L.append("")

    # ---- P3 ----
    v3 = p3["verdict"]
    L.append("## P3 选择原则（家族代表性的元学习雏形，LODO 回放）")
    L.append("")
    L.append(f"判定：R1 top-1 命中 {v3['r1_top1_hits']}/{v3['n_folds']} 折；"
             f"R1 regret < R0 regret（双折）："
             f"{'是' if v3['r1_regret_lt_r0_both_folds'] else '否'} -> "
             f"**{v3['tier']}**。（明示：2 个全引擎折 + 11 行 on-file 档案，"
             f"雏形级功效。）")
    L.append("")
    L.append("| train -> test | 真实最优 | R1 预测 | R1 hit | R1 regret | "
             "R0 预测 | R0 regret | random regret |")
    L.append("|---|---|---|---|---|---|---|---|")
    for f in p3["folds"]:
        L.append(f"| {f['train_unit']} -> {f['test_unit']} | "
                 f"{f['test_best_engine']} ({fmt(f['test_best_auroc'],4)}) | "
                 f"{f['R1']['predicted_engine']} | "
                 f"{'Y' if f['R1']['top1_hit'] else 'N'} | "
                 f"{fmt(f['R1']['regret'],4)} | "
                 f"{f['R0']['predicted_engine']} | "
                 f"{fmt(f['R0']['regret'],4)} | "
                 f"{fmt(f['random_reference_regret'],4)} |")
    L.append("")
    scf = p3["scf_descriptive"]
    L.append(f"- scf 描述性：scf_k562(diag) 就地 AUROC={fmt(scf['scf_auroc'],4)}；"
             f"K562 训练下 R1 选 {scf['R1_train_K562_picks']}、R0 选 "
             f"{scf['R0_train_K562_picks']}；RPE1 训练下 R1 选 "
             f"{scf['R1_train_RPE1_picks']}、R0 选 "
             f"{scf['R0_train_RPE1_picks']}。")
    L.append("")
    L.append("## 档案表（P3 输入，全 on-file）")
    L.append("")
    L.append("| key | engine | family | unit | n | balance | shared R² | "
             "own in-situ AUROC |")
    L.append("|---|---|---|---|---|---|---|---|")
    for r in p3["registry_rows"]:
        L.append(f"| {r['key']} | {r['engine']} | {r['family']} | "
                 f"{r['unit']} | {r['n']} | {fmt(r['label_balance'])} | "
                 f"{fmt(r['shared_r2'])} | {fmt(r['own_in_situ_auroc'],4)} |")
    L.append("")

    # ---- 汇总 ----
    L.append("## 汇总（三条独立）")
    L.append("")
    L.append(f"- P1: {v1['tier']}")
    L.append(f"- P2: {v2['tier']}")
    L.append(f"- P3: {v3['tier']}")
    L.append("")
    L.append("一句话原则（论文可用句按判定裁剪）：")
    L.append("- P1 过线：'同家族引擎的门错误冗余到可以只部署一台代表——"
             "对任一代表选择、任一预算档，单引擎门相对家族全集成门的"
             "效用损失 CI 均含 0。'")
    L.append("- P2 过线：'两家族门分数强烈分歧的扰动，家族平均误差确实"
             "更高——分歧度本身是可与两家族信号并列的第三信号。'")
    L.append("- P3 过线：'预算只容一台引擎时，家族归属 + 数据集属性可"
             "预测就地 AUROC 最高的引擎（LODO 回放）。'")
    L.append("- 任一未过线：对应句子不使用；阴性数字全部保留如上。")
    L.append("")
    L.append("与既有档案关系：G-B 定性 kappa -> G-B2 跨家族残差互补阴性 -> "
             "本步回答'家族结构能用来做什么'。三条判定独立，任一不过只禁"
             "对应原则句，不动 R2 结论（产品 = 符号层）。")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print("[SAVED]", path)


if __name__ == "__main__":
    main()
