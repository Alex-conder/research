# step22_task_utility_gc.py -- G-C task-level utility: gate-guided experimental
# budget allocation (VCC narrative: when may a virtual-cell prediction replace
# one real Perturb-seq experiment?)
#
# 预注册：prereg_expansion_gc.md（判据先于计算锁定）。纯 on-file 回放：
# 不训练任何基础模型；唯一拟合对象是家族既定规格的浅层 logistic 探针
# （step6 规格 OOF），与 step21 分析 A 同构造同数据。
#
# 回放：留出扰动集合上预算 B 个采信名额。总成本
#   C_arm(B) = sum_{采信} err + c*(n-B),  c = median(err)（臂间比较与 c 无关）
# 臂：门引导采信（probe_full/top20/top50, symbol, conformal_width=sigma 最窄）
#   vs 基线 random（精确期望）与 magnitude（delta_max 最高）。
# 统计：扰动级配对 bootstrap（固定分数不重拟合）B=1000, seed=20261030。
#
# 产物：outputs_rebuild/task_utility.{json,png} + task_utility_verdict.md
#       smoke 模式 -> outputs_rebuild/smoke/task_utility.smoke.*
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step22_task_utility_gc.py --smoke
import argparse
import json
import os
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

for _s in ("stdout", "stderr"):  # GBK 控制台纪律
    getattr(sys, _s).reconfigure(encoding="utf-8", errors="replace")
os.environ["OMP_NUM_THREADS"] = "1"
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

REPO = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(REPO, "outputs_rebuild")
SMOKE_DIR = os.path.join(OUT_DIR, "smoke")

BOOT_MASTER_SEED = 20261030
B_BOOT_FULL = 1000
B_BOOT_SMOKE = 200
OOF_SEED = 0
BS = [5, 10, 20, 50]
METRICS = ["full", "top20", "top50"]
GATES = ["probe_full", "probe_top20", "probe_top50", "symbol", "conformal_width"]
BASELINES = ["random", "magnitude"]
MIN_POS_DATASETS = 3   # >=3/5 注册数据集（锁定判据）
N_REGISTERED_DS = 5

from step21_ad_baseline import (  # noqa: E402
    SIGN_SOURCES, ad_scores, align, auroc_fast, error_metrics, fit_signs,
    load_bundle, load_product, oof_probe_scores, symbol_scores, y_from_e,
)

COMBOS = [
    {"key": "gears_k562_s0", "ds": "replogle_k562_essential", "out": "outputs",
     "seed": 0},
    {"key": "gears_k562_s1", "ds": "replogle_k562_essential", "out": "outputs",
     "seed": 1},
    {"key": "gears_k562_s2", "ds": "replogle_k562_essential", "out": "outputs",
     "seed": 2},
    {"key": "gears_rpe1_s0", "ds": "replogle_rpe1_essential",
     "out": "outputs_rpe1", "seed": 0},
    {"key": "gears_norman_s4", "ds": "norman",
     "out": "outputs_multidataset/norman", "seed": 4},
    {"key": "gears_adamson_s3", "ds": "adamson",
     "out": "outputs_multidataset/adamson", "seed": 3},
    {"key": "gears_adamson_s5", "ds": "adamson",
     "out": "outputs_multidataset/adamson", "seed": 5},
    {"key": "gears_dixit_s4", "ds": "dixit",
     "out": "outputs_multidataset/dixit", "seed": 4},
    {"key": "gears_dixit_s5", "ds": "dixit",
     "out": "outputs_multidataset/dixit", "seed": 5},
]

DS_LABEL = {"replogle_k562_essential": "K562 (GEARS, 3 seeds)",
            "replogle_rpe1_essential": "RPE1 (GEARS)",
            "norman": "Norman (GEARS s4)",
            "adamson": "Adamson (GEARS, 2 seeds)",
            "dixit": "Dixit (GEARS)"}


# ------------------------------------------------------------- 成本核
def topk_sum(values, scores, B):
    """按 scores 降序取 top-B（B 有效值 = min(B, n)）的 values 和。"""
    n = len(values)
    k = min(B, n)
    if k >= n:
        return float(values.sum())
    idx = np.argpartition(scores, -k)[-k:]
    return float(values[idx].sum())


def boot_diffs(err, s_gate, s_mag, B, rng, b_boot):
    """配对 bootstrap：cost_baseline - cost_gate 的 95% CI（逐档）。

    random 臂取精确期望 B_eff/n * sum(err[idx])；magnitude 臂按
    s_mag 在重采样集内重排序。退化档（B_eff = n）两臂恒等，差恰为 0。
    """
    n = len(err)
    k = min(B, n)
    if k >= n:
        z = {"diff_random": {"mean": 0.0, "ci_lo": 0.0, "ci_hi": 0.0,
                             "ci_excludes_zero": False, "favors": "none",
                             "degenerate": True},
             "diff_magnitude": {"mean": 0.0, "ci_lo": 0.0, "ci_hi": 0.0,
                                "ci_excludes_zero": False, "favors": "none",
                                "degenerate": True}}
        return z, True
    dr = np.empty(b_boot)
    dm = np.empty(b_boot)
    for b in range(b_boot):
        idx = rng.integers(0, n, n)
        e = err[idx]
        c_gate = topk_sum(e, s_gate[idx], B)
        dr[b] = k / n * e.sum() - c_gate
        dm[b] = topk_sum(e, s_mag[idx], B) - c_gate

    def _summ(d):
        lo, hi = np.percentile(d, [2.5, 97.5])
        return {"mean": float(d.mean()), "ci_lo": float(lo),
                "ci_hi": float(hi),
                "ci_excludes_zero": bool(lo > 0 or hi < 0),
                "favors": ("gate" if lo > 0 else "baseline" if hi < 0
                           else "undetermined"),
                "degenerate": False}
    return {"diff_random": _summ(dr), "diff_magnitude": _summ(dm)}, False


# ------------------------------------------------------------- 组合级分析
def run_combo(c, rng, b_boot, signs_cache):
    ds = c["ds"]
    t0 = time.time()
    bundle = load_bundle(ds)
    prod = load_product(ds, c["out"], c["seed"], bundle)
    perts = prod["perts"]
    n = len(perts)
    pred_delta = prod["pred_delta"]

    if ds not in signs_cache:
        srcs = SIGN_SOURCES[ds]
        Xs, ys = [], []
        for s_ds, s_out, s_seed in srcs:
            sb = load_bundle(s_ds)
            sp = load_product(s_ds, s_out, s_seed, sb)
            Xs.append(sp["X"])
            ys.append(sp["y_std"])
        signs_cache[ds] = fit_signs(np.vstack(Xs), np.concatenate(ys))
    s_sym = symbol_scores(prod["X"], signs_cache[ds])

    sigma = np.linalg.norm(pred_delta, axis=1) / np.sqrt(pred_delta.shape[1])

    errs = error_metrics(pred_delta, prod["true_delta"])
    y_lab = {m: y_from_e(errs[m]) for m in METRICS}
    s_probe = {m: oof_probe_scores(prod["X"], y_lab[m]) for m in METRICS}

    # 幅度基线（预测侧 delta_max，step6 特征同定义）
    s_mag = np.abs(pred_delta).max(axis=1)

    # 诊断门变体（不进判定）：AD-mean（step21 构造）
    ad_ref = align(bundle["var_genes"], bundle["train_delta_mean"],
                   prod["pred_genes"])
    s_ad = ad_scores(pred_delta, ad_ref)

    c_measure = float(np.median(errs["full"]))

    # 采信分数（越高越可靠 = 越该采信；探针/符号门为失败似然分数，取负；
    # conformal 采信 sigma 最窄；幅度基线为 delta_max 最高）
    trust_scores = {
        "probe_full": None if s_probe["full"] is None else -s_probe["full"],
        "probe_top20": None if s_probe["top20"] is None else -s_probe["top20"],
        "probe_top50": None if s_probe["top50"] is None else -s_probe["top50"],
        "symbol": -s_sym,
        "conformal_width": -sigma,
    }

    # AUROC 诊断分数（越高 = 更可能失败；conformal 用 +sigma）
    diag_scores = {
        "probe_full": s_probe["full"],
        "probe_top20": s_probe["top20"],
        "probe_top50": s_probe["top50"],
        "symbol": s_sym,
        "conformal_width": sigma,
        "ad_diagnostic": s_ad,
    }
    auroc_diag = {}
    for m in METRICS:
        if np.bincount(y_lab[m]).min() < 2:
            auroc_diag[m] = {"evaluable": False}
            continue
        row = {g: (auroc_fast(y_lab[m], s) if s is not None else None)
               for g, s in diag_scores.items()}
        auroc_diag[m] = {"evaluable": True, "auroc": row}

    evaluable_B = [B for B in BS if min(B, n) < n]
    out = {"dataset": ds, "seed": c["seed"], "n": int(n),
           "c_measure_median_err": c_measure,
           "measure_all_cost": c_measure * n,
           "evaluable_B": evaluable_B,
           "degenerate_all_B": len(evaluable_B) == 0,
           "auroc_diagnostic": auroc_diag, "per_B": {}, "per_metric_costs": {}}

    for B in BS:
        k = min(B, n)
        per_gate = {}
        for g in GATES:
            s = trust_scores[g]
            if s is None:
                per_gate[g] = {"non_evaluable": True}
                continue
            boots, deg = boot_diffs(errs["full"], s, s_mag, B, rng, b_boot)
            per_gate[g] = {
                "non_evaluable": False,
                "cost_gate": topk_sum(errs["full"], s, B),
                "cost_random": k / n * errs["full"].sum(),
                "cost_magnitude": topk_sum(errs["full"], s_mag, B),
                "savings_vs_measure_all": c_measure * k
                - topk_sum(errs["full"], s, B),
                "bootstrap": boots, "degenerate": deg}
        out["per_B"][str(B)] = {"B_eff": int(k), "gates": per_gate}

    # 敏感性：top20/top50 口径下各臂成本（点估计，不进判定）
    for m in ("top20", "top50"):
        rows = {}
        for g in GATES:
            s = trust_scores[g]
            rows[g] = {str(B): (topk_sum(errs[m], s, B) if s is not None
                                else None) for B in BS}
        rows["random"] = {str(B): min(B, n) / n * errs[m].sum() for B in BS}
        rows["magnitude"] = {str(B): topk_sum(errs[m], s_mag, B) for B in BS}
        out["per_metric_costs"][m] = rows

    print(f"[combo][{c['key']}] n={n} evaluable_B={evaluable_B} "
          f"({time.time()-t0:.0f}s)", flush=True)
    return out


# ------------------------------------------------------------- 判定（锁定）
def adjudicate(results):
    datasets = {}
    for key, r in results.items():
        datasets.setdefault(r["dataset"], []).append(key)
    verdict_gates = {}
    for g in GATES:
        ds_pos = {}
        for ds, keys in datasets.items():
            pos = True
            for k_ in keys:
                cell = results[k_]
                if cell["degenerate_all_B"]:
                    pos = False
                    continue
                for B in cell["evaluable_B"]:
                    pg = cell["per_B"][str(B)]["gates"].get(g, {})
                    if pg.get("non_evaluable"):
                        pos = False
                        break
                    for bl in BASELINES:
                        d = pg["bootstrap"]["diff_" + bl]
                        if d["degenerate"] or not d["ci_excludes_zero"] \
                                or d["favors"] != "gate":
                            pos = False
                            break
                    if not pos:
                        break
                if not pos:
                    break
            ds_pos[ds] = bool(pos)
        n_pos = sum(ds_pos.values())
        verdict_gates[g] = {
            "dataset_positive": ds_pos, "n_positive": int(n_pos),
            "confirmed": bool(n_pos >= MIN_POS_DATASETS),
            "rule": f"confirmed iff gate beats both baselines (95% paired "
                    f"bootstrap CI excludes 0 favoring gate at EVERY "
                    f"evaluable B) on every seed combo of the dataset, on "
                    f">={MIN_POS_DATASETS}/{N_REGISTERED_DS} registered "
                    f"datasets"}
    primary = verdict_gates["probe_full"]
    tier = "CONFIRMED" if primary["confirmed"] else "NOT_CONFIRMED"
    return {"tier": tier, "primary_gate": "probe_full",
            "per_gate": verdict_gates,
            "rule": "primary = probe_full; per-gate verdicts also reported; "
                    "negative result reported as-is (no utility is a "
                    "result); R2 product unchanged (symbol layer)"}


# ------------------------------------------------------------- 图
def make_figure(results, verdict, path, b_boot):
    order = ["replogle_k562_essential", "replogle_rpe1_essential", "norman",
             "adamson", "dixit"]
    by_ds = {}
    for key, r in results.items():
        by_ds.setdefault(r["dataset"], []).append((key, r))
    colors = {"probe_full": "#1f77b4", "probe_top20": "#17becf",
              "probe_top50": "#9edae5", "symbol": "#2ca02c",
              "conformal_width": "#9467bd", "random": "#7f7f7f",
              "magnitude": "#d62728"}
    fig, axes = plt.subplots(2, 3, figsize=(13.5, 7.6))
    axes = axes.ravel()
    for ax_i, ds in enumerate(order):
        ax = axes[ax_i]
        entries = by_ds.get(ds, [])
        if not entries:
            ax.axis("off")
            continue
        first = entries[0][1]
        if first["degenerate_all_B"]:
            ax.text(0.5, 0.5, "degenerate\n(n = %d < min B)" % first["n"],
                    ha="center", va="center", transform=ax.transAxes,
                    fontsize=10, color="gray")
            ax.set_title(DS_LABEL[ds], fontsize=10)
            ax.set_xticks(BS)
            continue
        for name in GATES + BASELINES:
            ys = []
            for B in BS:
                vals = []
                for _, r in entries:
                    cell = r["per_B"][str(B)]["gates"]
                    if name in BASELINES:
                        # 基线成本与门无关，存于每个门格内，取第一个可评测格
                        base_cell = next(
                            (cell[g] for g in GATES
                             if cell.get(g)
                             and not cell[g].get("non_evaluable")), None)
                        if base_cell is None:
                            vals.append(np.nan)
                        else:
                            vals.append(base_cell["cost_" + name])
                        continue
                    pg = cell.get(name)
                    if not pg or pg.get("non_evaluable"):
                        vals.append(np.nan)
                    else:
                        vals.append(pg["cost_gate"])
                ys.append(np.nanmean(vals))
            ax.plot(BS, ys, marker="o", ms=3, lw=1.4,
                    color=colors[name], label=name)
        ma = np.mean([r["measure_all_cost"] for _, r in entries])
        ax.axhline(ma, color="k", ls="--", lw=1.0)
        ax.text(BS[0], ma, "measure all", fontsize=6.5, va="bottom")
        ax.set_title(DS_LABEL[ds], fontsize=10)
        ax.set_xticks(BS)
        ax.grid(alpha=0.25)
    ax = axes[5]
    ax.axis("off")
    lines = ["Verdict (locked, prereg_expansion_gc.md)", ""]
    for g in GATES:
        v = verdict["per_gate"][g]
        tag = " (primary)" if g == "probe_full" else ""
        lines.append(f"{g}{tag}: {v['n_positive']}/"
                     f"{N_REGISTERED_DS} datasets -> "
                     f"{'CONFIRMED' if v['confirmed'] else 'not confirmed'}")
        for ds, p in v["dataset_positive"].items():
            lines.append(f"    {ds}: {'+' if p else '-'}")
    lines += ["", f"paired bootstrap B={b_boot}, seed={BOOT_MASTER_SEED}",
              "cost = sum trusted err + c*(n-B), c = median(err)"]
    ax.text(0.02, 0.98, "\n".join(lines), va="top", fontsize=8,
            transform=ax.transAxes, family="monospace")
    axes[0].set_ylabel("total cost (full-MAE units)", fontsize=9)
    axes[3].set_ylabel("total cost (full-MAE units)", fontsize=9)
    for ax in (axes[3], axes[4]):
        ax.set_xlabel("budget B (trusted perturbations)", fontsize=9)
    handles, labels = axes[0].get_legend_handles_labels()
    if handles:
        axes[4].legend(handles, labels, loc="upper left", fontsize=8)
    fig.suptitle("G-C task-level utility: gate-guided trust vs baselines "
                 "(cost of prediction errors + assay cost of the rest)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0.02, 1, 0.95])
    fig.savefig(path, dpi=300)
    print(f"[SAVED] {path}", flush=True)


# ------------------------------------------------------------- verdict md
def fmt_ci(d):
    if d.get("degenerate"):
        return "degenerate"
    return f"{d['mean']:+.3f} [{d['ci_lo']:+.3f},{d['ci_hi']:+.3f}] {d['favors']}"


def write_verdict_md(results, verdict, path, b_boot, smoke):
    L = ["# G-C task-level utility (step22, prereg_expansion_gc.md) - verdict",
         ""]
    L.append(f"- generated: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    L.append(f"- smoke: {smoke}")
    L.append("- replay: GEARS single-gene holdout products (step21 analysis A "
             "identical data/protocol); zero training, CPU, on-file")
    L.append("- cost: C(B) = sum of trusted prediction errors (full "
             "transcriptome delta MAE) + c*(n-B); c = per-combo median err "
             "(arm comparisons at fixed B are c-free)")
    L.append(f"- statistics: perturbation-level paired bootstrap, B={b_boot}, "
             f"master seed={BOOT_MASTER_SEED}; no model refits inside "
             f"bootstrap (score vectors fixed, preregistered)")
    L.append("- locked rule: combo positive vs a baseline iff the paired 95% "
             "CI of (cost_baseline - cost_gate) excludes 0 favoring the gate "
             "at EVERY evaluable B; dataset positive iff ALL its seed combos "
             "are positive vs BOTH baselines; confirmed iff >=3/5 registered "
             "datasets (dixit degenerate by construction -> auto-negative)")
    L.append("")
    for g in GATES:
        v = verdict["per_gate"][g]
        tag = " PRIMARY" if g == "probe_full" else ""
        L.append(f"## Gate arm: {g}{tag} -> "
                 f"{'CONFIRMED' if v['confirmed'] else 'NOT_CONFIRMED'} "
                 f"({v['n_positive']}/{N_REGISTERED_DS})")
        L.append("")
        L.append("| dataset | B | d(random-gate) [95% CI] | "
                 "d(magnitude-gate) [95% CI] |")
        L.append("|---|---|---|---|")
        for key, r in results.items():
            for B in BS:
                pg = r["per_B"][str(B)]["gates"].get(g, {})
                if not pg or pg.get("non_evaluable"):
                    L.append(f"| {key} | {B} | n/e | n/e |")
                    continue
                br = pg["bootstrap"]["diff_random"]
                bm = pg["bootstrap"]["diff_magnitude"]
                ds_tag = key if B == BS[0] else ""
                L.append(f"| {ds_tag} | {B} | {fmt_ci(br)} | {fmt_ci(bm)} |")
        L.append("")
        L.append("dataset_positive: "
                 + ", ".join(f"{d}={'Y' if p else 'N'}"
                             for d, p in v["dataset_positive"].items()))
        L.append("")
    L.append("## Verdict (locked)")
    L.append(f"- tier: **{verdict['tier']}** (primary gate = "
             f"{verdict['primary_gate']})")
    L.append(f"- rule: {verdict['rule']}")
    L.append("")
    L.append("## Where this goes in the paper (one paragraph)")
    L.append("Gating paper main text, 'task-level utility' panel: the gate is "
             "no longer only an error-scoring tool but a decision tool -- it "
             "decides which virtual-cell predictions may stand in for a real "
             "Perturb-seq measurement under an experimental budget (VCC "
             "'replace one real experiment' narrative). Per-dataset utility "
             "tables above; positive verdict licenses the claim 'the gate has "
             "positive utility under experimental budget constraints'; a "
             "negative verdict is reported as-is (a gate that scores errors "
             "but cannot save experiments), which would itself be a "
             "tier-level conclusion about reliability signals in this field.")
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"[SAVED] {path}", flush=True)


# ------------------------------------------------------------- 主流程
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    smoke = args.smoke
    b_boot = B_BOOT_SMOKE if smoke else B_BOOT_FULL
    out_dir = SMOKE_DIR if smoke else OUT_DIR
    os.makedirs(out_dir, exist_ok=True)

    combos = [c for c in COMBOS if c["key"] in (
        {"gears_k562_s0", "gears_norman_s4"} if smoke
        else {c2["key"] for c2 in COMBOS})]

    seed_seqs = np.random.SeedSequence(BOOT_MASTER_SEED).spawn(len(combos))
    signs_cache = {}
    results = {}
    t_all = time.time()
    for c, ss in zip(combos, seed_seqs):
        rng = np.random.default_rng(ss)
        res = run_combo(c, rng, b_boot, signs_cache)
        results[c["key"]] = res
    verdict = adjudicate(results)
    print(f"[adjudicate] tier={verdict['tier']} "
          f"{ {g: v['n_positive'] for g, v in verdict['per_gate'].items()} }",
          flush=True)

    res = {"task": "G-C task-level utility: gate-guided experimental budget "
                   "allocation (step22)",
           "prereg": "prereg_expansion_gc.md",
           "smoke": smoke,
           "spec": {"b_boot": b_boot, "boot_master_seed": BOOT_MASTER_SEED,
                    "oof_seed": OOF_SEED, "B_scan": BS,
                    "cost_metric_primary": "full delta MAE",
                    "c_measure": "per-combo median(full MAE)",
                    "gates": GATES, "baselines": BASELINES,
                    "conformal_width": "sigma_i = ||pred_delta_i||_2/sqrt(G) "
                                       "(G-D prereg_conformal_arm.md); trust "
                                       "= smallest sigma",
                    "magnitude_score": "delta_max = max_g |pred_delta| "
                                       "(prediction-side, step6 feature)",
                    "min_positive_datasets": MIN_POS_DATASETS,
                    "n_registered_datasets": N_REGISTERED_DS},
           "combos": results, "verdict": verdict,
           "timing_total_s": round(time.time() - t_all, 1)}
    stem = "task_utility.smoke" if smoke else "task_utility"
    out_json = os.path.join(out_dir, stem + ".json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print(f"[SAVED] {out_json}", flush=True)
    make_figure(results, verdict, os.path.join(out_dir, stem + ".png"),
                b_boot)
    if not smoke:
        write_verdict_md(results, verdict,
                         os.path.join(out_dir, "task_utility_verdict.md"),
                         b_boot, smoke)
    print(f"[DONE] verdict={verdict['tier']} "
          f"total={res['timing_total_s']}s", flush=True)


if __name__ == "__main__":
    main()
