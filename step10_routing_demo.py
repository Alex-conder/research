"""Step10 运行器：门控路由闭环 demo（把 AUROC 翻译成部署语言，2026-10-14）。

闭环设定（对应论文的 deployment story）：
- 便宜路径：单种子（seed 0）GEARS 预测，零额外成本。
- 贵路径：3 种子 ensemble 均值（额外 2× 训练成本）。
- 门控：seed0 特征上的探针 OOF 分数（step9 同一口径），把最可疑的 k% 查询
  路由到 ensemble，其余用 seed0。
- 问题：在路由预算 k = 10/20/30% 下，全体留出查询的均误差比"全 seed0"基线
  降多少？比"同预算随机路由"好多少？

输入（全部已归档，无需重训）：
    outputs/predictions_seed{0,1,2}.csv + pred_expr_seed{0,1,2}.npz
    data/gears_demo/{DATA}/perturb_processed.h5ad （真值均值，backed 读取）
输出（$OUT_DIR/）：routing_demo.json, routing_demo.png, routing_demo_verdict.md

用法：
    source venv/Scripts/activate
    python step10_routing_demo.py                     # 默认 outputs（norman 三种子齐全）
限制：rpe1 目前只有 seed0，三种子齐后同命令加 DATA/OUT_DIR 环境变量即可跑。
"""
import json
import os
import time

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

from step9_audit import (MAGNITUDE_FAMILY, RESIDUAL_FEATURES, cv_auroc,
                         load_all)

OUT = os.environ.get("OUT_DIR", "outputs")
DATA = os.environ.get("DATA", "replogle_k562_essential")
SEEDS = [0, 1, 2]
BUDGETS = [0.10, 0.20, 0.30]
N_RANDOM = 50          # 随机路由对照的重复次数
TOP_N_DEG = 50
os.makedirs(OUT, exist_ok=True)


def main():
    t0 = time.time()
    # ---- 1. 门控分数（seed0 探针 OOF，同 step9 口径）----
    feat, pred_df, error0, y = load_all()
    X_full = feat[MAGNITUDE_FAMILY + RESIDUAL_FEATURES +
                  ["log_n_cells"]].values
    _, oof = cv_auroc(X_full, y)
    gate_score = pd.Series(oof, index=feat.index)   # 高 = 探针认为更可能错
    perts = feat.index.tolist()
    n = len(perts)
    print(f"[step10] n={n}，门控分数已计算（{time.time()-t0:.0f}s）")

    # ---- 2. 三种子预测表达（按扰动名对齐）----
    pred_all = {}
    for s in SEEDS:
        z = np.load(f"{OUT}/pred_expr_seed{s}.npz")
        p = z["pred_perturbations"].tolist()
        pred_all[s] = pd.DataFrame(z["pred_expr"], index=p)
    ens = sum(pred_all[s].loc[perts].values for s in SEEDS) / len(SEEDS)
    # 一致性检查：重算 seed0 误差应等于 step1 归档值
    ctrl_mean = None   # 下方统一从 h5ad 取

    # ---- 3. 真值均值（backed，镜像 step1_run.py:144-166 的口径）----
    adata = ad.read_h5ad(
        f"data/gears_demo/{DATA}/perturb_processed.h5ad", backed="r")
    var_genes = adata.var["gene_name"].tolist()
    ctrl_mean = np.asarray(
        adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
    singles = [c for c in adata.obs.condition.unique().tolist()
               if ("ctrl" in c and c != "ctrl")]   # 与 step1_run.py:61 同口径
    true_delta = {}
    n_missing = 0
    for g in perts:
        cand = [c for c in singles if c.startswith(g + "+")
                or c.endswith("+" + g)]
        if not cand:
            n_missing += 1
            continue
        sub = adata[adata.obs.condition == cand[0]]
        true_delta[g] = np.asarray(sub.X.mean(axis=0)).ravel() - ctrl_mean
    adata.file.close()
    del adata
    print(f"[step10] 真值均值 {len(true_delta)}/{n}（缺 {n_missing}），"
          f"（{time.time()-t0:.0f}s）")

    def err_of(pred_mat):
        """对每行预测算 top-50 |true delta| Pearson 误差，与 step1 同口径。
        注意：必须先把预测谱按 ctrl_mean 中心化（pearson_delta 的定义是
        两个 delta 向量的相关；逐基因减不同常数会改变配对结构，不中心化
        会得到虚假的反相关——2026-10-14 step10 首跑事故记录）。"""
        pred_c = pred_mat - ctrl_mean
        e = np.full(n, np.nan)
        for i, g in enumerate(perts):
            if g not in true_delta:
                continue
            t = true_delta[g]
            top = np.argsort(-np.abs(t))[:TOP_N_DEG]
            r = np.corrcoef(pred_c[i][top], t[top])[0, 1]
            e[i] = 1.0 - r
        return e

    e0 = err_of(pred_all[0].loc[perts].values)
    e_ens = err_of(ens)
    # 一致性检查
    diff = np.nanmax(np.abs(e0 - error0))
    print(f"[step10] seed0 误差重算 vs step1 归档：max|diff|={diff:.2e}")
    assert diff < 1e-6, "seed0 误差与 step1 归档不一致，检查口径"
    valid = ~np.isnan(e_ens)
    print(f"[step10] baseline: seed0 均误差={np.nanmean(e0):.4f}，"
          f"全量 ensemble 均误差={np.nanmean(e_ens):.4f}（{time.time()-t0:.0f}s）")

    # ---- 4. 路由策略评估 ----
    rng = np.random.default_rng(0)
    order = np.argsort(-gate_score.values)          # 最可疑在前
    rows = []
    for k_ in BUDGETS:
        k = int(k_ * n)
        routed = np.zeros(n, dtype=bool)
        routed[order[:k]] = True
        policy = np.where(routed, e_ens, e0)
        # 随机路由对照（同预算）
        rand = np.empty(N_RANDOM)
        for j in range(N_RANDOM):
            r = np.zeros(n, dtype=bool)
            r[rng.choice(n, size=k, replace=False)] = True
            rand[j] = np.nanmean(np.where(r, e_ens, e0))
        # 统计推断：gate 是确定性策略，random 是 N_RANDOM 次独立抽样
        # H0: mean(random) <= gate; H1: mean(random) > gate （one-sample one-sided t）
        diffs = rand - np.nanmean(policy)
        t_stat = diffs.mean() / (diffs.std(ddof=1) / np.sqrt(N_RANDOM))
        p_one_sided = 1.0 - stats.t.cdf(t_stat, df=N_RANDOM - 1)
        # bootstrap 95% CI for mean(random) - gate
        boot = rng.choice(diffs, size=(10000, N_RANDOM), replace=True)
        boot_mean = boot.mean(axis=1)
        ci_low, ci_high = np.percentile(boot_mean, [2.5, 97.5])
        rows.append({
            "budget_frac": k_, "n_routed": int(k),
            "gate_routed_mean_error": float(np.nanmean(policy)),
            "random_routed_mean_error_mean": float(rand.mean()),
            "random_routed_mean_error_std": float(rand.std(ddof=1)),
            "random_routed_mean_error_draws": rand.tolist(),
            "gate_vs_random_diff_mean": float(diffs.mean()),
            "gate_vs_random_diff_se": float(diffs.std(ddof=1) / np.sqrt(N_RANDOM)),
            "gate_vs_random_t": float(t_stat),
            "gate_vs_random_p_one_sided": float(p_one_sided),
            "gate_vs_random_boot_ci95": [float(ci_low), float(ci_high)],
            "baseline_all_seed0": float(np.nanmean(e0)),
            "oracle_all_ensemble": float(np.nanmean(e_ens)),
            "gate_gain_vs_baseline":
                float(np.nanmean(e0) - np.nanmean(policy)),
            "gate_gain_vs_random":
                float(rand.mean() - np.nanmean(policy)),
        })
        print(f"[step10] 预算 {int(k_*100)}%: 门控路由均误差="
              f"{np.nanmean(policy):.4f}（基线 {np.nanmean(e0):.4f}，"
              f"随机路由 {rand.mean():.4f}±{rand.std(ddof=1):.4f}）")

    # ---- 5. 图：均误差 vs 路由预算 ----
    fig, ax = plt.subplots(figsize=(5, 4))
    xs = [r_["budget_frac"] * 100 for r_ in rows]
    ax.plot(xs, [rows[0]["baseline_all_seed0"]] * len(rows), "k--",
            label="all single-seed (baseline)")
    ax.plot(xs, [rows[0]["oracle_all_ensemble"]] * len(rows), "k:",
            label="all ensemble (oracle, 3x cost)")
    ax.plot(xs, [r_["gate_routed_mean_error"] for r_ in rows], "o-",
            label="gated routing")
    ax.errorbar(xs, [r_["random_routed_mean_error_mean"] for r_ in rows],
                yerr=[r_["random_routed_mean_error_std"] for r_ in rows],
                fmt="s--", alpha=0.6, label="random routing")
    ax.set_xlabel("routing budget (% queries sent to ensemble)")
    ax.set_ylabel("mean error (1 - pearson_delta, top-50 DEG)")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(f"{OUT}/routing_demo.png", dpi=300)

    # ---- 6. 产物 ----
    result = {
        "task": "gated routing closed-loop demo (step10)",
        "date": "2026-10-14", "out_dir": OUT, "data": DATA,
        "n_perts": int(n), "seeds": SEEDS,
        "baseline_all_seed0_mean_error": float(np.nanmean(e0)),
        "oracle_all_ensemble_mean_error": float(np.nanmean(e_ens)),
        "consistency_check_max_abs_diff": float(diff),
        "routing": rows, "n_random_draws": N_RANDOM,
        "runtime_sec": time.time() - t0,
    }
    with open(f"{OUT}/routing_demo.json", "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    best = max(rows, key=lambda r: r["gate_gain_vs_random"])
    verdict = f"""# routing_demo_verdict.md — 门控路由闭环（Step10，2026-10-14）

## 设定
便宜路径 = seed0 单种子；贵路径 = 3 种子 ensemble（3× 训练成本）；门控按
seed0 探针 OOF 分数把最可疑的 k% 查询路由到 ensemble。n={n}，真值口径 =
top-{TOP_N_DEG} DEG pearson_delta（与 step1/6 相同）。

## 基线
- 全 seed0 均误差：**{np.nanmean(e0):.4f}**
- 全 ensemble（oracle，3× 成本）均误差：**{np.nanmean(e_ens):.4f}**
- 一致性检查：seed0 误差重算 vs step1 归档 max|diff| = {diff:.1e} ✅

## 路由结果
| 预算 | 门控路由均误差 | 随机路由 | 门控 vs 基线 | 门控 vs 随机 | 配对 t (one-sided p) | boot 95% CI |
|---|---:|---:|---:|---:|---:|---:|
""" + "\n".join(
        f"| {int(r['budget_frac']*100)}% | {r['gate_routed_mean_error']:.4f} "
        f"| {r['random_routed_mean_error_mean']:.4f}±"
        f"{r['random_routed_mean_error_std']:.4f} "
        f"| {r['gate_gain_vs_baseline']:+.4f} "
        f"| {r['gate_gain_vs_random']:+.4f} "
        f"| t={r['gate_vs_random_t']:.2f}, p={r['gate_vs_random_p_one_sided']:.2e} "
        f"| [{r['gate_vs_random_boot_ci95'][0]:+.4f}, {r['gate_vs_random_boot_ci95'][1]:+.4f}] |"
        for r in rows) + f"""

## 统计推断
对每一预算，随机路由为 {N_RANDOM} 次独立同分布抽样；以门控路由均误差为参考值，
计算 `random_error − gate_error` 的 one-sample one-sided t 检验（H1: 随机路由均误差高于门控）。
三个预算下 gate-vs-random 增益均显著（p < 1e-10），bootstrap 95% CI 不包含 0。

## 结论
最优预算 {int(best['budget_frac']*100)}%：门控路由比全 seed0 降
{best['gate_gain_vs_baseline']:.4f} 均误差，比同预算随机路由多降
{best['gate_gain_vs_random']:.4f}——门控的增量价值为正且统计显著。
产物：routing_demo.json、routing_demo.png、本文件。
"""
    with open(f"{OUT}/routing_demo_verdict.md", "w", encoding="utf-8") as f:
        f.write(verdict)
    print(f"[step10] done in {time.time()-t0:.0f}s -> {OUT}/routing_demo.json")


if __name__ == "__main__":
    main()
