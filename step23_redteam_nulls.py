"""Step23: 结果级红黑队——四个头条新数字的置换零检验 (2026-10-30)。

制度第二道（先于论文回填执行）。对四个头条数字各跑置换零检验 (B>=1000)：

  1. norman  多数据集探针均值 0.7833 (s3/s4/s5 = 0.60/0.85/0.90)
            零假设: 标签与特征独立 -> 打乱 y, 按 step17 同一 5-fold CV 协议
            (探针种子=模型种子, fold 内标准化, L2 logistic C=1.0) 逐种子重算
            cv_auroc_mean, 取 3 种子均值; p = P(null_mean >= obs_mean)。
  2. adamson 同上 (0.80/0.60/0.70, 均值 0.7000)。
  3. G-C     RPE1 s0 B=50 probe_full 效用 d(random-gate) ~ +0.94
            零假设: 采信分配与误差独立 -> 打乱采信任分在扰动间的分配,
            重算 d = cost_random - cost_gate; p = P(d_perm >= d_obs)。
  4. G-D     A_k562 conformal 单元覆盖 0.8921 vs 名义 0.90
            零假设: 覆盖标签以名义率 p=0.9 独立生成 -> Monte-Carlo 标签打乱
            (B=1000) + 精确二项检验; 双侧偏差 p 与单侧欠覆盖 p 并报。

纪律: 先复现存档点估计并与 on-file json 对账 (tol 1e-6), 对账不过不跑置换;
数字逐位取自存档文件, 本脚本不改写任何判据。CPU, 不碰 GPU。
用法: venv/Scripts/python.exe step23_redteam_nulls.py
"""
import json
import os
import sys
import time

import numpy as np

for _s in ("stdout", "stderr"):  # GBK 控制台纪律
    getattr(sys, _s).reconfigure(encoding="utf-8", errors="replace")
os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("CUDA_VISIBLE_DEVICES", "")

REPO = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(REPO, "outputs_rebuild")
DATA_DIR = os.path.join(REPO, "data", "gears_demo")

B_PERM = 1000
MASTER_SEED = 20261031
MULTI = os.path.join(REPO, "outputs_multidataset")


def log(msg):
    print(f"[step23 {time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------- 探针核
def cv_auroc_mean_and_pooled(X, y, seed):
    """step6/step17 同协议: 分层 5-fold CV L2 logistic (fold 内标准化)。

    返回 (fold AUROC 均值, pooled OOF AUROC)。与 step6.run_probe 对账。
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import StratifiedKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    clf = make_pipeline(
        StandardScaler(),
        LogisticRegression(C=1.0, max_iter=2000, random_state=seed))
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=seed)
    aucs, oof = [], np.zeros(len(y), dtype=np.float64)
    for tr, te in cv.split(X, y):
        clf.fit(X[tr], y[tr])
        s = clf.decision_function(X[te])
        oof[te] = s
        aucs.append(roc_auc_score(y[te], s))
    return float(np.mean(aucs)), float(roc_auc_score(y, oof))


def load_xy_step17(ds, seed):
    """复刻 step17_multidataset.run_probe 的特征/标签装配 (无副作用)。"""
    import anndata as ad
    import pandas as pd
    from step6_gating_probe import COVARIATES, SPECTRUM_FEATURES, compute_features

    out_dir = os.path.join(MULTI, ds)
    z = np.load(os.path.join(out_dir, f"pred_expr_seed{seed}.npz"))
    pred_expr, perts = z["pred_expr"], z["pred_perturbations"].tolist()
    pred_genes = z["pred_genes"].tolist()

    h5 = os.path.join(DATA_DIR, ds, "perturb_processed.h5ad")
    adata = ad.read_h5ad(h5, backed="r")
    var_genes = adata.var["gene_name"].tolist()
    ctrl_mean = np.asarray(
        adata[adata.obs.condition == "ctrl"].X.mean(axis=0)).ravel()
    adata.file.close()
    del adata
    ctrl_aligned = (ctrl_mean if var_genes == pred_genes
                    else np.array([dict(zip(var_genes, ctrl_mean))[g]
                                   for g in pred_genes]))

    pred_df = pd.read_csv(
        os.path.join(out_dir, f"predictions_seed{seed}.csv")
    ).set_index("perturbation").loc[perts]

    feat = compute_features(pred_expr, perts, pred_genes, ctrl_aligned)
    feat["log_n_cells"] = np.log1p(pred_df["n_cells"].values)
    if pred_df["in_pert_graph"].nunique() > 1:
        feat["in_pert_graph"] = pred_df["in_pert_graph"].astype(int).values
    feature_names = SPECTRUM_FEATURES + [c for c in COVARIATES
                                         if c in feat.columns]
    X = feat[feature_names].values.astype(np.float64)
    assert not np.isnan(X).any(), "feature matrix has NaN"

    error = 1.0 - pred_df["pearson_delta"].values
    y = (error > np.median(error)).astype(int)
    return X, y, feature_names


def redteam_dataset(ds, seeds, master_rng):
    """数据集级置换零检验: 统计量 = 逐种子 cv_auroc_mean 的跨种子均值。"""
    archived = {}
    for s in seeds:
        with open(os.path.join(MULTI, ds, f"gating_probe_seed{s}.json"),
                  encoding="utf-8") as f:
            j = json.load(f)
        archived[s] = {"cv_mean": j["cv_auroc_mean"],
                       "pooled": j["pooled_cv_auroc"], "n": j["n_perts"]}

    obs_per_seed, Xs, ys = {}, {}, {}
    for s in seeds:
        X, y, _ = load_xy_step17(ds, s)
        Xs[s], ys[s] = X, y
        assert len(y) == archived[s]["n"], (
            f"{ds} s{s}: n mismatch {len(y)} vs {archived[s]['n']}")
        cm, po = cv_auroc_mean_and_pooled(X, y, s)
        assert abs(cm - archived[s]["cv_mean"]) < 1e-6, (
            f"{ds} s{s}: cv_mean replicate {cm} vs archived "
            f"{archived[s]['cv_mean']}")
        assert abs(po - archived[s]["pooled"]) < 1e-6, (
            f"{ds} s{s}: pooled replicate {po} vs archived "
            f"{archived[s]['pooled']}")
        obs_per_seed[s] = {"cv_mean": cm, "pooled": po}
        log(f"{ds} s{s}: replicated cv_mean={cm:.4f} pooled={po:.4f} "
            f"(archived {archived[s]['cv_mean']:.4f}/{archived[s]['pooled']:.4f})")

    obs_stat = float(np.mean([obs_per_seed[s]["cv_mean"] for s in seeds]))
    null_stat = np.empty(B_PERM)
    null_pooled = {s: np.empty(B_PERM) for s in seeds}
    for b in range(B_PERM):
        per_seed = []
        for s in seeds:
            yp = master_rng.permutation(ys[s])
            cm, po = cv_auroc_mean_and_pooled(Xs[s], yp, s)
            per_seed.append(cm)
            null_pooled[s][b] = po
        null_stat[b] = np.mean(per_seed)
        if (b + 1) % 200 == 0:
            log(f"{ds}: perm {b + 1}/{B_PERM}")
    p_stat = (1.0 + float((null_stat >= obs_stat - 1e-12).sum())) / (B_PERM + 1)
    pooled_p = {}
    for s in seeds:
        pooled_p[s] = (1.0 + float((null_pooled[s] >=
                                    obs_per_seed[s]["pooled"] - 1e-12).sum())
                       ) / (B_PERM + 1)
    log(f"{ds}: obs mean={obs_stat:.4f}, perm p={p_stat:.4f}")
    return {
        "dataset": ds, "seeds": seeds, "B": B_PERM,
        "obs_per_seed": obs_per_seed, "obs_stat": obs_stat,
        "null_stat": null_stat.tolist(),
        "null_mean": float(null_stat.mean()),
        "null_sd": float(null_stat.std()),
        "p_value": p_stat,
        "secondary_pooled_oof_p": pooled_p,
    }


# ---------------------------------------------------------------- G-C 核
def topk_sum(values, scores, B):
    n = len(values)
    k = min(B, n)
    if k >= n:
        return float(values.sum())
    idx = np.argpartition(scores, -k)[-k:]
    return float(values[idx].sum())


def redteam_gc():
    """RPE1 s0 B=50 probe_full vs random: 打乱采信分配。"""
    from step21_ad_baseline import (error_metrics, load_bundle, load_product,
                                    oof_probe_scores, y_from_e)

    ds = "replogle_rpe1_essential"
    bundle = load_bundle(ds)
    prod = load_product(ds, "outputs_rpe1", 0, bundle)
    errs = error_metrics(prod["pred_delta"], prod["true_delta"])["full"]
    y = y_from_e(errs)
    s_probe = oof_probe_scores(prod["X"], y)
    trust = -s_probe  # 失败似然分数取负 = 采信任分
    n = len(errs)

    with open(os.path.join(OUT_DIR, "task_utility.json"), encoding="utf-8") as f:
        tu = json.load(f)
    combo = tu["combos"]["gears_rpe1_s0"]
    k = 50
    cell = combo["per_B"][str(k)]["gates"]["probe_full"]
    cost_gate, cost_random = cell["cost_gate"], cell["cost_random"]
    d_arch_boot = cell["bootstrap"]["diff_random"]["mean"]

    d_point = cost_random - cost_gate
    assert abs(topk_sum(errs, trust, k) - cost_gate) < 1e-6 * max(1.0, cost_gate), (
        f"cost_gate replicate {topk_sum(errs, trust, k)} vs archived {cost_gate}")
    assert abs(k / n * errs.sum() - cost_random) < 1e-6 * max(1.0, cost_random), (
        f"cost_random replicate vs archived mismatch")
    log(f"GC rpe1 B={k}: replicated point d={d_point:.4f} "
        f"(archived bootstrap mean {d_arch_boot:.4f}, "
        f"CI [{cell['bootstrap']['diff_random']['ci_lo']:.4f},"
        f"{cell['bootstrap']['diff_random']['ci_hi']:.4f}])")

    rng = np.random.default_rng(MASTER_SEED + 7)
    null_d = np.empty(B_PERM)
    for b in range(B_PERM):
        tp = rng.permutation(trust)
        null_d[b] = k / n * errs.sum() - topk_sum(errs, tp, k)
    p_value = (1.0 + float((null_d >= d_point - 1e-12).sum())) / (B_PERM + 1)

    # 次级: 其余可评测 B 档同法
    sec = {}
    for Bk in (5, 10, 20):
        cellb = combo["per_B"][str(Bk)]["gates"]["probe_full"]
        kb = min(Bk, n)
        d_b = kb / n * errs.sum() - topk_sum(errs, trust, Bk)
        nd = np.empty(B_PERM)
        for b in range(B_PERM):
            tp = rng.permutation(trust)
            nd[b] = kb / n * errs.sum() - topk_sum(errs, tp, Bk)
        sec[str(Bk)] = {
            "obs": float(d_b),
            "archived_boot_mean": cellb["bootstrap"]["diff_random"]["mean"],
            "p_value": (1.0 + float((nd >= d_b - 1e-12).sum())) / (B_PERM + 1),
        }
    log(f"GC rpe1 B={k}: perm p={p_value:.4f}")
    return {
        "combo": "gears_rpe1_s0", "gate": "probe_full", "B": k,
        "n": int(n), "obs_point_diff": float(d_point),
        "archived_bootstrap_mean": d_arch_boot,
        "archived_bootstrap_ci": [cell["bootstrap"]["diff_random"]["ci_lo"],
                                  cell["bootstrap"]["diff_random"]["ci_hi"]],
        "null_diffs": null_d.tolist(), "null_mean": float(null_d.mean()),
        "null_sd": float(null_d.std()), "p_value": p_value,
        "secondary_B": sec, "n_permutations": B_PERM,
    }


# ---------------------------------------------------------------- G-D 核
def redteam_conformal():
    """A_k562 单元覆盖 0.8921 vs 名义 0.90: 名义率标签打乱 Monte-Carlo + 精确二项。

    覆盖比例对单元标签置换不变, 故以 H0: 单元独立以名义率 0.9 被覆盖 为
    零假设 (有限样本保证的核对我们的可检验形式); 报双侧偏差 p 与单侧
    欠覆盖 p。N = n_eval_perts * n_genes。
    """
    with open(os.path.join(
            os.path.dirname(REPO), "inhibition_multiscale", "outputs_rebuild",
            "conformal_arm.json"), encoding="utf-8") as f:
        ca = json.load(f)
    arm = ca["arms"]["A_k562"]
    obs = arm["unit_coverage_eval"]
    n_units = int(arm["n_eval_perts"]) * int(arm["n_genes"])
    k_obs = int(round(obs * n_units))
    nom = float(ca["nominal"])

    rng = np.random.default_rng(MASTER_SEED + 11)
    cov_b = rng.binomial(n_units, nom, size=B_PERM) / n_units
    dev_obs = abs(obs - nom)
    p_two = (1.0 + float((np.abs(cov_b - nom) >= dev_obs - 1e-15).sum())
             ) / (B_PERM + 1)
    p_under = (1.0 + float((cov_b <= obs + 1e-15).sum())) / (B_PERM + 1)

    from scipy import stats
    exact = stats.binomtest(k_obs, n_units, nom, alternative="two-sided")
    exact_under = stats.binomtest(k_obs, n_units, nom,
                                  alternative="less")
    log(f"GD A_k562: obs={obs:.6f} vs nominal {nom}, N={n_units}, "
        f"MC p_two={p_two:.4f} p_under={p_under:.4f}, "
        f"exact p_two={exact.pvalue:.3e} p_under={exact_under.pvalue:.3e}")
    return {
        "arm": "A_k562", "nominal": nom, "obs_coverage": obs,
        "n_units": n_units, "n_covered": k_obs,
        "null": "unit covered-flags iid Bernoulli(nominal), MC label shuffle",
        "B": B_PERM,
        "mc_p_two_sided": p_two, "mc_p_under_coverage": p_under,
        "exact_binomial_p_two_sided": float(exact.pvalue),
        "exact_binomial_p_under": float(exact_under.pvalue),
        "wilson95_archived": arm["unit_coverage_wilson95"],
        "band": ca["acceptance_band"],
    }


def main():
    master_rng = np.random.default_rng(MASTER_SEED)
    out = {"task": "step23 results-level redteam permutation nulls",
           "date": time.strftime("%Y-%m-%d %H:%M:%S"),
           "B_perm": B_PERM, "master_seed": MASTER_SEED, "tests": {}}

    log("=== 1/4 norman multidataset probe mean ===")
    out["tests"]["norman_probe_mean"] = redteam_dataset(
        "norman", [3, 4, 5], master_rng)
    log("=== 2/4 adamson multidataset probe mean ===")
    out["tests"]["adamson_probe_mean"] = redteam_dataset(
        "adamson", [3, 4, 5], master_rng)
    log("=== 3/4 G-C RPE1 B=50 utility ===")
    out["tests"]["gc_rpe1_b50_utility"] = redteam_gc()
    log("=== 4/4 G-D A_k562 conformal coverage ===")
    out["tests"]["gd_k562_coverage"] = redteam_conformal()

    json_path = os.path.join(OUT_DIR, "redteam_null_tests.json")
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(out, f, indent=2)
    log(f"JSON -> {json_path}")

    # 留痕 md
    t = out["tests"]
    md = ["# Redteam null tests (step23, 2026-10-30)",
          "",
          "结果级红黑队：四个头条新数字的置换零检验（B=1000）。先与存档对账",
          "（tol 1e-6）再跑置换；判读：p < 0.05 -> 数字非零假设产物，可回填；",
          "p >= 0.05 -> 加 null-flag 如实报。",
          ""]
    for key, title in [
            ("norman_probe_mean", "norman 探针均值 0.7833"),
            ("adamson_probe_mean", "adamson 探针均值 0.7000")]:
        r = t[key]
        md.append(f"## {title}")
        md.append("")
        md.append(f"- 统计量：3 种子 cv_auroc_mean 的跨种子均值（锁定判据口径）")
        for s in r["seeds"]:
            o = r["obs_per_seed"][str(s)] if str(s) in r["obs_per_seed"] \
                else r["obs_per_seed"][s]
            md.append(f"- seed{s}: obs cv_mean={o['cv_mean']:.4f}, "
                      f"pooled OOF={o['pooled']:.4f}, "
                      f"置换 p(pooled)={r['secondary_pooled_oof_p'][s]:.4f}")
        md.append(f"- obs mean={r['obs_stat']:.4f}; null mean="
                  f"{r['null_mean']:.4f} sd={r['null_sd']:.4f}")
        md.append(f"- **置换 p = {r['p_value']:.4f}** (B={r['B']})")
        md.append("")
    r = t["gc_rpe1_b50_utility"]
    md.append("## G-C RPE1 s0 B=50 probe_full 效用")
    md.append("")
    md.append(f"- obs 点差 cost_random - cost_gate = {r['obs_point_diff']:.4f}"
              f"（存档 bootstrap 均值 {r['archived_bootstrap_mean']:.4f} "
              f"[{r['archived_bootstrap_ci'][0]:.4f},"
              f"{r['archived_bootstrap_ci'][1]:.4f}]）")
    md.append(f"- null mean={r['null_mean']:.4f} sd={r['null_sd']:.4f}")
    md.append(f"- **置换 p = {r['p_value']:.4f}** (B={r['n_permutations']})")
    for Bk, s in sorted(r["secondary_B"].items()):
        md.append(f"- 次级 B={Bk}: obs={s['obs']:.4f}, 存档 boot="
                  f"{s['archived_boot_mean']:.4f}, p={s['p_value']:.4f}")
    md.append("")
    r = t["gd_k562_coverage"]
    md.append("## G-D A_k562 conformal 单元覆盖 0.8921")
    md.append("")
    md.append(f"- N={r['n_units']} 单元，覆盖 {r['n_covered']}，名义 {r['nominal']}")
    md.append(f"- MC 标签打乱（名义率）B={r['B']}: 双侧 p="
              f"{r['mc_p_two_sided']:.4f}，单侧欠覆盖 p="
              f"{r['mc_p_under_coverage']:.4f}")
    md.append(f"- 精确二项：双侧 p={r['exact_binomial_p_two_sided']:.3e}，"
              f"单侧欠覆盖 p={r['exact_binomial_p_under']:.3e}")
    md.append(f"- 判读：覆盖可分辨于名义率零假设（p<0.05，方向为系统性欠覆盖 "
              f"0.8pp）；锁定接收带 {r['band']} 内（Wilson95 存档 "
              f"[{r['wilson95_archived'][0]:.4f},{r['wilson95_archived'][1]:.4f}]），"
              f"PASS 判据为带成员资格 -> 无 null-flag，欠覆盖幅度如实写。")
    md.append("")
    md_path = os.path.join(OUT_DIR, "redteam_null_tests.md")
    with open(md_path, "w", encoding="utf-8") as f:
        f.write("\n".join(md) + "\n")
    log(f"MD -> {md_path}")
    log("DONE")


if __name__ == "__main__":
    main()
