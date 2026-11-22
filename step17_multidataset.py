"""Step17: G-A 多数据集扩展驱动（预注册 prereg_expansion_ga.md）。

管线与 K562/RPE1 完全同构：
  Phase A  逐数据集逐种子调 step1_run.py（GEARS 训练 3 epochs + 基因级留出
           SPLIT_SEED=0 + 预测落盘）
  Phase B  逐数据集逐种子就地门控探针（复用 step6_gating_probe 核心；
           探针种子 = 模型种子）
  Phase C  逐数据集 gating_verdict.md + 跨数据集 summary + 升层判定

断点手册：outputs_multidataset/STATUS.md（或 smoke 根目录下）逐单元增量更新；
中断后重跑同一命令自动跳过已完成单元（幂等续跑）。

用法：
  python step17_multidataset.py --smoke     # 各数据集 1 epoch、种子 3，输出
                                            # outputs_multidataset_smoke/（不污染正式产物）
  python step17_multidataset.py --full      # 正式：3 epochs、种子 3/4/5，
                                            # 输出 outputs_multidataset/ + summary

纪律：GBK 控制台——stdout 只打 ASCII；所有文件读写显式 encoding='utf-8'；
单进程（GPU 串行，子进程逐个等待）；不打印 emoji。
"""
import argparse
import json
import os
import subprocess
import sys
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data", "gears_demo")

MODEL_SEEDS = [3, 4, 5]
CRITERION = 0.6  # 项目内约定门控线（锁定，见 prereg_expansion_ga.md）


def log(msg):
    print(f"[step17 {time.strftime('%H:%M:%S')}] {msg}", flush=True)


def status_path(root):
    return os.path.join(root, "STATUS.md")


def status_append(root, line):
    with open(status_path(root), "a", encoding="utf-8") as f:
        f.write(line + "\n")


def status_init(root, mode, datasets, epochs):
    os.makedirs(root, exist_ok=True)
    if not os.path.exists(status_path(root)):
        with open(status_path(root), "w", encoding="utf-8") as f:
            f.write(f"# STATUS - step17 G-A multidataset ({mode})\n\n")
            f.write(f"- started: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"- pid: {os.getpid()}\n")
            f.write(f"- mode: {mode}, epochs={epochs}, datasets={datasets}\n")
            f.write("- resume: rerun the same command; finished units are skipped\n\n")


# ---------------- Phase A: train + holdout predictions ----------------
def run_training(dataset, seed, out_dir, epochs, root):
    """调 step1_run.py；返回 (ok, log_file)。"""
    log_file = os.path.join(out_dir, f"train_seed{seed}.log")
    env = dict(os.environ)
    env.update({
        "DATA": dataset,
        "OUT_DIR": out_dir,
        "SEED": str(seed),
        "SPLIT_SEED": "0",
        "EPOCHS": str(epochs),
        "BATCH": "32",
        "TEST_BATCH": "128",
        "PYTHONIOENCODING": "utf-8",
    })
    os.makedirs(out_dir, exist_ok=True)
    cmd = [sys.executable, os.path.join(HERE, "step1_run.py")]
    log(f"TRAIN {dataset} seed={seed} -> {log_file}")
    t0 = time.time()
    with open(log_file, "w", encoding="utf-8", errors="replace") as lf:
        lf.write(f"CMD: {' '.join(cmd)}\nENV: DATA={dataset} SEED={seed} "
                 f"SPLIT_SEED=0 EPOCHS={epochs}\n\n")
        lf.flush()
        proc = subprocess.Popen(cmd, cwd=HERE, env=env,
                                stdout=lf, stderr=subprocess.STDOUT)
        rc = proc.wait()
    dt = time.time() - t0
    marker = os.path.join(out_dir, f"predictions_seed{seed}.csv")
    ok = (rc == 0 and os.path.exists(marker))
    status_append(root,
                  f"- [{'x' if ok else ' '}] train {dataset} seed={seed} "
                  f"rc={rc} {dt:.0f}s ({time.strftime('%H:%M:%S')})")
    log(f"TRAIN {dataset} seed={seed} done ok={ok} rc={rc} {dt:.0f}s")
    return ok


# ---------------- Phase B: in-situ gating probe ----------------
def run_probe_loo(X, y, feature_names, seed=0):
    """留一法（LOO）探针：留出扰动数 <10 时 5-fold 分层 CV 无定义
    （修正案 1，prereg_expansion_ga.md；判定阈值不变）。其余纪律与
    run_probe 相同：L2 logistic、fold 内标准化、仅模型侧特征。
    返回键与 run_probe 对齐；auroc_std 无 fold 概念，置 None。
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score, roc_curve
    from sklearn.model_selection import LeaveOneOut
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    import numpy as np

    clf = make_pipeline(
        StandardScaler(),
        LogisticRegression(C=1.0, max_iter=2000, random_state=seed))
    loo = LeaveOneOut()
    scores = np.zeros(len(y), dtype=np.float64)
    for tr, te in loo.split(X):
        clf.fit(X[tr], y[tr])
        scores[te] = clf.decision_function(X[te])
    fpr, tpr, _ = roc_curve(y, scores)
    pooled = roc_auc_score(y, scores)

    single = {}
    for j, name in enumerate(feature_names):
        s = np.zeros(len(y), dtype=np.float64)
        for tr, te in loo.split(X):
            c = make_pipeline(
                StandardScaler(),
                LogisticRegression(C=1.0, max_iter=2000,
                                   random_state=seed))
            c.fit(X[tr][:, [j]], y[tr])
            s[te] = c.decision_function(X[te][:, [j]])
        single[name] = {"auroc_mean": float(roc_auc_score(y, s)),
                        "auroc_std": None}
    return {"auroc_mean": float(pooled), "auroc_std": None,
            "auroc_per_fold": None, "pooled_auroc": float(pooled),
            "single_feature": single, "roc_fpr": fpr, "roc_tpr": tpr}


def run_probe(dataset, seed, out_dir, root):
    """复用 step6_gating_probe 核心；探针种子=模型种子。返回 auroc 或 None。"""
    import io
    import numpy as np
    import pandas as pd
    import anndata as ad
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from step6_gating_probe import (COVARIATES, SPECTRUM_FEATURES,
                                    compute_features, run_probe)

    z = np.load(os.path.join(out_dir, f"pred_expr_seed{seed}.npz"))
    pred_expr, perts = z["pred_expr"], z["pred_perturbations"].tolist()
    pred_genes = z["pred_genes"].tolist()

    h5 = os.path.join(DATA_DIR, dataset, "perturb_processed.h5ad")
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
    excluded = {}
    if pred_df["in_pert_graph"].nunique() <= 1:
        excluded["in_pert_graph"] = "constant column, dropped"
    else:
        feat["in_pert_graph"] = pred_df["in_pert_graph"].astype(int).values
    excluded["strength"] = "truth-derived, excluded (leakage guard)"
    excluded["pearson_delta_all"] = "truth-derived (label source), excluded"
    feature_names = SPECTRUM_FEATURES + [c for c in COVARIATES
                                         if c in feat.columns]
    X = feat[feature_names].values.astype(np.float64)
    assert not np.isnan(X).any(), "feature matrix has NaN"

    error = 1.0 - pred_df["pearson_delta"].values
    med = float(np.median(error))
    y = (error > med).astype(int)

    if y.sum() < 2 or (len(y) - y.sum()) < 2:
        # 少数类 <2：任何 CV/LOO 都无从训练——探针不可评估（修正案 1）
        auroc = None
        with open(os.path.join(out_dir, f"gating_probe_seed{seed}.json"), "w",
                  encoding="utf-8") as f:
            json.dump({
                "task": "G-A multidataset gating probe",
                "dataset": dataset, "model_seed": seed, "probe_seed": seed,
                "n_perts": int(len(y)), "evaluable": False,
                "non_evaluable_reason": (
                    f"minority class <2 (pos={int(y.sum())}, "
                    f"neg={len(y)-int(y.sum())}); no CV/LOO can train "
                    "(amendment 1, prereg_expansion_ga.md)"),
                "label_definition": "y = (1 - pearson_delta) > median",
                "error_median": med, "n_positive": int(y.sum()),
                "cv_auroc_mean": None,
                "criterion": f"per-dataset mean over seeds >= {CRITERION}",
            }, f, indent=2, ensure_ascii=False)
        log(f"PROBE {dataset} seed={seed} NON-EVALUABLE "
            f"(pos={int(y.sum())}, neg={len(y)-int(y.sum())})")
        status_append(root,
                      f"- [x] probe {dataset} seed={seed} NON-EVALUABLE "
                      f"n={len(y)} ({time.strftime('%H:%M:%S')})")
        return None

    if len(y) < 10:
        res = run_probe_loo(X, y, feature_names, seed=seed)
        cv_scheme = "loo (5<n<10 fallback; amendment 1)"
        log(f"PROBE {dataset} seed={seed} using LOO (n={len(y)}<10)")
    else:
        res = run_probe(X, y, feature_names, seed=seed)
        cv_scheme = "stratified 5-fold CV"
    auroc = res["auroc_mean"]
    best_single = max(res["single_feature"],
                      key=lambda k: res["single_feature"][k]["auroc_mean"])
    log(f"PROBE {dataset} seed={seed} auroc={auroc:.4f} n={len(y)}")

    with open(os.path.join(out_dir, f"gating_probe_seed{seed}.json"), "w",
              encoding="utf-8") as f:
        json.dump({
            "task": "G-A multidataset gating probe",
            "dataset": dataset, "model_seed": seed, "probe_seed": seed,
            "n_perts": int(len(y)),
            "evaluable": True,
            "cv_scheme": cv_scheme,
            "low_power_warning": bool(len(y) < 10),
            "label_definition": "y = (1 - pearson_delta) > median",
            "features": feature_names,
            "excluded_features": excluded,
            "error_median": med, "n_positive": int(y.sum()),
            "cv_auroc_mean": auroc, "cv_auroc_std": res["auroc_std"],
            "cv_auroc_per_fold": res["auroc_per_fold"],
            "pooled_cv_auroc": res["pooled_auroc"],
            "best_single_feature": best_single,
            "single_feature_auroc": res["single_feature"],
            "criterion": f"per-dataset mean over seeds >= {CRITERION}",
        }, f, indent=2, ensure_ascii=False)

    fig, ax = plt.subplots(figsize=(6, 4.5))
    ax.plot(res["roc_fpr"], res["roc_tpr"],
            label=f"pooled ROC (AUROC={res['pooled_auroc']:.3f})")
    ax.plot([0, 1], [0, 1], "k--", lw=1, label="chance (0.5)")
    ax.set_xlabel("false positive rate")
    ax.set_ylabel("true positive rate")
    std_txt = (f"+- {res['auroc_std']:.3f}"
               if res["auroc_std"] is not None else "LOO")
    ax.set_title(f"{dataset} seed{seed} gating probe: "
                 f"AUROC={auroc:.3f} {std_txt}, n={len(y)}")
    ax.legend()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, f"gating_probe_cv_seed{seed}.png"),
                dpi=300)
    plt.close(fig)

    status_append(root,
                  f"- [x] probe {dataset} seed={seed} auroc={auroc:.4f} "
                  f"({time.strftime('%H:%M:%S')})")
    return auroc


# ---------------- Phase C: verdicts + summary ----------------
def fmt_a(v):
    return "non-evaluable" if v is None else f"{v:.4f}"


def dataset_verdict(dataset, out_dir, aurocs):
    """写 per-dataset gating_verdict.md；返回 dict 汇总。
    aurocs: {seed: float|None}；None = 探针不可评估（修正案 1），按未成立计数。
    """
    seeds = sorted(aurocs)
    vals = [aurocs[s] for s in seeds if aurocs[s] is not None]
    if vals:
        mean_a = sum(vals) / len(vals)
        passed = (len(vals) == len(MODEL_SEEDS)) and mean_a >= CRITERION
    else:
        mean_a = None
        passed = False
    lines = [f"# gating_verdict - {dataset} (G-A expansion)\n",
             f"- model/probe seeds: {seeds}",
             f"- n perts: see per-seed json",
             "- per-seed CV AUROC: " + ", ".join(
                 f"seed{s}={fmt_a(aurocs[s])}" for s in seeds),
             f"- mean AUROC = {fmt_a(mean_a)}",
             f"- criterion: mean >= {CRITERION} -> signal holds in this dataset",
             f"- verdict: {'PASS' if passed else 'FAIL / NOT ESTABLISHED'}\n"]
    with open(os.path.join(out_dir, "gating_verdict.md"), "w",
              encoding="utf-8") as f:
        f.write("\n".join(lines))
    return {"dataset": dataset, "per_seed": aurocs, "mean": mean_a,
            "pass": passed}


def write_summary(root, results, mode, epochs):
    n_pass = sum(1 for r in results if r["pass"])
    n_ds = len(results)
    tier_pass = n_pass >= 2
    lines = [f"# multidataset_summary ({mode}, epochs={epochs})\n",
             f"- generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
             f"- criterion: per-dataset mean of seed AUROCs >= {CRITERION}; "
             f"tier needs >=2/3 datasets\n",
             "## per-dataset\n"]
    for r in results:
        ps = ", ".join(f"seed{s}={fmt_a(v)}" for s, v in
                       sorted(r["per_seed"].items()))
        lines.append(f"- {r['dataset']}: mean={fmt_a(r['mean'])} "
                     f"({ps}) -> {'PASS' if r['pass'] else 'FAIL / NOT ESTABLISHED'}")
    lines += ["\n## tier verdict (locked)\n",
              f"- datasets passing: {n_pass}/{n_ds}",
              f"- tier verdict: {'PASS (>=2/3) -> cross-5-dataset claim '
                               if tier_pass else 'FAIL (<2/3) -> stratified '
                               'reporting'}",
              "", "## deviation log", "(filled after full execution)"]
    summary = os.path.join(root, "multidataset_summary.md")
    with open(summary, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")
    log(f"SUMMARY -> {summary} (tier_pass={tier_pass})")
    return tier_pass


# ---------------- main ----------------
class Tee(object):
    """同时写原 stdout 与日志文件的镜像流（GBK 安全：日志 UTF-8）。"""

    def __init__(self, orig, path):
        self.orig = orig
        self.f = open(path, "a", encoding="utf-8", buffering=1)

    def write(self, s):
        try:
            self.orig.write(s)
        except (UnicodeEncodeError, OSError):
            pass  # GBK 控制台遇到非 GBK 字符不炸
        self.f.write(s)
        return len(s)

    def flush(self):
        try:
            self.orig.flush()
        except Exception:
            pass
        self.f.flush()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--datasets", default="norman,adamson,dixit")
    args = ap.parse_args()
    datasets = [d.strip() for d in args.datasets.split(",") if d.strip()]

    if args.smoke and args.full:
        raise SystemExit("choose one of --smoke/--full")
    mode = "smoke" if args.smoke else "full"
    epochs = 1 if args.smoke else 3
    seeds = [3] if args.smoke else MODEL_SEEDS
    root = os.path.join(HERE, "outputs_multidataset_smoke" if args.smoke
                        else "outputs_multidataset")
    os.makedirs(root, exist_ok=True)
    # detached（pythonw）下控制台不存在：镜像打印到 <root>/<mode>_queue.log
    _qlog = os.path.join(root, f"{mode}_queue.log")
    sys.stdout = Tee(sys.__stdout__ if sys.__stdout__ else open(os.devnull, "w"), _qlog)
    sys.stderr = Tee(sys.__stderr__ if sys.__stderr__ else open(os.devnull, "w"), _qlog)

    status_init(root, mode, datasets, epochs)
    log(f"mode={mode} epochs={epochs} seeds={seeds} datasets={datasets}")
    log(f"pid={os.getpid()} root={root}")

    # sanity: data present
    for ds in datasets:
        if not os.path.isdir(os.path.join(DATA_DIR, ds)):
            status_append(root, f"- DATA MISSING: {ds}")
            raise SystemExit(f"data dir missing: {DATA_DIR}/{ds}")

    failures = []
    for ds in datasets:
        out_dir = os.path.join(root, ds)
        os.makedirs(out_dir, exist_ok=True)
        for seed in seeds:
            marker = os.path.join(out_dir, f"predictions_seed{seed}.csv")
            if not os.path.exists(marker):
                ok = run_training(ds, seed, out_dir, epochs, root)
                if not ok:
                    failures.append(f"train {ds} seed{seed}")
                    status_append(root, f"- FAIL train {ds} seed{seed}")
                    continue
            else:
                log(f"SKIP train {ds} seed{seed} (done)")
        for seed in seeds:
            pmarker = os.path.join(out_dir, f"gating_probe_seed{seed}.json")
            tmarker = os.path.join(out_dir, f"predictions_seed{seed}.csv")
            if os.path.exists(pmarker):
                log(f"SKIP probe {ds} seed{seed} (done)")
                continue
            if not os.path.exists(tmarker):
                failures.append(f"probe {ds} seed{seed} (no train output)")
                continue
            try:
                run_probe(ds, seed, out_dir, root)
            except Exception:
                failures.append(f"probe {ds} seed{seed}")
                status_append(root, f"- FAIL probe {ds} seed{seed}")
                with open(status_path(root), "a", encoding="utf-8") as f:
                    f.write("```\n" + traceback.format_exc() + "```\n")

    # Phase C (only uses completed probes; missing units recorded as failures)
    results = []
    for ds in datasets:
        out_dir = os.path.join(root, ds)
        aurocs = {}
        for seed in seeds:
            p = os.path.join(out_dir, f"gating_probe_seed{seed}.json")
            if os.path.exists(p):
                with open(p, encoding="utf-8") as f:
                    aurocs[seed] = json.load(f)["cv_auroc_mean"]
        if aurocs:
            results.append(dataset_verdict(ds, out_dir, aurocs))
    if results:
        write_summary(root, results, mode, epochs)
    if args.full and not failures:
        # 正式模式：summary 落到 outputs_rebuild/（预注册产物路径）
        dest = os.path.join(HERE, "outputs_rebuild", "multidataset_summary.md")
        src = os.path.join(root, "multidataset_summary.md")
        if os.path.exists(src):
            with open(src, encoding="utf-8") as f:
                content = f.read()
            with open(dest, "w", encoding="utf-8") as f:
                f.write(content)
            log(f"summary copied -> {dest}")
    status_append(root,
                  f"- ALL DONE {time.strftime('%Y-%m-%d %H:%M:%S')} "
                  f"failures={failures if failures else 'none'}")
    log(f"DONE failures={failures if failures else 'none'}")


if __name__ == "__main__":
    main()
