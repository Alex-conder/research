# step6f_crossmodel.py — 跨模型迁移验证：GEARS 训好的门控探针能否门控
# scGenePT / scFoundation 的预测？（分析脚本，服务器产物落地后直接运行）
#
# 背景与产物契约（源自 outputs_crossmodel/scgenept/RUNBOOK.md 与
# run_scgenept_k562.py，服务器端已锁定）：
#   服务器每模型每数据集产出 <out>_predictions.npz + _labels.csv + _meta.json：
#     pred_expr        (n_pert, n_matched) 匹配基因空间的预测均值（缺 condition 行为 NaN）
#     pred_genes       匹配基因名（与 ctrl_mean_matched 同序）
#     ctrl_mean_matched  ctrl 均值（匹配空间）
#     pred_perturbations held-out 扰动名（与本地 split_seed=0 同名单同序）
#     target_in_input_space  目标基因是否在模型输入空间（不在则预测≈ctrl，数据事实留痕）
#   _labels.csv 列：perturbation, pearson_delta, pearson_delta_all, n_cells,
#                 target_in_input_space（pearson_delta 口径与 step1_run.py 完全一致：
#                 top-50 |true_delta| 基因上的 Pearson）
#
# 跨模型特征可得性核查结论（2026-10-15，见 step6_gating_probe.py / step6_run.py）：
#   探针 9 特征 = 8 个预测谱统计（delta = pred_expr − ctrl_mean 的纯函数）
#   + log_n_cells。**全部不依赖 GEARS 内部量**（GEARS 原生 uncertainty 列在
#   step6 中已因全 NaN 排除；in_pert_graph 为常数已剔除）。因此只要模型输出
#   预测表达谱 + ctrl 均值 + n_cells，9 特征全部可算。如实记录的口径偏离：
#   (1) 特征在模型的匹配基因空间上计算（scGenePT K562 为 4813 基因，本地
#       GEARS 为 5000）——特征分布随基因空间略有偏移；
#   (2) 迁移探针系数为本地全量 in-sample 拟合（部署式用法，非 CV），
#       规格与预注册探针完全相同（StandardScaler + L2 logistic, C=1.0）。
#
# 分析：
#   A 迁移测试：本地训好的探针不重训，直接打分跨模型特征 → AUROC / P@20
#   B 就地参照：用目标模型自己的特征就地训同规格探针（5 个 CV 划分种子 OOF
#     均值）→ AUROC / P@20；同时给出两基线：|delta| L2 幅度（无训练）、
#     5 划分种子就地探针的集成方差
#   C 效应强度三分位分层：真值 strength（step1_run.py:106-115 口径，从本地
#     h5ad 算）三分位，各层迁移 AUROC（层 n<30 或单类则记 None）
#
# 产物：outputs_crossmodel/analysis/{model}_{dataset}.json
#       outputs_crossmodel/analysis/crossmodel_summary.png（≥1 个结果时出图，300dpi）
#
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step6f_crossmodel.py --dry-run
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step6f_crossmodel.py
import argparse
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.preprocessing import StandardScaler

from step6_gating_probe import COVARIATES, SPECTRUM_FEATURES, compute_features

REPO = os.path.dirname(os.path.abspath(__file__))
SEED = 0
IN_SITU_SEEDS = [0, 1, 2, 3, 4]
SMALL_N = 30

# scgpt 的 K562 产物 2026-10-16 在服务器跑完并落地
# （outputs_crossmodel/scgpt/），2026-10-17 起纳入分析。
# uce / geneformer 产物虽在，但为 embedding-only（in silico deletion 细胞
# 嵌入，无表达谱预测，_labels.csv pearson 列为 NaN，见各自 npz 的
# output_type），无法进入表达空间分析链，故不列入 MODELS。
MODELS = ["scgenept", "scfoundation", "scgpt"]
DATASETS = {
    "k562": {"data": "replogle_k562_essential", "out": "outputs",
             "tag": "Norman (K562)"},
    "rpe1": {"data": "replogle_rpe1_essential", "out": "outputs_rpe1",
             "tag": "RPE1"},
}
ANALYSIS_DIR = os.path.join(REPO, "outputs_crossmodel", "analysis")


# ---------------------------------------------------------------- 工具函数
def precision_at(scores, y, cov):
    """top-cov 分数覆盖带内，y=1（最差一半）的占比（与 step6b/6d 同）。"""
    order = np.argsort(-scores)
    k = max(1, int(round(cov * len(y))))
    return float(y[order[:k]].mean())


def oof_scores(X, y, split_seed):
    """同规格就地探针在某 CV 划分种子下的 OOF 概率。"""
    Xs = StandardScaler().fit_transform(X)
    cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=split_seed)
    return cross_val_predict(LogisticRegression(C=1.0, max_iter=2000),
                             Xs, y, cv=cv, method="predict_proba")[:, 1]


def load_local_xy(out_dir):
    """本地 GEARS 数据集的特征矩阵与标签（口径同 step6d.load_dataset）。"""
    from step6d_naive_baselines import load_dataset
    feat, feature_names, y, pred_df, _ = load_dataset(
        DATASETS["k562" if out_dir == "outputs" else "rpe1"]["data"], out_dir)
    return feat[feature_names].values.astype(np.float64), y, feature_names


def export_transfer_probes():
    """本地全量 in-sample 拟合迁移探针并落盘（部署式用法；规格同预注册）。"""
    saved = []
    for ds_key, ds in DATASETS.items():
        path = os.path.join(REPO, ds["out"], "gating_probe_transfer.json")
        X, y, feature_names = load_local_xy(ds["out"])
        pipe = StandardScaler().fit(X)
        Xs = pipe.transform(X)
        clf = LogisticRegression(C=1.0, max_iter=2000).fit(Xs, y)
        rec = {
            "task": "transfer probe: local GEARS probe fit in-sample "
                    "(deployment-style; same spec as pre-registered probe)",
            "dataset": ds["data"], "n": int(len(y)),
            "feature_names": feature_names,
            "scaler_mean": pipe.mean_.tolist(),
            "scaler_scale": pipe.scale_.tolist(),
            "coef": clf.coef_[0].tolist(), "intercept": float(clf.intercept_[0]),
            "fit": "StandardScaler + LogisticRegression(C=1.0, max_iter=2000), "
                   "in-sample on all held-out queries (no CV)",
        }
        with open(path, "w", encoding="utf-8") as f:
            json.dump(rec, f, ensure_ascii=False, indent=2)
        saved.append(path)
        print(f"[probe] exported transfer probe -> {path}")
    return saved


def load_transfer_probe(out_dir):
    path = os.path.join(REPO, out_dir, "gating_probe_transfer.json")
    if not os.path.exists(path):
        export_transfer_probes()
    with open(path, encoding="utf-8") as f:
        rec = json.load(f)
    return rec


def transfer_scores(rec, X):
    """迁移探针打分：(X - mean)/scale @ coef + intercept。"""
    Xs = (X - np.array(rec["scaler_mean"])) / np.array(rec["scaler_scale"])
    return Xs @ np.array(rec["coef"]) + rec["intercept"]


def compute_strength(data_name, perts):
    """真值效应强度（step1_run.py:106-115 口径）：留出组各细胞到 ctrl 均值
    的 L2 距离取平均。只用于分层诊断，绝不进特征（泄漏纪律）。"""
    import anndata as ad
    adata = ad.read_h5ad(
        f"{REPO}/data/gears_demo/{data_name}/perturb_processed.h5ad",
        backed="r")
    obs_cond = adata.obs["condition"].astype(str).values
    ctrl_mean = np.asarray(
        adata[obs_cond == "ctrl"].X.mean(axis=0)).ravel()
    cond_set = set(obs_cond)
    strength = {}
    for g in perts:
        cond = next((c for c in (f"ctrl+{g}", f"{g}+ctrl") if c in cond_set),
                    None)
        if cond is None:
            strength[g] = np.nan
            continue
        sub = adata[obs_cond == cond]
        d = sub.X.toarray() - ctrl_mean
        strength[g] = float(np.linalg.norm(d, axis=1).mean())
    adata.file.close()
    return np.array([strength[g] for g in perts])


# ---------------------------------------------------------------- 产物发现
def discover_products(model):
    """在 outputs_crossmodel/{model}/ 下找完整产物集，返回可用与缺失清单。"""
    mdir = os.path.join(REPO, "outputs_crossmodel", model)
    found = {}
    if not os.path.isdir(mdir):
        return found, [f"{mdir}/ 目录不存在"]
    missing = []
    for npz in sorted(glob.glob(os.path.join(mdir, "*_predictions.npz"))):
        base = npz[: -len("_predictions.npz")]
        stem = os.path.basename(base).lower()
        ds_key = ("k562" if "k562" in stem else
                  "rpe1" if "rpe1" in stem else None)
        if ds_key is None:  # 文件名判断不了时看 meta.json
            meta_f = base + "_meta.json"
            if os.path.exists(meta_f):
                meta = json.load(open(meta_f, encoding="utf-8"))
                txt = json.dumps(meta).lower()
                ds_key = ("k562" if "k562" in txt else
                          "rpe1" if "rpe1" in txt else None)
        labels = base + "_labels.csv"
        if ds_key is None:
            missing.append(f"{os.path.basename(npz)}: 无法判断数据集"
                           f"（文件名与 _meta.json 均无 k562/rpe1 线索）")
            continue
        if not os.path.exists(labels):
            missing.append(f"{os.path.basename(npz)}: 缺 "
                           f"{os.path.basename(labels)}")
            continue
        found[ds_key] = {"npz": npz, "labels": labels, "base": base}
    if not found:
        missing.append(f"{mdir}/ 下没有完整产物集（需要 "
                       f"*_predictions.npz + 同名 _labels.csv，"
                       f"文件名含 k562/rpe1 或 meta 可判数据集）")
    return found, missing


def load_crossmodel(npz_path, labels_path):
    """读服务器产物，对齐成 (feat_df, y, extra)。NaN 行（缺真值 condition）
    剔除并留痕。"""
    z = np.load(npz_path, allow_pickle=True)
    pred_expr = z["pred_expr"]
    perts = [str(g) for g in z["pred_perturbations"]]
    pred_genes = [str(g) for g in z["pred_genes"]]
    if "ctrl_mean_matched" in z.files:
        ctrl_mean = np.asarray(z["ctrl_mean_matched"], dtype=np.float64)
    elif "ctrl_mean_full" in z.files:  # scFoundation 概念验证批命名
        ctrl_mean = np.asarray(z["ctrl_mean_full"], dtype=np.float64)
    else:
        raise KeyError("npz 缺少 ctrl_mean（既无 ctrl_mean_matched 也无 ctrl_mean_full）")
    in_space = (z["target_in_input_space"].astype(bool)
                if "target_in_input_space" in z.files else None)
    n_genes = len(pred_genes)
    # 全 NaN 列（如 scFoundation 未映射基因）整列剔除，再查行
    nan_col = np.isnan(pred_expr).all(axis=0)
    if nan_col.any():
        pred_expr = pred_expr[:, ~nan_col]
        pred_genes = [g for g, keep in zip(pred_genes, ~nan_col) if keep]
        ctrl_mean = ctrl_mean[~nan_col]
        if in_space is not None:
            pass  # in_space 是逐扰动布尔，与基因无关
        print(f"  [load] 剔除全 NaN 基因列 {int(nan_col.sum())}，保留 {pred_expr.shape[1]} 基因")

    lab = pd.read_csv(labels_path)
    lab = lab.set_index("perturbation").loc[perts]

    err = 1.0 - lab["pearson_delta"].values
    ok = ~np.isnan(err)
    ok &= ~np.isnan(pred_expr).any(axis=1)
    dropped = int((~ok).sum())
    pred_expr, perts_arr = pred_expr[ok], np.array(perts)[ok]
    err = err[ok]
    y = (err > np.median(err)).astype(int)
    n_cells = lab["n_cells"].values[ok]

    feat = compute_features(pred_expr, list(perts_arr), pred_genes, ctrl_mean)
    feat["log_n_cells"] = np.log1p(n_cells)
    feature_names = SPECTRUM_FEATURES + \
        [c for c in COVARIATES if c in feat.columns]
    X = feat[feature_names].values.astype(np.float64)
    assert not np.isnan(X).any(), "特征矩阵含 NaN"

    extra = {"n_total": int(len(ok)), "n_used": int(len(y)),
             "n_dropped_no_truth": dropped,
             "n_genes_feature_space": int(n_genes),
             "in_input_space": (in_space[ok].tolist()
                                if in_space is not None else None)}
    return X, y, feature_names, list(perts_arr), extra


# ---------------------------------------------------------------- 主分析
def run_one(model, ds_key, prod, make_figure_data):
    ds = DATASETS[ds_key]
    X, y, feature_names, perts, extra = load_crossmodel(
        prod["npz"], prod["labels"])
    n = len(y)
    print(f"[{model}/{ds_key}] n={n}（总 {extra['n_total']}，剔除无真值 "
          f"{extra['n_dropped_no_truth']}），特征空间 "
          f"{extra['n_genes_feature_space']} 基因")

    # --- A 迁移探针（不重训）---
    rec = load_transfer_probe(ds["out"])
    assert rec["feature_names"] == feature_names, \
        f"特征列不一致: {rec['feature_names']} vs {feature_names}"
    s_transfer = transfer_scores(rec, X)

    # --- B 就地探针 + 两基线 ---
    in_situ = np.column_stack([oof_scores(X, y, s) for s in IN_SITU_SEEDS])
    s_in_situ_mean = in_situ.mean(axis=1)
    s_ensvar = in_situ.var(axis=1)
    j_l2 = feature_names.index("delta_l2")
    s_mag = X[:, j_l2]

    def metrics(scores):
        return {"auroc": float(roc_auc_score(y, scores)),
                "p_at_20": precision_at(scores, y, 0.20)}

    transfer = metrics(s_transfer)
    in_situ_m = {
        "auroc": float(np.mean([roc_auc_score(y, in_situ[:, s])
                                for s in range(len(IN_SITU_SEEDS))])),
        "p_at_20": float(np.mean([precision_at(in_situ[:, s], y, 0.20)
                                  for s in range(len(IN_SITU_SEEDS))])),
    }
    res = {
        "task": "cross-model transfer of gating probe (step6f)",
        "model": model, "dataset": ds["data"], "dataset_tag": ds["tag"],
        "n": n, "base_rate": 0.5,
        "label": "y = (1 - pearson_delta) > median (worst half), "
                 "same as step6",
        "feature_space_note":
            f"features computed on the model's matched gene space "
            f"({extra['n_genes_feature_space']} genes), not the local "
            f"GEARS 5000-gene space; all 9 features are pure functions of "
            f"pred_expr/ctrl_mean/n_cells, no GEARS-specific quantity",
        "transfer_probe_note":
            "local GEARS probe fit in-sample (deployment-style), identical "
            "spec; see gating_probe_transfer.json",
        "dropped_no_truth": extra["n_dropped_no_truth"],
        "target_in_input_space": extra["in_input_space"],
        "transfer": transfer,
        "in_situ": in_situ_m,
        "baselines": {
            "magnitude_l2": metrics(s_mag),
            "ensemble_variance": metrics(s_ensvar),
        },
    }

    # --- C 效应强度三分位分层（真值 strength，仅作诊断分层）---
    strength = compute_strength(ds["data"], perts)
    strata = {"n_nan_strength": int(np.isnan(strength).sum())}
    valid = ~np.isnan(strength)
    if valid.sum() >= 3:
        try:
            tl = pd.qcut(pd.Series(strength[valid]), 3,
                         labels=["low", "mid", "high"],
                         duplicates="drop").astype(str).values
        except ValueError:
            tl = None
        if tl is not None:
            for lvl in ["low", "mid", "high"]:
                if lvl not in tl:
                    continue
                m = np.zeros(len(y), dtype=bool)
                m[valid] = tl == lvl
                if m.sum() < SMALL_N or len(np.unique(y[m])) < 2:
                    strata[lvl] = {"n": int(m.sum()),
                                   "transfer_auroc": None}
                    continue
                strata[lvl] = {
                    "n": int(m.sum()),
                    "transfer_auroc": float(roc_auc_score(y[m], s_transfer[m])),
                    "worst_half_rate": float(y[m].mean())}
    res["strength_tertile_transfer"] = strata
    res["strength_note"] = ("strength = mean per-cell L2 to ctrl mean "
                            "(step1_run.py:106-115), truth-derived, "
                            "diagnostic strata only")

    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    out_json = os.path.join(ANALYSIS_DIR, f"{model}_{ds_key}.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    line = "  ".join(
        f"{k}={v:.3f}" for k, v in
        [("transfer", transfer["auroc"]), ("in_situ", in_situ_m["auroc"]),
         ("mag_l2", res["baselines"]["magnitude_l2"]["auroc"]),
         ("ensvar", res["baselines"]["ensemble_variance"]["auroc"])])
    print(f"[{model}/{ds_key}] AUROC: {line}")
    make_figure_data.append(
        {"model": model, "ds_key": ds_key, "tag": f"{model} {ds['tag']}",
         "transfer": transfer, "in_situ": in_situ_m,
         "mag": res["baselines"]["magnitude_l2"],
         "ensvar": res["baselines"]["ensemble_variance"]})


def make_summary_figure(figure_data):
    if not figure_data:
        return
    fig, ax = plt.subplots(figsize=(max(6.0, 2.2 * len(figure_data)), 3.8))
    methods = [("transfer", "GEARS probe (transfer)", "#1f77b4"),
               ("in_situ", "in-situ probe", "#2ca02c"),
               ("mag", "|delta| L2 (no train)", "#d62728"),
               ("ensvar", "ensemble variance", "#9467bd")]
    x = np.arange(len(figure_data))
    w = 0.2
    for mi, (key, label, color) in enumerate(methods):
        aucs = [d[key]["auroc"] for d in figure_data]
        p20s = [d[key]["p_at_20"] for d in figure_data]
        b1 = ax.bar(x + (mi - 1.5) * w - w / 4, aucs, w / 2, color=color,
                    edgecolor="gray", lw=0.4, label=f"{label} AUROC")
        b2 = ax.bar(x + (mi - 1.5) * w + w / 4, p20s, w / 2, color=color,
                    edgecolor="gray", lw=0.4, alpha=0.45,
                    label=f"{label} P@20")
        for b in list(b1) + list(b2):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.008,
                    f"{b.get_height():.2f}", ha="center", fontsize=6,
                    rotation=90)
    ax.axhline(0.5, color="red", ls=":", lw=1.2)
    ax.set_xticks(x)
    ax.set_xticklabels([d["tag"] for d in figure_data], fontsize=8)
    ax.set_ylim(0.0, 1.0)
    ax.set_ylabel("score (chance = 0.5)")
    ax.set_title("Cross-model transfer of the error gate")
    ax.text(0.99, 1.02, "solid = OOF AUROC, pale = P@20",
            transform=ax.transAxes, ha="right", fontsize=7, color="gray")
    ax.grid(axis="y", alpha=0.25)
    ax.legend(fontsize=7, ncol=2, loc="upper right")
    fig.tight_layout()
    path = os.path.join(ANALYSIS_DIR, "crossmodel_summary.png")
    fig.savefig(path, dpi=300)
    print("[SAVED]", path)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true",
                    help="只检查产物齐备性，缺什么报什么，不跑分析")
    ap.add_argument("--export-probes", action="store_true",
                    help="只导出本地迁移探针（gating_probe_transfer.json）")
    args = ap.parse_args()

    if args.export_probes:
        export_transfer_probes()
        return

    all_missing = []
    figure_data = []
    for model in MODELS:
        found, missing = discover_products(model)
        for m in missing:
            print(f"[MISSING] {m}")
        all_missing.extend(f"{model}: {m}" for m in missing)
        if args.dry_run:
            for ds_key, prod in found.items():
                print(f"[OK] {model}/{ds_key}: "
                      f"{os.path.basename(prod['npz'])} + _labels.csv 已就位")
            continue
        for ds_key, prod in found.items():
            run_one(model, ds_key, prod, figure_data)

    if args.dry_run:
        print(f"\n[DRY-RUN] 完整产物集 {sum(1 for m in MODELS if not any(x.startswith(m + ':') for x in all_missing))}/{len(MODELS)} 个模型齐备")
        return

    make_summary_figure(figure_data)
    if all_missing:
        print("\n[WARNING] 以下产物缺失，未分析：")
        for m in all_missing:
            print("  -", m)


if __name__ == "__main__":
    main()
