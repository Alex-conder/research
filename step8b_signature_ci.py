# step8b_signature_ci.py — 误差签名的 bootstrap 置信区间（森林图）
#
# 背景：outputs_crossmodel/analysis/error_signatures.json 已有 5 个误差签名
# （每签名 = 就地拟合 L2-logistic 探针的 9 个标准化系数，in-sample；
# scgpt_k562 由 step8e 于 2026-10-17 增补。uce/geneformer 为
# embedding-only 产物、无 pearson 标签，不进签名链）。
# 本脚本回答"这些系数有多稳"：对每个签名做 B=200 次扰动级 bootstrap
# （有放回重采样 → 重拟合同规格探针 StandardScaler + LogisticRegression
# (C=1.0, max_iter=2000)，in-sample），取每个系数的 95% percentile CI。
#
# 口径（严格复用，不自造加载逻辑）：
#   - GEARS(K562): step6d_naive_baselines.load_dataset("replogle_k562_essential",
#     "outputs")（seed 0 预测，step6 同款特征与标签）
#   - scgenept_k562 / scgenept_rpe1 / scf_k562:
#     step6f_crossmodel.load_crossmodel(npz, labels)
#   - 点估计核对：先在全量数据上重拟合，与 error_signatures.json 比对
#     （同一规格 + 同一数据 → 应数值一致；不一致即报错退出）。
#
# 随机源：BOOT_MASTER_SEED = 20261020（专属 seed；每签名由主 seed 派生
# 独立 Generator，互不复用同一随机流）。
#
# 判据留痕：CI 跨零的系数个数如实统计（n=45 的 scFoundation 必然很宽，
# 不隐藏、不截断）。
#
# 产物：outputs_crossmodel/analysis/signature_ci.json + signature_ci.png
#       （300 dpi，英文，5 模型 × 9 特征的系数点估计 + CI 森林图）
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step8b_signature_ci.py
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

REPO = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.join(REPO, "outputs_crossmodel", "analysis")

BOOT_MASTER_SEED = 20261020  # bootstrap 专属主 seed
B_BOOT = 200

FEATURE_ORDER = ["delta_l2", "delta_l1", "delta_max", "delta_top50_mean",
                 "n_de_025", "frac_de_025", "neg_frac", "delta_entropy",
                 "log_n_cells"]
SIG_ORDER = ["GEARS(K562)", "scgenept_k562", "scgenept_rpe1", "scf_k562",
             "scgpt_k562"]
SIG_LABELS = {"GEARS(K562)": "GEARS / K562",
              "scgenept_k562": "scGenePT / K562",
              "scgenept_rpe1": "scGenePT / RPE1",
              "scf_k562": "scFoundation / K562",
              "scgpt_k562": "scGPT / K562"}

SG_DIR = os.path.join(REPO, "outputs_crossmodel", "scgenept")
SCF_DIR = os.path.join(REPO, "outputs_crossmodel", "scfoundation")
SCGPT_DIR = os.path.join(REPO, "outputs_crossmodel", "scgpt")
UCE_DIR = os.path.join(REPO, "outputs_crossmodel", "uce")
GF_DIR = os.path.join(REPO, "outputs_crossmodel", "geneformer")
SOURCES = {
    "GEARS(K562)": ("gears", "replogle_k562_essential", "outputs"),
    "scgenept_k562": ("crossmodel", os.path.join(SG_DIR, "scgenept_k562_predictions.npz"),
                      os.path.join(SG_DIR, "scgenept_k562_labels.csv")),
    "scgenept_rpe1": ("crossmodel", os.path.join(SG_DIR, "scgenept_rpe1_predictions.npz"),
                      os.path.join(SG_DIR, "scgenept_rpe1_labels.csv")),
    "scf_k562": ("crossmodel", os.path.join(SCF_DIR, "scf_k562_predictions.npz"),
                 os.path.join(SCF_DIR, "scf_k562_labels.csv")),
    "scgpt_k562": ("crossmodel", os.path.join(SCGPT_DIR, "scgpt_k562_predictions.npz"),
                   os.path.join(SCGPT_DIR, "scgpt_k562_labels.csv")),
    "uce_k562": ("crossmodel", os.path.join(UCE_DIR, "uce_k562_predictions.npz"),
                 os.path.join(UCE_DIR, "uce_k562_labels.csv")),
    "geneformer_k562": ("crossmodel", os.path.join(GF_DIR, "geneformer_k562_predictions.npz"),
                        os.path.join(GF_DIR, "geneformer_k562_labels.csv")),
}


def load_xy(sig_name):
    """按签名名加载 (X, y, feature_names)，口径与 step6d/step6f 完全一致。"""
    kind = SOURCES[sig_name][0]
    if kind == "gears":
        from step6d_naive_baselines import load_dataset
        data_name, out_dir = SOURCES[sig_name][1], SOURCES[sig_name][2]
        feat, feature_names, y, _pred_df, _ = load_dataset(data_name, out_dir)
        X = feat[feature_names].values.astype(np.float64)
    else:
        from step6f_crossmodel import load_crossmodel
        X, y, feature_names, _perts, _extra = load_crossmodel(
            SOURCES[sig_name][1], SOURCES[sig_name][2])
    assert list(feature_names) == FEATURE_ORDER, \
        f"{sig_name} 特征列与签名顺序不一致: {feature_names}"
    return X, y, feature_names


def fit_coef(X, y):
    """同规格就地探针（in-sample）：StandardScaler + L2 logistic, C=1.0。"""
    Xs = StandardScaler().fit_transform(X)
    clf = LogisticRegression(C=1.0, max_iter=2000).fit(Xs, y)
    return clf.coef_[0]


def bootstrap_ci(X, y, rng, b_boot=B_BOOT):
    """扰动级 bootstrap：有放回重采样 → 重拟合 → 逐系数 95% percentile CI。"""
    n = len(y)
    coefs = np.empty((b_boot, X.shape[1]))
    for b in range(b_boot):
        idx = rng.integers(0, n, n)
        # 极端重采样可能单类（n 小概率非零），跳过并以新重采样替代
        while len(np.unique(y[idx])) < 2:
            idx = rng.integers(0, n, n)
        coefs[b] = fit_coef(X[idx], y[idx])
    lo = np.percentile(coefs, 2.5, axis=0)
    hi = np.percentile(coefs, 97.5, axis=0)
    return coefs, lo, hi


def main():
    with open(os.path.join(ANALYSIS_DIR, "error_signatures.json"),
              encoding="utf-8") as f:
        sigs = json.load(f)

    rng_master = np.random.default_rng(BOOT_MASTER_SEED)
    res = {"task": "bootstrap confidence intervals for error signatures "
                   "(step8b)",
           "method": "perturbation-level nonparametric bootstrap, B="
                     f"{B_BOOT}, refit same-spec in-situ probe "
                     "(StandardScaler + LogisticRegression C=1.0) on each "
                     "resample; 95% percentile CI per coefficient",
           "boot_master_seed": BOOT_MASTER_SEED,
           "label": "y = (1 - pearson_delta) > median (worst half), same as "
                    "step6 / step6f",
           "signatures": {}}

    for sig_name in SIG_ORDER:
        X, y, feature_names = load_xy(sig_name)
        n = len(y)

        # --- 点估计核对：全量重拟合必须与归档签名一致 ---
        point = fit_coef(X, y)
        archived = np.array([sigs[sig_name][f] for f in FEATURE_ORDER])
        max_dev = float(np.max(np.abs(point - archived)))
        assert max_dev < 1e-6, \
            f"{sig_name} 全量重拟合与归档签名不一致 (max dev {max_dev:.2e})"
        print(f"[{sig_name}] n={n}，点估计核对 max|dev|={max_dev:.1e} OK")

        rng = np.random.default_rng(rng_master.integers(0, 2**31))
        coefs, lo, hi = bootstrap_ci(X, y, rng)
        feats = {}
        n_cross = 0
        for j, feat in enumerate(FEATURE_ORDER):
            crosses = bool(lo[j] <= 0 <= hi[j])
            n_cross += crosses
            feats[feat] = {"point": float(point[j]),
                           "ci_lo": float(lo[j]), "ci_hi": float(hi[j]),
                           "ci_crosses_zero": crosses}
        res["signatures"][sig_name] = {
            "label": SIG_LABELS[sig_name], "n": int(n),
            "n_coefs_ci_crossing_zero": int(n_cross),
            "features": feats}
        print(f"[{sig_name}] CI 跨零系数 {n_cross}/9")

    os.makedirs(ANALYSIS_DIR, exist_ok=True)
    out_json = os.path.join(ANALYSIS_DIR, "signature_ci.json")
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=2)
    print("[SAVED]", out_json)
    make_figure(res, os.path.join(ANALYSIS_DIR, "signature_ci.png"))


def make_figure(res, path):
    n_sig = len(SIG_ORDER)
    fig, axes = plt.subplots(1, n_sig, figsize=(3.4 * n_sig, 4.4),
                             sharey=True)
    colors = {"GEARS(K562)": "#1f77b4", "scgenept_k562": "#d62728",
              "scgenept_rpe1": "#ff7f0e", "scf_k562": "#2ca02c",
              "scgpt_k562": "#9467bd"}
    for ax, sig_name in zip(axes, SIG_ORDER):
        sres = res["signatures"][sig_name]
        c = colors[sig_name]
        for j, feat in enumerate(FEATURE_ORDER):
            fr = sres["features"][feat]
            ypos = len(FEATURE_ORDER) - 1 - j
            ax.plot([fr["ci_lo"], fr["ci_hi"]], [ypos, ypos],
                    color=c, lw=1.4, solid_capstyle="round")
            ax.plot(fr["point"], ypos, "o", color=c, ms=4.5,
                    mec="black", mew=0.5)
        ax.axvline(0.0, color="gray", ls=":", lw=1.0)
        ax.set_title(f"{sres['label']}\n(n={sres['n']}, "
                     f"{sres['n_coefs_ci_crossing_zero']}/9 CI crosses 0)",
                     fontsize=9)
        ax.set_xlabel("standardized coefficient", fontsize=8)
        ax.grid(axis="x", alpha=0.25)
    axes[0].set_yticks(range(len(FEATURE_ORDER)))
    axes[0].set_yticklabels(FEATURE_ORDER[::-1], fontsize=8)
    fig.suptitle("Error signature coefficients with 95% bootstrap CIs "
                 f"(B={B_BOOT}, perturbation-level resampling)",
                 fontsize=11)
    fig.text(0.5, 0.005, "point = in-situ full-data fit; whiskers = 2.5-97.5 "
                          "percentile of bootstrap refits",
             ha="center", fontsize=7, color="gray")
    fig.tight_layout(rect=[0, 0.02, 1, 0.96])
    fig.savefig(path, dpi=300)
    print("[SAVED]", path)


if __name__ == "__main__":
    main()
