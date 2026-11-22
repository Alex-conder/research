# step8d_signature_space.py — 误差签名空间的 PCA 示意（2D 散点）
#
# 数据：outputs_crossmodel/analysis/error_signatures.json（5 个误差签名，
# 每个 = 9 维标准化系数向量）。流程：
#   1. 5×9 签名矩阵按列（特征）标准化（z-score across the 5 signatures）；
#   2. PCA(n_components=2) 降至 2D 画散点，点旁标注模型名；
#   3. 两两余弦相似度从 error_signatures.json 复算（不照抄任何归档数字），
#      与已流传数值比对留痕（assert 容差 0.01），图上标注关键几对；
#   4. 如实报告前两维方差解释率。
#
# 注意：5 个点在 9 维空间中的 2D 投影只是**示意**（第 1 主成分最多解释
# 5 个点的 4 个自由度），图注必须写明 "illustrative, 5 points in
# 9-dim space"。
#
# 产物：outputs_crossmodel/analysis/signature_space.png（300 dpi，英文）
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step8d_signature_space.py
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA

REPO = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.join(REPO, "outputs_crossmodel", "analysis")
SRC = os.path.join(ANALYSIS_DIR, "error_signatures.json")
OUT = os.path.join(ANALYSIS_DIR, "signature_space.png")

SIG_ORDER = ["GEARS(K562)", "scgenept_k562", "scgenept_rpe1", "scf_k562",
             "scgpt_k562"]
SIG_LABELS = {"GEARS(K562)": "GEARS / K562",
              "scgenept_k562": "scGenePT / K562",
              "scgenept_rpe1": "scGenePT / RPE1",
              "scf_k562": "scFoundation / K562",
              "scgpt_k562": "scGPT / K562"}
FEATURE_ORDER = ["delta_l2", "delta_l1", "delta_max", "delta_top50_mean",
                 "n_de_025", "frac_de_025", "neg_frac", "delta_entropy",
                 "log_n_cells"]
COLORS = {"GEARS(K562)": "#1f77b4", "scgenept_k562": "#d62728",
          "scgenept_rpe1": "#ff7f0e", "scf_k562": "#2ca02c",
          "scgpt_k562": "#9467bd"}

# 已流传的两两余弦值（仅作复算核对，不作数据源）
EXPECTED_COSINE = {
    ("GEARS(K562)", "scgenept_k562"): +0.58,
    ("GEARS(K562)", "scf_k562"): +0.07,
    ("scgenept_k562", "scf_k562"): -0.49,
    ("scgenept_k562", "scgenept_rpe1"): +0.01,
}


def cosine_table(M):
    """逐对余弦相似度（未标准化原始签名空间）。"""
    norms = np.linalg.norm(M, axis=1, keepdims=True)
    Mh = M / norms
    return Mh @ Mh.T


def main():
    with open(SRC, encoding="utf-8") as f:
        sigs = json.load(f)
    M = np.array([[sigs[s][feat] for feat in FEATURE_ORDER]
                  for s in SIG_ORDER])

    # --- 余弦相似度复算核对 ---
    cos = cosine_table(M)
    print("pairwise cosine similarity (recomputed from "
          "error_signatures.json):")
    checked = {}
    for i, a in enumerate(SIG_ORDER):
        for j, b in enumerate(SIG_ORDER):
            if j <= i:
                continue
            v = float(cos[i, j])
            print(f"  {a} <-> {b}: {v:+.3f}")
            for (ea, eb), exp in EXPECTED_COSINE.items():
                if {ea, eb} == {a, b}:
                    assert abs(v - exp) < 0.01, \
                        f"余弦核对失败 {a}/{b}: {v:.3f} vs {exp:+.2f}"
                    checked[(a, b)] = v
    assert len(checked) == len(EXPECTED_COSINE), "有核对项未匹配到签名对"

    # --- 按列标准化 → PCA ---
    Z = (M - M.mean(axis=0)) / M.std(axis=0, ddof=1)
    pca = PCA(n_components=2)
    P = pca.fit_transform(Z)
    evr = pca.explained_variance_ratio_
    print(f"PCA explained variance ratio: PC1={evr[0]:.3f}, "
          f"PC2={evr[1]:.3f}, cumulative={evr.sum():.3f}")

    # 产物侧记录（json 不落新文件，数值打印留痕；图内标注）
    fig, ax = plt.subplots(figsize=(6.4, 5.2))
    for i, s in enumerate(SIG_ORDER):
        ax.plot(P[i, 0], P[i, 1], "o", color=COLORS[s], ms=10,
                mec="black", mew=0.6)
        off = (8, 8) if s != "scgenept_rpe1" else (8, -14)
        ax.annotate(SIG_LABELS[s], (P[i, 0], P[i, 1]),
                    textcoords="offset points", xytext=off, fontsize=9)

    # 关键签名对的余弦标注（虚线连接 + 文本）
    ann = [(("GEARS(K562)", "scgenept_k562"), "top-right"),
           (("GEARS(K562)", "scf_k562"), "bottom"),
           (("scgenept_k562", "scf_k562"), "top-left")]
    for (a, b), _pos in ann:
        i, j = SIG_ORDER.index(a), SIG_ORDER.index(b)
        ax.plot([P[i, 0], P[j, 0]], [P[i, 1], P[j, 1]], color="gray",
                lw=0.7, ls=":", zorder=0)
        mx, my = (P[i, 0] + P[j, 0]) / 2, (P[i, 1] + P[j, 1]) / 2
        ax.text(mx, my, f"cos={checked.get((a, b), checked.get((b, a))):+.2f}",
                fontsize=7, color="dimgray", ha="center")

    ax.axhline(0.0, color="lightgray", lw=0.6)
    ax.axvline(0.0, color="lightgray", lw=0.6)
    ax.set_xlabel(f"PC1 ({evr[0] * 100:.1f}% of variance)")
    ax.set_ylabel(f"PC2 ({evr[1] * 100:.1f}% of variance)")
    ax.set_title("Error signatures in PCA space\n(illustrative, 5 points in "
                 "9-dim space)", fontsize=11)
    ax.grid(alpha=0.25)
    # 留边距防标注裁切
    xr = np.ptp(P[:, 0]) or 1.0
    yr = np.ptp(P[:, 1]) or 1.0
    ax.set_xlim(P[:, 0].min() - 0.25 * xr, P[:, 0].max() + 0.25 * xr)
    ax.set_ylim(P[:, 1].min() - 0.25 * yr, P[:, 1].max() + 0.25 * yr)
    fig.tight_layout()
    fig.savefig(OUT, dpi=300)
    print("[SAVED]", OUT)


if __name__ == "__main__":
    main()
