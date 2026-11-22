# step6g_signature_heatmap.py — 误差签名热力图（跨模型分析的核心增量图）
# 行 = 5 个 error signature（就地拟合探针的标准化系数，in-sample），
# 列 = 9 个探针特征；数值直接标注。
# 数据源：outputs_crossmodel/analysis/error_signatures.json
# 产出：outputs_crossmodel/analysis/signature_heatmap.png（300 dpi，英文）
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

REPO = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(REPO, "outputs_crossmodel", "analysis",
                   "error_signatures.json")
OUT = os.path.join(REPO, "outputs_crossmodel", "analysis",
                   "signature_heatmap.png")

SIG_ORDER = ["GEARS(K562)", "scgenept_k562", "scgenept_rpe1", "scf_k562",
             "scgpt_k562"]
SIG_LABELS = ["GEARS / K562 (this work)", "scGenePT / K562",
              "scGenePT / RPE1", "scFoundation / K562", "scGPT / K562"]
FEATURE_ORDER = ["delta_l2", "delta_l1", "delta_max", "delta_top50_mean",
                 "n_de_025", "frac_de_025", "neg_frac", "delta_entropy",
                 "log_n_cells"]

with open(SRC, encoding="utf-8") as f:
    sigs = json.load(f)
M = np.array([[sigs[s][feat] for feat in FEATURE_ORDER] for s in SIG_ORDER])

fig, ax = plt.subplots(figsize=(9.2, 5.2))
vmax = float(np.abs(M).max())
im = ax.imshow(M, cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="auto")
ax.set_xticks(range(len(FEATURE_ORDER)))
ax.set_xticklabels(FEATURE_ORDER, rotation=35, ha="right", fontsize=8)
ax.set_yticks(range(len(SIG_LABELS)))
ax.set_yticklabels(SIG_LABELS, fontsize=9)
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        ax.text(j, i, f"{M[i, j]:+.2f}", ha="center", va="center",
                fontsize=7.5,
                color="white" if abs(M[i, j]) > 0.6 * vmax else "black")
cb = fig.colorbar(im, ax=ax, shrink=0.85, pad=0.015)
cb.set_label("standardized in-situ probe coefficient", fontsize=8)
cb.ax.tick_params(labelsize=7)
ax.set_title("Error signatures across model$\\times$dataset pairs "
             "(in-sample coefficients, same 9 model-side features)",
             fontsize=10)
fig.tight_layout()
fig.savefig(OUT, dpi=300)
print("[SAVED]", OUT)
