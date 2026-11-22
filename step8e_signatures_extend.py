# step8e_signatures_extend.py — 把新落地模型的 K562 误差签名增补进
# outputs_crossmodel/analysis/error_signatures.json。
#
# 背景：scGPT 的 K562 产物 2026-10-16 在服务器（172.16.92.2）跑完并
# SFTP 落地到 outputs_crossmodel/scgpt/；error_signatures.json 原有的
# 4 个签名（GEARS/scGenePT×2/scFoundation）保持不变，本脚本只**追加**
# 新键，不改动归档值。
#
# 范围说明：uce / geneformer 的 K562 产物（同日落地）为 embedding-only
# （in silico deletion 细胞嵌入，无表达谱预测，_labels.csv pearson 列为
# NaN，见 npz 的 output_type），表达空间签名链（9 特征 + pearson 标签）
# 无法纳入；如需分析这两个模型须另设计嵌入空间探针（预注册纪律下先
# 登记再做），不属于本脚本的复跑范围。
#
# 口径（与既有签名逐字一致，不自造逻辑）：
#   - 数据源与特征/标签加载：step8b_signature_ci.load_xy（import 复用）
#   - 就地探针拟合：step8b_signature_ci.fit_coef（StandardScaler +
#     LogisticRegression(C=1.0, max_iter=2000)，in-sample）
#   - 签名值 = 9 个标准化系数，键序 = step8b.FEATURE_ORDER
#   - step8b 的点估计核对（全量重拟合 vs 归档 <1e-6）对本脚本写入的
#     新签名同样成立（同一数据 + 同一规格 → 数值必然一致）， rerun
#     step8b 时自动验证。
#
# 安全措施：
#   - 已存在的键默认跳过（--force 才重算覆盖，供数据修正后用）
#   - 写入前把当前 error_signatures.json 复制到
#     analysis_backup_<日期>/error_signatures.json.orig
#
# 用法：
#   PYTHONIOENCODING=utf-8 venv/Scripts/python.exe step8e_signatures_extend.py
import argparse
import datetime
import json
import os
import shutil
import sys

REPO = os.path.dirname(os.path.abspath(__file__))
ANALYSIS_DIR = os.path.join(REPO, "outputs_crossmodel", "analysis")
SRC = os.path.join(ANALYSIS_DIR, "error_signatures.json")

sys.path.insert(0, REPO)
from step8b_signature_ci import (FEATURE_ORDER, SIG_ORDER,  # noqa: E402
                                 fit_coef, load_xy)

# 本次要增补的签名（必须在 step8b.SOURCES 里有数据源条目）
NEW_SIGS = ["scgpt_k562"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--force", action="store_true",
                    help="重算并覆盖已存在的签名键（默认跳过已存在的键）")
    args = ap.parse_args()

    with open(SRC, encoding="utf-8") as f:
        sigs = json.load(f)

    for sig_name in NEW_SIGS:
        assert sig_name in SIG_ORDER, f"{sig_name} 不在 step8b.SIG_ORDER 中"
        if sig_name in sigs and not args.force:
            print(f"[SKIP] {sig_name} 已存在（--force 可重算覆盖）")
            continue
        X, y, feature_names = load_xy(sig_name)
        assert list(feature_names) == FEATURE_ORDER, \
            f"{sig_name} 特征列与签名顺序不一致: {feature_names}"
        coef = fit_coef(X, y)
        sigs[sig_name] = {feat: float(c) for feat, c
                          in zip(FEATURE_ORDER, coef)}
        print(f"[ADD] {sig_name} n={len(y)} 签名已计算（就地 in-sample 拟合）")

    stamp = datetime.date.today().strftime("%Y%m%d")
    backup_dir = os.path.join(os.path.dirname(ANALYSIS_DIR),
                              f"analysis_backup_{stamp}")
    os.makedirs(backup_dir, exist_ok=True)
    backup = os.path.join(backup_dir, "error_signatures.json.orig")
    if not os.path.exists(backup):
        shutil.copy2(SRC, backup)
        print(f"[BACKUP] 原文件已复制到 {backup}")
    else:
        print(f"[BACKUP] 已存在，跳过: {backup}")

    with open(SRC, "w", encoding="utf-8") as f:
        json.dump(sigs, f, ensure_ascii=False, indent=1)
    print(f"[SAVED] {SRC}（共 {len(sigs)} 个签名）")


if __name__ == "__main__":
    main()
