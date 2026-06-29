"""
第三优先级：外周-中枢配体-受体对接分析
发送者：GSE225948 PB 单核细胞（D02）
接收者：GSE233812 脑内小胶质细胞（D3）
重点配体-受体对：CCL2-CCR2、TNF-TNFR1、IL1B-IL1RAP、CXCL10-CXCR3、CSF1-CSF1R 等
"""
import os
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import seaborn as sns

OUT_DIR = "../04_reports/figures/peripheral_central_lr"
os.makedirs(OUT_DIR, exist_ok=True)

# 手动 curated 配体-受体对（小鼠基因名）
LR_PAIRS = [
    ("Ccl2", "Ccr2"),
    ("Ccl3", "Ccr1"),
    ("Ccl4", "Ccr5"),
    ("Ccl5", "Ccr5"),
    ("Ccl7", "Ccr2"),
    ("Cxcl10", "Cxcr3"),
    ("Cxcl9", "Cxcr3"),
    ("Tnf", "Tnfrsf1a"),
    ("Tnf", "Tnfrsf1b"),
    ("Il1b", "Il1r1"),
    ("Il1b", "Il1rap"),
    ("Il6", "Il6ra"),
    ("Il6", "Il6st"),
    ("Csf1", "Csf1r"),
    ("Csf2", "Csf2ra"),
    ("Ifng", "Ifngr1"),
    ("Lta", "Tnfrsf1a"),
    ("Ltb", "Ltbr"),
    ("Mif", "Cd74"),
    ("App", "Cd74"),
    ("Spp1", "Cd44"),
    ("Spp1", "Itgav"),
]


def load_sender_receiver():
    """加载发送者（PB Mo D02）和接收者（Brain Mg D3）。"""
    # 发送者
    ad225 = sc.read_h5ad('04_reports/figures/GSE225948_processed.h5ad')
    ad225.obs_names_make_unique()
    sender = ad225[
        (ad225.obs['tissue'] == 'PB') &
        (ad225.obs['time_point'] == 'D02') &
        (ad225.obs['parent_celltype'] == 'Mo')
    ].copy()
    print(f"Sender PB Mo D02: {sender.shape}")

    # 接收者
    ad812 = sc.read_h5ad('04_reports/figures/GSE233812_processed.h5ad')
    ad812.obs_names_make_unique()
    receiver = ad812[
        (ad812.obs['cell_type'] == 'Microglia') &
        (ad812.obs['time_point'] == 'D3')
    ].copy()
    print(f"Receiver Brain Mg D3: {receiver.shape}")

    return sender, receiver


def compute_lr_score(sender, receiver, ligand, receptor):
    """计算配体-受体对接分数：发送者配体平均表达 × 接收者受体平均表达。"""
    lig_expr = 0.0
    rec_expr = 0.0
    if ligand in sender.var_names:
        lig_expr = sender[:, ligand].X.toarray().flatten().mean()
    if receptor in receiver.var_names:
        rec_expr = receiver[:, receptor].X.toarray().flatten().mean()
    score = lig_expr * rec_expr
    return lig_expr, rec_expr, score


def compute_pair_fraction(sender, receiver, ligand, receptor):
    """计算配体+受体双阳细胞比例。"""
    lig_pos = 0
    rec_pos = 0
    if ligand in sender.var_names:
        lig_pos = (sender[:, ligand].X.toarray().flatten() > 0).mean()
    if receptor in receiver.var_names:
        rec_pos = (receiver[:, receptor].X.toarray().flatten() > 0).mean()
    return lig_pos, rec_pos


def main():
    sender, receiver = load_sender_receiver()

    # 计算所有 LR 对的分数
    records = []
    for lig, rec in LR_PAIRS:
        lig_expr, rec_expr, score = compute_lr_score(sender, receiver, lig, rec)
        lig_pos, rec_pos = compute_pair_fraction(sender, receiver, lig, rec)
        records.append({
            "ligand": lig,
            "receptor": rec,
            "ligand_expr_sender": lig_expr,
            "receptor_expr_receiver": rec_expr,
            "lr_score": score,
            "ligand_pos_fraction": lig_pos,
            "receptor_pos_fraction": rec_pos,
        })
    lr_df = pd.DataFrame(records).sort_values("lr_score", ascending=False)
    lr_df.to_csv(os.path.join(OUT_DIR, "peripheral_PB_Mo_D02_to_brain_Mg_D3_lr_scores.csv"), index=False)
    print("\n=== Top LR pairs (PB Mo D02 → Brain Mg D3) ===")
    print(lr_df.head(20).to_string())

    # 保存 Pirb 在接收者中的表达作为参考
    pirb_expr_receiver = receiver[:, "Pirb"].X.toarray().flatten()
    pirb_pos_fraction = (pirb_expr_receiver > 0).mean()
    print(f"\nPirb+ fraction in Brain Mg D3: {pirb_pos_fraction:.4f}")

    # 可视化 top LR pairs
    top_n = 15
    top = lr_df.head(top_n).copy()
    top["pair"] = top["ligand"] + " → " + top["receptor"]
    plt.figure(figsize=(10, 8))
    sns.barplot(data=top, y="pair", x="lr_score", palette="viridis")
    plt.title("Top peripheral-to-central ligand-receptor pairs\n(PB monocytes D02 → Brain microglia D3)")
    plt.xlabel("LR score (ligand expr × receptor expr)")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "top_lr_pairs.png"), dpi=300)
    plt.close()

    # 热图：发送者配体 vs 接收者受体表达
    plot_df = lr_df.head(top_n).copy()
    plt.figure(figsize=(8, 8))
    heatmap_data = plot_df[["ligand_expr_sender", "receptor_expr_receiver"]].T
    heatmap_data.columns = plot_df["ligand"] + "→" + plot_df["receptor"]
    sns.heatmap(heatmap_data, annot=True, fmt=".3f", cmap="YlOrRd")
    plt.title("Ligand expression in sender vs receptor expression in receiver")
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, "lr_heatmap.png"), dpi=300)
    plt.close()

    print(f"\n[SAVED] all results to {OUT_DIR}")


if __name__ == "__main__":
    main()
