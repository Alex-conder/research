"""
TTA: Test-Time Adaptation for cross-cell-line gating.
Use control profiles to standardize features per cell line, then transfer K562 gate.
"""
import os, json
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt

OUTDIR = 'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rebuild/tta_mcf7'
os.makedirs(OUTDIR, exist_ok=True)


def compute_features(pred, ctrl):
    delta = pred - ctrl
    feats = pd.DataFrame({
        'neg_frac': (delta < 0).mean(axis=1),
        'delta_max': np.abs(delta).max(axis=1),
        'delta_l2': np.linalg.norm(delta, axis=1),
        'delta_l1': np.abs(delta).sum(axis=1),
        'delta_top50_mean': np.partition(np.abs(delta), -50, axis=1)[:, -50:].mean(axis=1),
        'delta_entropy': -np.sum(np.abs(delta) / (np.abs(delta).sum(axis=1, keepdims=True) + 1e-9) *
                                np.log(np.abs(delta) / (np.abs(delta).sum(axis=1, keepdims=True) + 1e-9) + 1e-9), axis=1),
    })
    return feats


def main():
    # ---- K562 source domain ----
    k562_gt, k562_ctrl, k562_gene_symbols = load_h5ad(
        'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/data/gears_demo/replogle_k562_essential/perturb_processed.h5ad'
    )
    k562_perts, k562_pred, k562_genes = load_gears_npz(
        'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs/pred_expr_seed0.npz'
    )
    # Align K562 predictions to h5ad gene symbols
    k562_true, k562_ctrl_aligned = align_to_gears(k562_gt, k562_ctrl, k562_gene_symbols, k562_genes)
    k562_pred_aligned = k562_pred  # already in gears gene space
    k562_feats = compute_features(k562_pred_aligned, k562_ctrl_aligned)
    k562_error = compute_error(k562_pred_aligned, k562_true, k562_ctrl_aligned)
    k562_y = k562_error > np.median(k562_error)

    # Compute K562 control-derived scaler: use control mean as a single pseudo-sample
    ctrl_pred = np.tile(k562_ctrl_aligned, (1, 1))  # (1, n_genes)
    ctrl_feats = compute_features(ctrl_pred, k562_ctrl_aligned)
    # For meaningful scaler, add small Gaussian noise around control to simulate control cells
    np.random.seed(42)
    ctrl_samples = k562_ctrl_aligned + np.random.normal(0, 0.05, size=(50, len(k562_ctrl_aligned)))
    ctrl_feats_samples = compute_features(ctrl_samples, k562_ctrl_aligned)
    scaler_k562 = StandardScaler().fit(ctrl_feats_samples)

    # ---- MCF7 target domain ----
    mcf7_preds = np.load('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/gse281048_testbed/outputs/gears_predictions.npz', allow_pickle=True)
    mcf7_pred = mcf7_preds['predictions']
    mcf7_genes = mcf7_preds['gene_names']
    mcf7_ctrl = mcf7_preds['ctrl_mean']
    mcf7_gene_list = mcf7_preds['genes']  # 189 perturbed genes

    # Align MCF7 to K562 GEARS gene space
    common_genes = sorted(set(k562_genes) & set(mcf7_genes))
    k_idx = [list(k562_genes).index(g) for g in common_genes]
    m_idx = [list(mcf7_genes).index(g) for g in common_genes]
    mcf7_pred_c = mcf7_pred[:, m_idx]
    mcf7_ctrl_c = mcf7_ctrl[m_idx]
    k562_ctrl_c = k562_ctrl_aligned[k_idx]

    # MCF7 control-derived scaler
    np.random.seed(43)
    mcf7_ctrl_samples = mcf7_ctrl_c + np.random.normal(0, 0.05, size=(50, len(mcf7_ctrl_c)))
    mcf7_ctrl_feats = compute_features(mcf7_ctrl_samples, mcf7_ctrl_c)
    scaler_mcf7 = StandardScaler().fit(mcf7_ctrl_feats)

    mcf7_feats = compute_features(mcf7_pred_c, mcf7_ctrl_c)

    # ---- TTA: standardize MCF7 using MCF7 control scaler ----
    mcf7_feats_tta = pd.DataFrame(scaler_mcf7.transform(mcf7_feats), columns=mcf7_feats.columns)
    k562_feats_tta = pd.DataFrame(scaler_k562.transform(k562_feats), columns=k562_feats.columns)

    # ---- Train gate on K562 TTA ----
    # Use sign-layer (neg_frac) and magnitude (delta_max)
    best_auroc = 0
    best_t = None
    best_feat = None
    for feat_name in ['neg_frac', 'delta_max', 'delta_l2']:
        score = k562_feats_tta[feat_name].values
        for t in np.linspace(score.min(), score.max(), 100):
            pred_y = score >= t
            if len(np.unique(pred_y)) < 2:
                continue
            auroc = roc_auc_score(k562_y, pred_y.astype(int))
            # also try reverse
            auroc_rev = roc_auc_score(k562_y, (~pred_y).astype(int))
            a = max(auroc, auroc_rev)
            if a > best_auroc:
                best_auroc = a
                best_t = t
                best_feat = feat_name
                best_direction = 1 if auroc >= auroc_rev else -1

    print(f'Best K562 TTA gate: {best_feat}, threshold={best_t:.3f}, direction={best_direction}, AUROC={best_auroc:.3f}')

    # ---- Evaluate on MCF7 ----
    # Need MCF7 error labels from panel_candidates.csv
    pan = pd.read_csv('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/gse281048_testbed/outputs/panel_candidates.csv')
    gene_to_error = dict(zip(pan['gene'], pan['error']))
    gene_to_pearson = dict(zip(pan['gene'], pan['median_pearson_r']))

    mcf7_error = np.array([gene_to_error.get(g, np.nan) for g in mcf7_gene_list])
    mcf7_pearson = np.array([gene_to_pearson.get(g, np.nan) for g in mcf7_gene_list])
    valid = ~np.isnan(mcf7_error)

    mcf7_y = mcf7_error[valid] > np.median(mcf7_error[valid])
    mcf7_score = mcf7_feats_tta[best_feat].values[valid]
    if best_direction == -1:
        mcf7_score = -mcf7_score
    mcf7_pred = (mcf7_score >= best_t).astype(int)
    mcf7_auroc = roc_auc_score(mcf7_y, mcf7_pred) if len(np.unique(mcf7_pred)) > 1 else 0.5

    # Baseline: no TTA (raw features, K562 threshold)
    best_t_raw = None
    best_auroc_raw = 0
    score_raw = k562_feats[best_feat].values
    for t in np.linspace(score_raw.min(), score_raw.max(), 100):
        pred_y = score_raw >= t
        if len(np.unique(pred_y)) < 2:
            continue
        auroc = roc_auc_score(k562_y, pred_y.astype(int))
        auroc_rev = roc_auc_score(k562_y, (~pred_y).astype(int))
        a = max(auroc, auroc_rev)
        if a > best_auroc_raw:
            best_auroc_raw = a
            best_t_raw = t
            dir_raw = 1 if auroc >= auroc_rev else -1

    mcf7_score_raw = mcf7_feats[best_feat].values[valid]
    if dir_raw == -1:
        mcf7_score_raw = -mcf7_score_raw
    mcf7_pred_raw = (mcf7_score_raw >= best_t_raw).astype(int)
    mcf7_auroc_raw = roc_auc_score(mcf7_y, mcf7_pred_raw) if len(np.unique(mcf7_pred_raw)) > 1 else 0.5

    # In-situ MCF7 calibration (cheating baseline)
    best_auroc_insitu = 0
    for t in np.linspace(mcf7_score.min(), mcf7_score.max(), 100):
        pred_y = mcf7_score >= t
        if len(np.unique(pred_y)) < 2:
            continue
        auroc = roc_auc_score(mcf7_y, pred_y.astype(int))
        auroc_rev = roc_auc_score(mcf7_y, (~pred_y).astype(int))
        a = max(auroc, auroc_rev)
        if a > best_auroc_insitu:
            best_auroc_insitu = a

    print(f'MCF7 AUROC: TTA={mcf7_auroc:.3f}, raw transfer={mcf7_auroc_raw:.3f}, in-situ={best_auroc_insitu:.3f}')

    # Save
    results = {
        'best_k562_feature': best_feat,
        'best_k562_threshold': float(best_t),
        'best_k562_direction': best_direction,
        'k562_tta_auroc': float(best_auroc),
        'mcf7_tta_auroc': float(mcf7_auroc),
        'mcf7_raw_transfer_auroc': float(mcf7_auroc_raw),
        'mcf7_insitu_auroc': float(best_auroc_insitu),
        'n_mcf7_valid': int(valid.sum()),
    }
    with open(f'{OUTDIR}/tta_results.json', 'w') as f:
        json.dump(results, f, indent=2)

    # Plot
    fig, ax = plt.subplots(figsize=(6, 4))
    methods = ['raw transfer', 'TTA', 'in-situ']
    vals = [mcf7_auroc_raw, mcf7_auroc, best_auroc_insitu]
    colors = ['red', 'green', 'blue']
    ax.bar(methods, vals, color=colors)
    ax.axhline(0.5, color='k', linestyle='--')
    ax.set_ylabel('MCF7 AUROC')
    ax.set_ylim(0.4, 0.8)
    for i, v in enumerate(vals):
        ax.text(i, v + 0.01, f'{v:.3f}', ha='center')
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/tta_mcf7_auroc.png', dpi=300)
    plt.close()

    # Feature distribution before/after TTA
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    for ax, df, title in zip(axes, [k562_feats, mcf7_feats, mcf7_feats_tta], ['K562 raw', 'MCF7 raw', 'MCF7 TTA']):
        ax.hist(df[best_feat], bins=30, alpha=0.7)
        ax.set_xlabel(best_feat)
        ax.set_ylabel('Count')
        ax.set_title(title)
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/feature_distribution_tta.png', dpi=300)
    plt.close()

    with open(f'{OUTDIR}/tta_summary.md', 'w') as f:
        f.write('# TTA for Cross-Cell-Line Gating\n\n')
        f.write(f'- Best K562 feature: {best_feat}, threshold={best_t:.3f}, AUROC={best_auroc:.3f}\n')
        f.write(f'- MCF7 raw transfer AUROC: {mcf7_auroc_raw:.3f}\n')
        f.write(f'- MCF7 TTA AUROC: {mcf7_auroc:.3f}\n')
        f.write(f'- MCF7 in-situ AUROC: {best_auroc_insitu:.3f}\n')
        f.write('\n## Narrative\n')
        if mcf7_auroc > mcf7_auroc_raw + 0.05:
            f.write('Control-based TTA substantially recovers cross-cell-line gate performance.\n')
        else:
            f.write('TTA provides only marginal improvement over raw transfer; MCF7 cytokine perturbations remain structurally different from CRISPRi data.\n')

    print('TTA complete. Outputs:', OUTDIR)


def load_h5ad(path):
    import scanpy as sc
    ad = sc.read_h5ad(path)
    perts = ad.obs['condition'].astype(str).values
    pert_names = np.array([p.split('+')[0] if '+' in p else p for p in perts])
    genes = np.array(ad.var['gene_name'].astype(str))
    expr = ad.X.toarray() if hasattr(ad.X, 'toarray') else np.array(ad.X)
    gt = {p: expr[pert_names == p].mean(axis=0) for p in np.unique(pert_names)}
    ctrl = expr[perts == 'ctrl'].mean(axis=0)
    return gt, ctrl, genes


def load_gears_npz(path):
    d = np.load(path, allow_pickle=True)
    return d['pred_perturbations'], d['pred_expr'], d['pred_genes']


def align_to_gears(gt, ctrl, gt_genes, gears_genes):
    gt_idx = {g: i for i, g in enumerate(gt_genes)}
    true = np.array([gt[p][gt_idx[g]] if g in gt_idx else np.nan for p in sorted(gt.keys()) for g in gears_genes])
    true = true.reshape(len(gt), len(gears_genes))
    ctrl_aligned = np.array([ctrl[gt_idx[g]] if g in gt_idx else np.nan for g in gears_genes])
    return true, ctrl_aligned


def compute_error(pred, true, ctrl):
    delta_pred = pred - ctrl
    delta_true = true - ctrl
    corrs = []
    for i in range(pred.shape[0]):
        top = np.argsort(np.abs(delta_true[i]))[-50:]
        r = np.corrcoef(delta_pred[i, top], delta_true[i, top])[0, 1]
        corrs.append(0 if np.isnan(r) else r)
    return 1 - np.array(corrs)


if __name__ == '__main__':
    main()
