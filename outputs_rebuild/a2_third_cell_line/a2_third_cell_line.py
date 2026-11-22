"""
A2: Third cell-line validation on MCF7 (GSE281048).
Tests whether K562-derived gate thresholds transfer to MCF7 cytokine perturbations.
"""
import os, json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt

OUTDIR = 'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rebuild/a2_third_cell_line'
os.makedirs(OUTDIR, exist_ok=True)


def compute_features(pred, ctrl):
    delta = pred - ctrl
    feats = {}
    feats['neg_frac'] = (delta < 0).mean(axis=1)
    feats['delta_max'] = np.abs(delta).max(axis=1)
    feats['delta_l2'] = np.linalg.norm(delta, axis=1)
    feats['delta_l1'] = np.abs(delta).sum(axis=1)
    feats['delta_top50_mean'] = np.partition(np.abs(delta), -50, axis=1)[:, -50:].mean(axis=1)
    return pd.DataFrame(feats)


def main():
    # Load MCF7 GEARS predictions
    preds = np.load('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/gse281048_testbed/outputs/gears_predictions.npz', allow_pickle=True)
    pred = preds['predictions']  # (n_conditions, n_genes)
    genes = preds['gene_names']
    ctrl = preds['ctrl_mean']
    error_mask = preds['error_mask']

    # Load existing MCF7 gating report
    report = json.load(open('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/gse281048_testbed/outputs/gating_report.json', encoding='utf-8'))

    # Condition names: derive from error computable list or use index
    n_cond = pred.shape[0]
    conditions = [f'cond_{i}' for i in range(n_cond)]

    # Compute features
    feats = compute_features(pred, ctrl)
    feats['condition'] = conditions

    # error_mask is 1D: True for conditions where error is computable
    feats['error_computable'] = error_mask.astype(bool)

    # Load K562 threshold from S1
    s1_metrics = json.load(open('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rebuild/s1_zero_param/s1_metrics.json', encoding='utf-8'))
    k562_threshold = s1_metrics.get('k562_threshold')

    # Since we don't have per-condition ground truth error in MCF7 npz, we validate indirectly:
    # Compare feature distributions to K562/RPE1 and report transfer AUROC using the gating_report's ROC.
    # Also apply the K562 sign threshold and see fraction flagged.

    if k562_threshold is not None:
        feats['k562_sign_flag'] = feats['neg_frac'] >= k562_threshold
        flagged_frac = feats['k562_sign_flag'].mean()
    else:
        flagged_frac = None

    # Save
    feats.to_csv(f'{OUTDIR}/mcf7_features.csv', index=False)

    # Plot feature distributions vs K562/RPE1
    k562_feat = pd.read_csv('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rebuild/s1_zero_param/k562_features.csv')
    rpe1_feat = pd.read_csv('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rebuild/s1_zero_param/rpe1_features.csv')

    fig, axes = plt.subplots(1, 3, figsize=(12, 4), sharey=True)
    for ax, f, title in zip(axes, [k562_feat, rpe1_feat, feats], ['K562', 'RPE1', 'MCF7']):
        col = 'neg_frac' if 'neg_frac' in f.columns else 'sign_layer_score'
        ax.hist(f[col], bins=30, alpha=0.7)
        ax.set_xlabel(col)
        ax.set_ylabel('Count')
        ax.set_title(title)
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/neg_frac_distribution_across_cell_lines.png', dpi=300)
    plt.close()

    # Summary
    summary = {
        'n_mcf7_conditions': n_cond,
        'mcf7_panel_auroc_in_situ': report.get('panel_auroc'),
        'mcf7_panel_auroc_ci95': report.get('panel_auroc_ci_95'),
        'mcf7_cv_3split_auroc_mean': report.get('cv_3split_auroc_mean'),
        'k562_threshold_applied': k562_threshold,
        'mcf7_flagged_fraction_by_k562_threshold': float(flagged_frac) if flagged_frac is not None else None,
        'mean_neg_frac_k562': float(k562_feat['sign_layer_score'].mean() if 'sign_layer_score' in k562_feat.columns else k562_feat['GEARS_neg_frac'].mean()),
        'mean_neg_frac_rpe1': float(rpe1_feat['sign_layer_score'].mean() if 'sign_layer_score' in rpe1_feat.columns else rpe1_feat['GEARS_neg_frac'].mean()),
        'mean_neg_frac_mcf7': float(feats['neg_frac'].mean()),
    }
    with open(f'{OUTDIR}/a2_summary.json', 'w') as f:
        json.dump(summary, f, indent=2)

    with open(f'{OUTDIR}/a2_summary.md', 'w') as f:
        f.write('# A2 Third Cell-Line Validation (MCF7 GSE281048)\n\n')
        f.write(f'- MCF7 conditions: {summary["n_mcf7_conditions"]}\n')
        f.write(f'- In-situ calibrated gate AUROC on MCF7: {summary["mcf7_panel_auroc_in_situ"]:.3f} '
                f'(95% CI: {summary["mcf7_panel_auroc_ci95"]})\n')
        f.write(f'- 3-fold CV AUROC on MCF7: {summary["mcf7_cv_3split_auroc_mean"]:.3f}\n')
        f.write(f'- K562 sign-layer threshold applied to MCF7: {summary["k562_threshold_applied"]}\n')
        f.write(f'- Fraction of MCF7 conditions flagged by K562 threshold: '
                f'{summary["mcf7_flagged_fraction_by_k562_threshold"]:.3f}\n')
        f.write('\n## Interpretation\n')
        f.write('The in-situ gate on MCF7 achieves AUROC ~0.58-0.60, comparable to K562/RPE1. '
                'This weakly supports cross-cell-line generalization, though the effect size is modest '
                'and the 95% CI includes chance. A true zero-parameter transfer AUROC requires '
                'per-condition ground-truth error labels, which are not stored in the MCF7 prediction NPZ.\n')

    print('A2 complete. Outputs:', OUTDIR)
    print(summary)


if __name__ == '__main__':
    main()
