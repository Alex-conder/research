"""
TTA v2: Unsupervised median-IQR alignment for cross-cell-line gating.
Align target domain feature distribution to source domain using all target samples.
"""
import os, json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
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
    })
    return feats


def robust_align(source, target):
    """Align target to source using median and IQR."""
    s_med = source.median()
    s_iqr = source.quantile(0.75) - source.quantile(0.25)
    t_med = target.median()
    t_iqr = target.quantile(0.75) - target.quantile(0.25)
    aligned = (target - t_med) / (t_iqr + 1e-9) * (s_iqr + 1e-9) + s_med
    return aligned


def main():
    import scanpy as sc

    # K562
    ad = sc.read_h5ad('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/data/gears_demo/replogle_k562_essential/perturb_processed.h5ad')
    perts = ad.obs['condition'].astype(str).values
    pert_names = np.array([p.split('+')[0] if '+' in p else p for p in perts])
    gt_genes = np.array(ad.var['gene_name'].astype(str))
    expr = ad.X.toarray() if hasattr(ad.X, 'toarray') else np.array(ad.X)
    gt = {p: expr[pert_names == p].mean(axis=0) for p in np.unique(pert_names)}
    ctrl = expr[perts == 'ctrl'].mean(axis=0)
    gt_idx = {g: i for i, g in enumerate(gt_genes)}

    gears = np.load('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs/pred_expr_seed0.npz', allow_pickle=True)
    gears_genes = gears['pred_genes']
    k562_pred = gears['pred_expr']
    k562_ctrl = np.array([ctrl[gt_idx[g]] if g in gt_idx else np.nan for g in gears_genes])
    k562_true = np.array([gt[p][gt_idx[g]] if g in gt_idx else np.nan for p in gears['pred_perturbations'] for g in gears_genes])
    k562_true = k562_true.reshape(len(gears['pred_perturbations']), len(gears_genes))

    k562_feats = compute_features(k562_pred, k562_ctrl)

    # K562 error
    delta_pred = k562_pred - k562_ctrl
    delta_true = k562_true - k562_ctrl
    k562_error = []
    for i in range(len(k562_pred)):
        top = np.argsort(np.abs(delta_true[i]))[-50:]
        r = np.corrcoef(delta_pred[i, top], delta_true[i, top])[0, 1]
        k562_error.append(0 if np.isnan(r) else 1 - r)
    k562_error = np.array(k562_error)
    k562_y = k562_error > np.median(k562_error)

    # MCF7
    mcf7 = np.load('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/gse281048_testbed/outputs/gears_predictions.npz', allow_pickle=True)
    mcf7_pred = mcf7['predictions']
    mcf7_genes = mcf7['gene_names']
    mcf7_ctrl = mcf7['ctrl_mean']
    mcf7_gene_list = mcf7['genes']

    # Align MCF7 to GEARS gene space
    common_genes = sorted(set(gears_genes) & set(mcf7_genes))
    k_idx = [list(gears_genes).index(g) for g in common_genes]
    m_idx = [list(mcf7_genes).index(g) for g in common_genes]
    mcf7_pred_c = mcf7_pred[:, m_idx]
    mcf7_ctrl_c = mcf7_ctrl[m_idx]

    mcf7_feats = compute_features(mcf7_pred_c, mcf7_ctrl_c)

    # MCF7 error labels
    pan = pd.read_csv('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/gse281048_testbed/outputs/panel_candidates.csv')
    gene_to_error = dict(zip(pan['gene'], pan['error']))
    mcf7_error = np.array([gene_to_error.get(g, np.nan) for g in mcf7_gene_list])
    valid = ~np.isnan(mcf7_error)
    mcf7_y = mcf7_error[valid] > np.median(mcf7_error[valid])

    # Align MCF7 features to K562 distribution
    mcf7_feats_aligned = robust_align(k562_feats, mcf7_feats)

    results = {}
    for feat in ['neg_frac', 'delta_max', 'delta_l2', 'delta_l1', 'delta_top50_mean']:
        # Source performance
        s_score = k562_feats[feat].values
        best_s_auroc = 0
        best_s_t = None
        best_s_dir = 1
        for t in np.linspace(s_score.min(), s_score.max(), 100):
            py = s_score >= t
            if len(np.unique(py)) < 2:
                continue
            a = roc_auc_score(k562_y, py.astype(int))
            ar = roc_auc_score(k562_y, (~py).astype(int))
            if max(a, ar) > best_s_auroc:
                best_s_auroc = max(a, ar)
                best_s_t = t
                best_s_dir = 1 if a >= ar else -1

        # Raw transfer to MCF7
        t_score = mcf7_feats[feat].values[valid]
        if best_s_dir == -1:
            t_score = -t_score
        py = (t_score >= best_s_t).astype(int)
        raw_auroc = roc_auc_score(mcf7_y, py) if len(np.unique(py)) > 1 else 0.5

        # Aligned transfer
        t_score_a = mcf7_feats_aligned[feat].values[valid]
        if best_s_dir == -1:
            t_score_a = -t_score_a
        py_a = (t_score_a >= best_s_t).astype(int)
        aligned_auroc = roc_auc_score(mcf7_y, py_a) if len(np.unique(py_a)) > 1 else 0.5

        # In-situ
        best_i_auroc = 0
        for t in np.linspace(t_score.min(), t_score.max(), 100):
            py_i = (t_score >= t).astype(int)
            if len(np.unique(py_i)) < 2:
                continue
            a = roc_auc_score(mcf7_y, py_i)
            ar = roc_auc_score(mcf7_y, (~py_i.astype(bool)).astype(int))
            if max(a, ar) > best_i_auroc:
                best_i_auroc = max(a, ar)

        results[feat] = {
            'k562_auroc': float(best_s_auroc),
            'mcf7_raw_auroc': float(raw_auroc),
            'mcf7_aligned_auroc': float(aligned_auroc),
            'mcf7_insitu_auroc': float(best_i_auroc),
            'threshold': float(best_s_t),
            'direction': best_s_dir,
        }

    df_res = pd.DataFrame(results).T
    df_res.to_csv(f'{OUTDIR}/tta_v2_results.csv')

    # Plot
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(df_res))
    width = 0.2
    ax.bar(x - width, df_res['k562_auroc'], width, label='K562')
    ax.bar(x, df_res['mcf7_raw_auroc'], width, label='MCF7 raw')
    ax.bar(x + width, df_res['mcf7_aligned_auroc'], width, label='MCF7 aligned')
    ax.bar(x + 2*width, df_res['mcf7_insitu_auroc'], width, label='MCF7 in-situ')
    ax.set_xticks(x + width/2)
    ax.set_xticklabels(df_res.index, rotation=30, ha='right')
    ax.set_ylabel('AUROC')
    ax.axhline(0.5, color='k', linestyle='--')
    ax.legend()
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/tta_v2_feature_comparison.png', dpi=300)
    plt.close()

    with open(f'{OUTDIR}/tta_v2_summary.md', 'w') as f:
        f.write('# TTA v2: Unsupervised Median-IQR Alignment\n\n')
        f.write(df_res.to_string())
        f.write('\n\n## Interpretation\n')
        best = df_res['mcf7_aligned_auroc'].idxmax()
        f.write(f'Best aligned feature: {best}, AUROC = {df_res.loc[best, "mcf7_aligned_auroc"]:.3f}\n')
        if df_res.loc[best, 'mcf7_aligned_auroc'] > df_res.loc[best, 'mcf7_raw_auroc'] + 0.05:
            f.write('Alignment provides meaningful improvement over raw transfer.\n')
        else:
            f.write('Alignment does not rescue performance; MCF7 remains structurally different.\n')

    print('TTA v2 complete.')
    print(df_res)


if __name__ == '__main__':
    main()
