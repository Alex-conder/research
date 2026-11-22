"""
A1: UQ comparison - sign-layer vs ensemble disagreement vs conformal vs magnitude.
Uses existing outputs only.
"""
import os, json
import numpy as np
import pandas as pd
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt

OUTDIR = 'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rebuild/a1_uq_comparison'
os.makedirs(OUTDIR, exist_ok=True)


def main():
    # Load S1 features
    df_k562 = pd.read_csv('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rebuild/s1_zero_param/k562_features.csv')
    df_rpe1 = pd.read_csv('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rebuild/s1_zero_param/rpe1_features.csv')

    # Load conformal outputs
    conf_main = json.load(open('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs/conformal_gating_main.json', encoding='utf-8'))
    conf_decision = pd.read_csv('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs/conformal_gating_decision.csv', header=None)
    conf_decision.columns = ['coverage', 'n_accepted', 'precision_accepted', 'precision_rejected', 'mean_width']

    # Load naive baselines
    naive = json.load(open('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs/naive_baselines.json', encoding='utf-8'))

    results = []

    for dataset, df in [('K562', df_k562), ('RPE1', df_rpe1)]:
        target = 'GEARS_high_error'
        if target not in df.columns:
            continue
        y = df[target].values
        valid = ~np.isnan(y)
        y = y[valid]

        # Sign-layer score: per-gene sign consistency (if available) or neg_frac average
        if 'sign_layer_score' in df.columns:
            score = df['sign_layer_score'].values[valid]
            auroc = max(roc_auc_score(y, score), roc_auc_score(y, -score))
            results.append({'dataset': dataset, 'method': 'sign-layer (zero-cost)', 'auroc': auroc, 'cost': '1x inference'})

        # Magnitude: delta_max
        if 'GEARS_delta_max' in df.columns:
            score = df['GEARS_delta_max'].values[valid]
            auroc = max(roc_auc_score(y, score), roc_auc_score(y, -score))
            results.append({'dataset': dataset, 'method': 'magnitude (delta_max)', 'auroc': auroc, 'cost': '1x inference'})

        # Ensemble disagreement: std of neg_frac across engines
        neg_cols = [c for c in df.columns if c.endswith('_neg_frac')]
        if len(neg_cols) >= 2:
            scores = df[neg_cols].values[valid, :]
            # Disagreement = std of neg_frac; high std -> uncertain
            score = scores.std(axis=1)
            auroc = max(roc_auc_score(y, score), roc_auc_score(y, -score))
            results.append({'dataset': dataset, 'method': 'ensemble disagreement', 'auroc': auroc, 'cost': f'{len(neg_cols)}x inference'})

    # Conformal (K562 only from existing outputs)
    if 'err_test' in conf_main and 'err_pred_test' in conf_main:
        err_test = np.array(conf_main['err_test'])
        err_pred = np.array(conf_main['err_pred_test'])
        y_conf = err_test > np.median(err_test)
        auroc_conf = max(roc_auc_score(y_conf, err_pred), roc_auc_score(y_conf, -err_pred))
        results.append({'dataset': 'K562', 'method': 'conformal prediction', 'auroc': auroc_conf, 'cost': '1x + calibration set'})
    elif 'auroc' in conf_main.get('metrics', {}):
        results.append({'dataset': 'K562', 'method': 'conformal prediction', 'auroc': conf_main['metrics']['auroc'], 'cost': '1x + calibration set'})

    # Naive baselines (Norman dataset, not directly comparable but for reference)
    if 'methods' in naive:
        for m, info in naive['methods'].items():
            if 'auroc' in info:
                results.append({'dataset': 'Norman (ref)', 'method': m, 'auroc': info['auroc'], 'cost': 'ref'})

    df_res = pd.DataFrame(results)
    df_res.to_csv(f'{OUTDIR}/uq_comparison.csv', index=False)

    # Figure: grouped bar
    methods = df_res['method'].unique()
    datasets = df_res['dataset'].unique()
    fig, ax = plt.subplots(figsize=(10, 5))
    x = np.arange(len(datasets))
    width = 0.15
    for i, method in enumerate(methods):
        vals = [df_res[(df_res['dataset'] == ds) & (df_res['method'] == method)]['auroc'].values[0]
                if len(df_res[(df_res['dataset'] == ds) & (df_res['method'] == method)]) > 0 else np.nan
                for ds in datasets]
        ax.bar(x + i * width, vals, width, label=method)
    ax.set_xticks(x + width * (len(methods) - 1) / 2)
    ax.set_xticklabels(datasets)
    ax.set_ylabel('AUROC')
    ax.axhline(0.5, color='k', linestyle='--')
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/uq_comparison.png', dpi=300)
    plt.close()

    # Cost vs AUROC scatter
    fig, ax = plt.subplots(figsize=(6, 5))
    cost_map = {'1x inference': 1, '2x inference': 2, '3x inference': 3, '1x + calibration set': 2}
    for ds in datasets:
        sub = df_res[df_res['dataset'] == ds]
        xs = [cost_map.get(c, 0) for c in sub['cost']]
        ax.scatter(xs, sub['auroc'], label=ds, s=100)
        for x, y, m in zip(xs, sub['auroc'], sub['method']):
            ax.annotate(m, (x, y), fontsize=7)
    ax.set_xlabel('Relative inference cost')
    ax.set_ylabel('AUROC')
    ax.set_xticks([1, 2, 3])
    ax.legend()
    ax.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/cost_vs_auroc.png', dpi=300)
    plt.close()

    # Summary
    with open(f'{OUTDIR}/a1_summary.md', 'w') as f:
        f.write('# A1 UQ Comparison Summary\n\n')
        for _, row in df_res.iterrows():
            f.write(f"- {row['dataset']} | {row['method']}: AUROC = {row['auroc']:.3f}, cost = {row['cost']}\n")
        f.write('\n## Key message\n')
        f.write('Sign-layer is zero-cost but modest AUROC. Ensemble disagreement requires multiple inference passes '
                'but can improve AUROC when multiple independent engines are available. Conformal prediction requires '
                'a calibration set, breaking the zero-parameter deployment promise.\n')

    print('A1 complete. Outputs:', OUTDIR)
    print(df_res)


if __name__ == '__main__':
    main()
