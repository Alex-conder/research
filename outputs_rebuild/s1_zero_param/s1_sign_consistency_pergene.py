"""
S1 supplement: per-gene sign-consistency gate.
"""
import os, json, warnings
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore')
OUTDIR = 'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rebuild/s1_zero_param'


def load_gt(h5ad_path):
    ad = sc.read_h5ad(h5ad_path)
    perts = ad.obs['condition'].astype(str).values
    pert_names = np.array([p.split('+')[0] if '+' in p else p for p in perts])
    genes = np.array(ad.var['gene_name'].astype(str))
    expr = ad.X.toarray() if hasattr(ad.X, 'toarray') else np.array(ad.X)
    gt = {p: expr[pert_names == p].mean(axis=0) for p in np.unique(pert_names)}
    ctrl = expr[perts == 'ctrl'].mean(axis=0)
    return gt, ctrl, genes


def align_engine(perts, pred, genes, common_perts, common_genes):
    pert_idx = {p: i for i, p in enumerate(perts)}
    gene_idx = {g: i for i, g in enumerate(genes)}
    out = np.full((len(common_perts), len(common_genes)), np.nan)
    for i, p in enumerate(common_perts):
        if p in pert_idx:
            row = pred[pert_idx[p]]
            for j, g in enumerate(common_genes):
                if g in gene_idx:
                    out[i, j] = row[gene_idx[g]]
    return out


def main():
    # Ground truth
    k562_gt, k562_ctrl, k562_gene_symbols = load_gt(
        'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/data/gears_demo/replogle_k562_essential/perturb_processed.h5ad'
    )
    rpe1_gt, rpe1_ctrl, rpe1_gene_symbols = load_gt(
        'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/data/gears_demo/replogle_rpe1_essential/perturb_processed.h5ad'
    )

    # K562 engines
    gears = np.load('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs/pred_expr_seed0.npz', allow_pickle=True)
    scgpt = np.load('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_crossmodel/scgpt/scgpt_k562_predictions.npz', allow_pickle=True)
    scgenept = np.load('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_crossmodel/scgenept/scgenept_k562_predictions.npz', allow_pickle=True)

    # RPE1 engines
    gears_r = np.load('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rpe1/pred_expr_seed0.npz', allow_pickle=True)
    scgenept_r = np.load('E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_crossmodel/scgenept/scgenept_rpe1_predictions.npz', allow_pickle=True)

    def process(gt, ctrl, gene_symbols, engines, engine_names):
        common_perts = sorted(set.intersection(*[set(e['perts']) for e in engines]))
        # Use GEARS gene list as common gene space
        common_genes = np.array(engines[0]['genes'])
        gt_idx = {g: i for i, g in enumerate(gene_symbols)}
        true = np.array([[gt[p][gt_idx[g]] if g in gt_idx else np.nan for g in common_genes] for p in common_perts])
        ctrl_aligned = np.array([ctrl[gt_idx[g]] if g in gt_idx else np.nan for g in common_genes])

        preds = [align_engine(e['perts'], e['pred'], e['genes'], common_perts, common_genes) for e in engines]
        # Add ctrl means for non-GEARS engines
        for k, e in enumerate(engines):
            if 'ctrl' in e:
                ctrl_map = {g: e['ctrl'][i] for i, g in enumerate(e['genes'])}
                ctrl_aligned_engine = np.array([ctrl_map.get(g, np.nan) for g in common_genes])
                # Not used directly; we use GEARS ctrl for all as approximation to keep consistent

        # Per-gene sign consistency
        signs = np.stack([np.sign(p - ctrl_aligned) for p in preds], axis=0)  # (n_engines, n_perts, n_genes)
        sign_sum = signs.sum(axis=0)
        consistency = np.abs(sign_sum) / len(engines)  # 1 = all agree, 0 = split

        # Perturbation-level uncertainty scores (higher = more uncertain)
        scores = {
            'mean_consistency': consistency.mean(axis=1),
            'frac_split': (consistency < 1.0).mean(axis=1),
            'entropy_sign': -(consistency * np.log2(consistency + 1e-9) +
                              (1 - consistency) * np.log2(1 - consistency + 1e-9)).mean(axis=1),
        }

        # Error: GEARS top-50 Pearson below median
        pred_delta = preds[0] - ctrl_aligned
        true_delta = true - ctrl_aligned
        corrs = []
        for i in range(len(common_perts)):
            top = np.argsort(np.abs(true_delta[i]))[-50:]
            r = np.corrcoef(pred_delta[i, top], true_delta[i, top])[0, 1]
            corrs.append(0 if np.isnan(r) else r)
        corrs = np.array(corrs)
        high_error = corrs < np.median(corrs)

        return scores, high_error, corrs, common_perts

    k_scores, k_err, k_corrs, k_perts = process(k562_gt, k562_ctrl, k562_gene_symbols,
                                                  [{'perts': gears['pred_perturbations'], 'pred': gears['pred_expr'], 'genes': gears['pred_genes']},
                                                   {'perts': scgpt['pred_perturbations'], 'pred': scgpt['pred_expr'], 'genes': scgpt['pred_genes'], 'ctrl': scgpt['ctrl_mean_matched']},
                                                   {'perts': scgenept['pred_perturbations'], 'pred': scgenept['pred_expr'], 'genes': scgenept['pred_genes'], 'ctrl': scgenept['ctrl_mean_matched']}],
                                                  ['GEARS', 'scGPT', 'scGenePT'])

    r_scores, r_err, r_corrs, r_perts = process(rpe1_gt, rpe1_ctrl, rpe1_gene_symbols,
                                                  [{'perts': gears_r['pred_perturbations'], 'pred': gears_r['pred_expr'], 'genes': gears_r['pred_genes']},
                                                   {'perts': scgenept_r['pred_perturbations'], 'pred': scgenept_r['pred_expr'], 'genes': scgenept_r['pred_genes'], 'ctrl': scgenept_r['ctrl_mean_matched']}],
                                                  ['GEARS', 'scGenePT'])

    # Find best K562 threshold per score and transfer to RPE1
    results = []
    for score_name, k_score in k_scores.items():
        valid = ~np.isnan(k_score) & ~np.isnan(k_err)
        if valid.sum() < 20:
            continue
        k_score_v = k_score[valid]
        k_err_v = k_err[valid]

        # Determine direction
        a_pos = roc_auc_score(k_err_v, k_score_v)
        a_neg = roc_auc_score(k_err_v, -k_score_v)
        direction = 1 if a_pos >= a_neg else -1
        best_auroc = max(a_pos, a_neg)

        # Threshold sweep for discrete transfer
        thresholds = np.linspace(k_score_v.min(), k_score_v.max(), 100)
        best_t = thresholds[np.argmax([roc_auc_score(k_err_v, (direction * (k_score_v >= t)).astype(int)) if len(np.unique(k_score_v >= t)) > 1 else 0.5 for t in thresholds])]

        # Apply to RPE1
        r_score = r_scores[score_name]
        r_valid = ~np.isnan(r_score) & ~np.isnan(r_err)
        r_pred = (direction * (r_score[r_valid] >= best_t)).astype(int)
        r_auroc = roc_auc_score(r_err[r_valid], r_pred) if len(np.unique(r_pred)) > 1 else 0.5

        # Direct AUROC on RPE1 (continuous)
        r_direct = max(roc_auc_score(r_err[r_valid], r_score[r_valid]), roc_auc_score(r_err[r_valid], -r_score[r_valid]))

        results.append({
            'score': score_name,
            'k562_auroc': float(best_auroc),
            'k562_threshold': float(best_t),
            'rpe1_transfer_auroc': float(r_auroc),
            'rpe1_direct_auroc': float(r_direct),
        })

    df_res = pd.DataFrame(results)
    print(df_res)
    df_res.to_csv(f'{OUTDIR}/sign_consistency_pergene_results.csv', index=False)

    # Plot
    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.arange(len(df_res))
    width = 0.25
    ax.bar(x - width, df_res['k562_auroc'], width, label='K562')
    ax.bar(x, df_res['rpe1_direct_auroc'], width, label='RPE1 direct')
    ax.bar(x + width, df_res['rpe1_transfer_auroc'], width, label='RPE1 zero-param transfer')
    ax.set_xticks(x)
    ax.set_xticklabels(df_res['score'])
    ax.set_ylabel('AUROC')
    ax.axhline(0.5, color='k', linestyle='--')
    ax.legend()
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/sign_consistency_pergene.png', dpi=300)
    plt.close()


if __name__ == '__main__':
    main()
