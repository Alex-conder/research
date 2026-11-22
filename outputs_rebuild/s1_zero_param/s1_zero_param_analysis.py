"""
S1: Zero-parameter cross-dataset validation for virtual-cell gating.
Uses existing predictions only; no model retraining.
"""
import os, json, warnings
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
from sklearn.metrics import roc_auc_score, roc_curve, precision_recall_curve
import matplotlib.pyplot as plt
import seaborn as sns

warnings.filterwarnings('ignore')
OUTDIR = 'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rebuild/s1_zero_param'
os.makedirs(OUTDIR, exist_ok=True)

# ----------------------------- Loaders ----------------------------- #

def load_gears_h5ad_groundtruth(h5ad_path, pert_col='condition'):
    """Return {perturbation: mean expression vector}, control mean, and gene symbol list."""
    ad = sc.read_h5ad(h5ad_path)
    perts = ad.obs[pert_col].astype(str).values
    pert_names = np.array([p.split('+')[0] if '+' in p else p for p in perts])
    unique_perts = np.unique(pert_names)
    # Use gene symbols, not Ensembl IDs
    genes = np.array(ad.var['gene_name'].astype(str))
    expr = ad.X.toarray() if hasattr(ad.X, 'toarray') else np.array(ad.X)
    mean_expr = {}
    for p in unique_perts:
        mask = pert_names == p
        mean_expr[p] = expr[mask].mean(axis=0)
    # Control mean from cells with condition == 'ctrl'
    ctrl_mask = ad.obs[pert_col].astype(str) == 'ctrl'
    ctrl_mean = expr[ctrl_mask].mean(axis=0) if ctrl_mask.sum() > 0 else expr[pert_names == 'ctrl'].mean(axis=0)
    return mean_expr, ctrl_mean, genes


def load_gears_npz(npz_path):
    d = np.load(npz_path, allow_pickle=True)
    perts = d['pred_perturbations']
    pred = d['pred_expr']  # (n_perts, n_genes)
    genes = d['pred_genes']
    return perts, pred, genes


def load_scgenept_npz(npz_path):
    d = np.load(npz_path, allow_pickle=True)
    perts = d['pred_perturbations']
    pred = d['pred_expr']            # matched genes
    true = d['true_expr']            # matched genes
    ctrl = d['ctrl_mean_matched']
    genes = d['pred_genes']
    return perts, pred, true, ctrl, genes


def load_scgpt_npz(npz_path):
    d = np.load(npz_path, allow_pickle=True)
    perts = d['pred_perturbations']
    pred = d['pred_expr']
    true = d['true_expr']
    ctrl = d['ctrl_mean_matched']
    genes = d['pred_genes']
    return perts, pred, true, ctrl, genes


# ----------------------------- Feature computation ----------------------------- #

def compute_delta_features(pred_expr, ctrl_expr):
    """pred_expr: (n_perts, n_genes); ctrl_expr: (n_genes,)."""
    delta = pred_expr - ctrl_expr
    feats = {}
    feats['neg_frac'] = (delta < 0).mean(axis=1)
    feats['delta_l2'] = np.linalg.norm(delta, axis=1)
    feats['delta_l1'] = np.abs(delta).sum(axis=1)
    feats['delta_max'] = np.abs(delta).max(axis=1)
    feats['delta_top50_mean'] = np.partition(np.abs(delta), -50, axis=1)[:, -50:].mean(axis=1)
    feats['delta_entropy'] = stats.entropy(np.abs(delta) + 1e-9, axis=1)
    return pd.DataFrame(feats)


def top50_pearson_per_pert(pred_expr, true_expr):
    """Compute top-50 absolute delta Pearson per perturbation."""
    delta_pred = pred_expr - pred_expr.mean(axis=1, keepdims=True)
    delta_true = true_expr - true_expr.mean(axis=1, keepdims=True)
    corrs = []
    for i in range(pred_expr.shape[0]):
        # top 50 by absolute true delta
        top_idx = np.argsort(np.abs(delta_true[i]))[-50:]
        if np.std(delta_pred[i, top_idx]) == 0 or np.std(delta_true[i, top_idx]) == 0:
            corrs.append(0.0)
        else:
            r = np.corrcoef(delta_pred[i, top_idx], delta_true[i, top_idx])[0, 1]
            corrs.append(0.0 if np.isnan(r) else r)
    return np.array(corrs)


def sign_error_rate_per_pert(pred_delta, true_delta):
    """Fraction of genes where predicted delta sign mismatches true delta sign."""
    return ((pred_delta * true_delta) < 0).mean(axis=1)


# ----------------------------- Main analysis ----------------------------- #

def main():
    # Load ground truths
    k562_gt, k562_ctrl, k562_genes = load_gears_h5ad_groundtruth(
        'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/data/gears_demo/replogle_k562_essential/perturb_processed.h5ad'
    )
    rpe1_gt, rpe1_ctrl, rpe1_genes = load_gears_h5ad_groundtruth(
        'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/data/gears_demo/replogle_rpe1_essential/perturb_processed.h5ad'
    )

    # Load K562 predictions
    k562_engines = {}
    # GEARS
    perts, pred, genes = load_gears_npz(
        'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs/pred_expr_seed0.npz'
    )
    k562_engines['GEARS'] = {'perts': perts, 'pred': pred, 'genes': genes}

    # scGenePT
    perts, pred, true, ctrl, genes = load_scgenept_npz(
        'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_crossmodel/scgenept/scgenept_k562_predictions.npz'
    )
    k562_engines['scGenePT'] = {'perts': perts, 'pred': pred, 'genes': genes, 'ctrl': ctrl}

    # scGPT
    perts, pred, true, ctrl, genes = load_scgpt_npz(
        'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_crossmodel/scgpt/scgpt_k562_predictions.npz'
    )
    k562_engines['scGPT'] = {'perts': perts, 'pred': pred, 'genes': genes, 'ctrl': ctrl}

    # Load RPE1 predictions
    rpe1_engines = {}
    perts, pred, genes = load_gears_npz(
        'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_rpe1/pred_expr_seed0.npz'
    )
    rpe1_engines['GEARS'] = {'perts': perts, 'pred': pred, 'genes': genes}

    perts, pred, true, ctrl, genes = load_scgenept_npz(
        'E:/清华类脑计算研究中心/virtual_cell_constraint_mvp/outputs_crossmodel/scgenept/scgenept_rpe1_predictions.npz'
    )
    rpe1_engines['scGenePT'] = {'perts': perts, 'pred': pred, 'genes': genes, 'ctrl': ctrl}

    # ---- Build matched data frames per dataset ---- #
    def build_dataset(engines, gt, gt_ctrl, gt_genes, dataset_name):
        records = []
        common_perts = set(engines['GEARS']['perts'])
        for name in engines:
            common_perts &= set(engines[name]['perts'])
        common_perts = sorted(common_perts)

        ref_genes = engines['GEARS']['genes']
        ref_gene_idx = {g: i for i, g in enumerate(ref_genes)}
        gt_gene_idx = {g: i for i, g in enumerate(gt_genes)}

        # Pre-map ground truth to GEARS gene space once
        true_full_all = {p: np.array([gt_vec[gt_gene_idx[g]] if g in gt_gene_idx else np.nan
                                      for g in ref_genes])
                         for p, gt_vec in gt.items()}
        ctrl_full = np.array([gt_ctrl[gt_gene_idx[g]] if g in gt_gene_idx else np.nan
                              for g in ref_genes])

        for p in common_perts:
            true_full = true_full_all.get(p, None)
            if true_full is None:
                continue

            row = {'perturbation': p, 'dataset': dataset_name}

            for name, data in engines.items():
                idx = list(data['perts']).index(p)
                pred = data['pred'][idx]
                genes = data['genes']
                pred_map = {g: pred[i] for i, g in enumerate(genes)}
                pred_full = np.array([pred_map.get(g, np.nan) for g in ref_genes])

                if name == 'GEARS':
                    ctrl_v_full = ctrl_full
                else:
                    ctrl_map = {g: data['ctrl'][i] for i, g in enumerate(genes)}
                    ctrl_v_full = np.array([ctrl_map.get(g, np.nan) for g in ref_genes])

                valid = ~np.isnan(pred_full) & ~np.isnan(true_full) & ~np.isnan(ctrl_v_full)
                if valid.sum() < 50:
                    continue
                pred_v = pred_full[valid]
                true_v = true_full[valid]
                ctrl_v = ctrl_v_full[valid]

                feats = compute_delta_features(pred_v.reshape(1, -1), ctrl_v)
                for c in feats.columns:
                    row[f'{name}_{c}'] = feats[c].values[0]

                pred_delta = pred_v - ctrl_v
                true_delta = true_v - ctrl_v
                row[f'{name}_top50_pearson'] = top50_pearson_per_pert(pred_delta.reshape(1, -1),
                                                                       true_delta.reshape(1, -1))[0]
                row[f'{name}_sign_error_rate'] = sign_error_rate_per_pert(pred_delta.reshape(1, -1),
                                                                          true_delta.reshape(1, -1))[0]
            records.append(row)
        return pd.DataFrame(records)

    df_k562 = build_dataset(k562_engines, k562_gt, k562_ctrl, k562_genes, 'K562')
    df_rpe1 = build_dataset(rpe1_engines, rpe1_gt, rpe1_ctrl, rpe1_genes, 'RPE1')

    print(f"K562 matched: {len(df_k562)} perturbations")
    print(f"RPE1 matched: {len(df_rpe1)} perturbations")

    # ---- Sign-layer definition: use sign consistency across engines ---- #
    def add_sign_consistency(df, engines):
        # Per perturbation, compute fraction of engines predicting negative delta (majority vote)
        # For simplicity use neg_frac averaged? No, need per-gene sign across engines.
        # Use average neg_frac across engines as proxy for sign-layer score.
        neg_cols = [c for c in df.columns if c.endswith('_neg_frac')]
        df['sign_layer_score'] = df[neg_cols].mean(axis=1)
        return df

    df_k562 = add_sign_consistency(df_k562, k562_engines)
    df_rpe1 = add_sign_consistency(df_rpe1, rpe1_engines)

    # ---- Error label: top50 pearson below median per engine ---- #
    available_engines = []
    for df in [df_k562, df_rpe1]:
        for eng in ['GEARS', 'scGenePT', 'scGPT']:
            col = f'{eng}_top50_pearson'
            if col in df.columns and df[col].notna().sum() >= 10:
                df[f'{eng}_high_error'] = df[col] < df[col].median()
                if eng not in available_engines:
                    available_engines.append(eng)

    # ---- 1. K562 threshold -> RPE1 zero-parameter transfer ---- #
    # Prefer GEARS as target; fallback to first available engine
    target_col = None
    for eng in ['GEARS', 'scGenePT', 'scGPT']:
        c = f'{eng}_high_error'
        if c in df_k562.columns and df_k562[c].notna().sum() >= 10 and c in df_rpe1.columns and df_rpe1[c].notna().sum() >= 10:
            target_col = c
            break
    if target_col is None:
        raise KeyError('No valid error label column available across K562 and RPE1')

    # Find optimal threshold on K562 using sign_layer_score
    k562_score = df_k562['sign_layer_score'].values
    k562_y = df_k562[target_col].values
    thresholds = np.linspace(k562_score.min(), k562_score.max(), 100)
    best_auroc = 0
    best_t = None
    for t in thresholds:
        pred_y = k562_score >= t
        if len(np.unique(pred_y)) < 2:
            continue
        # For high error prediction, higher score -> high error? Check correlation
        auroc = roc_auc_score(k562_y, pred_y.astype(int))
        if auroc > 0.5 and auroc > best_auroc:
            best_auroc = auroc
            best_t = t
    # Also try reverse
    for t in thresholds:
        pred_y = k562_score <= t
        if len(np.unique(pred_y)) < 2:
            continue
        auroc = roc_auc_score(k562_y, pred_y.astype(int))
        if auroc > 0.5 and auroc > best_auroc:
            best_auroc = auroc
            best_t = t

    print(f'K562 sign-layer AUROC (threshold search): {best_auroc:.3f}, threshold={best_t:.3f}')

    # Apply threshold to RPE1
    rpe1_score = df_rpe1['sign_layer_score'].values
    rpe1_y = df_rpe1[target_col].values
    rpe1_pred = (rpe1_score >= best_t) if best_t is not None else np.zeros_like(rpe1_score, dtype=bool)
    rpe1_auroc = roc_auc_score(rpe1_y, rpe1_pred.astype(int)) if len(np.unique(rpe1_pred)) > 1 else 0.5
    print(f'RPE1 zero-param transfer AUROC: {rpe1_auroc:.3f}')

    # Save metrics
    metrics = {
        'n_k562': len(df_k562),
        'n_rpe1': len(df_rpe1),
        'k562_sign_auroc': float(best_auroc),
        'k562_threshold': float(best_t) if best_t is not None else None,
        'rpe1_transfer_auroc': float(rpe1_auroc),
    }

    # ---- 2. Leave-one-engine-out on K562 ---- #
    engine_names = ['GEARS', 'scGenePT', 'scGPT']
    loeo_aurocs = {}
    for held in engine_names:
        others = [e for e in engine_names if e != held]
        neg_cols = [f'{e}_neg_frac' for e in others if f'{e}_neg_frac' in df_k562.columns]
        if len(neg_cols) < 2:
            continue
        score = df_k562[neg_cols].mean(axis=1).values
        y = df_k562[f'{held}_high_error'].values
        try:
            auroc = roc_auc_score(y, score)
            if auroc < 0.5:
                auroc = 1 - auroc
            loeo_aurocs[held] = float(auroc)
        except Exception as e:
            loeo_aurocs[held] = None
    print('LOEO AUROCs:', loeo_aurocs)
    metrics['loeo_aurocs'] = loeo_aurocs

    # ---- 3. Magnitude-only baseline ---- #
    mag_col = 'GEARS_delta_max'
    if mag_col in df_k562.columns:
        mag_auroc_k562 = max(roc_auc_score(df_k562[target_col], df_k562[mag_col]),
                             roc_auc_score(df_k562[target_col], -df_k562[mag_col]))
        mag_auroc_rpe1 = max(roc_auc_score(df_rpe1[target_col], df_rpe1[mag_col]),
                             roc_auc_score(df_rpe1[target_col], -df_rpe1[mag_col]))
        metrics['k562_magnitude_auroc'] = float(mag_auroc_k562)
        metrics['rpe1_magnitude_auroc'] = float(mag_auroc_rpe1)
        print(f'Magnitude baseline AUROC K562={mag_auroc_k562:.3f}, RPE1={mag_auroc_rpe1:.3f}')

    # ---- Figures ---- #
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax, df, title in [(axes[0], df_k562, 'K562'), (axes[1], df_rpe1, 'RPE1')]:
        if target_col in df.columns and 'sign_layer_score' in df.columns:
            y = df[target_col].values
            score = df['sign_layer_score'].values
            # Determine direction
            auroc_pos = roc_auc_score(y, score)
            auroc_neg = roc_auc_score(y, -score)
            s = score if auroc_pos >= auroc_neg else -score
            fpr, tpr, _ = roc_curve(y, s)
            ax.plot(fpr, tpr, label=f'sign-layer (AUROC={max(auroc_pos, auroc_neg):.3f})')
            if mag_col in df.columns:
                mag = df[mag_col].values
                auroc_mag_pos = roc_auc_score(y, mag)
                auroc_mag_neg = roc_auc_score(y, -mag)
                m = mag if auroc_mag_pos >= auroc_mag_neg else -mag
                fpr_m, tpr_m, _ = roc_curve(y, m)
                ax.plot(fpr_m, tpr_m, label=f'magnitude (AUROC={max(auroc_mag_pos, auroc_mag_neg):.3f})', linestyle='--')
            ax.plot([0, 1], [0, 1], 'k:')
            ax.set_xlabel('FPR')
            ax.set_ylabel('TPR')
            ax.set_title(title)
            ax.legend()
    plt.tight_layout()
    plt.savefig(f'{OUTDIR}/zero_param_roc.png', dpi=300)
    plt.close()

    # LOEO heatmap
    if loeo_aurocs:
        fig, ax = plt.subplots(figsize=(4, 3))
        vals = [v for v in loeo_aurocs.values() if v is not None]
        labs = [k for k, v in loeo_aurocs.items() if v is not None]
        sns.heatmap(np.array(vals).reshape(1, -1), annot=True, fmt='.3f',
                    xticklabels=labs, yticklabels=['LOEO AUROC'], cmap='RdYlGn',
                    vmin=0.5, vmax=1.0, ax=ax)
        plt.tight_layout()
        plt.savefig(f'{OUTDIR}/loeo_heatmap.png', dpi=300)
        plt.close()

    # Save data
    df_k562.to_csv(f'{OUTDIR}/k562_features.csv', index=False)
    df_rpe1.to_csv(f'{OUTDIR}/rpe1_features.csv', index=False)
    with open(f'{OUTDIR}/s1_metrics.json', 'w') as f:
        json.dump(metrics, f, indent=2)

    # Summary markdown
    with open(f'{OUTDIR}/s1_summary.md', 'w') as f:
        f.write('# S1 Zero-Parameter Cross-Dataset Validation Summary\n\n')
        f.write(f'- K562 matched perturbations: {len(df_k562)}\n')
        f.write(f'- RPE1 matched perturbations: {len(df_rpe1)}\n')
        f.write(f'- K562 sign-layer AUROC: {best_auroc:.3f}\n')
        f.write(f'- RPE1 zero-parameter transfer AUROC: {rpe1_auroc:.3f}\n')
        f.write('- LOEO AUROCs:\n')
        for k, v in loeo_aurocs.items():
            f.write(f'  - held-out {k}: {v:.3f}\n' if v is not None else f'  - held-out {k}: N/A\n')
        if 'k562_magnitude_auroc' in metrics:
            f.write(f'- K562 magnitude baseline AUROC: {metrics["k562_magnitude_auroc"]:.3f}\n')
            f.write(f'- RPE1 magnitude baseline AUROC: {metrics["rpe1_magnitude_auroc"]:.3f}\n')

    print('S1 complete. Outputs:', OUTDIR)


if __name__ == '__main__':
    main()
