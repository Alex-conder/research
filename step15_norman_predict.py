# -*- coding: utf-8 -*-
"""Predict Norman 2019 single-gene perturbations with the main K562 GEARS seed0 model."""
import os
import time
import numpy as np
import pandas as pd
import anndata
import scipy.sparse as sp
import torch
from gears import GEARS, PertData

BASE = os.path.dirname(os.path.abspath(__file__))
ADATA = os.path.join(BASE, "data/gears_demo/norman/perturb_processed.h5ad")
WEIGHTS = os.path.join(BASE, "outputs/gears_model_seed0/model.pt")
OUT_DIR = os.path.join(BASE, "outputs_external")
os.makedirs(OUT_DIR, exist_ok=True)
DEVICE = "cuda" if torch.cuda.is_available() else "cpu"


def main():
    t0 = time.time()
    adata = anndata.read_h5ad(ADATA)
    print(f"loaded {adata.shape} ({time.time()-t0:.1f}s)", flush=True)
    # The K562 model expects exactly 5000 genes in the same order. Norman
    # measures only a subset, so unmeasured genes are zero-filled for model input.
    k562 = np.load(os.path.join(BASE, "outputs/pred_expr_seed0.npz"), allow_pickle=True)
    k562_genes = k562["pred_genes"].astype(str).tolist()
    norman_var = adata.var["gene_name"].astype(str).tolist()
    pos = {g: i for i, g in enumerate(norman_var)}
    cols = []
    for g in k562_genes:
        if g in pos:
            cols.append(adata.X[:, pos[g]])
        else:
            cols.append(sp.csr_matrix((adata.n_obs, 1), dtype=adata.X.dtype))
    X = sp.hstack(cols, format="csr")
    var = pd.DataFrame({"gene_name": k562_genes}, index=k562_genes)
    adata = anndata.AnnData(X=X, obs=adata.obs.copy(), var=var)
    adata.uns["non_zeros_gene_idx"] = {g: i for i, g in enumerate(k562_genes)}
    print(f"aligned to K562 gene space {adata.shape} (zero-filled {len(k562_genes)-len(pos)} genes)", flush=True)
    conditions = adata.obs["condition"].astype(str).to_numpy()
    single_conditions = []
    for c in np.unique(conditions):
        parts = c.split("+")
        if len(parts) == 2 and "ctrl" in parts and c != "ctrl":
            single_conditions.append(c)
    genes = sorted({c.split("+")[0] if c.split("+")[0] != "ctrl" else c.split("+")[1]
                    for c in single_conditions})
    print(f"single conditions={len(single_conditions)}, unique genes={len(genes)}", flush=True)

    pert_data = PertData(os.path.join(BASE, "data/gears_demo"))
    pert_data.new_data_process("norman_for_gears", adata=adata, skip_calc_de=True)
    pert_data.prepare_split(split="no_test", seed=0)
    if "val" not in pert_data.set2conditions:
        pert_data.set2conditions["val"] = pert_data.set2conditions.get("train", ["ctrl"])
    pert_data.get_dataloader(batch_size=32, test_batch_size=128)
    model = GEARS(pert_data, device=DEVICE)
    model.model_initialize(hidden_size=64)
    graph = set(pert_data.pert_names)
    pred_genes = sorted([g for g in genes if g in graph])
    print(f"predictable genes={len(pred_genes)}", flush=True)
    gene_names = pert_data.adata.var["gene_name"].astype(str).tolist()
    ctrl = adata.X[conditions == "ctrl"].mean(axis=0)
    ctrl_mean = np.asarray(ctrl).ravel().astype(np.float32)

    for seed in [0, 1, 2]:
        weights = os.path.join(BASE, f"outputs/gears_model_seed{seed}/model.pt")
        sd = torch.load(weights, map_location=DEVICE, weights_only=True)
        model.model.load_state_dict(sd)
        model.best_model.load_state_dict(sd)
        # GEARS caches predictions in saved_pred; clear it between model seeds.
        model.saved_pred = {}
        model.saved_logvar_sum = {}
        pred = model.predict([(g,) for g in pred_genes])
        pred_arr = np.stack([pred[g].astype(np.float32) for g in pred_genes], axis=0)
        np.savez_compressed(
            os.path.join(OUT_DIR, f"norman_predictions_seed{seed}.npz"),
            genes=np.array(pred_genes),
            gene_names=np.array(gene_names),
            predictions=pred_arr,
            ctrl_mean=ctrl_mean,
        )
        print(f"saved seed{seed} predictions {pred_arr.shape} ({time.time()-t0:.1f}s)", flush=True)


if __name__ == "__main__":
    main()
