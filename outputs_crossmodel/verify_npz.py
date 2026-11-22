# -*- coding: utf-8 -*-
"""验证本机 outputs_crossmodel/<model>/<model>_k562_predictions.npz：
- 键齐全（与 scgenept 同构骨架）
- shape 正确（217 扰动；pred_expr 列为匹配基因数）
- 无全 NaN 行/列
- 217 扰动与 scgenept_k562_labels.csv 顺序一致
- labels csv 列齐全

用法: C:/Python314/python.exe verify_npz.py <model>   (model = scgpt|geneformer|uce)
"""
import os
import sys

import numpy as np
import pandas as pd

BASE = os.path.dirname(os.path.abspath(__file__))
SCGENEPT = os.path.join(BASE, "scgenept", "scgenept_k562_labels.csv")

REQUIRED_KEYS = [
    "pred_perturbations", "target_in_input_space",
    "ctrl_mean_full", "full_gene_space",
]
EXPR_KEYS = ["pred_expr", "pred_expr_full", "pred_genes",
             "pred_gene_indices_in_h5ad", "true_expr", "ctrl_mean_matched"]
EMB_KEYS = ["pred_emb", "pred_delta_emb", "ctrl_emb"]


def main():
    model = sys.argv[1]
    d = os.path.join(BASE, model)
    z = np.load(os.path.join(d, f"{model}_k562_predictions.npz"),
                allow_pickle=True)
    keys = set(z.files)
    print(f"[{model}] keys: {sorted(keys)}")

    ref = pd.read_csv(SCGENEPT)
    assert len(ref) == 217
    ref_perts = ref.perturbation.astype(str).tolist()

    perts = [str(p) for p in z["pred_perturbations"]]
    assert perts == ref_perts, "perturbation 顺序与 scgenept 标签不一致!"
    print(f"[{model}] 217 扰动顺序与 scgenept labels 一致 OK")

    assert z["full_gene_space"].shape == (5000,)
    assert z["ctrl_mean_full"].shape == (5000,)
    tis = z["target_in_input_space"]
    assert tis.shape == (217,)
    print(f"[{model}] full_gene_space/ctrl_mean_full (5000,) OK; "
          f"target_in_input_space: {int(tis.sum())} True")

    is_expr = "pred_expr" in keys
    if is_expr:
        for k in EXPR_KEYS:
            assert k in keys, f"缺键 {k}"
        pe = z["pred_expr"]
        assert pe.shape[0] == 217 and pe.shape[1] == len(z["pred_genes"]), \
            f"pred_expr shape {pe.shape}"
        assert pe.shape[1] == 4813, f"matched genes {pe.shape[1]} != 4813"
        assert z["pred_expr_full"].shape == (217, 5000)
        assert z["true_expr"].shape == pe.shape
        # 无全 NaN 行/列
        row_nan = np.isnan(pe).all(axis=1).sum()
        col_nan = np.isnan(pe).all(axis=0).sum()
        assert row_nan == 0, f"{row_nan} 个全 NaN 行"
        assert col_nan == 0, f"{col_nan} 个全 NaN 列"
        finite = np.isfinite(pe).mean()
        print(f"[{model}] pred_expr {pe.shape} OK; 全NaN行=0 列=0; "
              f"finite {finite*100:.1f}%")
        # 预测与 ctrl 不应逐元素完全相同
        ctrl = z["ctrl_mean_matched"]
        same = np.allclose(pe, ctrl[None, :], equal_nan=True)
        print(f"[{model}] pred==ctrl 全同: {same} (mask协议下预期 False)")
    else:
        for k in EMB_KEYS:
            assert k in keys, f"缺键 {k}"
        emb = z["pred_emb"]
        assert emb.shape[0] == 217, f"pred_emb {emb.shape}"
        assert np.isfinite(emb).all(), "pred_emb 含 NaN"
        assert np.isfinite(z["ctrl_emb"]).all()
        dn = np.linalg.norm(z["pred_delta_emb"], axis=1)
        print(f"[{model}] pred_emb {emb.shape} OK; delta_emb_norm "
              f"median={np.median(dn):.4f} max={dn.max():.4f} "
              f"零值(未删除)={(dn < 1e-8).sum()}")
        assert "output_type" in keys
        print(f"[{model}] output_type: {str(z['output_type'])[:120]}")

    lab = pd.read_csv(os.path.join(d, f"{model}_k562_labels.csv"))
    assert len(lab) == 217
    assert lab.perturbation.astype(str).tolist() == ref_perts
    need_cols = {"perturbation", "pearson_delta", "pearson_delta_all",
                 "n_cells", "target_in_input_space"}
    assert need_cols <= set(lab.columns), f"labels 缺列: {need_cols - set(lab.columns)}"
    if is_expr:
        v = lab.dropna(subset=["pearson_delta"])
        print(f"[{model}] labels OK; pearson_delta n={len(v)} "
              f"mean={v.pearson_delta.mean():.4f} "
              f"median={v.pearson_delta.median():.4f}")
    else:
        assert lab.pearson_delta.isna().all(), "embedding 模型 pearson 应为 NaN"
        print(f"[{model}] labels OK; pearson 全 NaN (by design)")
    print(f"[{model}] ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
