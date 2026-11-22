"""Synthetic end-to-end tests for vcgate (no external data needed)."""
import numpy as np
import pandas as pd
import pytest

from vcgate import (compute_features, gate_auroc, median_split_labels,
                    oof_gate_scores, cv_gate_auroc, SPECTRUM_FEATURES)
from vcgate.cli import main as cli_main

rng = np.random.default_rng(42)
N_PERT, N_GENE = 120, 300


@pytest.fixture()
def synthetic():
    """构造有真信号的合成数据：预测幅度小 -> 误差大（模拟真实门控信号）。"""
    perts = [f"GENE{i}" for i in range(N_PERT)]
    genes = [f"g{i}" for i in range(N_GENE)]
    ctrl_mean = rng.normal(2.0, 0.5, size=N_GENE).astype(np.float32)
    amp = rng.uniform(0.02, 1.5, size=N_PERT).astype(np.float32)  # 扰动幅度
    pred_expr = ctrl_mean[None, :] + amp[:, None] * rng.normal(
        0, 1, size=(N_PERT, N_GENE)).astype(np.float32)
    # 误差与幅度负相关（小幅度 -> Pearson 病态 -> 大误差）
    error = 0.8 - 0.4 * amp + rng.normal(0, 0.05, size=N_PERT)
    return dict(perts=perts, genes=genes, ctrl_mean=ctrl_mean,
                pred_expr=pred_expr, error=error,
                n_cells=rng.integers(20, 300, size=N_PERT))


def test_features_shape_and_columns(synthetic):
    feat = compute_features(synthetic["pred_expr"], synthetic["perts"],
                            synthetic["genes"], synthetic["ctrl_mean"])
    assert feat.shape == (N_PERT, len(SPECTRUM_FEATURES))
    assert list(feat.columns) == SPECTRUM_FEATURES
    assert feat.index.name == "perturbation"
    assert not feat.isna().any().any()


def test_delta_max_correct(synthetic):
    feat = compute_features(synthetic["pred_expr"], synthetic["perts"],
                            synthetic["genes"], synthetic["ctrl_mean"])
    manual = np.abs(synthetic["pred_expr"] - synthetic["ctrl_mean"]).max(axis=1)
    assert np.allclose(feat["delta_max"].values, manual, atol=1e-5)


def test_gate_recovers_signal(synthetic):
    feat = compute_features(synthetic["pred_expr"], synthetic["perts"],
                            synthetic["genes"], synthetic["ctrl_mean"])
    y = median_split_labels(synthetic["error"])
    scores = oof_gate_scores(feat.values, y)
    assert scores.shape == (N_PERT,)
    assert gate_auroc(y, scores) > 0.7      # 合成信号应被检出
    mean, std = cv_gate_auroc(feat.values, y)
    assert 0.0 < mean <= 1.0 and std >= 0.0


def test_median_split_definition():
    err = np.array([0.9, 0.1, 0.6, 0.5, 0.2])
    y = median_split_labels(err)
    assert y.sum() == 2                     # 严格大于中位数才算 1
    assert y.tolist() == [1, 0, 1, 0, 0]


def test_cli_end_to_end(synthetic, tmp_path):
    pred_expr_path = tmp_path / "pred_expr.npz"
    np.savez_compressed(
        pred_expr_path,
        pred_expr=synthetic["pred_expr"],
        pred_perturbations=np.array(synthetic["perts"]),
        pred_genes=np.array(synthetic["genes"]))
    ctrl_path = tmp_path / "ctrl.npy"
    np.save(ctrl_path, synthetic["ctrl_mean"])
    pred_csv = tmp_path / "predictions.csv"
    pd.DataFrame({
        "perturbation": synthetic["perts"],
        "pearson_delta": 1.0 - synthetic["error"],
        "n_cells": synthetic["n_cells"],
    }).to_csv(pred_csv, index=False)
    out_csv = tmp_path / "gate_scores.csv"

    rc = cli_main(["oof", "--pred-expr", str(pred_expr_path),
                   "--predictions", str(pred_csv),
                   "--ctrl-mean", str(ctrl_path), "--out", str(out_csv)])
    assert rc == 0
    out = pd.read_csv(out_csv)
    assert len(out) == N_PERT
    assert {"perturbation", "gate_score", "error",
            "y_worst_half"} <= set(out.columns)
