#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Plant LC-MS metabolomics analysis pipeline
====================================================
Reads SPSS (.sav) abundance matrices, performs:
  1. Adaptive loading / QC / imputation / normalization
  2. Local or PubChem annotation
  3. Univariate + OPLS-DA multivariate differential analysis
  4. Volcano plot + OPLS-DA score plot
  5. Execution timeline log

Install dependencies:
    pip install -r requirements.txt

Example usage:
    python run_analysis.py --input ./data/raw_data.sav --lib ./library/plant_library.csv --output ./output
"""

import argparse
import logging
import re
import sys
import time
import warnings
from datetime import datetime
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore")
sns.set_style("whitegrid")
plt.rcParams["figure.dpi"] = 300

# ---------------------------------------------------------------------------
# Logging setup
# ---------------------------------------------------------------------------
def setup_logging(output_dir: Path):
    output_dir.mkdir(parents=True, exist_ok=True)
    log_path = output_dir / "metabolomics_analysis.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler(log_path, mode="a", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )
    logging.info(f"Logging initialized: {log_path}")


def log_step(step_name: str, start: float, end: float):
    logging.info(f"{step_name}: {end - start:.2f} s")


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------
def parse_args():
    parser = argparse.ArgumentParser(description="Plant LC-MS metabolomics pipeline")
    parser.add_argument("--input", default="./data/raw_data.sav", help="Input SPSS .sav file")
    parser.add_argument("--lib", default="./library/plant_library.csv", help="Local compound library CSV")
    parser.add_argument("--output", default="./output", help="Output directory")
    parser.add_argument("--group-col", default="Group", help="Sample group column name")
    parser.add_argument("--ppm", type=float, default=5.0, help="Mass tolerance in ppm")
    parser.add_argument("--rt-tol", type=float, default=0.15, help="RT tolerance in minutes")
    parser.add_argument("--fc", type=float, default=1.0, help="|Log2FC| threshold")
    parser.add_argument("--pval", type=float, default=0.05, help="Adjusted p-value threshold")
    parser.add_argument("--vip", type=float, default=1.0, help="VIP threshold")
    parser.add_argument("--isotope-mz-tol", type=float, default=0.005, help="m/z tolerance for isotope pairs (Da)")
    parser.add_argument("--isotope-rt-tol", type=float, default=0.05, help="RT tolerance for isotope pairs (min)")
    return parser.parse_args()


# ---------------------------------------------------------------------------
# Metadata extraction helpers
# ---------------------------------------------------------------------------
def _find_column(columns, candidates):
    for c in columns:
        if c.lower() in candidates:
            return c
    return None


def _parse_mz_rt_from_label(label: str):
    if not isinstance(label, str):
        return np.nan, np.nan
    mz_match = re.search(r"m[/\\]?z[:\s=]*([0-9]+\.?[0-9]*)", label, re.I)
    rt_match = re.search(r"rt[:\s=]*([0-9]+\.?[0-9]*)", label, re.I)
    rt_match2 = re.search(r"retention[_\s]?time[:\s=]*([0-9]+\.?[0-9]*)", label, re.I)
    mz = float(mz_match.group(1)) if mz_match else np.nan
    rt = float(rt_match.group(1)) if rt_match else (float(rt_match2.group(1)) if rt_match2 else np.nan)
    return mz, rt


def _parse_mz_rt_from_name(name: str):
    if not isinstance(name, str):
        return np.nan, np.nan
    # Format like M123.45T6.78 or M123.45_T6.78
    m = re.search(r"M(\d+\.?\d*)[\-_\s]?T(\d+\.?\d*)", name, re.I)
    if m:
        return float(m.group(1)), float(m.group(2))
    # Otherwise extract all numeric tokens from the string
    nums = [float(x) for x in re.findall(r"\d+\.?\d*", name)]
    if len(nums) >= 2:
        # Heuristic: larger value is m/z, smaller is retention time
        nums = sorted(nums)
        if nums[-1] > 50 and nums[0] < 50:
            return nums[-1], nums[0]
        return nums[0], nums[1]
    if len(nums) == 1:
        return (nums[0], np.nan) if nums[0] > 50 else (np.nan, nums[0])
    return np.nan, np.nan


def extract_feature_metadata(feature_cols, meta=None):
    """Build a DataFrame of Mz / RT for each feature column."""
    meta_df = pd.DataFrame(index=feature_cols, columns=["Mz", "RT"], dtype=float)

    # Try column labels from SPSS metadata first
    if meta is not None and hasattr(meta, "column_labels") and meta.column_labels:
        labels = meta.column_labels
        names = meta.column_names if hasattr(meta, "column_names") else feature_cols
        label_map = dict(zip(names, labels)) if len(names) == len(labels) else {}
        for feat in feature_cols:
            label = label_map.get(feat, "")
            mz, rt = _parse_mz_rt_from_label(label)
            if not (pd.isna(mz) and pd.isna(rt)):
                meta_df.at[feat, "Mz"] = mz
                meta_df.at[feat, "RT"] = rt

    # Fill missing from column names
    for feat in feature_cols:
        if pd.isna(meta_df.at[feat, "Mz"]) and pd.isna(meta_df.at[feat, "RT"]):
            mz, rt = _parse_mz_rt_from_name(feat)
            meta_df.at[feat, "Mz"] = mz
            meta_df.at[feat, "RT"] = rt

    return meta_df


# ---------------------------------------------------------------------------
# Isotope peak removal
# ---------------------------------------------------------------------------
def remove_isotope_peaks(intensity_df: pd.DataFrame, meta_df: pd.DataFrame,
                         mz_tol: float = 0.005, rt_tol: float = 0.05):
    """Remove isotope features (m/z ~ +1.003 Da, similar RT) keeping the most intense."""
    if intensity_df.shape[1] == 0:
        return intensity_df, meta_df

    mean_intensity = intensity_df.mean(axis=0)
    feat_order = mean_intensity.sort_values(ascending=False).index.tolist()
    keep = set()
    removed = set()

    # Work on sorted-by-intensity list so the first kept is the most intense
    for feat in feat_order:
        if feat in removed:
            continue
        keep.add(feat)
        mz_i = meta_df.at[feat, "Mz"]
        rt_i = meta_df.at[feat, "RT"]
        if pd.isna(mz_i):
            continue
        for other in feat_order:
            if other == feat or other in keep or other in removed:
                continue
            mz_j = meta_df.at[other, "Mz"]
            rt_j = meta_df.at[other, "RT"]
            if pd.isna(mz_j):
                continue
            if abs(abs(mz_j - mz_i) - 1.003) <= mz_tol and abs(rt_j - rt_i) <= rt_tol:
                removed.add(other)

    keep = sorted(keep, key=lambda f: intensity_df.columns.get_loc(f))
    logging.info(f"Isotope removal: {intensity_df.shape[1]} -> {len(keep)} features")
    return intensity_df[keep].copy(), meta_df.loc[keep].copy()


# ---------------------------------------------------------------------------
# Module 1: Load & preprocess
# ---------------------------------------------------------------------------
def module_load_preprocess(input_path: Path, output_dir: Path, group_col: str):
    start = time.time()
    logging.info(f"[A] Start data loading/preprocessing at {datetime.fromtimestamp(start)}")

    try:
        import pyreadstat
        df, meta = pyreadstat.read_sav(str(input_path))
        logging.info(f"Loaded {input_path}: shape={df.shape}")
    except Exception as e:
        raise RuntimeError(f"Failed to read SPSS file: {e}")

    # Locate group column (case-insensitive)
    gcol = _find_column(df.columns, {group_col.lower()})
    if gcol is None:
        raise ValueError(f"Group column '{group_col}' not found in {list(df.columns[:20])}...")
    groups = df[gcol].astype(str)
    logging.info(f"Groups: {groups.value_counts().to_dict()}")

    # Drop group column and obvious metadata columns from abundance matrix
    feature_df = df.drop(columns=[gcol])
    meta_cols = [c for c in feature_df.columns if c.lower() in {"mz", "m/z", "mass", "rt", "retention_time"}]
    if meta_cols:
        logging.info(f"Dropping non-abundance metadata columns: {meta_cols}")
        feature_df = feature_df.drop(columns=meta_cols)

    # Ensure numeric abundance values
    feature_df = feature_df.apply(pd.to_numeric, errors="coerce")

    # Extract m/z and RT for each feature
    meta_df = extract_feature_metadata(list(feature_df.columns), meta)

    # Remove isotope peaks using raw intensities
    intensity_df, meta_df = remove_isotope_peaks(
        feature_df, meta_df,
        mz_tol=0.005, rt_tol=0.05
    )

    # Missing-value filtering and adaptive imputation
    missing_rate = intensity_df.isna().mean()
    drop_cols = missing_rate[missing_rate > 0.30].index.tolist()
    if drop_cols:
        logging.info(f"Dropping {len(drop_cols)} features with >30% missing values")
        intensity_df = intensity_df.drop(columns=drop_cols)
        meta_df = meta_df.loc[intensity_df.columns]
        missing_rate = missing_rate[intensity_df.columns]

    from sklearn.impute import KNNImputer

    low_missing = missing_rate[(missing_rate > 0) & (missing_rate <= 0.05)].index.tolist()
    mid_missing = missing_rate[(missing_rate > 0.05) & (missing_rate <= 0.30)].index.tolist()

    if low_missing:
        intensity_df[low_missing] = intensity_df[low_missing].fillna(intensity_df[low_missing].median())
        logging.info(f"Median-imputed {len(low_missing)} features (<5% missing)")
    if mid_missing:
        k = min(5, intensity_df.shape[0] - 1) or 1
        imputer = KNNImputer(n_neighbors=k)
        imputed = imputer.fit_transform(intensity_df[mid_missing])
        intensity_df[mid_missing] = pd.DataFrame(
            imputed, index=intensity_df.index, columns=mid_missing
        )
        logging.info(f"KNN-imputed {len(mid_missing)} features (5-30% missing, k={k})")

    # Log2 transform
    xmin = intensity_df.min().min()
    pseudo = 1.0 if xmin >= 0 else 1.0 - xmin
    log2_df = np.log2(intensity_df + pseudo)
    logging.info(f"Log2-transformed abundance matrix (pseudo={pseudo:.4f})")

    # Adaptive normalization
    from scipy.stats import shapiro
    from sklearn.preprocessing import StandardScaler, QuantileTransformer

    flat = log2_df.values.flatten()
    if flat.size > 5000:
        flat = np.random.choice(flat, size=5000, replace=False)
    try:
        _, shapiro_p = shapiro(flat)
    except Exception:
        shapiro_p = 0.0
    logging.info(f"Shapiro-Wilk p={shapiro_p:.4g}")

    if shapiro_p > 0.05:
        scaler = StandardScaler()
        norm_arr = scaler.fit_transform(log2_df)
        norm_method = "Z-score (StandardScaler)"
    else:
        scaler = QuantileTransformer(output_distribution="normal", random_state=42)
        norm_arr = scaler.fit_transform(log2_df)
        norm_method = "Quantile normalization (QuantileTransformer)"

    norm_df = pd.DataFrame(norm_arr, index=log2_df.index, columns=log2_df.columns)
    logging.info(f"Normalization method: {norm_method}")

    # Save preprocessed data
    out_csv = output_dir / "preprocessed_data.csv"
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    norm_with_group = norm_df.copy()
    norm_with_group.insert(0, "Group", groups)
    norm_with_group.to_csv(out_csv, encoding="utf-8-sig")

    meta_csv = output_dir / "feature_metadata.csv"
    meta_df.to_csv(meta_csv, encoding="utf-8-sig")

    end = time.time()
    log_step("Module 1 (load & preprocess)", start, end)
    return {
        "raw_intensity": intensity_df,
        "log2": log2_df,
        "normalized": norm_df,
        "groups": groups,
        "meta": meta_df,
        "duration": end - start,
    }


# ---------------------------------------------------------------------------
# Module 2: Annotation
# ---------------------------------------------------------------------------
def module_annotation(meta_df: pd.DataFrame, lib_path: Path, ppm_tol: float, rt_tol: float):
    start = time.time()
    logging.info(f"[B] Start annotation at {datetime.fromtimestamp(start)}")

    annotations = pd.Series(index=meta_df.index, dtype=object).fillna("")

    if lib_path.exists():
        logging.info(f"Using local library: {lib_path}")
        try:
            lib = pd.read_csv(lib_path)
            name_col = _find_column(lib.columns, {"name"})
            mz_col = _find_column(lib.columns, {"mz", "m/z", "mass"})
            rt_col = _find_column(lib.columns, {"rt", "retention_time"})
            if name_col is None or mz_col is None or rt_col is None:
                raise ValueError("Library must contain Name, Mz and RT columns")

            lib_mz = lib[mz_col].astype(float).values
            lib_rt = lib[rt_col].astype(float).values
            lib_names = lib[name_col].astype(str).values

            for feat in meta_df.index:
                mz = meta_df.at[feat, "Mz"]
                rt = meta_df.at[feat, "RT"]
                if pd.isna(mz) or pd.isna(rt):
                    continue
                ppm = 1e6 * np.abs(lib_mz - mz) / np.maximum(lib_mz, 1e-9)
                rt_diff = np.abs(lib_rt - rt)
                mask = (ppm <= ppm_tol) & (rt_diff <= rt_tol)
                if mask.any():
                    best = np.argmin(ppm[mask])
                    idx = np.where(mask)[0][best]
                    annotations.at[feat] = lib_names[idx]
        except Exception as e:
            logging.error(f"Local library matching failed: {e}")
    else:
        logging.info("No local library found; querying PubChem by monoisotopic mass ...")
        try:
            import pubchempy as pcp
            import urllib.request
            import json

            unique_mzs = meta_df["Mz"].dropna().unique()
            mz_to_name = {}
            for mz in unique_mzs:
                try:
                    delta = mz * ppm_tol / 1e6
                    low, high = mz - delta, mz + delta
                    url = (
                        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
                        f"monoisotopic_mass/range/{low:.6f}/{high:.6f}/cids/JSON"
                    )
                    with urllib.request.urlopen(url, timeout=15) as resp:
                        data = json.loads(resp.read().decode())
                    cids = data.get("IdentifierList", {}).get("CID", [])[:50]
                    if not cids:
                        continue
                    props = pcp.get_properties(
                        ["IUPACName", "Title"], cids, namespace="cid"
                    )
                    names = [
                        d.get("IUPACName") or d.get("Title")
                        for d in props
                        if (d.get("IUPACName") or d.get("Title"))
                    ]
                    mz_to_name[mz] = "; ".join(names[:5]) if names else ""
                except Exception as inner:
                    logging.warning(f"PubChem query failed for m/z {mz}: {inner}")
                time.sleep(0.25)
            annotations = meta_df["Mz"].map(mz_to_name).fillna("")
        except Exception as e:
            logging.error(f"PubChem annotation failed: {e}")

    n_annotated = (annotations != "").sum()
    logging.info(f"Annotated {n_annotated} / {len(annotations)} features")

    end = time.time()
    log_step("Module 2 (annotation)", start, end)
    return {"annotations": annotations, "duration": end - start}


# ---------------------------------------------------------------------------
# OPLS-DA helper
# ---------------------------------------------------------------------------
def opls_da(X: np.ndarray, y: np.ndarray, n_orth: int = 1, n_components: int = 1):
    from sklearn.cross_decomposition import PLSRegression

    X = np.asarray(X, dtype=float)
    y = np.asarray(y, dtype=float)
    yc = y - y.mean()

    # Predictive weight from covariance with y
    w = X.T @ yc
    norm = np.linalg.norm(w)
    if norm > 0:
        w = w / norm
    t = X @ w

    # Remove orthogonal variation
    X_osc = X.copy()
    t_orth_list = []
    for _ in range(n_orth):
        p = X_osc.T @ t / (t @ t)
        X_osc = X_osc - np.outer(t, p)
        t_orth_list.append(t.copy())

    # Fit PLS on orthogonal-filtered data
    pls = PLSRegression(n_components=n_components)
    pls.fit(X_osc, y)

    # VIP (single predictive component)
    n_features = X.shape[1]
    vip = np.sqrt(n_features) * np.abs(pls.x_weights_[:, 0])

    return {
        "pls": pls,
        "vip": vip,
        "t_pred": pls.x_scores_[:, 0],
        "t_orth": t_orth_list[0] if t_orth_list else np.zeros_like(pls.x_scores_[:, 0]),
        "w": pls.x_weights_[:, 0],
    }


# ---------------------------------------------------------------------------
# Module 3: Differential analysis
# ---------------------------------------------------------------------------
def module_differential(log2_df: pd.DataFrame, norm_df: pd.DataFrame,
                        groups: pd.Series, fc_thr: float, pval_thr: float, vip_thr: float):
    start = time.time()
    logging.info(f"[C] Start differential analysis at {datetime.fromtimestamp(start)}")

    group_vals = groups.unique()
    if len(group_vals) != 2:
        raise ValueError(f"Exactly two groups required, found {group_vals}")

    # Identify control group heuristically
    control = None
    for g in group_vals:
        if re.search(r"control|ctrl|con|blank|untreated|vehicle|mock", g, re.I):
            control = g
            break
    if control is None:
        control = group_vals[0]
    treatment = [g for g in group_vals if g != control][0]
    logging.info(f"Control={control}, Treatment={treatment}")

    mask_ctrl = groups == control
    mask_trt = groups == treatment

    # Log2 fold change
    log2fc = log2_df.loc[mask_trt].mean(axis=0) - log2_df.loc[mask_ctrl].mean(axis=0)

    # Univariate tests
    from scipy.stats import levene, mannwhitneyu, ttest_ind

    def per_feature(col):
        a = col.loc[mask_ctrl].dropna().values
        b = col.loc[mask_trt].dropna().values
        if len(a) < 1 or len(b) < 1:
            return 1.0
        # 单样本分组时无法进行方差齐性检验，直接用 Mann-Whitney U（最小 n=1 仍可计算）
        if len(a) < 2 or len(b) < 2:
            try:
                _, p = mannwhitneyu(a, b, alternative="two-sided")
            except Exception:
                p = 1.0
            return p
        try:
            _, p_lev = levene(a, b)
        except Exception:
            p_lev = 0.0
        if p_lev > 0.05:
            try:
                _, p = ttest_ind(a, b, equal_var=True)
            except Exception:
                p = 1.0
        else:
            try:
                _, p = mannwhitneyu(a, b, alternative="two-sided")
            except Exception:
                p = 1.0
        return p

    p_raw = log2_df.apply(per_feature)
    p_raw = p_raw.fillna(1.0).clip(lower=1e-300)

    # FDR correction
    from statsmodels.stats.multitest import multipletests
    _, p_adj, _, _ = multipletests(p_raw.values, method="fdr_bh")
    p_adj = pd.Series(p_adj, index=p_raw.index).clip(lower=1e-300)

    # OPLS-DA
    y = (groups == treatment).astype(int).values
    opls = opls_da(norm_df.values, y)
    vip = pd.Series(opls["vip"], index=norm_df.columns)

    significant = (np.abs(log2fc) >= fc_thr) & (p_adj < pval_thr) & (vip > vip_thr)
    logging.info(f"Significant hits: {significant.sum()}")

    end = time.time()
    log_step("Module 3 (differential analysis)", start, end)
    return {
        "log2fc": log2fc,
        "p_raw": p_raw,
        "p_adj": p_adj,
        "vip": vip,
        "significant": significant,
        "control": control,
        "treatment": treatment,
        "opls_scores": (opls["t_pred"], opls["t_orth"]),
        "duration": end - start,
    }


# ---------------------------------------------------------------------------
# Module 4: Visualization & Excel export
# ---------------------------------------------------------------------------
def module_visualization(results: pd.DataFrame, significant: pd.Series,
                         groups: pd.Series, opls_scores: tuple,
                         control: str, treatment: str, output_dir: Path,
                         fc_thr: float, pval_thr: float, vip_thr: float):
    start = time.time()
    logging.info(f"[D] Start visualization/output at {datetime.fromtimestamp(start)}")

    plot_dir = output_dir / "plots"
    plot_dir.mkdir(parents=True, exist_ok=True)

    log2fc = results["Log2FC"]
    p_adj = results["P_adj"].clip(lower=1e-300)
    neg_log_p = -np.log10(p_adj)

    # Volcano plot
    fig, ax = plt.subplots(figsize=(10, 8))
    sig_up = significant & (log2fc > 0)
    sig_down = significant & (log2fc < 0)
    ns = ~significant

    ax.scatter(log2fc[ns], neg_log_p[ns], c="#999999", alpha=0.5, s=30, label="Not significant")
    ax.scatter(log2fc[sig_up], neg_log_p[sig_up], c="#E41A1C", alpha=0.8, s=50, label="Up-regulated")
    ax.scatter(log2fc[sig_down], neg_log_p[sig_down], c="#377EB8", alpha=0.8, s=50, label="Down-regulated")

    ax.axhline(-np.log10(pval_thr), color="black", linestyle="--", linewidth=0.8)
    ax.axvline(fc_thr, color="black", linestyle="--", linewidth=0.8)
    ax.axvline(-fc_thr, color="black", linestyle="--", linewidth=0.8)

    # Annotate top hits by p-value
    top_hits = results.loc[significant].nsmallest(20, "P_adj")
    for idx, row in top_hits.iterrows():
        name = row["Annotation"] if row["Annotation"] else idx
        ax.annotate(name, (row["Log2FC"], -np.log10(row["P_adj"])),
                    fontsize=7, ha="left", va="bottom")

    ax.set_xlabel("Log2 Fold Change", fontsize=12)
    ax.set_ylabel("-Log10 Adjusted P-value", fontsize=12)
    ax.set_title("Volcano Plot of Differential Metabolites", fontsize=14)
    ax.legend(loc="best")
    sns.despine()
    volcano_path = plot_dir / "volcano_plot.png"
    fig.savefig(volcano_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logging.info(f"Saved volcano plot: {volcano_path}")

    # OPLS-DA score plot
    t_pred, t_orth = opls_scores
    fig, ax = plt.subplots(figsize=(8, 8))
    palette = {control: "#377EB8", treatment: "#E41A1C"}
    for grp in [control, treatment]:
        mask = groups == grp
        ax.scatter(t_pred[mask], t_orth[mask], c=palette[grp], label=grp, alpha=0.8, s=80, edgecolor="k")
    ax.set_xlabel("Predictive score (t_pred)", fontsize=12)
    ax.set_ylabel("Orthogonal score (t_orth)", fontsize=12)
    ax.set_title("OPLS-DA Score Plot", fontsize=14)
    ax.legend(loc="best", title="Group")
    sns.despine()
    score_path = plot_dir / "oplsda_score_plot.png"
    fig.savefig(score_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    logging.info(f"Saved OPLS-DA score plot: {score_path}")

    # Excel outputs
    all_path = output_dir / "all_features_annotated.xlsx"
    sig_path = output_dir / "significant_hits.xlsx"
    results.to_excel(all_path, index=False, engine="openpyxl")
    results.loc[significant].to_excel(sig_path, index=False, engine="openpyxl")
    logging.info(f"Saved Excel: {all_path}, {sig_path}")

    end = time.time()
    log_step("Module 4 (visualization/output)", start, end)
    return {"duration": end - start}


# ---------------------------------------------------------------------------
# Module 5: Timeline summary
# ---------------------------------------------------------------------------
def write_timeline(durations: dict, output_dir: Path):
    start_total = durations["total_start"]
    total_elapsed = time.time() - start_total
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    timeline_path = output_dir / f"execution_timeline_{timestamp}.txt"

    lines = [
        "Metabolomics Analysis Execution Timeline",
        f"Generated: {datetime.now().isoformat()}",
        "=" * 50,
        f"Total runtime (s): {total_elapsed:.2f}",
        f"Total runtime (min): {total_elapsed / 60:.2f}",
        "-" * 50,
        f"Module 1 - Load & preprocess (s): {durations.get('module1', 0):.2f}",
        f"Module 2 - Annotation (s): {durations.get('module2', 0):.2f}",
        f"Module 3 - Differential analysis (s): {durations.get('module3', 0):.2f}",
        f"Module 4 - Visualization & export (s): {durations.get('module4', 0):.2f}",
        "=" * 50,
    ]
    timeline_path.write_text("\n".join(lines), encoding="utf-8")
    logging.info(f"Saved timeline: {timeline_path}")

    if total_elapsed > 1800:
        logging.warning("[WARNING] 计算耗时过长，建议检查数据维度或启用PCA降维预处理")


# ---------------------------------------------------------------------------
# Main entry
# ---------------------------------------------------------------------------
def main():
    args = parse_args()
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    setup_logging(output_dir)

    total_start = time.time()
    durations = {"total_start": total_start}

    logging.info("=" * 60)
    logging.info("Plant LC-MS metabolomics pipeline started")
    logging.info(f"Input: {args.input}")
    logging.info(f"Library: {args.lib}")
    logging.info(f"Output: {args.output}")
    logging.info("=" * 60)

    # Module 1
    try:
        prep = module_load_preprocess(Path(args.input), output_dir, args.group_col)
        durations["module1"] = prep["duration"]
    except Exception as e:
        logging.exception(f"Module 1 failed: {e}")
        return 1

    # Module 2
    try:
        annot = module_annotation(prep["meta"], Path(args.lib), args.ppm, args.rt_tol)
        durations["module2"] = annot["duration"]
    except Exception as e:
        logging.exception(f"Module 2 failed: {e}")
        annot = {"annotations": pd.Series(index=prep["meta"].index, dtype=object).fillna(""),
                 "duration": 0}

    # Module 3
    try:
        diff = module_differential(
            prep["log2"], prep["normalized"], prep["groups"],
            fc_thr=args.fc, pval_thr=args.pval, vip_thr=args.vip
        )
        durations["module3"] = diff["duration"]
    except Exception as e:
        logging.exception(f"Module 3 failed: {e}")
        return 1

    # Build results table
    results = pd.DataFrame({
        "Feature": prep["meta"].index,
        "Mz": prep["meta"]["Mz"].values,
        "RT": prep["meta"]["RT"].values,
        "Annotation": annot["annotations"].reindex(prep["meta"].index).fillna("").values,
        "Log2FC": diff["log2fc"].reindex(prep["meta"].index).values,
        "P_raw": diff["p_raw"].reindex(prep["meta"].index).values,
        "P_adj": diff["p_adj"].reindex(prep["meta"].index).values,
        "VIP": diff["vip"].reindex(prep["meta"].index).values,
    })
    significant = diff["significant"].reindex(prep["meta"].index).fillna(False)
    results.set_index("Feature", inplace=True)

    # Module 4
    try:
        vis = module_visualization(
            results, significant, prep["groups"], diff["opls_scores"],
            diff["control"], diff["treatment"], output_dir,
            fc_thr=args.fc, pval_thr=args.pval, vip_thr=args.vip
        )
        durations["module4"] = vis["duration"]
    except Exception as e:
        logging.exception(f"Module 4 failed: {e}")
        return 1

    # Module 5
    try:
        write_timeline(durations, output_dir)
    except Exception as e:
        logging.exception(f"Module 5 failed: {e}")

    logging.info("Pipeline completed successfully.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
