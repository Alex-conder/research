#!/usr/bin/env python3
"""
Build a feature abundance matrix from Waters mzML files (profile or centroided).

Steps:
  1. For each file, read function=1 MS1 spectra.
  2. Apply vectorized local-maxima peak picking on profile data.
  3. Group peaks into chromatographic features (mz, rt, apex intensity).
  4. Align features across samples by mz/RT tolerance.
  5. Output sample x feature matrix as CSV and SPSS .sav.

Output:
  feature_matrix_{POS|NEG}.csv
  feature_matrix_{POS|NEG}.sav
  feature_metadata_{POS|NEG}.csv
"""
import argparse
import base64
import logging
import re
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pyreadstat


def decode_array(binary_data_array_elem, ns):
    """Decode a single mzML binaryDataArray to numpy array."""
    precision = np.float64
    compressed = False
    for cv in binary_data_array_elem.iter('{http://psi.hupo.org/ms/mzml}cvParam'):
        acc = cv.get('accession')
        if acc == 'MS:1000521':
            precision = np.float32
        elif acc == 'MS:1000523':
            precision = np.float64
        elif acc == 'MS:1000574':
            compressed = True
        elif acc == 'MS:1000576':
            compressed = False
    binary = binary_data_array_elem.find('.//mzml:binary', ns)
    if binary is None or binary.text is None:
        return np.array([])
    raw = base64.b64decode(binary.text)
    if compressed:
        raw = zlib.decompress(raw)
    return np.frombuffer(raw, dtype=precision)


def get_spectrum_arrays(spectrum_elem, ns):
    """Return (mz, intensity) arrays from a spectrum element."""
    mz_arr = None
    int_arr = None
    for bda in spectrum_elem.iter('{http://psi.hupo.org/ms/mzml}binaryDataArray'):
        for cv in bda.iter('{http://psi.hupo.org/ms/mzml}cvParam'):
            acc = cv.get('accession')
            if acc == 'MS:1000514' and mz_arr is None:
                mz_arr = decode_array(bda, ns)
            elif acc == 'MS:1000515' and int_arr is None:
                int_arr = decode_array(bda, ns)
    return mz_arr, int_arr


def extract_features_from_file(mzml_path: Path, mz_tol: float = 0.01, rt_tol: float = 0.1,
                                min_scans: int = 5, intensity_percentile: float = 0.05,
                                local_max_threshold_ratio: float = 0.001):
    """Extract chromatographic features (mz, rt, apex intensity) from one mzML."""
    import xml.etree.ElementTree as ET

    ns = {'mzml': 'http://psi.hupo.org/ms/mzml'}
    mzs, rts, ints = [], [], []
    count = 0

    for event, elem in ET.iterparse(str(mzml_path), events=('end',)):
        if not elem.tag.endswith('spectrum'):
            continue
        native_id = elem.get('id', '')
        if 'function=1' not in native_id:
            elem.clear()
            continue

        rt_elem = elem.find('.//mzml:scanList/mzml:scan/mzml:cvParam[@accession="MS:1000016"]', ns)
        rt = float(rt_elem.get('value')) if rt_elem is not None else 0.0

        mz_arr, int_arr = get_spectrum_arrays(elem, ns)
        elem.clear()

        if mz_arr is None or int_arr is None or len(mz_arr) < 3:
            continue
        count += 1

        # Vectorized local maxima + adaptive intensity threshold
        max_int = int_arr.max()
        threshold = max(100.0, max_int * local_max_threshold_ratio)
        lm = (
            (int_arr[1:-1] > int_arr[:-2]) &
            (int_arr[1:-1] > int_arr[2:]) &
            (int_arr[1:-1] > threshold)
        )
        idx = np.where(lm)[0] + 1
        if len(idx) == 0:
            continue

        mzs.append(mz_arr[idx])
        ints.append(int_arr[idx])
        rts.append(np.full(len(idx), rt))

    if not mzs:
        return pd.DataFrame(columns=["Mz", "RT", "Intensity"])

    all_mz = np.concatenate(mzs)
    all_rt = np.concatenate(rts)
    all_int = np.concatenate(ints)

    # Adaptive intensity threshold on picked peaks
    non_zero = all_int[all_int > 0]
    if len(non_zero) == 0:
        return pd.DataFrame(columns=["Mz", "RT", "Intensity"])
    threshold = max(np.percentile(non_zero, intensity_percentile * 100), 100.0)
    mask = all_int >= threshold
    all_mz = all_mz[mask]
    all_rt = all_rt[mask]
    all_int = all_int[mask]

    if len(all_mz) == 0:
        return pd.DataFrame(columns=["Mz", "RT", "Intensity"])

    # Sort by mz then rt
    order = np.lexsort((all_rt, all_mz))
    all_mz = all_mz[order]
    all_rt = all_rt[order]
    all_int = all_int[order]

    # Group by mz tolerance then rt gap
    features = []
    mz_diff = np.diff(all_mz)
    mz_splits = np.where(mz_diff > mz_tol)[0] + 1
    start = 0
    for end in np.append(mz_splits, len(all_mz)):
        seg_mz = all_mz[start:end]
        seg_rt = all_rt[start:end]
        seg_int = all_int[start:end]
        if len(seg_mz) == 0:
            start = end
            continue

        o = np.argsort(seg_rt)
        seg_rt = seg_rt[o]
        seg_int = seg_int[o]
        seg_mz = seg_mz[o]

        rt_diff = np.diff(seg_rt)
        rt_splits = np.where(rt_diff > rt_tol)[0] + 1
        s = 0
        for e in np.append(rt_splits, len(seg_mz)):
            if e - s < min_scans:
                s = e
                continue
            cluster_mz = seg_mz[s:e]
            cluster_rt = seg_rt[s:e]
            cluster_int = seg_int[s:e]
            apex = np.argmax(cluster_int)
            features.append({
                "Mz": float(np.median(cluster_mz)),
                "RT": float(cluster_rt[apex]),
                "Intensity": float(cluster_int[apex]),
                "N": int(e - s),
            })
            s = e
        start = end

    return pd.DataFrame(features)


def align_features(feature_dfs: dict, ppm_tol: float = 5.0, rt_tol: float = 0.2,
                   min_samples: int = 2) -> pd.DataFrame:
    """
    Align feature lists across samples.
    feature_dfs: {sample_name: DataFrame(Mz, RT, Intensity)}
    Returns long DataFrame with Sample, FeatureID, Mz, RT, Intensity.
    """
    all_feats = []
    for sample, df in feature_dfs.items():
        if df.empty:
            continue
        df = df.copy()
        df["Sample"] = sample
        all_feats.append(df[["Sample", "Mz", "RT", "Intensity"]])
    if not all_feats:
        return pd.DataFrame(columns=["Sample", "FeatureID", "Mz", "RT", "Intensity"])
    all_df = pd.concat(all_feats, ignore_index=True)
    all_df = all_df.sort_values(["Mz", "RT"]).reset_index(drop=True)

    feature_ids = np.full(len(all_df), -1, dtype=int)
    mz_arr = all_df["Mz"].values
    rt_arr = all_df["RT"].values

    # Group by relative mz tolerance (ppm)
    rel_tol = mz_arr[:-1] * ppm_tol / 1e6
    mz_diff = np.diff(mz_arr)
    mz_splits = np.where(mz_diff > rel_tol)[0] + 1

    fid = 0
    start = 0
    for end in np.append(mz_splits, len(mz_arr)):
        seg = all_df.iloc[start:end].copy()
        if len(seg) == 0:
            start = end
            continue
        seg = seg.sort_values("RT")
        rt_values = seg["RT"].values
        rt_diff2 = np.diff(rt_values)
        rt_splits2 = np.where(rt_diff2 > rt_tol)[0] + 1
        s = 0
        for e in np.append(rt_splits2, len(seg)):
            cluster = seg.iloc[s:e]
            if cluster["Sample"].nunique() >= min_samples:
                idx = cluster.index
                feature_ids[idx] = fid
                fid += 1
            s = e
        start = end

    all_df["FeatureID"] = feature_ids
    all_df = all_df[all_df["FeatureID"] >= 0]
    return all_df


def build_matrix(long_df: pd.DataFrame, metadata: pd.DataFrame) -> pd.DataFrame:
    """Pivot long table to sample x feature matrix, fill gaps with 0."""
    mat = long_df.pivot_table(index="Sample", columns="FeatureID",
                               values="Intensity", aggfunc="max").fillna(0)
    feat_meta = metadata.set_index("FeatureID")
    mat.columns = [f"M{feat_meta.loc[c, 'Mz']:.4f}_T{feat_meta.loc[c, 'RT']:.4f}" for c in mat.columns]
    return mat


def assign_group(name: str) -> str:
    upper = name.upper()
    if "BLANK" in upper or "NED" in upper:
        return "Control"
    if any(x in upper for x in ["S1", "S2", "S3"]):
        return "Treatment"
    return "QC"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", default="./data/raw_mzml", help="Directory of mzML files")
    parser.add_argument("--output-dir", default="./output", help="Output directory")
    parser.add_argument("--mz-tol", type=float, default=0.01, help="m/z tolerance for within-file grouping (Da)")
    parser.add_argument("--rt-tol", type=float, default=0.1, help="RT gap tolerance for within-file grouping (min)")
    parser.add_argument("--align-ppm", type=float, default=5.0, help="m/z tolerance for cross-sample alignment (ppm)")
    parser.add_argument("--align-rt", type=float, default=0.2, help="RT tolerance for cross-sample alignment (min)")
    parser.add_argument("--min-samples", type=int, default=2, help="Min samples a feature must appear in")
    parser.add_argument("--min-scans", type=int, default=5, help="Min scans to define a chromatographic feature")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler("build_matrix.log", mode="a", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )

    input_dir = Path(args.input_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    mzml_files = sorted(input_dir.glob("*.mzML"))
    logging.info(f"Found {len(mzml_files)} mzML files")

    pos_files = [f for f in mzml_files if "POS" in f.name.upper()]
    neg_files = [f for f in mzml_files if "NEG" in f.name.upper() or "NED" in f.name.upper()]

    def process_mode(files, mode):
        feature_dfs = {}
        for f in files:
            logging.info(f"[{mode}] Extracting features from {f.name}")
            t0 = time.time()
            feats = extract_features_from_file(
                f, args.mz_tol, args.rt_tol, args.min_scans
            )
            logging.info(f"[{mode}] {f.name}: {len(feats)} features in {time.time()-t0:.1f}s")
            sample_name = f.stem
            feature_dfs[sample_name] = feats

        if not feature_dfs:
            logging.warning(f"[{mode}] No features found")
            return

        logging.info(f"[{mode}] Aligning features across {len(feature_dfs)} samples")
        t0 = time.time()
        long_df = align_features(feature_dfs, args.align_ppm, args.align_rt, args.min_samples)
        logging.info(f"[{mode}] Aligned to {long_df['FeatureID'].nunique()} consensus features in {time.time()-t0:.1f}s")

        if long_df.empty:
            logging.warning(f"[{mode}] No consensus features")
            return

        meta = long_df.groupby("FeatureID").agg(Mz=("Mz", "median"), RT=("RT", "median")).reset_index()
        mat = build_matrix(long_df, meta)
        mat.insert(0, "Group", [assign_group(s) for s in mat.index])

        csv_path = outdir / f"feature_matrix_{mode}.csv"
        mat.to_csv(csv_path, encoding="utf-8-sig")
        logging.info(f"[{mode}] Saved {csv_path}")

        sav_path = outdir / f"feature_matrix_{mode}.sav"
        pyreadstat.write_sav(mat, str(sav_path))
        logging.info(f"[{mode}] Saved {sav_path}")

        meta_path = outdir / f"feature_metadata_{mode}.csv"
        meta.to_csv(meta_path, index=False, encoding="utf-8-sig")
        logging.info(f"[{mode}] Saved {meta_path}")

    process_mode(pos_files, "POS")
    process_mode(neg_files, "NEG")

    logging.info("Matrix build finished.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
