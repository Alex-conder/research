#!/usr/bin/env python3
"""
Annotate significant hits via PubChem (monoisotopic mass search, +/- 5 ppm).
Only queries unique m/z values from significant_hits.xlsx to save time.
"""
import argparse
import json
import logging
import sys
import time
import urllib.request
from pathlib import Path

import pandas as pd
import pubchempy as pcp


def query_pubchem_names(mz: float, ppm_tol: float = 5.0, top_n: int = 5) -> str:
    delta = mz * ppm_tol / 1e6
    low, high = mz - delta, mz + delta
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
        f"monoisotopic_mass/range/{low:.6f}/{high:.6f}/cids/JSON"
    )
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            data = json.loads(resp.read().decode())
        cids = data.get("IdentifierList", {}).get("CID", [])[:50]
        if not cids:
            return ""
        props = pcp.get_properties(["IUPACName", "Title"], cids, namespace="cid")
        names = [d.get("IUPACName") or d.get("Title") for d in props if (d.get("IUPACName") or d.get("Title"))]
        return "; ".join(names[:top_n])
    except Exception as e:
        logging.warning(f"PubChem query failed for m/z {mz}: {e}")
        return ""


def annotate_excel(xlsx_path: Path, ppm_tol: float = 5.0) -> Path:
    df = pd.read_excel(xlsx_path)
    if "Mz" not in df.columns:
        logging.error(f"{xlsx_path} has no 'Mz' column")
        return None
    unique_mzs = df["Mz"].dropna().unique()
    mz_to_name = {}
    for mz in unique_mzs:
        mz_to_name[mz] = query_pubchem_names(float(mz), ppm_tol)
        time.sleep(0.25)
    df["PubChem_Name"] = df["Mz"].map(mz_to_name)
    out_path = xlsx_path.with_name(xlsx_path.stem + "_annotated.xlsx")
    df.to_excel(out_path, index=False, engine="openpyxl")
    logging.info(f"Saved annotated hits: {out_path}")
    return out_path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--hits", required=True, help="significant_hits.xlsx path")
    parser.add_argument("--ppm", type=float, default=5.0, help="Mass tolerance in ppm")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)],
    )

    annotate_excel(Path(args.hits), args.ppm)
    return 0


if __name__ == "__main__":
    sys.exit(main())
