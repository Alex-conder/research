#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
adaptive_annotation.py - 自适应阈值扩展注释

对显著特征按质量从严格到宽松（默认 5/10/20/50 ppm）迭代查询 KEGG、PubChem、ChEBI，
自动停止当新增 KEGG ID 趋于平稳；输出扩展后的 cell_library_adaptive.csv 与迭代日志。

用法示例：
    python adaptive_annotation.py \
        --pos-hits output/pos_analysis/significant_hits.xlsx \
        --neg-hits output/neg_analysis_qcctrl/significant_hits.xlsx \
        --output-library library/cell_library_adaptive.csv \
        --output-dir output/adaptive
"""
import argparse
import hashlib
import json
import logging
import re
import sqlite3
import sys
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from tenacity import retry, retry_if_exception_type, stop_after_attempt, wait_exponential

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[logging.StreamHandler(sys.stdout)],
)

# ---------------------------------------------------------------------------
# 加合物中性质量转换：neutral = observed_mz + shift
# ---------------------------------------------------------------------------
ADDUCT_SHIFTS = {
    "POS": {
        "[M+H]+": -1.00727646688,
        "[M+Na]+": -22.989218,
        "[M+K]+": -38.963158,
        "[M+NH4]+": -18.033823,
        "[M-H2O+H]+": -19.017841,
    },
    "NEG": {
        "[M-H]-": 1.00727646688,
        "[M+Cl]-": -34.968853,
        "[M+HCOO]-": -44.998214,
        "[M+CH3COO]-": -59.013851,
    },
}

DEFAULT_PPM_SEQUENCE = [5.0, 10.0, 20.0, 50.0]
REQUEST_INTERVAL = 0.2  # 秒，控制 API 请求速率

# ---------------------------------------------------------------------------
# API URLs
# ---------------------------------------------------------------------------
KEGG_FIND_URL = "https://rest.kegg.jp/find/compound/{low:.6f}-{high:.6f}/exact_mass"
KEGG_GET_URL = "https://rest.kegg.jp/get/{entry}"
KEGG_LINK_PW_URL = "https://rest.kegg.jp/link/pathway/{entry}"
KEGG_CONV_PUBCHEM_URL = "https://rest.kegg.jp/conv/pubchem/{entries}"
KEGG_CONV_CHEBI_URL = "https://rest.kegg.jp/conv/chebi/{entries}"

PUBCHEM_MASS_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
    "monoisotopic_mass/range/{low:.6f}/{high:.6f}/cids/JSON"
)
PUBCHEM_PROP_URL = (
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
    "cid/{cids}/property/IUPACName,Title,CanonicalSMILES,IsomericSMILES,MonoisotopicMass/JSON"
)

CHEBI_ADV_URL = "https://www.ebi.ac.uk/chebi/backend/api/public/advanced_search/"

# ---------------------------------------------------------------------------
# SQLite 缓存
# ---------------------------------------------------------------------------
class APICache:
    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self.conn = sqlite3.connect(str(path))
        self._init_tables()

    def _init_tables(self):
        self.conn.execute(
            "CREATE TABLE IF NOT EXISTS api_cache (key TEXT PRIMARY KEY, value TEXT, ts REAL)"
        )
        self.conn.execute(
            """CREATE TABLE IF NOT EXISTS kegg_cpd (
                cpd_id TEXT PRIMARY KEY,
                name TEXT,
                exact_mass REAL,
                pubchem_cid TEXT,
                chebi_id TEXT,
                pathways TEXT
            )"""
        )
        self.conn.commit()

    @staticmethod
    def _key(source: str, **kwargs) -> str:
        s = json.dumps(kwargs, sort_keys=True, default=str)
        return f"{source}:{hashlib.md5(s.encode()).hexdigest()}"

    def get(self, source: str, **kwargs):
        key = self._key(source, **kwargs)
        row = self.conn.execute(
            "SELECT value FROM api_cache WHERE key=?", (key,)
        ).fetchone()
        if row:
            return json.loads(row[0])
        return None

    def set(self, source: str, value, **kwargs):
        key = self._key(source, **kwargs)
        self.conn.execute(
            "INSERT OR REPLACE INTO api_cache (key, value, ts) VALUES (?, ?, ?)",
            (key, json.dumps(value, default=str), time.time()),
        )
        self.conn.commit()

    def get_kegg_cpd(self, cpd_id: str):
        row = self.conn.execute(
            "SELECT name, exact_mass, pubchem_cid, chebi_id, pathways FROM kegg_cpd WHERE cpd_id=?",
            (cpd_id,),
        ).fetchone()
        if row:
            return {
                "name": row[0],
                "exact_mass": row[1],
                "pubchem_cid": row[2],
                "chebi_id": row[3],
                "pathways": json.loads(row[4]) if row[4] else [],
            }
        return None

    def set_kegg_cpd(self, cpd_id: str, name: str, exact_mass: float,
                     pubchem_cid: str, chebi_id: str, pathways: list):
        self.conn.execute(
            "INSERT OR REPLACE INTO kegg_cpd VALUES (?, ?, ?, ?, ?, ?)",
            (cpd_id, name, exact_mass, pubchem_cid, chebi_id, json.dumps(pathways)),
        )
        self.conn.commit()


# ---------------------------------------------------------------------------
# HTTP 工具
# ---------------------------------------------------------------------------
@retry(
    retry=retry_if_exception_type((requests.exceptions.RequestException,)),
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=15),
    reraise=True,
)
def http_get(url: str, timeout: int = 60):
    logging.debug("GET %s", url)
    resp = requests.get(url, timeout=timeout)
    resp.raise_for_status()
    return resp


@retry(
    retry=retry_if_exception_type((requests.exceptions.RequestException,)),
    stop=stop_after_attempt(3),
    wait=wait_exponential(multiplier=1, min=1, max=15),
    reraise=True,
)
def http_post(url: str, json_data: dict, timeout: int = 60):
    logging.debug("POST %s", url)
    resp = requests.post(url, json=json_data, timeout=timeout)
    resp.raise_for_status()
    return resp


# ---------------------------------------------------------------------------
# KEGG 查询
# ---------------------------------------------------------------------------
def query_kegg_by_mass_range(low: float, high: float, cache: APICache):
    cached = cache.get("kegg_find", low=low, high=high)
    if cached is not None:
        return cached
    url = KEGG_FIND_URL.format(low=low, high=high)
    try:
        resp = http_get(url, timeout=120)
    except requests.HTTPError as e:
        if e.response is not None and e.response.status_code == 404:
            cache.set("kegg_find", [], low=low, high=high)
            return []
        raise
    results = []
    for line in resp.text.strip().splitlines():
        parts = line.split("\t")
        if len(parts) >= 2:
            cpd_id = parts[0].replace("cpd:", "").strip()
            try:
                exact_mass = float(parts[1].strip())
            except ValueError:
                continue
            results.append({"cpd_id": cpd_id, "exact_mass": exact_mass})
    cache.set("kegg_find", results, low=low, high=high)
    time.sleep(REQUEST_INTERVAL)
    return results


def query_kegg_by_mass(neutral_mass: float, ppm_tol: float, cache: APICache):
    tol = neutral_mass * ppm_tol / 1e6
    return query_kegg_by_mass_range(neutral_mass - tol, neutral_mass + tol, cache)


def fetch_kegg_cpd_info(cpd_ids: list, cache: APICache):
    """批量拉取 KEGG compound 名称、质量、外部链接、通路。 chunk=10"""
    cpd_ids = [c for c in cpd_ids if c and not cache.get_kegg_cpd(c)]
    if not cpd_ids:
        return
    for i in range(0, len(cpd_ids), 10):
        chunk = cpd_ids[i : i + 10]
        entries = "+".join(chunk)
        # get
        try:
            text = http_get(KEGG_GET_URL.format(entry=entries), timeout=60).text
        except requests.HTTPError:
            continue
        # parse entries
        current_id = None
        data = {}
        for line in text.splitlines():
            if line.startswith("ENTRY"):
                m = re.search(r"C\d+", line)
                if m:
                    current_id = m.group(0)
                    data[current_id] = {"name": "", "exact_mass": None, "pubchem_cid": "", "chebi_id": ""}
            elif current_id and line.startswith("NAME"):
                name = line.replace("NAME", "", 1).strip().rstrip(";").split(";")[0].strip()
                if name:
                    data[current_id]["name"] = name
            elif current_id and line.startswith("EXACT_MASS"):
                m = re.search(r"\d+\.?\d*", line)
                if m:
                    data[current_id]["exact_mass"] = float(m.group(0))
            elif current_id and line.startswith("DBLINKS"):
                rest = line.replace("DBLINKS", "", 1).strip()
                if "PubChem:" in rest:
                    data[current_id]["pubchem_cid"] = rest.split("PubChem:")[1].split()[0].strip()
                if "ChEBI:" in rest:
                    data[current_id]["chebi_id"] = rest.split("ChEBI:")[1].split()[0].strip()
            elif current_id and line.startswith("///"):
                current_id = None
        # link pathway
        try:
            pw_text = http_get(KEGG_LINK_PW_URL.format(entry=entries), timeout=60).text
        except requests.HTTPError:
            pw_text = ""
        pw_map = defaultdict(list)
        for line in pw_text.strip().splitlines():
            parts = line.split("\t")
            if len(parts) == 2:
                cpd = parts[0].replace("cpd:", "").strip()
                pw = parts[1].replace("path:", "").strip()
                pw_map[cpd].append(pw)
        # conv pubchem / chebi
        pc_map = {}
        cb_map = {}
        try:
            pc_text = http_get(KEGG_CONV_PUBCHEM_URL.format(entries=entries), timeout=60).text
            for line in pc_text.strip().splitlines():
                parts = line.split("\t")
                if len(parts) == 2:
                    pc_map[parts[0].replace("cpd:", "")] = parts[1].replace("pubchem:", "")
        except requests.HTTPError:
            pass
        try:
            cb_text = http_get(KEGG_CONV_CHEBI_URL.format(entries=entries), timeout=60).text
            for line in cb_text.strip().splitlines():
                parts = line.split("\t")
                if len(parts) == 2:
                    cb_map[parts[0].replace("cpd:", "")] = parts[1].replace("chebi:", "")
        except requests.HTTPError:
            pass
        for cid, rec in data.items():
            if not rec["pubchem_cid"] and cid in pc_map:
                rec["pubchem_cid"] = pc_map[cid]
            if not rec["chebi_id"] and cid in cb_map:
                rec["chebi_id"] = cb_map[cid]
            cache.set_kegg_cpd(
                cid,
                rec["name"],
                rec["exact_mass"],
                rec["pubchem_cid"],
                rec["chebi_id"],
                pw_map.get(cid, []),
            )
        time.sleep(REQUEST_INTERVAL)


# ---------------------------------------------------------------------------
# PubChem 查询
# ---------------------------------------------------------------------------
def query_pubchem_by_mass_range(low: float, high: float, cache: APICache, max_cids: int = 200):
    cached = cache.get("pubchem_mass", low=low, high=high)
    if cached is not None:
        return cached
    url = PUBCHEM_MASS_URL.format(low=low, high=high)
    try:
        resp = http_get(url, timeout=120)
    except requests.HTTPError as e:
        if e.response is not None and e.response.status_code in (404, 500):
            cache.set("pubchem_mass", [], low=low, high=high)
            return []
        raise
    data = resp.json()
    cids = data.get("IdentifierList", {}).get("CID", [])[:max_cids]
    props = []
    if cids:
        for j in range(0, len(cids), 100):
            chunk = cids[j : j + 100]
            prop_url = PUBCHEM_PROP_URL.format(cids=",".join(map(str, chunk)))
            prop_resp = http_get(prop_url, timeout=120)
            props.extend(prop_resp.json().get("PropertyTable", {}).get("Properties", []))
    results = []
    for p in props:
        mass = p.get("MonoisotopicMass")
        if mass is None:
            continue
        results.append({
            "cid": str(p.get("CID")),
            "name": p.get("IUPACName") or p.get("Title") or "",
            "exact_mass": float(mass),
            "canonical_smiles": p.get("CanonicalSMILES", ""),
            "isomeric_smiles": p.get("IsomericSMILES", ""),
        })
    cache.set("pubchem_mass", results, low=low, high=high)
    time.sleep(REQUEST_INTERVAL)
    return results


def query_pubchem_by_mass(neutral_mass: float, ppm_tol: float, cache: APICache):
    tol = neutral_mass * ppm_tol / 1e6
    return query_pubchem_by_mass_range(neutral_mass - tol, neutral_mass + tol, cache)


# ---------------------------------------------------------------------------
# ChEBI 查询
# ---------------------------------------------------------------------------
def query_chebi_by_mass_range(low: float, high: float, cache: APICache, size: int = 200):
    cached = cache.get("chebi_mass", low=low, high=high, size=size)
    if cached is not None:
        return cached
    body = {
        "monoisotopicmass_specification": {
            "and_specification": [{"init_range": low, "final_range": high}]
        }
    }
    params = {"three_star_only": "false", "size": size}
    try:
        resp = http_post(CHEBI_ADV_URL + "?" + "&".join(f"{k}={v}" for k, v in params.items()),
                         json_data=body, timeout=120)
    except requests.HTTPError as e:
        if e.response is not None and e.response.status_code in (404, 400):
            cache.set("chebi_mass", [], low=low, high=high, size=size)
            return []
        raise
    data = resp.json()
    results = []
    for hit in data.get("results", []):
        src = hit.get("_source", {})
        mass = src.get("monoisotopicmass")
        if mass is None:
            continue
        results.append({
            "chebi_id": src.get("chebi_accession", ""),
            "name": src.get("ascii_name") or src.get("name", ""),
            "exact_mass": float(mass),
            "smiles": src.get("smiles", ""),
            "formula": src.get("formula", ""),
        })
    cache.set("chebi_mass", results, low=low, high=high, size=size)
    time.sleep(REQUEST_INTERVAL)
    return results


def query_chebi_by_mass(neutral_mass: float, ppm_tol: float, cache: APICache, size: int = 20):
    tol = neutral_mass * ppm_tol / 1e6
    return query_chebi_by_mass_range(neutral_mass - tol, neutral_mass + tol, cache, size=size)


# ---------------------------------------------------------------------------
# 显著特征读取
# ---------------------------------------------------------------------------
def load_hits(pos_path: Path, neg_path: Path):
    dfs = []
    for path, mode in [(pos_path, "POS"), (neg_path, "NEG")]:
        if not path.exists():
            logging.warning("%s 不存在，跳过", path)
            continue
        df = pd.read_excel(path)
        df["Mode"] = mode
        for col in ["Mz", "RT", "P_Value", "Fold_Change", "Adjusted_P_Value"]:
            if col not in df.columns:
                logging.warning("%s 缺少列 %s", path, col)
        dfs.append(df)
    return pd.concat(dfs, ignore_index=True) if dfs else pd.DataFrame()


# ---------------------------------------------------------------------------
# 候选集合构建
# ---------------------------------------------------------------------------
def _mass_bins(low: float, high: float, bin_size: float = 20.0):
    """生成 [low, high] 的左闭右开质量 bin。"""
    import math
    start = math.floor(low / bin_size) * bin_size
    while start < high:
        yield max(0.0, start), start + bin_size
        start += bin_size


def _query_db_bins(query_fn, low: float, high: float, cache: APICache,
                   bin_size: float = 20.0, id_key: str = "id", **kwargs):
    """按 bin 分片查询数据库并去重。"""
    seen = {}
    for b_low, b_high in _mass_bins(low, high, bin_size):
        try:
            hits = query_fn(b_low, b_high, cache, **kwargs)
            for h in hits:
                key = h.get(id_key)
                if key and key not in seen:
                    seen[key] = h
        except Exception as e:
            logging.warning("bin 查询失败 [%.4f-%.4f]: %s", b_low, b_high, e)
    return list(seen.values())


def _vector_match(hit_masses, neutral_masses, max_ppm):
    """返回 mask: 命中 x 特征 的布尔矩阵。"""
    if len(hit_masses) == 0 or len(neutral_masses) == 0:
        return np.zeros((len(hit_masses), len(neutral_masses)), dtype=bool)
    hm = np.asarray(hit_masses, dtype=float).reshape(-1, 1)
    nm = np.asarray(neutral_masses, dtype=float).reshape(1, -1)
    ppm = (hm - nm) / nm * 1e6
    return np.abs(ppm) <= max_ppm


def build_candidates(features_df: pd.DataFrame, max_ppm: float, cache: APICache,
                     adducts_pos=None, adducts_neg=None):
    """按 (mode, adduct) 全局查询 KEGG 与 ChEBI，再用矩阵向量化匹配每个特征候选。"""
    adducts_pos = adducts_pos or list(ADDUCT_SHIFTS["POS"].keys())
    adducts_neg = adducts_neg or list(ADDUCT_SHIFTS["NEG"].keys())

    grouped = features_df.groupby(["Mz", "RT", "Mode"], as_index=False).first()
    records = []

    for mode, adduct_map in (("POS", adducts_pos), ("NEG", adducts_neg)):
        sub = grouped[grouped["Mode"] == mode].reset_index(drop=True)
        if sub.empty:
            continue
        for adduct in adduct_map:
            neutral_masses = sub["Mz"].astype(float) + ADDUCT_SHIFTS[mode][adduct]
            valid = neutral_masses > 0
            if not valid.any():
                continue
            nm_min, nm_max = neutral_masses[valid].min(), neutral_masses[valid].max()
            tol_min, tol_max = nm_min * max_ppm / 1e6, nm_max * max_ppm / 1e6
            low, high = max(0.0, nm_min - tol_min), nm_max + tol_max
            logging.info("%s %s: 范围 %.4f - %.4f (%.1f ppm), %d 个特征",
                         mode, adduct, low, high, max_ppm, len(sub))

            # KEGG
            kegg_hits = []
            try:
                kegg_hits = query_kegg_by_mass_range(low, high, cache)
                logging.info("  KEGG 命中 %d", len(kegg_hits))
            except Exception as e:
                logging.warning("KEGG 查询失败 %s %s: %s", mode, adduct, e)
            # ChEBI（全局查询，size 限制）
            chebi_hits = []
            try:
                chebi_hits = query_chebi_by_mass_range(low, high, cache, size=500)
                logging.info("  ChEBI 命中 %d", len(chebi_hits))
            except Exception as e:
                logging.warning("ChEBI 查询失败 %s %s: %s", mode, adduct, e)

            # 向量化匹配
            all_hits = [
                (kegg_hits, "KEGG", {"Name": "", "SMILES": "", "Formula": "",
                                      "PubChem_CID": "", "ChEBI_ID": ""}),
                (chebi_hits, "ChEBI", {}),
            ]
            for hits, source, defaults in all_hits:
                if not hits:
                    continue
                masses = [h["exact_mass"] for h in hits]
                mask = _vector_match(masses, neutral_masses.values, max_ppm)
                for feat_idx in range(mask.shape[1]):
                    hit_indices = np.where(mask[:, feat_idx])[0]
                    if len(hit_indices) == 0:
                        continue
                    row = sub.iloc[feat_idx]
                    mz = float(row["Mz"])
                    rt = float(row["RT"]) if pd.notna(row["RT"]) else 0.0
                    neutral = float(neutral_masses.iloc[feat_idx])
                    for hi in hit_indices:
                        h = hits[hi]
                        ppm = (h["exact_mass"] - neutral) / neutral * 1e6
                        rec = {
                            "Feature_Mz": mz, "Feature_RT": rt, "Mode": mode, "Adduct": adduct,
                            "Neutral_Mass": neutral, "Source": source, "ID": h.get("cpd_id") or h.get("chebi_id"),
                            "Name": h.get("name", ""), "Exact_Mass": h["exact_mass"], "PPM_Diff": ppm,
                            "SMILES": h.get("smiles", ""),
                            "Formula": h.get("formula", ""),
                            "PubChem_CID": "",
                            "ChEBI_ID": h.get("chebi_id", ""),
                        }
                        if source == "KEGG":
                            rec["ChEBI_ID"] = ""
                        records.append(rec)
            logging.info("  %s %s 候选生成完成", mode, adduct)
    return pd.DataFrame(records)


# ---------------------------------------------------------------------------
# 迭代选择最佳候选
# ---------------------------------------------------------------------------
def select_best_per_feature(candidates_df: pd.DataFrame, ppm: float, cache: APICache):
    """给定 ppm 阈值，为每个特征选一个最佳候选。"""
    if candidates_df.empty:
        return pd.DataFrame()
    sub = candidates_df[candidates_df["PPM_Diff"].abs() <= ppm].copy()
    if sub.empty:
        return pd.DataFrame()
    # 优先级 KEGG > ChEBI > PubChem
    src_rank = {"KEGG": 0, "ChEBI": 1, "PubChem": 2}
    sub["src_rank"] = sub["Source"].map(src_rank).fillna(9)
    sub["abs_ppm"] = sub["PPM_Diff"].abs()
    best = sub.sort_values(["Feature_Mz", "Feature_RT", "Mode", "src_rank", "abs_ppm"]).groupby(
        ["Feature_Mz", "Feature_RT", "Mode"], as_index=False
    ).first()
    # 拉取 KEGG 信息
    kegg_ids = best[best["Source"] == "KEGG"]["ID"].dropna().unique().tolist()
    if kegg_ids:
        fetch_kegg_cpd_info(kegg_ids, cache)
        for _, row in best.iterrows():
            if row["Source"] != "KEGG":
                continue
            info = cache.get_kegg_cpd(row["ID"])
            if info:
                best.loc[_, "Name"] = info["name"] or best.loc[_, "Name"]
                best.loc[_, "Exact_Mass"] = info["exact_mass"] if info["exact_mass"] else best.loc[_, "Exact_Mass"]
                best.loc[_, "PubChem_CID"] = info["pubchem_cid"] or best.loc[_, "PubChem_CID"]
                best.loc[_, "ChEBI_ID"] = info["chebi_id"] or best.loc[_, "ChEBI_ID"]
    return best.drop(columns=["src_rank", "abs_ppm"], errors="ignore")


def run_adaptive_iterations(candidates_df: pd.DataFrame, cache: APICache,
                            ppm_sequence=None, convergence_ratio=0.05):
    ppm_sequence = ppm_sequence or DEFAULT_PPM_SEQUENCE
    log = []
    best_df = pd.DataFrame()
    prev_unique = 0
    chosen_ppm = ppm_sequence[0]
    for ppm in ppm_sequence:
        lib = select_best_per_feature(candidates_df, ppm, cache)
        unique_kegg = lib[lib["Source"] == "KEGG"]["ID"].dropna().nunique()
        unique_total = lib["ID"].dropna().nunique()
        new_kegg = unique_kegg - prev_unique
        ratio = new_kegg / prev_unique if prev_unique else float("inf")
        log.append({
            "ppm": ppm,
            "features_annotated": len(lib),
            "unique_annotations": unique_total,
            "unique_kegg_ids": unique_kegg,
            "new_kegg_ids": new_kegg,
            "new_ratio": ratio if ratio != float("inf") else None,
        })
        logging.info(
            "ppm=%.1f 注释=%d unique_KEGG=%d new=%d ratio=%.3f",
            ppm, len(lib), unique_kegg, new_kegg, ratio if ratio != float("inf") else -1,
        )
        # 只在确实有结果或更优时更新
        if len(lib) >= len(best_df):
            best_df = lib
            chosen_ppm = ppm
        # 收敛停止：新增 KEGG ID 比例低于阈值且已经有过结果
        if prev_unique > 0 and ratio < convergence_ratio and unique_kegg >= 5:
            logging.info("收敛于 ppm=%.1f (新增比例 %.2f%% < %.2f%%)", ppm, ratio * 100, convergence_ratio * 100)
            break
        prev_unique = unique_kegg
    return best_df, pd.DataFrame(log), chosen_ppm


# ---------------------------------------------------------------------------
# 输出格式转换：转为 cell_library 列名
# ---------------------------------------------------------------------------
def library_to_cell_format(lib_df: pd.DataFrame, mode_col="Mode"):
    """将 lib_df 转为 prepare_pathway.py 可识别的 cell_library 格式。"""
    if lib_df.empty:
        return pd.DataFrame(columns=[
            "Mz", "RT", "Mode", "Adduct", "Neutral_Mass", "Library_Match",
            "KEGG_ID", "Compound_Name", "Exact_Mass", "PPM_Error",
            "Formula", "SMILES", "PubChem_CID", "ChEBI_ID", "Pathways"
        ])
    out = pd.DataFrame({
        "Mz": lib_df["Feature_Mz"],
        "RT": lib_df["Feature_RT"],
        "Mode": lib_df[mode_col],
        "Adduct": lib_df["Adduct"],
        "Neutral_Mass": lib_df["Neutral_Mass"],
        "Library_Match": lib_df["Source"],
        "KEGG_ID": lib_df.apply(lambda r: r["ID"] if r["Source"] == "KEGG" else "", axis=1),
        "Compound_Name": lib_df["Name"],
        "Exact_Mass": lib_df["Exact_Mass"],
        "PPM_Error": lib_df["PPM_Diff"],
        "Formula": lib_df["Formula"],
        "SMILES": lib_df["SMILES"],
        "PubChem_CID": lib_df["PubChem_CID"],
        "ChEBI_ID": lib_df["ChEBI_ID"],
        "Pathways": "",
    })
    return out


# ---------------------------------------------------------------------------
# 主入口
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="自适应阈值扩展代谢物注释")
    parser.add_argument("--pos-hits", default="output/pos_analysis/significant_hits.xlsx", type=Path)
    parser.add_argument("--neg-hits", default="output/neg_analysis_qcctrl/significant_hits.xlsx", type=Path)
    parser.add_argument("--output-dir", default="output/adaptive", type=Path)
    parser.add_argument("--output-library", default="library/cell_library_adaptive.csv", type=Path)
    parser.add_argument("--cache-path", default="cache/adaptive_cache.sqlite", type=Path)
    parser.add_argument("--ppm-seq", nargs="+", type=float, default=None)
    parser.add_argument("--max-ppm", type=float, default=50.0,
                        help="用于构建候选池的最大质量容差（默认 50 ppm）")
    parser.add_argument("--convergence-ratio", type=float, default=0.05,
                        help="连续两轮新增 KEGG ID 比例低于此值时停止")
    parser.add_argument("--use-existing-candidates", action="store_true",
                        help="若存在 candidates CSV 则直接加载，跳过 API 查询")
    parser.add_argument("--max-features", type=int, default=None,
                        help="仅用于调试：限制处理前 N 个特征")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache = APICache(args.cache_path)

    hits = load_hits(args.pos_hits, args.neg_hits)
    if hits.empty:
        logging.error("未读取到任何显著特征，退出")
        sys.exit(1)
    logging.info("共读取显著特征 %d 行，唯一 m/z %d", len(hits), hits["Mz"].nunique())

    if args.max_features:
        hits = hits.head(args.max_features)

    candidates_path = args.output_dir / "all_candidates.csv"
    if args.use_existing_candidates and candidates_path.exists():
        logging.info("加载已有候选表 %s", candidates_path)
        candidates = pd.read_csv(candidates_path)
    else:
        logging.info("构建候选池（最大 ppm=%.1f）...", args.max_ppm)
        candidates = build_candidates(hits, args.max_ppm, cache)
        candidates.to_csv(candidates_path, index=False)
        logging.info("候选池已保存：%d 条", len(candidates))

    # 迭代收敛
    ppm_seq = args.ppm_seq or DEFAULT_PPM_SEQUENCE
    ppm_seq = [p for p in ppm_seq if p <= args.max_ppm]
    logging.info("开始迭代：ppm 序列 %s", ppm_seq)
    best_lib, log_df, chosen_ppm = run_adaptive_iterations(
        candidates, cache, ppm_sequence=ppm_seq, convergence_ratio=args.convergence_ratio
    )

    cell_lib = library_to_cell_format(best_lib)
    # 补通路信息
    for idx, row in cell_lib.iterrows():
        if row["KEGG_ID"]:
            info = cache.get_kegg_cpd(row["KEGG_ID"])
            if info and info["pathways"]:
                cell_lib.at[idx, "Pathways"] = ";".join(info["pathways"])

    # 保存
    cell_lib.to_csv(args.output_library, index=False)
    log_df.to_csv(args.output_dir / "iteration_log.csv", index=False)
    best_lib.to_csv(args.output_dir / "best_annotation.csv", index=False)

    # 基线摘要
    summary = {
        "input_features": int(len(hits)),
        "unique_mz": int(hits["Mz"].nunique()),
        "chosen_ppm": float(chosen_ppm),
        "annotated_features": int(len(cell_lib)),
        "unique_annotations": int(cell_lib["KEGG_ID"].replace("", pd.NA).dropna().nunique() +
                                   cell_lib["PubChem_CID"].replace("", pd.NA).dropna().nunique() +
                                   cell_lib["ChEBI_ID"].replace("", pd.NA).dropna().nunique()),
        "unique_kegg_ids": int(cell_lib["KEGG_ID"].replace("", pd.NA).dropna().nunique()),
        "iteration_log": log_df.to_dict(orient="records"),
    }
    with open(args.output_dir / "adaptive_summary.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2, ensure_ascii=False)
    logging.info(
        "完成：最终选用 ppm=%.1f，注释 %d 个特征，%d 个 unique KEGG ID",
        chosen_ppm, len(cell_lib), summary["unique_kegg_ids"],
    )


if __name__ == "__main__":
    main()
