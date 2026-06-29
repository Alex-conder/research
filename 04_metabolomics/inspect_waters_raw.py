#!/usr/bin/env python3
"""Inspect Waters .raw directory headers and produce a summary table."""
import re
from pathlib import Path
from datetime import datetime
import pandas as pd

ROOT = Path("D:/mASS SPEC")
raw_dirs = sorted([p for p in ROOT.iterdir() if p.is_dir() and p.suffix.lower() == ".raw"])


def parse_header_text(path: Path):
    text = path.read_text(encoding="utf-8", errors="ignore")
    info = {}
    patterns = {
        "Acquired Name": r"\$\$ Acquired Name:\s*(.*)",
        "Acquired Date": r"\$\$ Acquired Date:\s*(.*)",
        "Acquired Time": r"\$\$ Acquired Time:\s*(.*)",
        "Instrument": r"\$\$ Instrument:\s*(.*)",
        "Sample Description": r"\$\$ Sample Description:\s*(.*)",
        "Job Code": r"\$\$ Job Code:\s*(.*)",
        "MS Method": r"\$\$ MS Method:\s*(.*)",
        "Inlet Method": r"\$\$ Inlet Method:\s*(.*)",
        "Bottle Number": r"\$\$ Bottle Number:\s*(.*)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, text)
        info[key] = m.group(1).strip() if m else ""
    return info


def parse_extern_inf(path: Path):
    text = path.read_text(encoding="utf-8", errors="ignore")
    info = {}
    patterns = {
        "Polarity": r"Polarity\s+(.*)",
        "Capillary kV": r"Capillary \(kV\)\s+(.*)",
        "Source Temp C": r"Source Temperature \(.*?C\)\s+(.*)",
        "Desolvation Temp C": r"Desolvation Temperature \(.*?C\)\s+(.*)",
        "Cone Gas Flow": r"Cone Gas Flow \(L/Hr\)\s+(.*)",
        "Desolvation Gas Flow": r"Desolvation Gas Flow \(L/Hr\)\s+(.*)",
        "Resolution": r"Resolution\s+(.*)",
        "IMS Gas Flow": r"IMS Gas Flow \(mL/min\)\s+(.*)",
        "Acquisition Device": r"Acquisition Device\s+(.*)",
        "Acquisition Algorithm": r"Acquisition Algorithm\s+(.*)",
    }
    for key, pat in patterns.items():
        m = re.search(pat, text)
        info[key] = m.group(1).strip() if m else ""
    return info


def count_functions(raw_dir: Path):
    return len(list(raw_dir.glob("_FUNC*.DAT")))


def dir_size_gb(raw_dir: Path):
    total = sum(f.stat().st_size for f in raw_dir.rglob("*") if f.is_file())
    return total / (1024 ** 3)


records = []
for d in raw_dirs:
    header = parse_header_text(d / "_HEADER.TXT") if (d / "_HEADER.TXT").exists() else {}
    extern = parse_extern_inf(d / "_extern.inf") if (d / "_extern.inf").exists() else {}
    n_func = count_functions(d)
    size_gb = dir_size_gb(d)
    polarity = extern.get("Polarity", "")
    mode = "POS" if polarity == "ES+" else ("NEG" if polarity == "ES-" else "?")
    if mode == "?":
        mode = "POS" if "POS" in d.name.upper() else ("NEG" if "NEG" in d.name.upper() else "?")
    rec_type = "Blank" if "BLANK" in d.name.upper() else ("Std" if "STD" in d.name.upper() else "Sample")
    remark = ""
    if n_func == 0:
        remark = "Empty/incomplete (0 functions)"
    elif n_func < 4:
        remark = f"Possibly truncated ({n_func} functions)"
    records.append({
        "File": d.name,
        "Mode": mode,
        "Type": rec_type,
        "Acquired": f"{header.get('Acquired Date','')} {header.get('Acquired Time','')}".strip(),
        "Instrument": header.get("Instrument", ""),
        "Polarity": polarity,
        "Functions": n_func,
        "Size_GB": round(size_gb, 2),
        "SampleDesc": header.get("Sample Description", ""),
        "MS_Method": Path(header.get("MS Method", "")).name,
        "Capillary_kV": extern.get("Capillary kV", ""),
        "SourceTemp_C": extern.get("Source Temp C", ""),
        "DesolvTemp_C": extern.get("Desolvation Temp C", ""),
        "Resolution": extern.get("Resolution", ""),
        "IMS_GasFlow": extern.get("IMS Gas Flow", ""),
        "Remark": remark,
    })

df = pd.DataFrame(records)
out_path = Path("./output/raw_files_summary.csv")
out_path.parent.mkdir(parents=True, exist_ok=True)
df.to_csv(out_path, index=False, encoding="utf-8-sig")
print(df.to_string(index=False))
print(f"\nSaved summary to {out_path}")
