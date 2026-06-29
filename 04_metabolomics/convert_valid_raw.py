#!/usr/bin/env python3
"""
Batch convert Waters .raw files to mzML using msconvert,
with recommended settings for SYNAPT G2-Si HDMS data.

Recommended for metabolomics:
  --combineIonMobilitySpectra : collapse ion-mobility drift bins
  --filter "scanEvent 1"      : keep only the low-energy MS1 function

Files with 0 functions or missing headers are skipped automatically.
"""
import argparse
import logging
import subprocess
import sys
from pathlib import Path

# Default msconvert path (update if your version differs)
DEFAULT_MSCONVERT = r"path/to/msconvert.exe"  # Update to your local ProteoWizard msconvert path


def count_functions(raw_dir: Path) -> int:
    return len(list(raw_dir.glob("_FUNC*.DAT")))


def has_header(raw_dir: Path) -> bool:
    return (raw_dir / "_HEADER.TXT").exists()


def convert_file(raw_dir: Path, outdir: Path, msconvert: Path, extra_args: list) -> bool:
    outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        str(msconvert),
        str(raw_dir),
        "--mzML",
        "--combineIonMobilitySpectra",
        "-o", str(outdir),
    ] + extra_args
    logging.info(f"Running: {' '.join(cmd)}")
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            logging.error(f"msconvert failed for {raw_dir.name}\n{result.stderr}")
            return False
        logging.info(f"Success: {raw_dir.name}")
        if result.stdout:
            logging.debug(result.stdout)
        return True
    except subprocess.TimeoutExpired:
        logging.error(f"Timeout converting {raw_dir.name}")
        return False
    except Exception as e:
        logging.error(f"Exception converting {raw_dir.name}: {e}")
        return False


def main():
    parser = argparse.ArgumentParser(description="Batch convert Waters RAW to mzML")
    parser.add_argument("--input-dir", default="D:/mASS SPEC", help="Directory containing .raw folders")
    parser.add_argument("--output-dir", default="D:/mASS SPEC/mzML_scan1", help="Output directory")
    parser.add_argument("--msconvert", default=DEFAULT_MSCONVERT, help="Path to msconvert.exe")
    parser.add_argument("--scan-event", type=int, default=1, help="Scan event (function) to keep")
    parser.add_argument("--min-functions", type=int, default=2, help="Skip files with fewer functions")
    parser.add_argument("--skip", nargs="*", default=[], help="Additional filenames to skip")
    parser.add_argument("--extra", nargs="*", default=[], help="Extra msconvert arguments")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler("msconvert_batch.log", mode="a", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )

    msconvert = Path(args.msconvert)
    if not msconvert.exists():
        logging.error(f"msconvert not found: {msconvert}")
        return 1

    input_dir = Path(args.input_dir)
    outdir = Path(args.output_dir)
    raw_dirs = sorted([p for p in input_dir.iterdir() if p.is_dir() and p.suffix.lower() == ".raw"])

    extra = [f'--filter', f'scanEvent {args.scan_event}'] + args.extra
    skip_set = set(args.skip)

    for d in raw_dirs:
        if d.name in skip_set:
            logging.info(f"Skipping {d.name} (user specified)")
            continue
        n_func = count_functions(d)
        if n_func == 0:
            logging.warning(f"Skipping {d.name}: 0 functions (empty/corrupt)")
            continue
        if n_func < args.min_functions:
            logging.warning(f"Skipping {d.name}: only {n_func} functions (incomplete)")
            continue
        if not has_header(d):
            logging.warning(f"Skipping {d.name}: missing _HEADER.TXT (incomplete)")
            continue
        convert_file(d, outdir, msconvert, extra)

    logging.info("Batch conversion finished.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
