#!/usr/bin/env python3
"""Re-process existing mzML files to centroided MS1-only mzML."""
import argparse
import logging
import subprocess
import sys
from pathlib import Path

DEFAULT_MSCONVERT = r"path/to/msconvert.exe"  # Update to your local ProteoWizard msconvert path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", default="./data/raw_mzml")
    parser.add_argument("--output-dir", default="./data/mzml_centroided")
    parser.add_argument("--msconvert", default=DEFAULT_MSCONVERT)
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[
            logging.FileHandler("centroid_batch.log", mode="a", encoding="utf-8"),
            logging.StreamHandler(sys.stdout),
        ],
    )

    msconvert = Path(args.msconvert)
    input_dir = Path(args.input_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    mzml_files = sorted(input_dir.glob("*.mzML"))
    logging.info(f"Found {len(mzml_files)} mzML files")

    for f in mzml_files:
        cmd = [
            str(msconvert),
            str(f),
            "--mzML",
            "--filter", "peakPicking true 1-",
            "--filter", "msLevel 1",
            "-o", str(outdir),
        ]
        logging.info(f"Processing {f.name}")
        try:
            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            if result.returncode != 0:
                logging.error(f"Failed {f.name}: {result.stderr}")
            else:
                logging.info(f"Finished {f.name}")
        except subprocess.TimeoutExpired:
            logging.error(f"Timeout {f.name}")
        except Exception as e:
            logging.error(f"Exception {f.name}: {e}")

    logging.info("Batch centroiding finished.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
