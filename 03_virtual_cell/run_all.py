"""
一键运行所有核心 notebook。

用法：
    python run_all.py

要求：
    - 已安装 jupyter 和 nbconvert
    - 已激活虚拟环境
"""

import os
import sys
import subprocess

NOTEBOOKS = [
    "notebooks/01_hh_advanced.ipynb",
    "notebooks/02_fba_genome_scale.ipynb",
    "notebooks/03_scvi_validation.ipynb",
    "notebooks/04_integration_interface.ipynb",
]


def run_notebook(path):
    print(f"\n{'='*60}")
    print(f"Running: {path}")
    print('=' * 60)
    cmd = [
        sys.executable, "-m", "nbconvert",
        "--to", "notebook",
        "--execute",
        "--inplace",
        "--ExecutePreprocessor.timeout=600",
        path,
    ]
    result = subprocess.run(cmd, cwd=os.path.dirname(os.path.abspath(__file__)))
    if result.returncode != 0:
        print(f"[ERROR] Failed: {path}")
        sys.exit(result.returncode)
    print(f"[OK] Completed: {path}")


def main():
    for nb in NOTEBOOKS:
        if not os.path.exists(nb):
            print(f"[ERROR] Notebook not found: {nb}")
            sys.exit(1)
    
    for nb in NOTEBOOKS:
        run_notebook(nb)
    
    print("\n" + "=" * 60)
    print("All notebooks completed successfully.")
    print("=" * 60)


if __name__ == "__main__":
    main()
