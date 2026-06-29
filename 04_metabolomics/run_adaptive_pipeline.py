#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_adaptive_pipeline.py - 自适应注释与通路分析全流程

依次执行：
  1. adaptive_annotation.py（KEGG + ChEBI 自适应阈值注释）
  2. consensus_annotation.py（多库共识打分）
  3. prepare_pathway.py（生成通路输入/映射）
  4. kegg_ora.py（本地 KEGG ORA）
  5. mummichog_local.py（m/z-to-pathway 兜底）
  6. generate_adaptive_report.py（收敛报告）

MetaboAnalystR ORA 为可选项（--run-metaboanalyst），因为该步骤依赖 R 环境
且常因化合物名称映射不足而失败。
"""
import argparse
import os
import shutil
import subprocess
import sys
from pathlib import Path, PureWindowsPath

PYTHON = Path("venv/Scripts/python.exe").resolve()
RSCRIPT = shutil.which("Rscript") or "C:/Program Files/R/R-4.5.2/bin/x64/Rscript.exe"
# 兼容 Git Bash / MSYS2 路径
if "/bin/Rscript" in RSCRIPT and "Rscript.exe" not in RSCRIPT:
    RSCRIPT = RSCRIPT.replace("/bin/Rscript", "/bin/x64/Rscript.exe")


def run(cmd, check=True, **kwargs):
    print(f"\n>>> {' '.join(str(c) for c in cmd)}")
    return subprocess.run(cmd, check=check, **kwargs)


def run_or_warn(cmd):
    """运行命令，失败时打印警告但不中断流水线。"""
    try:
        run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[WARN] 步骤失败（exit={e.returncode}），继续后续步骤...", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(description="Run adaptive annotation pipeline")
    parser.add_argument("--max-ppm", type=float, default=20.0)
    parser.add_argument("--ppm-seq", nargs="+", type=float, default=[5.0, 10.0, 20.0])
    parser.add_argument("--skip-adaptive", action="store_true")
    parser.add_argument("--skip-consensus", action="store_true")
    parser.add_argument("--skip-mummichog", action="store_true")
    parser.add_argument("--run-metaboanalyst", action="store_true",
                        help="尝试运行 MetaboAnalystR ORA（可能因名称映射失败）")
    args = parser.parse_args()

    if not PYTHON.exists():
        print(f"未找到 venv Python: {PYTHON}", file=sys.stderr)
        sys.exit(1)

    # 1. adaptive annotation
    if not args.skip_adaptive:
        run([
            str(PYTHON), "adaptive_annotation.py",
            "--max-ppm", str(args.max_ppm),
            "--ppm-seq"] + [str(p) for p in args.ppm_seq]
        )

    # 2. consensus annotation
    if not args.skip_consensus:
        run([
            str(PYTHON), "consensus_annotation.py",
            "--candidates", "output/adaptive/all_candidates.csv",
            "--output-dir", "output/consensus",
            "--output-library", "library/cell_library_consensus.csv",
            "--max-ppm", str(args.max_ppm),
        ])

    # 3. prepare pathway
    run([
        str(PYTHON), "prepare_pathway.py",
        "--lib", "library/cell_library_consensus.csv",
        "--pos-hits", "output/pos_analysis/significant_hits.xlsx",
        "--neg-hits", "output/neg_analysis_qcctrl/significant_hits.xlsx",
        "--output-dir", "output/pathway_consensus",
    ])

    # 4. KEGG ORA
    run([
        str(PYTHON), "kegg_ora.py",
        "--mapping", "output/pathway_consensus/pathway_mapping.csv",
        "--summary", "output/pathway_consensus/pathway_summary.csv",
        "--output", "output/pathway_consensus/pathway_enrichment.csv",
    ])

    # 5. MetaboAnalystR ORA（可选，失败不中断）
    if args.run_metaboanalyst:
        run_or_warn([
            RSCRIPT, "run_metaboanalyst_local.R",
            "consensus",
            "output/pathway_consensus",
        ])

    # 6. Mummichog
    if not args.skip_mummichog:
        run([
            str(PYTHON), "mummichog_local.py",
            "--ppm", str(args.max_ppm),
            "--output-dir", "output/mummichog",
        ])

    # 7. Report
    run([
        str(PYTHON), "generate_adaptive_report.py",
        "--output-dir", "output",
        "--out-md", "output/adaptive_convergence_report.md",
    ])

    print("\nPipeline completed.")


if __name__ == "__main__":
    main()
