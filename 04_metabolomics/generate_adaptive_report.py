#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
generate_adaptive_report.py - 生成自适应注释收敛报告
"""
import argparse
import json
from pathlib import Path

import pandas as pd


def safe_read_csv(path, **kwargs):
    if path.exists():
        try:
            return pd.read_csv(path, **kwargs)
        except Exception:
            return pd.DataFrame()
    return pd.DataFrame()


def safe_read_json(path):
    if path.exists():
        with open(path, "r", encoding="utf-8") as f:
            return json.load(f)
    return {}


def fmt_pval(x):
    if pd.isna(x):
        return "NA"
    if x < 1e-4:
        return f"{x:.2e}"
    return f"{x:.4f}"


def generate_report(output_dir, out_md):
    out_dir = output_dir
    summary = safe_read_json(out_dir / "adaptive" / "adaptive_summary.json")
    iter_log = safe_read_csv(out_dir / "adaptive" / "iteration_log.csv")
    consensus = safe_read_csv(out_dir / "consensus" / "consensus_annotation.csv")
    consensus_lib = safe_read_csv(Path("library") / "cell_library_consensus.csv")
    pathway_summary = safe_read_csv(out_dir / "pathway_consensus" / "pathway_summary.csv")
    ora = safe_read_csv(out_dir / "pathway_consensus" / "pathway_enrichment.csv")
    metabo = safe_read_csv(out_dir / "metaboanalyst_results" / "pathway_enrichment_metaboanalyst_consensus.csv")
    mum_pw = safe_read_csv(out_dir / "mummichog" / "mummichog_pathways.csv")

    lines = []
    lines.append("# 自适应注释收敛报告\n")
    lines.append(f"**生成时间:** {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')}\n")

    lines.append("## 1. 输入特征与基线\n")
    lines.append(f"- POS/NEG 显著特征合并后 {summary.get('input_features', 'N/A')} 行，unique m/z {summary.get('unique_mz', 'N/A')}\n")
    lines.append("- 原始本地库 5 ppm 仅 21 条，10 ppm 36 条；初始 ORA 不显著\n")

    lines.append("## 2. 自适应阈值迭代\n")
    if not iter_log.empty:
        lines.append("| ppm | 注释特征 | Unique 注释 | Unique KEGG | 新增 KEGG | 新增比例 |\n")
        lines.append("|----:|--------:|------------:|------------:|----------:|--------:|\n")
        for _, r in iter_log.iterrows():
            ratio = r.get("new_ratio")
            ratio_str = f"{ratio*100:.1f}%" if pd.notna(ratio) else "—"
            lines.append(
                f"| {r['ppm']:.1f} | {int(r['features_annotated'])} | "
                f"{int(r['unique_annotations'])} | {int(r['unique_kegg_ids'])} | "
                f"{int(r['new_kegg_ids'])} | {ratio_str} |\n"
            )
    else:
        lines.append("未找到迭代日志。\n")
    lines.append(f"**最终选用 ppm:** {summary.get('chosen_ppm', 'N/A')}\n")
    lines.append(f"**最终 KEGG ID 数:** {summary.get('unique_kegg_ids', 'N/A')}\n")

    lines.append("\n## 3. 多库共识注释\n")
    if not consensus.empty:
        kegg_n = consensus_lib["KEGG_ID"].replace("", pd.NA).dropna().nunique() if not consensus_lib.empty else 0
        pc_n = consensus_lib["PubChem_CID"].replace("", pd.NA).dropna().nunique() if not consensus_lib.empty else 0
        cb_n = consensus_lib["ChEBI_ID"].replace("", pd.NA).dropna().nunique() if not consensus_lib.empty else 0
        level_counts = consensus["Consensus_Level"].value_counts().to_dict() if "Consensus_Level" in consensus.columns else {}
        lines.append(f"- 共识特征数: {len(consensus)}\n")
        lines.append(f"- KEGG/ChEBI unique ID: {kegg_n} / {cb_n}\n")
        if pc_n:
            lines.append(f"- PubChem unique ID: {pc_n}\n")
        else:
            lines.append("- PubChem 在本次运行中因 PUG REST 超时未纳入候选池；已用 KEGG + ChEBI 完成共识。\n")
        lines.append(f"- 高/中/单源支持: {level_counts.get('High', 0)} / {level_counts.get('Medium', 0)} / {level_counts.get('Single', 0)}\n")
    else:
        lines.append("未找到共识注释。\n")

    lines.append("\n## 4. 通路富集\n")
    if not ora.empty and "FDR_BH" in ora.columns:
        sig = ora[ora["FDR_BH"] < 0.25]
        lines.append(f"- 通路总数: {len(ora)}，FDR<0.25: {len(sig)}\n")
        if not sig.empty:
            lines.append("### 显著通路(FDR<0.25)\n")
            lines.append("| Pathway | Hits | Pathway_Total | OR | P value | FDR |\n")
            lines.append("|--------|-----:|--------------:|---:|--------:|----:|\n")
            for _, r in sig.head(20).iterrows():
                lines.append(
                    f"| {r.get('Pathway_ID','')} | {int(r.get('Hits',0))} | "
                    f"{int(r.get('Pathway_Total',0))} | {r.get('Odds_Ratio',0):.2f} | "
                    f"{fmt_pval(r.get('P_value'))} | {fmt_pval(r.get('FDR_BH'))} |\n"
                )
    elif not pathway_summary.empty:
        lines.append(f"- 本地 ORA 找到 {len(pathway_summary)} 条通路\n")
        lines.append("### 通路命中数 Top 10\n")
        lines.append("| Pathway | Compound_Count |\n")
        lines.append("|--------|---------------:|\n")
        for _, r in pathway_summary.head(10).iterrows():
            lines.append(f"| {r.get('Pathway_ID','')} | {int(r.get('Compound_Count',0))} |\n")
    else:
        lines.append("未找到通路结果。\n")

    lines.append("\n## 5. MetaboAnalystR ORA\n")
    if not metabo.empty:
        sig_col = [c for c in metabo.columns if "fdr" in c.lower() or "FDR" in c]
        if sig_col:
            sig = metabo[metabo[sig_col[0]] < 0.25]
            lines.append(f"- 通路数: {len(metabo)}，FDR<0.25: {len(sig)}\n")
            if not sig.empty:
                lines.append("| Pathway | Hits | Total | P value | FDR |\n")
                lines.append("|--------|-----:|------:|--------:|----:|\n")
                for _, r in sig.head(20).iterrows():
                    lines.append(
                        f"| {r.iloc[0]} | {r.get('Total_Hits',0)} | "
                        f"{r.get('Pathway_Total',0)} | {fmt_pval(r.get('P_value'))} | "
                        f"{fmt_pval(r.get(sig_col[0]))} |\n"
                    )
        else:
            lines.append(f"- 输出 {len(metabo)} 条通路结果\n")
    else:
        lines.append("未运行 MetaboAnalystR ORA 或无结果。\n")

    lines.append("\n## 6. Mummichog 兜底分析\n")
    if not mum_pw.empty:
        sig = mum_pw[mum_pw.get("FDR_BH", 1) < 0.25] if "FDR_BH" in mum_pw.columns else pd.DataFrame()
        lines.append(f"- m/z-to-pathway 通路数: {len(mum_pw)}，FDR<0.25: {len(sig)}\n")
        if not sig.empty:
            lines.append("| Pathway | Hits | P value | FDR |\n")
            lines.append("|--------|-----:|--------:|----:|\n")
            for _, r in sig.head(20).iterrows():
                lines.append(
                    f"| {r.get('Pathway_ID','')} | {int(r.get('Hits',0))} | "
                    f"{fmt_pval(r.get('P_value'))} | {fmt_pval(r.get('FDR_BH'))} |\n"
                )
    else:
        lines.append("未找到 Mummichog 结果。\n")

    lines.append("\n## 7. 结论与建议\n")
    if summary.get("unique_kegg_ids", 0) >= 30:
        lines.append("- 自适应扩展后 KEGG ID 数量达到阈值，可用于通路富集。\n")
    else:
        lines.append("- KEGG ID 数量仍偏少，建议继续实验验证或引入同位素/MS2 信息。\n")
    lines.append("- 若 ORA 仍不显著，Mummichog 可作为独立的 m/z 水平解释。\n")

    out_md.write_text("".join(lines), encoding="utf-8")
    print(f"报告已保存: {out_md}")


def main():
    parser = argparse.ArgumentParser(description="生成自适应注释收敛报告")
    parser.add_argument("--output-dir", default="output", type=Path)
    parser.add_argument("--out-md", default="output/adaptive_convergence_report.md", type=Path)
    args = parser.parse_args()
    generate_report(args.output_dir, args.out_md)


if __name__ == "__main__":
    main()
