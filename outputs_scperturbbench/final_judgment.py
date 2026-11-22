# -*- coding: utf-8 -*-
"""scPerturBench 口径对齐重跑：终版对照判定（主会话亲算，agent 六次超时后接管）。
锚点（STATUS 处置：仓库可核对值为主，锁定值并列）：norman top5000=0.5582, top100=0.7148;
rpe1 top5000=0.3479, top100=0.4642。判据（预注册锁定）：|我方3种子均值-论文值|<=0.05=一致。
"""
import json
import pandas as pd

base = r"E:\清华类脑计算研究中心\virtual_cell_constraint_mvp\outputs_scperturbbench"
df = pd.read_csv(base + r"\per_seed_results.csv")

ANCHORS = {
    ("norman", 5000): ("norman top-5000", 0.5582, 0.7395),
    ("norman", 100): ("norman top-100", 0.7148, 0.7395),
    ("rpe1", 5000): ("rpe1 top-5000", 0.3479, 0.4923),
    ("rpe1", 100): ("rpe1 top-100", 0.4642, 0.4923),
}

rows = []
for (ds, ndeg), (label, repo_val, locked_val) in ANCHORS.items():
    sub = df[(df.dataset == ds) & (df.numDEG == ndeg)]
    per_seed = sub.groupby("seed")["cor"].mean()
    ours_mean = float(per_seed.mean())
    diff = ours_mean - repo_val
    rows.append({
        "对照点": label,
        "仓库锚点": repo_val,
        "预注册锁定值(并列参考)": locked_val,
        "我方逐种子均值": {int(k): round(float(v), 4) for k, v in per_seed.items()},
        "我方3种子均值": round(ours_mean, 4),
        "差值": round(diff, 4),
        "判定(|差|<=0.05)": "一致" if abs(diff) <= 0.05 else "不一致",
    })

all_inconsistent = all(r["判定(|差|<=0.05)"] == "不一致" for r in rows)
all_consistent = all(r["判定(|差|<=0.05)"] == "一致" for r in rows)
verdict = "一致" if all_consistent else ("不一致但可归因" if all_inconsistent else "混合")
out = {
    "date": "2026-10-27",
    "tolerance": 0.05,
    "caveat_D1": "epochs=3（预注册资源锁定）vs 论文 myGears=15——最可能归因；判定须带此 caveat",
    "verdict": verdict,
    "points": rows,
    "n_groups": 6,
}
with open(base + r"\comparison_table.json", "w", encoding="utf-8") as f:
    json.dump(out, f, ensure_ascii=False, indent=1)
for r in rows:
    print(r["对照点"], "| repo", r["仓库锚点"], "| ours", r["我方3种子均值"],
          "| diff", r["差值"], "|", r["判定(|差|<=0.05)"])
print("VERDICT:", verdict)
