# P0 裁决报告：External stress-tests 存档-稿件矛盾（final_audit_20261030.md 条目 1.6 / M1）

- 裁决日期：2026-10-30
- 裁决人：agent-63（只读调查 + 本机 CPU 重跑验证；未触碰在跑 GPU 任务，全程 `CUDA_VISIBLE_DEVICES=""`）
- 结论一句话：**稿件数字是规范运行（10-30 远程规范运行，step16 锁定脚本 + 锁定种子 + 10-30 重生成 MCF7 预测输入），本地 10-29 的 `mcf7_protocol.json` 是被旧输入取代的 preliminary 运行，其"显著 +0.0509"结论不采信；norman_protocol.json 从未丢失，在远程工作目录中完整存在并已拉回入档。**

---

## 1. 矛盾双方

| 项 | 正文（main.tex） | 本地存档（10-29 23:25） |
|---|---|---|
| MCF7 paired target-z−sign | +0.0075 [-0.0264,+0.0412]（不显著） | +0.0509 [0.0025,0.0988]（显著） |
| MCF7 sign AUROC | 0.465 [0.353,0.582] | 0.486 [0.384,0.597] |
| MCF7 target-z AUROC | 0.515 [0.393,0.626] | 0.626 [0.515,0.738] |
| MCF7 in-situ / few-shot | 0.599 [0.438,0.740] / 0.576 [0.415,0.709] | 0.630 [0.464,0.761] / 0.624 [0.451,0.758] |
| Norman 七臂 + 配对 | 正文 14 个数字 | 本地无任何 json 存档（仅 png） |

## 2. 出处链（调查所得）

时间线（均为 2026-10，+0800）：

1. **10-21 00:28** 本地 `gse281048_testbed/outputs/gears_predictions.npz` 生成（旧版预测流水线，输入 md5 `681c281d…`）。
2. **10-29 09:57** 本地 step14 → `mcf7_validation.json`（同样基于旧 10-21 预测）。
3. **10-29 23:25** 本地 step16 首跑（run#1，旧 10-21 预测输入）→ 本地现存 `mcf7_protocol.json/.md/_scores.csv`（+0.0509 显著）。
4. **10-29 23:26** `prereg_external_norman_v1.md` 明示 "recorded before running step16 norman **on the new predictions**"——即当时已计划用重新生成的预测重跑。
5. **10-30 08:59–09:33** 远程服务器（xmu_ssh 助手所指 172.16.92.2，`/opt/research/gating_external/`）执行规范运行：
   - `gears_predict_all.py`（09:09 版）用 A_k562 模型权重重新生成 MCF7 预测（npz 09:19:05，md5 `d6fd2a0b…`，**与本地 10-21 npz 不同**）；
   - step15 生成 Norman 三 seed 预测（09:24–09:32）；
   - step16 `mcf7` 09:24:53 完成、`norman` 09:33:00 完成 → 远程 `outputs_external/{mcf7,norman}_protocol.{json,md,csv,png}` 全套存档。
   - 远程 step16 与本地 step16 **逐字节一致**（diff 为空）。
6. **10-30 11:16:08** 本地拉取两张 png（SFTP 批量 get，本地 mtime=拉取时刻）；**json/md/csv 未拉回** → 形成"norman json 缺失、mcf7 json 停留 10-29 旧值"的存档断裂。
7. 稿件图 `fig:external_protocols`（main.tex 行 733–734）引用的正是这两张 png；本地 png 与远程 png **md5 逐字节一致**（`76b463fb…` / `c0d238af…`）。v3.1 备份（10-30 15:06）已含全部稿件数字 → 数字于 10-30 上午至 15:06 之间据远程运行结果写入正文。
8. 该远程目录同时是 10-30 其余全部外部分析的实跑现场（variants/warmstart/reliability/quality_threshold/external_independent/sibling_transfer 等，日志齐全），与本地 `outputs_external` 后来拉回的存档一一对应——**远程运行就是论文的规范计算现场**。

## 3. 重跑验证（裁决执行）

以锁定脚本（本地 step16，与远程逐字节一致；常量 B=1000、boot_seed=20261029、OOF seed 0、few-shot seed 20261030）+ **拉回远程的规范输入**（MCF7 重生成 npz/gene_errors + Norman 三 seed npz），在本机 venv 全 CPU 重跑（沙盒 `outputs_external/_p0_rerun/sandbox/`，不改动本地原始输入）：

| 验证 | 结果 |
|---|---|
| MCF7 重跑 vs 远程 canonical json | 171 个数值叶全部一致，**最大差 6.2e-9**；paired +0.0075 [-0.0264,+0.0412] 逐位复现 |
| Norman 重跑 vs 远程 canonical json | 193 个数值叶全部一致，**最大差 8.7e-9**；paired +0.0305 [0.0022,0.0594] 逐位复现 |
| MCF7 重跑/远程 canonical vs 稿件 | 正文 684–706、摘要、Discussion、部署表所引全部数字**逐位一致** |
| 反向对照：本地 step16 + 本地 10-21 旧输入 | 精确复现 10-29 存档（点估计 4 位小数全同，逐基因分数差 ≤1e-13）→ 确认 10-29 存档=旧输入 preliminary 运行 |

输入差异根因：两份 MCF7 预测 npz 均为 A_k562 模型权重产物（**model.pt 本地↔远程 md5 相同** `35d24990…`），但预测流水线版本不同（10-21 本地旧版 vs 10-30 远程 09:09 版 `gears_predict_all.py`），且误差文件 `gears_gene_errors.csv` 随之不同（md5 不同）→ 特征与 median-split 标签均变 → 全部臂数字位移。**非转录笔误、非种子漂移，是输入版本漂移。**

## 4. 裁决

- **以稿件（= 远程 canonical = 本机重跑）为准。** 依据：锁定脚本可复现（本地重跑与 canonical 差 ≤1e-8）；稿件图与远程 png 字节一致；10-30 全部关联分析同现场同输入，内部自洽；prereg 指明以"new predictions"运行为准。
- **10-29 本地存档判定为 superseded**，其 +0.0509 显著结论不再代表论文分析；已按审计修复建议②"找回并入档"处理（见第 6 节）。

## 5. 稿件改动清单（main.tex，供主代理执行；P0 结论：正文无需改写"不显著"叙事）

必改（1 行）：

1. **行 1348**（`tab:deployment` 第 3 行）：`Few-shot AUROC $0.695$ on Norman and $0.578$ on MCF7` → `$0.576$`（即审计 M2；canonical few-shot 0.576 [0.415,0.709]，正文行 702 已是 0.576，部署表笔误）。

可选（2 行，求逐位严格时）：

2. **行 701**：δmax 与 source 合并写作 `0.523 [0.402,0.631]`；canonical 中 source 为 [0.402,0.631]、δmax 为 [0.402,0.629]，若严格要求逐臂 CI 归属可将该句改为分别标注。
3. **行 294 附近 Methods 段**："K562 models (seeds 0/1/2) provide predictions" 对 Norman 准确（3 seed），但 MCF7 实为单一 A_k562 消融模型预测（canonical json "seed predictions=1"，正文也未报 MCF7 ensemble disagreement，不矛盾）；可斟酌补半句限定，避免审稿人按 seeds 0/1/2 理解 MCF7。

P0 牵连面核查结论：行 684–706、摘要行 104–106、Introduction 行 167–168、Discussion 行 1319–1324、部署表行 1346–1349 的其余数字（Norman 七臂、MCF7 全臂、+0.0305/+0.0075 两处配对、0.465/0.705/0.599 等）**与 canonical 逐位一致，无需改动**；"MCF7 不显著 / 负向 stress test"叙事成立。

## 6. 存档侧已执行动作（本调查所做，全部可回退）

- `outputs_external/mcf7_protocol.{json,md}` + `_scores.csv`：10-29 run#1 原样改名为 `mcf7_protocol_localrun_20261029_superseded.*`（保留，未删）；正式名替换为远程 canonical 版（= 本机重跑版，差 ≤1e-8）。png 本就与 canonical 字节一致，未动。
- `outputs_external/norman_protocol.{json,md}` + `_scores.csv`：自远程拉回并入档（补齐缺失；png 已在库且字节一致）。
- 证据目录（供复核）：`outputs_external/_p0_remote/`（远程 canonical 四件套 + 远程 step16 副本）、`outputs_external/_p0_rerun/`（沙盒重跑脚本、规范输入、重跑全部产物、`norman_rerun.log`）。
- 远程侧仅执行只读操作（ls/cat/md5sum），未做任何修改。

## 7. 遗留事项（非 P0，供数据治理跟进）

1. `outputs_external/mcf7_validation.json/.md`（step14，10-29 09:57）同样基于 10-21 旧预测；稿件未引用其数字（审计未列差异），如需外部臂存档全集合输入一致，可用规范输入重跑 step14（本次未跑）。
2. 建议把"远程 `/opt/research/gating_external/` 为规范计算现场、本地为部分镜像"写入数据治理说明，避免再发生"只拉回 png 未拉回 json"。
3. 审计另列的两项存档侧非稿件问题（quality_threshold 0.0645 vs rpe1_stratified 0.0667 并存等）不在本次裁决范围，维持审计原建议。
