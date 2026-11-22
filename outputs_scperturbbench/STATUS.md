# STATUS.md — scPerturBench 口径重跑 · 断点状态

> 纪律文件：每个会话开始先读本文件再动手；每次会话结束更新"当前状态"段。
> 预注册：`paper_reproduction/prereg_scperturbbench.md`（已锁定）。

## 断点盘点（2026-10-24，会话三开始时）

**已完成：**
- 论文仓库已克隆到 `paper_reproduction/repos/scPerturBench/`（main 分支
  tarball，commit `63fd9b95071b894c4f39c8f7eb2bd42ca9a70a77`，2026-10-22，
  已用 jsdelivr 按该 commit 取 README.md 逐字节比对验证 tarball 一致性；
  GitHub 直连被墙，走 gh-proxy.com 镜像）。
- 关键代码已读：`Perturbation_generalization/Genetic/myGears.py`（训练+split）、
  `../calPerformance_genetic.py`（评测）。
- 论文结果表结构已摸清：`Results/Genetic_perturbation/` 下
  single/combo × top5000/top100 四张 CSV，行为逐扰动（op 带 _1/_2/_3 种子
  后缀，3 种子聚合），主指标列 `cor`（pearson_distance）。
- 机器状态：无 python 进程、GPU 空闲（46MiB/8GB）、无孤儿训练。

**未完成（前两个会话 6 小时零训练产出的原因，如实记录）：**
1. 会话一被打断时只完成了上面"已完成"项，尚未写任何运行代码。
2. 会话二继续读结果表（single CSV 的 GEARS 均值表），同样未写运行代码、
   未起训练。**根因：探索期过长，没有先落 STATUS/先起 smoke。**
3. 对照数字提取未收尾：Norman 在 combo CSV（未算均值）；RPE1 single
   top5000 GEARS 均值 = 0.3479、top100 = 0.4642（与预注册锁定的
   "RPE1 最高 0.4923" 的对应关系待确认——0.4923 可能来自 top100 或
   其他变体行，跑前必须逐位核对）。

## 阻塞诊断（会话三，15 分钟内）

**诊断结论（决策支点 (a) 成立，阻塞均为机械性可修）：**
1. **数据**：论文预处理文件（`filter_hvg5000_logNor.h5ad`）不在 repo 内；
   Figshare API 本机 403（被墙）；本机已有 GEARS 官方处理版 Norman（91205×5045）
   与 replogle_rpe1_essential（162733×5000），condition 格式与 myGears 兼容——
   **走本机复用**（预注册允许，记偏离）。
2. **依赖**：`myGears.py` import 的 `myUtil` 未随仓库发布（仅有内容等价的
   `myUtil1.py`）→ 不 import 论文脚本，改为逐行复刻其 split/评测语义；
   `calPerformance` 需要的 pertpy 1.4.0 **可装**（dry-run 不触碰 pandas<3 /
   torch 等锁定依赖），已安装。pertpy pearson_distance 定义确认：
   `perf = 1 - pearsonr(组均值)`，论文结果表 `cor` 列为 r=1-perf
   （依据方法排名方向判定，备注进预注册）。
3. **Windows 兼容**：myGears 硬编码 `/home/project/...` 路径与 `cuda:1`——
   机械适配，不影响口径。
4. **Podman 备选评估**：镜像 40GB（Zenodo 10 分卷/Onedrive），含全部依赖与
   数据，md5 清单在 repo；Windows 需 Podman Desktop + WSL2，WSL2 默认不直通
   本机 8GB GPU（需 WSL 内装 CUDA，训练意义丧失）。**结论：仅作为数据/依赖
   兜底，本机主路径可行时不启用**；若将来走服务器复现，Podman 是首选。
5. **对照锚点核对（材料级发现）**：仓库 Results 表可逐位核对的 GEARS 均值为
   Norman top5000 = **0.5582**（与预注册锁定一致）；但 Norman "top-50" 锁定值
   0.7395 与仓库 top100 表 0.7148 不符（名称与值都对不上）；RPE1 锁定
   "最高 0.4923" 与仓库 single 表 top5000=0.3479 / top100=0.4642 均不符。
   0.7395/0.4923 在仓库中仅作为无关表的巧合值出现。**处置**：判定以仓库可
   核对锚点为主（0.5582 / 0.7148 / 0.4642），锁定值并列报告（已写入预注册备注）。

## 决策记录（30 分钟时点）

**选 (a)**：阻塞可修 → 修。epochs=3（按预注册资源锁定与上级指令；论文
myGears 为 15 epochs，记为第一偏离，也是最可能影响对照的因素）。

## 当前状态

**2026-10-24 22:00（会话四结束状态）——6 组已全部排入队列并起活。**

- 脚本：
  - 训练+评测一体 `paper_reproduction/code/run_spb_train.py`
    （myGears split 复刻 + calPerformance pearson_distance 复刻 + pertpy；
    逐种子把 result.h5ad 与 cor 行落盘 outputs_scperturbbench/）
  - 串行调度器 `paper_reproduction/code/orchestrate_spb.py`
    （subprocess 逐组排队，单组失败不阻断后续；pid 记录在 progress.json）
- 发射（脱离 agent 生命周期）：
  `cmd //c start //b //min pythonw.exe orchestrator_spb.py`；
  pythonw 无控制台，agent 被 kill 不影响训练。
- 进程：pythonw 调度器 + 工作进程在跑（tasklist 查 pythonw.exe）；
  当前在跑 **norman_seed3**（new_data_process 阶段，进度 95%，RAM ~4.8GB）。
- 日志：`outputs_scperturbbench/train_{norman,rpe1}_seed{3,4,5}.log`（逐组）、
  `orchestrator.log`（调度）、`progress.json`（增量状态）。
- 排队顺序（串行，8GB 显存单进程纪律）：
  norman_s3 -> norman_s4 -> norman_s5 -> rpe1_s3 -> rpe1_s4 -> rpe1_s5。
- 预计完成：norman ~50min/组、rpe1 ~70min/组 + 数据预处理 ~20min，
  全部约 6 小时，即 **2026-10-25 04:00 前后**。
- 事故记录：首次发射后 GEARS 走 dataverse 慢速下载 gene2go（12kB/s）——
  已杀进程、预置本机 pkl 副本（md5 校验一致）到 work 数据目录，
  脚本 ensure_processed 已加自动预置逻辑，重跑不复发。
- 评测对照脚本 `compare_spb.py` 待写（下个会话：训练完成后逐对照点判定）。

**下个会话动作**：读 progress.json → 若 6 组全 done：写 compare_spb.py 出对照表、
填预注册偏离记录、写 comparison_verdict.md、更新 deliverables；若有 failed：
查对应 train 日志归因，按预注册降级条款处置。


## 当前状态（2026-10-27 收档）

**任务完成。** 6 组训练全齐，终版判定落盘：`comparison_table.json`（4 对照点）、
`comparison_verdict.md`（四节）、`final_judgment.py`（可复跑，0.4s）。
判定：不一致（4/4 超 ±0.05 锁定容差）；Norman 两点归因于 3-vs-15 epochs；
RPE1 top-5000 正向超额 +0.150 未归因（数据预处理文件差异为候选，硬约束挂账）。
deliverables 对照文档与预注册偏离（D1-D5）已同步。无未决动作。
