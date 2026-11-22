# 预注册：门控论文升层 G-A —— 多数据集扩展（工作量硬增量）

- 日期：2026-10-30
- 状态：**已锁定**（判据锁定于任何本扩展计算运行之前；smoke 只验证可执行性，不改判据）
- 起因：主 claim「门控信号跨数据集成立」目前只有 K562（AUROC 0.6556, n=217）与
  RPE1（0.7965, n=306）两个数据集就地支撑。升层 G-A 要求
  **检测有效性跨 5 数据集**（K562 / RPE1 / Norman / Adamson / Dixit），
  本预注册锁定后三个新数据集的就地检验协议与判定档线。
- 纪律：与 prereg_rebuild_v2.md 同——两档制，事先锁死，跑出来哪档写哪档；
  部分数据集失败本身就是发现，如实归因，不许和稀泥。

## 数据（判据锁定后不动）

- GEARS 内置公开扰动数据集三个，本机 GEARS 管线直接加载（`PertData.load`，
  同 replogle 口径）：
  - **Norman**（combinatorial；本机已有 `data/gears_demo/norman/`，
    91,205 cells / 5,045 genes / 284 conditions = 152 单扰动 + 131 组合 + ctrl）
  - **Adamson**（本机需下载，dataverse file 6154417）
  - **Dixit**（本机需下载，dataverse file 6154416）
- 下载纪律沿用 setup_log.md E9 教训：`curl -sSL --retry 8
  --retry-all-errors`（不带 -C），完成后 `testzip()` 验证，再解压；
  `mkdir data/gears_demo/<name>` 使 GEARS 跳过其不支持续传的下载封装。

## 管线（与 K562/RPE1 完全同构，代码复用不重写）

1. **训练 + 留出**（`step1_run.py`，env：`DATA=<name>
   OUT_DIR=outputs_multidataset/<name> SPLIT_SEED=0`）：
   - GEARS 训练 **3 epochs**、hidden 64、batch 32，GPU（本机 venv, torch cu128）
   - **基因级留出**：20% 单扰动基因（HOLDOUT_FRAC=0.2），划分种子固定 SPLIT_SEED=0
     （与主 runs 同考场纪律；模型种子只驱动初始化）
   - 产物：`predictions_seed{S}.csv`、`pred_expr_seed{S}.npz`、`gears_model_seed{S}`
2. **逐种子探针**（模型种子 **3/4/5**；探针实现复用 `step6_gating_probe.py`）：
   - 特征（模型侧，泄漏纪律与 K562 完全相同）：预测谱 8 项
     （delta_l2/delta_l1/delta_max/delta_top50_mean/n_de_025/frac_de_025/
     neg_frac/delta_entropy）+ log_n_cells；
     **strength、pearson_delta_all 排除**（真值派生）；in_pert_graph 常数自动剔除留痕
   - 标签（定义性选择，不动）：y = (1 − pearson_delta) > 全体留出扰动误差中位数
   - 探针：L2 logistic（C=1.0）+ fold 内标准化，分层 5-fold CV
     （shuffle），**探针种子 = 模型种子**（3/4/5 各自一套 CV 划分）
3. **执行**：detached 发射（pythonw + Start-Process 启动器，脱离会话存活），
   **GPU 串行队列**（单进程，逐数据集逐种子），预计 overnight；
   GBK 控制台纪律：不打印 emoji，日志 UTF-8 落盘。

## 判据（锁定）

- **数据集级**：某数据集 3 个种子探针 CV AUROC 的**均值 ≥ 0.6** →
  「门控信号在该数据集成立」（0.6 线沿用项目内约定，与 K562/RPE1 主 runs
  同线；逐种子值与极差一并报告，最小值留痕但不另设否决线）。
- **升层级**：3 个新数据集中 **≥2 个成立** → 「跨 5 数据集检测有效」进入论文主
  claim（K562/RPE1 两票 + 本扩展至少两票，共 ≥4/5）。
  **1 或 0 个成立** → 主 claim 保持 2 数据集口径，本扩展**如实报分层结论**：
  成立的数据集（信号普适的证据）与不成立的数据集（失败归因：n 太小 /
  扰动类型差异 / 标签中位数退化等，逐条分析）各自成文——部分数据集失败
  本身就是发现，不许藏。
- 无论哪档，逐数据集逐种子的全部数字进 `outputs_multidataset/` 与
  `outputs_rebuild/multidataset_summary.md`，阴性留档。

## 产物

- `outputs_multidataset/{norman,adamson,dixit}/`：每数据集
  `predictions_seed{3,4,5}.csv`、`pred_expr_seed{3,4,5}.npz`、
  `gating_probe_seed{3,4,5}.json`、`gating_probe_cv_seed{3,4,5}.png`、
  `gating_verdict.md`
- `outputs_rebuild/multidataset_summary.md`：跨数据集汇总 + 升层判定
- `outputs_multidataset/STATUS.md`：**断点手册**——逐单元（数据集 × 种子）增量
  落盘，每单元完成即更新；中断后重跑同一命令自动跳过已完成单元
  （幂等续跑，防超时丢进度）

## smoke 先行

- 三数据集各 1 epoch、种子 3 单跑（输出至 `outputs_multidataset_smoke/`，
  不污染正式产物），仅验证链路可执行 + 探针能出数；smoke 数字不进门控判定。

## smoke 结果（2026-10-30，正式全量执行前）

三数据集 smoke（各 1 epoch、种子 3）全部跑通链路，暴露一个口径问题（修正案 1）：

- norman：train 440s，留出 30 基因 / 扰动图内可预测 20；探针 n=20 可跑，1-epoch
  AUROC=0.5000（逐 fold 0.75/0.5/0.5/0.5/0.25，小样本噪声；最好单特征
  delta_max=0.80）——链路 OK，数字不作数。
- adamson：train 436s，留出 17 / 可预测 16；探针 n=16 可跑，1-epoch
  AUROC=0.5000（逐 fold 1/0/0/1/0.5；最好单特征 delta_l2=0.80）——链路 OK。
- dixit：train 260s 正常，留出 3（19 单扰动 × 20%）→ **5-fold 分层 CV 无定义**
  （sklearn ValueError: n_splits=5 > n_samples=3）→ 修正案 1。

## 修正案 1（2026-10-30，依 smoke 暴露的口径问题；判定阈值与两档制不变）

- **A1（探针 CV 口径的小样本回退与不可评估声明）**：留出协议保持锁定不动
  （20% 基因级留出、SPLIT_SEED=0）。smoke 实测 dixit 留出扰动数 n=3，
  5-fold 分层 CV 无定义（sklearn ValueError）；进一步验证发现 n=3 时标签
  正类仅 1 个，LOO 同样无从训练（留出正类后训练集单类）。裁定分两级，阈值不变：
  1. **少数类 < 2**（dixit 即此：pos=1/neg=2）→ 该数据集探针**不可评估**：
     AUROC 记 null，JSON 记 `evaluable=false` + 归因；数据集级按
     **未成立**计数，升层级锁定规则照算（3 集中 ≥2 成立才升层）——
     即 dixit 不可评估时须 norman、adamson 双双成立才可升层，等效更严。
  2. **少数类 ≥2 且 n < 10**（本次未触发，通用回退）→ 探针用**留一法（LOO）**：
     逐点以其余 n−1 个扰动训练同一 L2 logistic（fold 内标准化等泄漏纪律
     不变），汇总 LOO 决策分得 pooled AUROC；标注低效力。
  norman（n=20）与 adamson（n=16）维持 5-fold 分层 CV。
  归因原则：dixit 不可评估本身是发现（GEARS 内置 dixit 子集仅 19 单扰动，
  基因级留出协议下样本量退化），进论文方法学讨论，不藏。

## 偏离记录（执行后填）

- D1（发射机制，执行前置）：任务指定的 `pythonw + start //b` 在本机 Git Bash
  下因 MSYS 参数转换 + cmd 对 GBK bat 的路径解析双重问题不可用（中文项目路径
  写进 bat 即被 codepage 破坏）。改用**等价 detached 机制**：PowerShell
  `Start-Process pythonw.exe -WindowStyle Hidden -WorkingDirectory <项目>`，
  脱离会话存活语义相同（无窗口、独立进程树）；驱动自身把 stdout/stderr 镜像到
  `outputs_multidataset*/<mode>_queue.log`（pythonw 无控制台）。下载同理：
  ASCII bat + Start-Process -WorkingDirectory。
  不改变任何计算口径与判据，纯发射通道替换，如实留痕。
- D2（smoke 暴露的 dixit 口径）：见修正案 1。非偏离，修正案本身留痕。
- 其余执行偏离：（全量完成后填）

## 全量执行记录（2026-10-30 填）

- **发射**：2026-10-30 10:53:45，detached（PowerShell Start-Process pythonw，
  驱动 PID 66248；STATUS 与 queue 日志自落盘，断点续跑 = 重跑同一命令）。
  预计完成：训练 9 单元 ≈ 2.5–3 h（norman/adamson 各 ~55 min、
  dixit ~40 min，smoke 已付 pyg 缓存），加探针/汇总 ≈ 5 min，
  **ETA 2026-10-30 13:45–14:15**。
- **结果**：（跑完后填——逐数据集 3 种子 AUROC、数据集级判定、升层级判定）
