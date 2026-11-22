# Smoke 报告：重构预注册 R1 / R2（prereg_rebuild_v2.md）

> **纪律声明（先读）**：本报告所有 smoke 数字**不带入锁定判据**。锁定判据
> （R1 的 ≥4/12、R2 的 ≥8/12 两档制）只对全量锁定执行生效；smoke 的唯一目的
> 是验证两条代码路径可端到端执行、产出 schema 正确、全量耗时外推合理。
> smoke 组合仅为 GEARS-K562 + Geneformer-K562（一表达谱类 + 一 embedding 类）。

- 日期：2026-10-24
- 代码：`step10_signature_resolution.py`（R1）、`step10_two_layer_gate.py`（R2）
- smoke 产物：`outputs_rebuild/smoke/signature_resolution.smoke.{json,png}`、
  `outputs_rebuild/smoke/two_layer_gate.smoke.{json,png}`
- 运行环境：本机 venv（`venv/Scripts/python.exe`），CPU 单进程
  （`OMP_NUM_THREADS=1`），GBK 控制台纪律（stdout/stderr reconfigure utf-8
  errors=replace，无 emoji）

## 1. 可执行性结论

| 实验 | smoke 规格 | 实测墙钟 | 结论 |
|---|---|---|---|
| R1 签名解析度 | B=50，n∈{15,60,217}，单种子，2 组合 | 3.6–5.0 s（≤10 min 限制） | **可执行** |
| R2 两层门对决 | B=100，coverage∈{1.0,0.75,0.5}，2 组合 | 9.4–11.0 s（≤10 min 限制） | **可执行** |

两条路径均一次跑通（R2 首轮发现并修复一个臂间配对键名 bug——
`symbol_minus_two_layer` 与预期间 `two_layer_minus_symbol` 的命名错位，
已改为显式有序配对）。产出 JSON schema 与锁定产物
（`outputs_rebuild/signature_resolution.json`、`two_layer_gate.json`）同构：
判定字段齐备（`verdict.tier` 在 smoke 模式固定为
`smoke_only_not_adjudicated`），smoke 与全量只差 `spec` 数值与输出目录。

## 2. Smoke 数字（仅留痕，不判定）

### R1（B=50，单种子；稳定特征 = 95% bootstrap CI 不含零）

| 组合 | n | 特征数 | 全量 n 稳定特征 | 曲线（n : 稳定特征均值） |
|---|---|---|---|---|
| GEARS / K562 | 217 | 9 | 2（delta_max, delta_entropy） | 15:1.0，60:3.0，217:2.0 |
| Geneformer / K562 | 217 | 4（嵌入几何） | 1（delta_emb_norm） | 15:0.0，60:1.0，217:1.0 |

### R2（B=100；带内 precision = step6b 同口径；risk = 1 − precision）

| 组合 | 臂 | P@cov 1.00 / 0.75 / 0.50 | two_layer − symbol（跨 coverage 均值配对差，95% CI） |
|---|---|---|---|
| GEARS / K562 | in_situ / symbol / two_layer | .498/.571/.602；.498/.552/.593；.498/.534/.556 | −0.0019 [−0.0185, +0.0154]，含 0 |
| Geneformer / K562 | 同左 | .498/.472/.481；.498/.485/.519；.498/.497/.500 | +0.0016 [−0.0232, +0.0316]，含 0 |

诊断列（非锁定）：spearman(symbol, two_layer) = 0.81 / 0.47；最大秩位移
131 / 162（n=217）——见风险 R-1。

## 3. 全量耗时外推（按 smoke 实测外推；本机另跑 GEARS 训练，留 2× 余量）

### R1 全量（B=500 × n 档 {15,30,60,120,全量} × 3 种子 × 12 组合）

- 实测：就地探针单次重拟合 ≈ 1.1–1.9 ms（n=217，9 或 4 特征）
- 每组合重拟合数：500（全量）+ 4 子采样档 × 3 种子 × 500 = **6500 次**
  （scf_k562 n=45 时网格收缩为 {15,30,全量} → 3500 次）
- 每组合计算 ≈ 13 s；12 组合 ≈ 2.5–3 min；加载（2 个 GEARS h5ad ≈ 5 s/个、
  scGenePT/scFoundation npz、6 份 wave csv）≈ 1 min
- **预计全量墙钟 ≈ 4–6 min（保守 ≤15 min）**

### R2 全量（B=1000 × 11 coverage 档 {1.0→0.5 步 0.05} × 12 组合）

- 实测：每次配对 bootstrap 重采样（3 臂全部重算：5 折就地探针 + 符号打分 +
  5 折 isotonic）≈ 8.9–12.6 ms
- B=1000 → 每组合 ≈ 9–13 s；12 组合 ≈ 2–3 min；加载（表达谱类组合各需
  目标 + 同引擎兄弟数据集，GEARS h5ad 重复加载）≈ 0.5–1 min
- **预计全量墙钟 ≈ 3–6 min（保守 ≤15 min）**

两实验全量合计 ≈ 10–30 min，单进程 CPU 无压力，不需要并行化（保持对
GEARS 训练的 CPU 克制）。

## 4. 风险清单（全量执行前需主代理裁定的按序排列）

- **R-1（最高优先，R2 结构性）**：锁定文本的两层臂 = 符号层 + isotonic 就地
  校准（5-fold OOF）。isotonic 单调，同一符号分经逐折校准后，全局排序变化
  只来自**跨折非单调性**（各折校准函数不同），即折划分噪声而非真实增益
  信号。smoke 实测两臂差均值 |Δ|≈0.002、CI 含 0，且符号分与两层分秩相关
  仅 0.47–0.81、最大秩位移达 131–162/217。若按锁定判据逐组合比 (c) vs (b)
  的秩口径 precision 差，结果近似掷硬币，两档判定可能被噪声主导。裁定选项：
  (a) 确认"两层门"的第二层是否实为就地探针分数的 isotonic 校准（与符号层
  组合而非校准符号层本身）；(b) 或改用校准敏感指标（Brier/log-loss）。
  **全量跑之前必须裁定，否则 R2 判定无效。**
- **R-2**：scFoundation 无 RPE1 产物 → 最多 11/12 组合；且 scf_k562 无同引擎
  兄弟数据集，**符号臂不可得**（代码已置 `symbol arm unavailable` 并跳过
  b/c 臂）。锁定判据写的是 /12 分母，需裁定分母口径（11？scf_k562 出局后
  R2 实际可评组合 = 11）。代码不静默改分母，verdict 块带
  `denominator_note`。
- **R-3（R1 特征口径）**：embedding 类特征取 pred_emb_norm / delta_emb_norm /
  log_n_cells / target_in_input_space（4 个）；`true_emb_norm` 为真值侧量，
  按泄漏纪律排除。锁定文本"emb_norm 等"未指明 pred/true 侧，本实现选择已
  在 JSON `features_note` 留痕；若全量前改主意须先改预注册。
- **R-4（操作化口径，低风险但需知情）**：R1"解析首个稳定签名所需扰动数"实现为
  曲线上首个均值稳定特征数 ≥1（另报 ≥2 档）的**网格点** n，不做插值；
  R2 每组合判定实现为"跨 11 个 coverage 档的配对差均值"的 CI 是否含 0
  （配对 bootstrap 同重采样跨档配对）。两者均为锁定文本的合理操作化，
  非文本明写，留痕待认。
- **R-5（加载耗时一次性）**：GEARS 两数据集的 h5ad 读取 ≈ 5 s/个，R2 中
  目标+兄弟各读一次（无缓存，K562/RPE1 两组合会各读两次 h5ad）。
  全量总加载 <2 min，不优化；如需复跑提速可加 lru_cache，但会改变
  留痕口径，先不动。
- **R-6**：scf_k562 n=45，R1 全量 bootstrap（B=500）在该组合上 CI 必然极宽，
  判定时"≥2 稳定特征"大概率不过——这是数据事实，如实呈现，不隐藏。

## 5. 复跑命令

```
PYTHONIOENCODING=utf-8 OMP_NUM_THREADS=1 venv/Scripts/python.exe step10_signature_resolution.py --smoke
PYTHONIOENCODING=utf-8 OMP_NUM_THREADS=1 venv/Scripts/python.exe step10_two_layer_gate.py --smoke
# 全量（待主代理确认后）：去掉 --smoke，产物落 outputs_rebuild/
```

随机源（专属，已锁在代码常量）：R1 bootstrap 20261024、R1 子采样 20261025、
R2 bootstrap 20261026；OOF 划分种子 0（step6 同口径）。
