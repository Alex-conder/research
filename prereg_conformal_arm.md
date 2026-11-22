# 预注册：G-D conformal 弃权臂——扰动响应回归的 split-conformal 残差带（gating 升层）

- 日期：2026-10-30
- 状态：**已锁定**（全部判据锁定于任何 conformal 计算运行之前；本文件第 1 节空位复核先于实算完成并留痕）
- 起因：`meeting_prep/同行扫描_虚拟细胞可靠性_20261030.md` §3 第 5 条与 §5 第 1 条——分类侧已有 conformal 细胞注释，但检索未见"扰动响应回归 + conformal"工作，判为**空位**（非结论）；投稿前必须再复核一轮（本文件 §1 执行记录）。
- 纪律：判定标准事先锁死，跑出来哪档写哪档；失败也归档；如实报全部数字；GBK 控制台不打 emoji/非 ASCII。

## 1. 空位复核（投稿声明前的强制动作，先于实算）

- 复核日期：**2026-10-30**（与同行扫描同日，投稿前再查）。
- 检索渠道：Google Scholar / web 检索，检索词（原文记录）：
  1. `conformal prediction perturbation response gene expression regression`
  2. `conformalized quantile regression single-cell perturbation prediction uncertainty`
  3. `conformal "perturb-seq" OR "gene perturbation" prediction interval coverage selective prediction`
- 复核结论（锁定写入）：**窄口径空位仍成立，但邻接工作已出现，定位从"空白"降级为"早期未覆盖区"**：
  - 未见对**扰动响应回归/表达谱预测器**（GEARS 类 profile predictor）做 split-conformal 残差带并报告名义 vs 经验覆盖的同行工作；
  - 邻接 1（最接近，必须引用）：Asiaee et al., *Partial Causal Structure Learning for Valid Selective Conformal Inference under Interventions*（arXiv 2603.02204，v1 2026-03-02，v2 2026-07-09）——选择性 conformal 在干预（基因扰动实验）下选校准集，Replogle K562 CRISPRi 数据上评估；任务是"为标量目标选可交换校准干预"，**不是**对 profile 回归预测器做残差带；
  - 邻接 2：分类侧 conformal 细胞类型注释（Bioinformatics 2025, btaf521）与图约束 conformal 注释（arXiv 2410.23786）——分类非回归；
  - 邻接 3：`miguelrosell/adaptive-conformal-omics`（GitHub 2026-02）——单细胞 CQR，但为**合成数据**作品集项目，非真实扰动数据集、非扰动响应回归。
- 定位语句（锁定）：**"对扰动响应回归预测器给出有限样本覆盖保证的 conformal 弃权包装器"在检索可见范围内尚无直接同行；最近邻为干预下选择性 conformal（Asiaee et al. 2026），其不覆盖 profile 回归臂。**

## 2. 实验对象（已有产物，零训练）

| 臂 | 来源 | n 扰动 | 说明 |
|---|---|---|---|
| A_k562 | `outputs_ablation/A_k562/` | 217 | K562 分布内 |
| B_rpe1 | `outputs_ablation/B_rpe1/` | 306 | RPE1 分布内 |
| D_zeroshot | `outputs_ablation/D_zeroshot/` | 见产物 | K562→RPE1 零样本（OOD 主检验场） |
| norman_s4 | `outputs_multidataset/norman/` seed4 | 20 | 单基因留出，小样本 |
| adamson_s3 / s5 | `outputs_multidataset/adamson/` seed3/5 | 16/16 | 单基因留出，小样本 |

- 真值口径（锁定）：A/B/D 臂直接用 npz 内 `true_expr`（每扰动条件均值，与 A_k562_meta eval_protocol 同口径）；norman/adamson 由 `data/gears_demo/{ds}/perturb_processed.h5ad` 重算条件均值（`condition` 去 "ctrl" 成分后精确匹配扰动基因名），ctrl 取全部 ctrl 细胞均值。**重算真值必须与 `predictions_seed{seed}.csv` 的 `pearson_delta` 对账**（|Δr| 中位数量级检查），对账不过则该臂标 INVALID 不出数字。
- 误差口径（锁定）：主标签 `error = 1 - pearson_delta`（top-50 |true delta| 基因上的 Pearson，step6 探针同约定）；分层用同一 error。

## 3. Conformal 协议（锁定）

- **单元**：(扰动, 基因) 对上的 delta 残差 `r_{i,g} = true_expr[i,g] - pred_expr[i,g]`（ctrl 相消；log-normalized 空间）。
- **校准/评估划分（锁定）**：每臂内按扰动名确定性散列 `sha256("gating-GD-v1|" + 臂名 + "|" + 扰动名) % 2`，0→校准、1→评估，50/50。**必须按扰动划分**（同扰动内基因残差相关，按单元划分会破坏交换性——写入协议并遵守）。
- **主法：自适应 split-conformal**（标准化残差带）：
  - 扰动级尺度 `sigma_i = ||delta_pred_i||_2 / sqrt(G)`（纯模型侧特征，无真值泄漏）；
  - 符合性分数 `s_{i,g} = |r_{i,g}| / sigma_i`（在校准半估计）；
  - 名义 90% 区间：评估扰动 i 的基因 g 上半宽 `q90 * sigma_i`，
    `q90 = ceil((n_cal_units + 1) * 0.9) / n_cal_units` 阶的校准分数经验分位（有限样本修正）。
  - 次法（留痕对照）：非自适应常带宽 split-conformal（`sigma_i ≡ 1`），仅报覆盖与宽度，不进判定。
- **保证的层数（如实声明）**：交换性在"评估扰动"层成立（扰动间 i.i.d. 假设），单元级边际覆盖是有限样本保证；同扰动内基因残差相关 → 另报**扰动级覆盖**（每扰动被覆盖基因比例，在评估半上取中位/最小），它不携带保证，只作诊断。

## 4. 指标（全部预锁）

1. **名义 90% vs 经验覆盖**（评估半，单元级边际）；二项扰动级视角：被覆盖比例 ≥90% 的评估扰动占比。
2. **区间宽度**：评估半平均半宽、按 `sigma_i` 分位点报告。
3. **按误差分层的覆盖均衡性**：评估扰动按 error 三分层（低/中/高），各层报单元覆盖 + 平均半宽。**关键诊断：高误差（OOD 样）层覆盖是否仍接近名义**。
4. **Risk-coverage 并排臂**（与符号门/探针对照的价值所在）：
   - 弃权分数（越大越先弃）：`conformal_width` = `sigma_i`（常数 q 下单调于 sigma）；
   - 对照分数（均只用校准半拟合、评估半出分，与 conformal 同预算）：`probe_cal`（step6 九谱特征 + log_n_cells 的 L2 logistic，校准半拟合决策函数）、`sign_gate_cal`（同拟合系数符号向量点评估半 z 分数）；
   - 曲线：按分数降序逐个弃权，x=保留率，y=保留集 mean error（=1-pearson_delta）；报 AURC、保留 80%/50% 时 error 与相对全量降幅、随机基线。
5. **判定（锁）**：主判定只看一条——评估半单元级经验覆盖 ∈ **[85%, 95%]**（名义 90%）→ **conformal 校准成立（有限样本保证兑现）**；<85% 或 >95% → 如实报偏差方向与幅度，并给出诊断（分层表 + 宽度表），不修饰。小样本臂（norman/adamson）覆盖区间天然宽，照实报并报二项 CI，不调阈值。
- 与探针对照的定位语句（锁定）：探针给**无保证的排序**（CV AUROC），conformal 给**有保证的弃权**（覆盖率可核验）——两者并排即 G-D 臂的贡献。

## 5. 产物与复跑

- `inhibition_multiscale/outputs_rebuild/conformal_arm.json`（全数字）、`conformal_arm.png`（覆盖/宽度/分层 + risk-coverage 并排）、`conformal_arm_verdict.md`（判定 + 一句话定位）、`conformal_arm.py`（可复跑脚本）。
- 运行环境：本机 venv CPU，不碰在跑训练。
