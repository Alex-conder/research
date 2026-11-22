# 多模型波（scGPT / Geneformer / UCE）ah-server 准备状态盘点

> 盘点日期：2026-10-24。方式：ssh ah-server 只读盘点（docker ps / pip list / du / cat），
> 未安装任何东西、未下载任何大文件，动物房四个生产容器（backend / backend_staging /
> mysql / ocr）零接触。本文件是 `virtual_cell_constraint_mvp` 门控探针框架多模型波的
> 接入准备清单。

## 0. 结论（TL;DR）

**准备状态远超任务假设：三个引擎不但环境/权重/脚本全部就位，wave 推理已经全部跑完。**

- 三容器（research_scgpt / research_geneformer / research_uce，2026-10-16 建）健在，
  pip 环境经 `pip list` 实测核对，与各自 RUNBOOK 记录一致。
- 权重全部经 hf-mirror 拉齐（scGPT 205MB、Geneformer 417MB、UCE 8.5GB）。
- wave 产物（per-cell 嵌入，供门控探针用）6 份齐全，2026-10-19 至 10-21 完成：
  - K562：217 扰动 × 3 引擎；RPE1：306 扰动 × 3 引擎。
  - 位置：`ah-server:/opt/research/outputs_crossmodel/{scgpt,geneformer,uce}/`
  - 每个引擎 `_labels.csv` 已带 `pred_emb_norm / true_emb_norm / delta_emb_norm / err_emb_dist`
    探针列；`pearson_delta` 列全 NaN 为设计如此（输出是嵌入不是表达谱，见 §5）。
- **剩余准备项只有一件**：6 份产物（合计约 330MB）尚未拉回本地
  `outputs_crossmodel/wave/`（本地 wave/ 目前只有 3 个推理脚本和 fetch_b64.py）。
  拉回后用 `verify_npz.py` 口径校验即可进入探针集成。
- 磁盘/环境**无硬阻塞**：系统盘 694G 可用、/home 55T 可用、load 0、内存 840G 可用。

## 1. ah-server 现状盘点（只读）

### 1.1 容器快照（docker ps -a，2026-10-24）

| 容器 | 镜像 | 状态 | 资源配额 |
|---|---|---|---|
| research_scgenept | python:3.10-slim-bookworm | Up 8 天 | cpus 32 / 128g（GEARS/scGenePT K562 线） |
| research_scgenept_rpe1 | research:scgenept-py310-v2 | Up 8 天 | （RPE1 线） |
| research_scfoundation | research_scf_env:latest | Up 8 天 | 10.7GB 镜像 |
| research_gears | python:3.10-slim-bookworm | Up 7 天 | |
| **research_scgpt** | research:scgenept-py310-v2 | Up 8 天 | cpus 64，cpuset socket0（0-63,128-191），128g |
| **research_geneformer** | python:3.10-slim-bookworm | Up 8 天 | cpus 64，cpuset socket1（64-127,192-255），128g |
| **research_uce** | python:3.10-slim-bookworm | Up 8 天 | cpus 48，floating 避让，128g |
| animal_house_backend / _staging / mysql / ocr | 自有镜像 | Up | **生产容器，全程未触碰** |

全部自建容器 `--restart no`、显式配额，符合 RUNBOOK 纪律。

### 1.2 空间与负载

- `/opt/research` 共 238G：external 209G、ablation 13G、models 11G（UCE 占 8.5G）、data 5.4G、
  outputs_crossmodel 312M。
- 系统盘（/opt 所在）1TB 已用 330G，**可用 694G**；/home 62T 已用 7.3T，可用 55T。
- load 0.00 / 256 线程；机主 llama 服务常驻 ~32GB，881GB 内存可用 840G。

### 1.3 RUNBOOK 中与装新环境相关的坑（汇总三份引擎 RUNBOOK + 主 RUNBOOK）

1. **Python 必须 3.10**：scgpt 0.2.1 要求 <3.11（3.11 只能装 0.2.2+，未验证）。三容器均 3.10。
2. **torch 必须先装 CPU 版**：`pip install torch==2.1.2+cpu --index-url https://download.pytorch.org/whl/cpu`，
   否则默认拉 ~2.5GB CUDA 版。三容器实测均为 torch 2.1.2+cpu。
3. **numpy 钉死 1.26.4**：scgpt 依赖链会把 numpy 顶到 2.x，torch 2.1.2 是 numpy1 ABI，混用直接炸。
4. **hf-mirror 是权重下载主通道**（服务器直连 S3 仅 ~200KB/s）：
   - scGenePT / scFoundation / scGPT-human / Geneformer 均走 `hf-mirror.com/.../resolve/main/...`，
     实测 8.4MB/s，秒到分钟级。
   - **UCE 例外**：官方 figshare 在本网络一律 HTTP 202 空响应（区域限制），改用非官方镜像
     `hf-mirror.com/lza1/uce_hf`（model_files 原样拷贝官方文件）。
5. **github.com 不稳定**（HTTP/2 断流）：repo 一律本机下 tarball 后 SFTP 上传，commit 写进
   `GIT_COMMIT`（uce_repo/UCE-main 即此法）。
6. **大文件走 SFTP**：本机↔服务器 ~4-8MB/s；>5GB 的服务器直拉不划算。
7. **docker exec 路径是容器内路径**（/research/... = 宿主机 /opt/research/...）。
8. **pkill -f 会匹配自身命令行把 SSH 会话杀死**：杀任务用
   `docker top <容器> | grep <脚本> | awk '{print $2}' | xargs -r kill`。
9. **BLAS 线程必须压到容器配额内**：与主任务并行时用 OMP/MKL/OPENBLAS_NUM_THREADS 显式设，
   否则 numpy C 扩展在 pids-limit 下导入即崩（`PyCapsule_Import` 假象报错）。
10. **cpuset 分工**：scGPT 钉 socket0、Geneformer 钉 socket1、UCE floating；实测高负载下
    32 线程比 64 线程单扰动更快（scGenePT/scGPT 同款现象，NUMA/超线程争用）。

## 2. 三引擎 × 接入成本矩阵

时长口径说明：任务要求按 scGenePT K562 8.4h/217 扰动外推。wave 协议与首轮 scGenePT 协议
不完全相同（wave 对每扰动出 ctrl 池 32 + true 细胞 ≤32 的 per-cell 嵌入，细胞数约翻倍，
且保留 per-cell 输出），因此"外推估计"仅供参考，**实测数字已存在，直接采用实测**。

### 2.1 scGPT whole-human

- **容器方案（实际）**：复用 `research:scgenept-py310-v2` 镜像新建 research_scgpt
  （scGenePT 修复后的 env：torch 2.1.2+cpu / scgpt 0.2.1 / anndata 0.10.9 / numpy 1.26.4），
  cpus 64 钉 socket0。**禁止在容器里升级 torch/numpy**（ABI 事故史）。
- **依赖清单（pip list 实测）**：torch 2.1.2+cpu、torchtext 0.16.2+cpu、scgpt 0.2.1、
  numpy 1.26.4、scipy 1.13.1、pandas 2.2.3、anndata 0.10.9、scanpy 1.10.4、
  scvi-tools 0.20.3、torch_geometric 2.5.3、local-attention 1.11.2、hyper-connections 0.4.11。
- **权重获取（已完成）**：`hf-mirror.com/perturblab/scgpt-human` →
  best_model.pt 205MB + vocab.json + args.json（存 `/opt/research/models/scgpt_whole_human/`，
  sha256 6cb5d451…已记入 meta）。官方 whole-human 镜像，60697 token / 12 层 / 51.3M 参数。
- **推理时长**：外推估计 ~8h/数据集（与 scGenePT 同构、参数量更小）；**实测**：
  K562 62812s@32线程 ≈ **17.4h**（290s/扰动），RPE1 59409s ≈ **16.5h**。
- **风险**：
  - checkpoint 是 fast_transformer 保存的，CPU 路径需 `Wqkv → in_proj_` 重命名 + strict=False
    （scgpt 0.2.1 的键前缀是 `.self_attn.Wqkv.`，脚本已修）；
  - 60697-vocab MLM 头极大，wave 脚本已改用 `_encode + _get_cell_emb_from_layer` 跳过
    （与全 forward 逐元素验证一致）；
  - 首轮 expression 口径（mask 协议）下预测谱整体≈ctrl，pearson_delta 主要反映扰动基因
    单点——wave 已改 embedding 口径规避，但若上游要表达谱需另补协议。

### 2.2 Geneformer V2-104M

- **容器方案（实际）**：python:3.10-slim-bookworm 新建 research_geneformer，cpus 64 钉
  socket1，128g。**必须先 `pip install --upgrade pip`** 再装 torch（旧 pip 会触发
  typing-extensions sdist 构建失败）。
- **依赖清单（pip list 实测）**：torch 2.1.2+cpu、transformers 4.44.2、numpy 1.26.4、
  scipy 1.13.1、pandas 2.2.3、anndata 0.10.9、safetensors 0.8.0、huggingface_hub 0.36.2。
- **权重获取（已完成）**：`hf-mirror.com/ctheodoris/Geneformer` → Geneformer-V2-104M/
  （config.json + generation_config.json + model.safetensors 417MB）+ 2 个 token 字典 pkl
  （存 `/opt/research/models/geneformer/`，sha256 fff5cba2…）。Vocab 20275，hidden 768，12 层。
- **推理时长**：外推估计 2-4h/数据集（序列短 ≤2042、无大输出头）；**实测**：
  K562 8918s@48线程 ≈ **2.5h**（41s/扰动），RPE1 16611s ≈ **4.6h**。**三引擎中最快**。
- **风险**：
  - 特殊 token 是 `<pad>=0 <mask>=1 <cls>=2 <eos>=3`，**没有 `<sep>`**（踩过）；
  - pooler 未初始化警告无害（细胞嵌入 = 基因 token last_hidden_state 均值，不用 pooler）；
  - 217 扰动中仅 73 个可真实 in-silico 删除（134 目标不在 5000 输入基因内，另 6 个无
    Ensembl/token 映射），其余为退化基线（labels 有 target_deleted 列留痕）；
  - **输出只有 cell embedding，没有表达谱**——pearson_delta 口径不适用（见 §5）。

### 2.3 UCE 33L

- **容器方案（实际）**：python:3.10-slim-bookworm 新建 research_uce，cpus 48 floating，
  128g；代码 `/opt/research/uce_repo/UCE-main`（GitHub tarball SFTP 上传）。
- **依赖清单（pip list 实测）**：torch 2.1.2+cpu、numpy 1.26.4、pandas 2.2.3、
  scipy 1.13.1、anndata 0.10.9、scanpy 1.10.4、tqdm 4.70.1（依赖最简）。
- **权重获取（已完成）**：官方 figshare 被 HTTP 202 卡死 → `hf-mirror.com/lza1/uce_hf`：
  33l_8ep_1024t_1280.torch 5.69GB + all_tokens.torch 2.98GB + human protein embeddings
  411MB + species_chrom/offsets（合计 8.5GB，存 `/opt/research/models/uce/model_files/`，
  sha256 3f458726…）。
- **推理时长**：外推估计 8-16h/数据集（33 层 × d1280，前向最重）；**实测**：
  K562 40066s@48线程 ≈ **11.1h**（185s/扰动），RPE1 68513s ≈ **19h**。
- **风险**：
  - **权重最大（8.5GB）+ 前向最重**，是成本最大项；
  - `data_utils.py` 模块级读 CWD 相对路径 `./model_files/...` → 必须在 /research 下放
    model_files 目录（只有表头 14 字节的小 csv 即可）；
  - 官方采样用无 seed 的 np.random → wave 脚本已固定 per-cell seed 保证可复现；
  - checkpoint 已含训练好的 pe_embedding（145469×5120），与 all_tokens 二选一由官方
    加载逻辑决定，无需额外处理。

## 3. 产物清单（服务器端，待拉回本地）

| 产物 | 行数 | npz 大小 | 完成时间(UTC) | 输出类型 |
|---|---|---|---|---|
| scgpt_wave_k562 / rpe1 | 217 / 306 | 27MB / 38MB | 10-19 / 10-20 | 512-d CLS cell_emb |
| geneformer_wave_k562 / rpe1 | 217 / 306 | 41MB / 57MB | 10-18 / 10-20 | 768-d mean-pool cell_emb |
| uce_wave_k562 / rpe1 | 217 / 306 | 68MB / 95MB | 10-19 / 10-20 | 1280-d CLS cell_emb |

每个产物三件套齐全（`_predictions.npz` + `_labels.csv` + `_meta.json`），meta 含
权重 sha256、命令行、线程数、秒级耗时。labels 含探针列 pred/true/delta_emb_norm 与
err_emb_dist。合计约 330MB，SFTP 拉回约 1-2 分钟。

## 4. 推荐执行顺序

**当前状态（环境/权重/推理全就绪）下的剩余顺序：**

1. SFTP 拉回 6 份 wave 产物到 `outputs_crossmodel/wave/{scgpt,geneformer,uce}/`，
   按 `verify_npz.py` 口径校验（键齐、行序与 held-out 名单一致、值有限）；
2. 门控探针集成：以 labels 的 emb_norm / err_emb_dist 列对接
   virtual_cell_constraint_mvp 探针框架（scGenePT/scFoundation 首轮产物同目录可比）；
3. 若上游分析必须要表达谱 pearson_delta：scGPT 首轮 expression 产物
   （scgpt_k562_predictions.npz，本地已有）可作一家之言；Geneformer/UCE 需补
   表达头协议（成本各约一轮 wave 时长）或改用嵌入距离口径。

**若需重跑/扩数据集（新 h5ad）的参考顺序**（按单数据集耗时升序，先快后慢验证管线）：

1. Geneformer（~2.5h，依赖最简，最快暴露管线问题）；
2. UCE（~11h，权重已就位，注意 8.5GB 加载内存峰值）；
3. scGPT（~17h，注意 cpuset 与 MLM 头跳过）。
   三容器配额合计 176 线程 < 256，可并行但建议错峰（实测高负载下单扰动耗时上升）。

## 5. 口径提醒（与任务假设的差异，需主代理确认）

任务背景假设三引擎"均为表达谱类、误差口径 pearson_delta"。**实际 wave 产物全部是
per-cell cell embedding**，三引擎 labels 的 pearson 列全 NaN 是设计如此，探针列是
emb_norm / err_emb_dist。这与门控探针框架直接兼容，但与"表达谱 pearson_delta"口径
不同——scGPT 另有首轮 expression 产物（mask 协议，局限已在 RUNBOOK 记录），
Geneformer/UCE 无表达谱输出。若 pearson_delta 是硬需求，需补协议（见 §4.3）。

## 6. 硬阻塞清单

无。磁盘（694G / 55T）、内存（840G 可用）、容器（3/3 在跑）、权重（3/3 就位）、
算力余量（load 0）全部满足；唯一动作项是 §4.1 的 330MB 产物回拉。
