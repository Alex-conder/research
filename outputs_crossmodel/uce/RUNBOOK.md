# RUNBOOK — UCE 33-layer 跨模型泛化验证（Replogle K562）

> 服务器：172.16.92.2（Rocky Linux 9.8，256 线程 EPYC 7702 / 881GB RAM / 无 GPU）
> 工作区：`/opt/research`；硬约束：绝不碰 `animal_house_*` 容器。

## 0. 产物约定

`/research/results/uce_k562_predictions.npz`：
`pred_perturbations (217,)` / `pred_emb (217, 1280)`（L2 归一化 CLS 嵌入）/
`pred_delta_emb` / `ctrl_emb (1280,)` / `ctrl_emb_cells (32, 1280)` /
`target_in_input_space` / `target_deleted` / `ctrl_mean_full` /
`uce_gene_space`（UCE human 映射子集）/ `full_gene_space` / `output_type`
（**输出为 cell embedding，非表达谱**；labels pearson 列 NaN 为设计如此）。

## 1. 容器与环境

```bash
docker run -d --name research_uce --restart no --cpus 48 --memory 128g \
  --pids-limit 512 -v /opt/research:/research -w /research \
  python:3.10-slim-bookworm sleep infinity
# 容器内：
pip install --upgrade pip
pip install torch==2.1.2+cpu --index-url https://download.pytorch.org/whl/cpu
pip install "numpy==1.26.4" "pandas<2.3" "scipy<1.14" "anndata<0.11" "scanpy<1.11" tqdm
```

- UCE 代码：`/opt/research/uce_repo/UCE-main`（GitHub main tarball，本机下载 SFTP，
  服务器直连 github 不稳定）。
- 无 cpuset（floating 48 核）：scGPT 占 socket0、Geneformer 占 socket1，
  UCE 浮动避让；scgenept 任务结束后可 `docker update --cpuset-cpus` 钉扎。

## 2. 权重（坑：figshare 在本网络被 202 卡死，改用 HF 镜像）

官方权重在 figshare（ndownloader/files/42706558 等），但服务器与本机访问
figshare.com 一律返回 HTTP 202 空响应（疑似区域限制，重试无效）。
**改用 HF 镜像 `lza1/uce_hf`**（非官方 HF 实现，model_files 原样拷贝官方文件）：

```bash
cd /opt/research/models/uce/model_files
B=https://hf-mirror.com/lza1/uce_hf/resolve/main/model_files
curl -sL -o species_chrom.csv $B/species_chrom.csv                 # 4MB
curl -sL -o species_offsets.pkl $B/species_offsets.pkl             # 139B
curl -sL -o protein_embeddings/Homo_sapiens.GRCh38.gene_symbol_to_embedding_ESM2.pt \
  $B/protein_embeddings/Homo_sapiens.GRCh38.gene_symbol_to_embedding_ESM2.pt   # 411MB
curl -sL -o 33l_8ep_1024t_1280.torch $B/33l_8ep_1024t_1280.torch   # 5.69GB
curl -sL -o all_tokens.torch $B/all_tokens.torch                   # 2.98GB
```

- 33L 模型含训练好的 `pe_embedding`（145469×5120），与 all_tokens.torch 二选一
  由官方加载逻辑决定（checkpoint 已含 145469 行时保留 checkpoint 的）。
- 只需要 human 的 protein embedding 文件（基因符号→PE 行号），其他物种不下载。

## 3. 推理协议（官方 eval 流程复刻 + 确定性 seed + in silico deletion）

- 基因子集：h5ad gene_name.upper() 必须在 human PE 字典 AND species_chrom.csv；
- 细胞句：权重 ∝ 表达（h5ad 为 log-normalized；log1p(expm1(x))=x，故 ∝ X 与
  ∝ log1p(counts) 采样概率一致，协议偏离已记录）；每细胞有放回采 1024 个基因
  token（官方 np.random 无 seed → 本脚本每细胞固定 seed 保证可复现）；
  染色体乱序 → 起点排序 → 包裹染色体左右 token → CLS 头 → PAD 1536；
- 扰动：采样前目标基因权重置 0（=从细胞句删除）；
- 嵌入 = decoder CLS 输出 L2 归一化；pred_emb[pi] = 池(32)均值。
- 33L 前向很重（33 层 × d1280 × seq~1000），32 细胞 × 217 扰动，CPU 小时级。

## 4. 运行

```bash
# 冒烟（3 扰动）：
docker exec research_uce python -u /research/scripts/run_uce_k562.py \
  --limit 3 --threads 32 --out /research/results/uce_smoke
# 全量：
setsid nohup docker exec -e SCT_THREADS=32 -e OMP_NUM_THREADS=32 research_uce \
  python -u /research/scripts/run_uce_k562.py --threads 32 \
  --out /research/results/uce_k562 \
  > /opt/research/results/uce_k562.log 2>&1 < /dev/null &
```

- 断点：每扰动写 `uce_k562_checkpoint.npz`；基线 ctrl 嵌入不缓存（仅 90s 量级，
  缓存逻辑未加，重启会重算一次基线）。
- 杀任务：`docker top research_uce | grep run_uce | awk '{print $2}' | xargs -r kill`。

- 冒烟（3 扰动）：PRPF3/COX17 deleted=True（各 ~300-330s），POP7 退化基线；
  产物写盘正常（10.5min 端到端）。
- 基线 32 细胞 358s（~11s/细胞，32 线程，机器 load 130 高负载）；全量已启动：
  `setsid nohup docker exec -e SCT_THREADS=32 research_uce python -u
  /research/scripts/run_uce_k562.py --threads 32 --out /research/results/uce_k562`
  （日志 `results/uce_k562.log`）。高负载 ETA ~18-20h；其他任务结束后加速。
- **坑**：`data_utils.py` 模块级读 `./model_files/new_species_protein_embeddings.csv`
  （CWD 相对路径！容器 CWD=/research → 需放 `/opt/research/model_files/` 下，
  不是模型目录）。文件只有表头 14 字节。
- **坑**：宿主机 ssh 里 `/research` 不是挂载点，容器内 `/research` = 宿主机
  `/opt/research`（创建文件时注意路径）。

## 5. 收尾

```bash
docker stop research_uce
```
