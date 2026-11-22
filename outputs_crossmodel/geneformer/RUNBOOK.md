# RUNBOOK — Geneformer V2-104M 跨模型泛化验证（Replogle K562）

> 服务器：172.16.92.2（Rocky Linux 9.8，256 线程 EPYC 7702 / 881GB RAM / 无 GPU）
> 工作区：`/opt/research`；硬约束：绝不碰 `animal_house_*` 容器。

## 0. 产物约定

`/research/results/geneformer_k562_predictions.npz`：
`pred_perturbations (217,)` / `pred_emb (217, 768)`（扰动后细胞嵌入）/
`pred_delta_emb (217, 768)` / `ctrl_emb (768,)` / `ctrl_emb_cells (32, 768)` /
`ctrl_pool_indices` / `target_in_input_space (217,)` / `target_deleted (217,)` /
`ctrl_mean_full (5000,)` / `full_gene_space (5000,)` / `output_type`（明确标注
**输出为 cell embedding，非表达谱**；`_labels.csv` 的 pearson 列全 NaN 是设计如此）。

## 1. 容器与环境

```bash
docker run -d --name research_geneformer --restart no --cpus 64 --memory 128g \
  --pids-limit 512 --cpuset-cpus "64-127,192-255" \
  -v /opt/research:/research -w /research python:3.10-slim-bookworm sleep infinity
# 容器内（必须先进 pip，旧 pip 装 torch 会触发 typing-extensions sdist 构建失败）：
pip install --upgrade pip
pip install torch==2.1.2+cpu --index-url https://download.pytorch.org/whl/cpu
pip install "numpy==1.26.4" "pandas<2.3" "scipy<1.14" "transformers==4.44.2" "anndata<0.11" tqdm
```

## 2. 权重（hf-mirror 直下载）

```bash
cd /opt/research/models/geneformer
curl -sL -o token_dictionary_gc104M.pkl https://hf-mirror.com/ctheodoris/Geneformer/resolve/main/geneformer/token_dictionary_gc104M.pkl
curl -sL -o gene_name_id_dict_gc104M.pkl  https://hf-mirror.com/ctheodoris/Geneformer/resolve/main/geneformer/gene_name_id_dict_gc104M.pkl
cd Geneformer-V2-104M   # mkdir -p  beforehand
for f in config.json generation_config.json model.safetensors; do
  curl -sL -o $f https://hf-mirror.com/ctheodoris/Geneformer/resolve/main/Geneformer-V2-104M/$f
done
```

- V2-104M：BertForMaskedLM(hidden 768, 12 层, vocab 20275)，safetensors 417MB。
- **坑**：特殊 token 是 `<pad>=0 <mask>=1 <cls>=2 <eos>=3`（没有 `<sep>`！）。
- pooler 权重未初始化警告无害（不用 pooler；细胞嵌入 = 基因 token last_hidden_state 均值）。
- K562 基因名 → Ensembl：`gene_name_id_dict_gc104M.pkl`（63675 个符号）。
  217 扰动中 134 个目标不在 5000 输入基因内；输入空间内的 79 个里 6 个无
  Ensembl/token 映射 → 可删除 73 个。

## 3. 推理协议（in silico deletion，如实标注输出类型）

- ctrl 池（pool=32, seed=0，与 scGenePT/scGPT 同池同 seed）；每细胞按表达降序
  取非零基因（截到 2042）构成 cell sentence；
- 基线嵌入 = BERT 前向，基因 token（去 CLS/EOS）last_hidden_state 均值；
- 扰动 = 从句中删除该基因 token，重算嵌入；`pred_emb[pi]` = 池均值；
- 无表达谱输出 → pearson 评估不适用，labels 的 pearson 列留 NaN。

## 4. 运行

```bash
setsid nohup docker exec -e SCT_THREADS=32 -e OMP_NUM_THREADS=32 research_geneformer \
  python -u /research/scripts/run_geneformer_k562.py --threads 32 --batch 16 \
  --out /research/results/geneformer_k562 \
  > /opt/research/results/geneformer_k562.log 2>&1 < /dev/null &
# 断点：每扰动写 geneformer_k562_checkpoint.npz（含 ctrl_emb 缓存 .ctrl_emb.npy）
# 杀任务：docker top research_geneformer | grep run_geneformer | awk '{print $2}' | xargs -r kill
```

- 冒烟（3 扰动）：每 pert ~65-85s（32 线程，高负载）；"sys 未导入"崩溃一次
  （meta 段 `sys.argv`），已修，重新上传后全量启动。
- **全量完成（2026-10-16）**：117.9 min 跑完 217 扰动（73 个真实删除 ×~100s +
  144 个退化 0s）；产物已下载本机 `outputs_crossmodel/geneformer/` 并通过
  `verify_npz.py geneformer` 全部校验（键齐、217 顺序一致、pred_emb (217,768)
  有限值、138 个零 delta、pearson 全 NaN by design）。
- 注意：`target_in_input_space` 本模型定义为"目标 ∈ h5ad var 基因"（83 个 True），
  比 scgenept 的 vocab 口径（79）略宽；真实可删除 73 个（79 减去 6 个无
  Ensembl/token 映射）。

## 5. 收尾

```bash
docker stop research_geneformer   # restart=no
```
