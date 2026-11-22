# scGenePT 跨模型泛化验证 — 服务器侧环境搭建脚本
# 在宿主机执行：bash /opt/research/scripts/setup_env.sh > /opt/research/scripts/setup.log 2>&1
set -x

# 1. 容器（显式资源限制；restart=no；不碰 animal_house_*）
docker rm -f research_scgenept 2>/dev/null
docker run -d --name research_scgenept \
  --restart no --cpus 32 --memory 128g --pids-limit 512 \
  -v /opt/research:/research -w /research \
  python:3.10-slim-bookworm sleep infinity

# 2. Python 依赖（CPU 版 torch 先行，防 pip 拉 CUDA 版；numpy 钉 1.26.4 与 repo 一致）
docker exec research_scgenept pip install --no-cache-dir --upgrade pip
docker exec research_scgenept pip install --no-cache-dir \
  torch==2.1.2+cpu torchtext==0.16.2+cpu \
  --index-url https://download.pytorch.org/whl/cpu
docker exec research_scgenept pip install --no-cache-dir \
  "numpy==1.26.4" "pandas<2.3" "scipy<1.14" "scikit-learn<1.6"
docker exec research_scgenept pip install --no-cache-dir \
  scgpt==0.2.1 gears torch_geometric==2.5.3
docker exec research_scgenept pip install --no-cache-dir \
  "anndata<0.11" "scanpy<1.11" matplotlib tqdm numba h5py pyarrow

docker exec research_scgenept python -c \
  "import torch,scgpt,gears,scanpy,anndata,torchtext;print('torch',torch.__version__,'cuda',torch.cuda.is_available());print('ok')"

# 3. scGenePT 代码（固定 commit）
cd /opt/research/scgenept_repo
git clone https://github.com/czi-ai/scGenePT.git repo
cd repo && git rev-parse HEAD && cd /opt/research

# 4. 模型权重与基因嵌入（S3 公开桶，服务器直连可达）
mkdir -p /opt/research/models/finetuned/scgenept_ncbi+uniprot/norman
curl -sL --retry 3 -o "/opt/research/models/finetuned/scgenept_ncbi+uniprot/norman/best_model_gpt3.5_ada_rnd_seed_42.pt" \
  "https://czi-scgenept-public.s3.amazonaws.com/models/finetuned/scgenept_ncbi%2Buniprot/norman/best_model_gpt3.5_ada_rnd_seed_42.pt"
curl -sL --retry 3 -o /opt/research/models/gene_embeddings/NCBI+UniProt_embeddings-gpt3.5-ada.pkl \
  "https://czi-scgenept-public.s3.amazonaws.com/models/gene_embeddings/NCBI%2BUniProt_embeddings-gpt3.5-ada.pkl"
# vocab.json 已在 repo 内 models/pretrained/scgpt/vocab.json
ls -la /opt/research/models/finetuned/scgenept_ncbi+uniprot/norman/ /opt/research/models/gene_embeddings/
echo SETUP_DONE
