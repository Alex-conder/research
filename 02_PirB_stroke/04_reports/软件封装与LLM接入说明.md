# Pirb-Stroke 分析流程软件封装与 LLM 接入说明

## 一、封装目标

将本次"脑缺血后 Pirb 阳性细胞单细胞图谱"完整分析流程（数据加载、质控、差异表达、跨数据集整合、空间分析、报告生成）封装为可复用的软件模板，接入 LLM 接口后可稳定输出可信化分析报告。

## 二、软件模板位置

本地目录：

```
D:\Pirb_stroke_project\05_software_template\
```

已一并打包至桌面交付包：

```
C:\Users\ASUS\Desktop\Pirb_stroke_project_最终交付包.zip\13_软件封装模板\
```

## 三、模板目录结构

```
05_software_template/
├── README.md                    # 项目说明
├── AGENTS.md                    # 编码与 LLM 使用规范
├── config.yaml                  # 可配置文件
├── requirements.txt             # Python 依赖
├── src/                         # 核心源代码
│   ├── pipeline.py              # 主流程
│   ├── data_loader.py           # 数据加载
│   ├── qc_module.py             # 质控
│   ├── de_analysis.py           # 差异表达
│   ├── integration.py           # 跨数据集整合
│   ├── spatial_analysis.py      # Visium 空间分析
│   ├── report_generator.py      # 报告生成
│   └── llm_interface.py         # LLM API 封装
├── templates/                   # 报告模板（MD + DOCX）
│   ├── nc_report_template_chinese.md
│   ├── original_report_template.md
│   ├── cover_letter_template.md
│   └── response_to_reviewers_template.md
├── prompts/                     # LLM 提示词模板
│   ├── abstract_generation.txt
│   ├── discussion_generation.txt
│   ├── mechanism_inference.txt
│   ├── figure_legends_generation.txt
│   └── limitations_suggestions.txt
├── tests/                       # 单元测试
└── docs/                        # 架构与部署文档
    ├── architecture.md
    ├── deployment_guide.md
    └── trustworthiness_checklist.md
```

## 四、接入 LLM 接口步骤

### 1. 安装依赖

```bash
cd D:/Pirb_stroke_project/05_software_template
python -m venv venv
venv\Scripts\activate  # Windows
pip install -r requirements.txt
```

### 2. 配置 LLM

编辑 `config.yaml`：

```yaml
llm:
  provider: "openai"      # 可选：openai, anthropic, moonshot, qwen, http
  model: "gpt-4o"
  api_key_env: "LLM_API_KEY"
  base_url: null
  temperature: 0.3
  max_tokens: 2048
  append_review_notice: true
```

设置环境变量：

```bash
set LLM_API_KEY=your-api-key  # Windows CMD
# 或 $env:LLM_API_KEY="your-api-key"  # PowerShell
```

### 3. 接入国产/私有化模型示例

**Moonshot / Kimi**

```yaml
llm:
  provider: "http"
  model: "moonshot-v1-8k"
  base_url: "https://api.moonshot.cn/v1/chat/completions"
  api_key_env: "MOONSHOT_API_KEY"
```

**通义千问**

```yaml
llm:
  provider: "http"
  model: "qwen-max"
  base_url: "https://dashscope.aliyuncs.com/api/v1/services/aigc/text-generation/generation"
  api_key_env: "DASHSCOPE_API_KEY"
```

**私有化部署（兼容 OpenAI 接口）**

```yaml
llm:
  provider: "openai"
  model: "your-local-model"
  base_url: "http://localhost:8000/v1"
  api_key_env: "LLM_API_KEY"
```

### 4. 运行完整流程

```bash
python -m src.pipeline --config config.yaml
```

## 五、可信化设计

1. **数据隔离**：LLM 仅接收脱敏后的统计摘要、差异基因列表、通路名称，不接触原始表达矩阵。
2. **复核声明**：所有 LLM 生成文本自动附加"AI 生成草稿，需专家复核"。
3. **数值校验**：模板中的定量占位符（如 Pirb+ 比例、p 值）必须来自真实分析结果。
4. **证据强度标注**：机制推断、讨论段落要求标注置信度（高/中/低）。
5. **检查清单**：`docs/trustworthiness_checklist.md` 用于发布前逐项确认。

## 六、可稳定输出的报告类型

- Nature Communications 英文版
- Nature Communications 中文版
- 中文阶段性报告（原报告格式）
- Cover Letter
- Response to Reviewers
- Supplementary Tables Excel
- Graphical Abstract（PNG/SVG/PDF）
- 提升点对比 PDF 表格

## 七、扩展建议

1. **前端界面**：可用 Gradio/Streamlit 封装为 Web 应用，用户上传配置即可运行。
2. **工作流编排**：可用 Snakemake/Nextflow 管理分析依赖。
3. **数据库持久化**：将分析结果存入 SQLite/PostgreSQL，支持历史版本追溯。
4. **模型微调**：基于本领域文献微调开源 LLM，提升摘要和机制推断的专业度。
5. **CI/CD**：GitHub Actions 自动运行单元测试和示例流程。

## 八、核心代码入口

```
src/pipeline.py      # 主流程
src/llm_interface.py # LLM 封装
src/report_generator.py # 报告渲染
```

## 九、注意事项

- 当前模板为 v1.0 框架，部分模块（如空间分析）为伪代码框架，需根据实际数据结构补全。
- 报告 DOCX 生成采用简单 Markdown→段落映射，复杂排版建议后续改用 python-docx 精细控制。
- LLM 生成内容不可替代专家审核，发布前必须人工复核。
