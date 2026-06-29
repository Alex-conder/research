# LLM-Assisted Research Automation

## Overview

This directory documents my workflow for using Large Language Model (LLM) APIs as a **rapid-prototyping accelerator** in computational research. The key principle is: **LLMs accelerate code generation and documentation drafting, while scientific design, validation, and interpretation remain under my full control.**

## Automated Modules

| Module | Purpose | Example Output |
|---|---|---|
| Technical Chain Automation | Generate analysis pipeline skeletons from high-level intent | GEO download → QC → DE → pathway scripts |
| Boundary Definition | Auto-generate parameter config files and validation checks | `config.yaml`, QC thresholds, logging |
| Report Drafting | Convert statistical summaries into structured reports | Nature Communications-style manuscript drafts |
| Code Refactoring | Modularize scripts into reusable Python packages | `src/` modules for docking, scVI, FBA |

## Principles

1. **Human-in-the-loop**: All LLM-generated code is reviewed, tested, and validated before use.
2. **Modularity**: Reusable components are extracted into functions/classes.
3. **Reproducibility**: Every pipeline includes versioned dependencies and execution logs.
4. **Transparency**: LLM-assisted sections are clearly marked in reports.

## Tools Used

- OpenAI API / 通义千问 (Qwen) / Moonshot / DeepSeek
- Python `openai` SDK and custom HTTP wrappers
- Prompt templates for code generation, debugging, and documentation
