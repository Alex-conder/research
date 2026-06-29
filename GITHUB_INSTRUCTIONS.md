# GitHub 发布指南（已选方案 B：新建 research-portfolio）

> **已确认方案**：新建公开仓库 `research-portfolio` 作为学术作品集，不破坏原有 `Kaelis-archive`。
>
> **原因**：`research-portfolio` 名称更专业、风险更低，且便于与工程项目 `Kaelis-archive` / `Kaelis-main` 区分。

---

## 第一步：在 GitHub 上创建 `research-portfolio` 仓库

1. 打开浏览器，登录 GitHub：https://github.com
2. 点击右上角 `+` → `New repository`
3. 填写信息：
   - **Repository name**: `research-portfolio`
   - **Description**: `Academic research portfolio: molecular docking, multi-omics, virtual cell modeling, and metabolomics`
   - **Visibility**: `Public`
   - **Initialize this repository with**: **全部不勾选**（因为本地已有 README）
4. 点击 `Create repository`

---

## 第二步：将本地作品集推送到新仓库

打开 Git Bash 或 Windows Terminal，执行以下命令：

```bash
# 1. 进入作品集目录
cd D:\Kaelis-portfolio

# 2. 初始化 Git
git init

# 3. 添加所有文件
git add .

# 4. 提交
git commit -m "feat: initial academic research portfolio

- 6UJN sequential docking with Amber MD/MM-GBSA
- PirB stroke multi-omics (scRNA-seq, snRNA-seq, Visium)
- HH-FBA-scVI multi-scale virtual cell modeling
- LC-MS non-targeted metabolomics pipeline
- LLM-assisted research automation workflow"

# 5. 关联远程仓库
git branch -M main
git remote add origin https://github.com/Alex-conder/research-portfolio.git

# 6. 推送
git push -u origin main
```

---

## 第三步：验证新仓库

浏览器打开：https://github.com/Alex-conder/research-portfolio

确认以下事项：
- [ ] 首页显示 `Huang Teng — Academic Research Portfolio`
- [ ] 5 个项目目录（01_6UJN_docking 至 05_llm_automation）清晰可见
- [ ] 每个项目目录下都有 README.md
- [ ] 没有 `.env`、API key、密码等敏感信息泄露
- [ ] 仓库大小合理（当前约 150MB，GitHub 单文件限制 100MB）

---

## 第四步：更新套磁信中的 GitHub 链接

已将所有套磁信文档中的链接替换为新的仓库地址。请使用以下文件：

- `d:/vscode workplace/硕士套磁信_黄腾_四导师独立风格版_终稿.md`
- `d:/vscode workplace/硕士套磁信_黄腾_四导师终稿_research_portfolio.docx`

新链接统一为：

```
https://github.com/Alex-conder/research-portfolio
```

**如果你之前已经发送过任何邮件**，无需追回，后续邮件使用新链接即可。

---

## 第五步：处理原有 `Kaelis-archive`（三选一）

### 选择 1：保留不变（最简单）
- 不操作 `Kaelis-archive`。
- 风险：如果教授误点旧链接，会看到 AI Agent OS 项目，但不影响新链接。

### 选择 2：添加跳转提示（推荐）
1. 打开 `https://github.com/Alex-conder/Kaelis-archive`
2. 编辑 `README.md`，在**最顶部**添加：
   ```markdown
   > 📁 My academic research portfolio has moved to [Alex-conder/research-portfolio](https://github.com/Alex-conder/research-portfolio).
   ```
3. 提交更改。

### 选择 3：归档或设为私有
- **归档**：Settings → Danger Zone → Archive this repository
  - 归档后仓库只读，外人仍可查看但会显示 "This repository has been archived".
- **设为私有**：Settings → Danger Zone → Change repository visibility → Private
  - 设为私有后外人无法访问，但会消耗一个私有仓库名额（GitHub 免费账户有数量限制）。

> **建议**：如果 `Kaelis-archive` 还有参考价值，选择 **选择 2（添加跳转提示）**。

---

## 第六步：将 `Kaelis-main` 设为私有

### 情况 1：`Kaelis-main` 已经是 GitHub 仓库

1. 打开 `https://github.com/Alex-conder/Kaelis-main`
2. 点击页面上方的 **Settings**
3. 左侧菜单拉到最下方，点击 **Danger Zone** 下的 **Change repository visibility**
4. 选择 **Change to private**
5. 输入仓库名确认 → 点击 **I understand, change repository visibility**

### 情况 2：`Kaelis-main` 只在本地

1. 在 GitHub 上创建新仓库：`Alex-conder/Kaelis-main`
2. Visibility 选择 **`Private`**
3. 执行：
   ```bash
   cd D:\Kaelis-main
   git init
   git add .
   git commit -m "feat: initial private Kaelis-main project"
   git branch -M main
   git remote add origin https://github.com/Alex-conder/Kaelis-main.git
   git push -u origin main
   ```

---

## 安全提醒

1. **不要上传成绩单、简历到 GitHub**：这些通过邮件附件发送。
2. **不要上传 API key 或 `.env` 文件**：本作品集已配置 `.gitignore`，但请再次检查。
3. **不要上传原始测序/质谱数据**：`.h5ad`、`.mzML`、FASTQ 等已排除。

---

## 完成后的最终状态

| 仓库 | 可见性 | 用途 |
|---|---|---|
| `research-portfolio` | Public | 学术作品集，套磁信中展示 |
| `Kaelis-archive` | Public（可选：添加跳转提示） | 原有 AI Agent OS 项目 |
| `Kaelis-main` | Private | 完整 Kaelis 项目 |
