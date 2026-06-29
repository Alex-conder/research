# 强制推送前备份检查清单 + 新建仓库方案

> 本指南提供两套方案：
> 1. **方案 A**：强制推送覆盖 `Kaelis-archive`，但需先备份。
> 2. **方案 B**：新建 `research` 仓库，更安全，已最终采用。
>
> **推荐方案 B**：避免破坏 `Kaelis-archive` 历史，降低操作风险。

---

## 方案 A：强制推送覆盖 Kaelis-archive（高风险）

### 适用场景
- 你确定 `Kaelis-archive` 当前内容可以丢弃或已迁移。
- 你希望保留原来的仓库名和 Stars/Forks（如果有）。
- 你不介意丢失原仓库的 commit 历史。

### 强制推送前备份检查清单

- [ ] **1. 确认本地是否有 `Kaelis-archive` 的完整副本**
  ```bash
  # 检查本地目录
  ls D:\Kaelis-archive
  # 或检查是否有其他位置的备份
  ```

- [ ] **2. 下载 GitHub 上 `Kaelis-archive` 的完整备份 ZIP**
  1. 浏览器打开 `https://github.com/Alex-conder/Kaelis-archive`
  2. 点击绿色 `<> Code` 按钮
  3. 选择 `Download ZIP`
  4. 保存到 `D:\Backups\Kaelis-archive-backup-2026-06-28.zip`

- [ ] **3. 确认 `Kaelis-main` 已设为私有或已备份**
  - 如果 `Kaelis-main` 和 `Kaelis-archive` 有重叠代码，确保 `Kaelis-main` 不会因为 `Kaelis-archive` 被覆盖而丢失。

- [ ] **4. 检查原 `Kaelis-archive` 是否有重要 Issue / Pull Request / Wiki**
  - 打开 `https://github.com/Alex-conder/Kaelis-archive/issues`
  - 打开 `https://github.com/Alex-conder/Kaelis-archive/pulls`
  - 如有重要内容，截图或导出保存。

- [ ] **5. 确认原 `Kaelis-archive` 没有未同步到本地的提交**
  ```bash
  cd D:\Kaelis-archive
  git fetch origin
  git status
  ```
  如果有未推送的本地提交，请先备份。

- [ ] **6. 确认 `.env`、API key、密码等敏感信息未出现在 `Kaelis-portfolio` 中**
  ```bash
  cd D:\Kaelis-portfolio
  find . -name ".env" -o -name "*.key" -o -name "*secret*"
  ```

- [ ] **7. 确认大文件已排除**
  ```bash
  cd D:\Kaelis-portfolio
  find . -type f -size +50M
  ```
  应无输出。

- [ ] **8. 准备好强制推送命令**
  ```bash
  cd D:\Kaelis-portfolio
  git init
  git add .
  git commit -m "feat: add academic research portfolio (docking, multi-omics, virtual cell, metabolomics)"
  git branch -M main
  git remote add origin https://github.com/Alex-conder/Kaelis-archive.git
  git push -u origin main --force
  ```

- [ ] **9. 强制推送后立即验证**
  - 浏览器打开 `https://github.com/Alex-conder/Kaelis-archive`
  - 确认首页为学术作品集 README
  - 确认 5 个项目目录存在

---

## 方案 B：新建 `research` 仓库（推荐）

### 适用场景
- 你不想破坏 `Kaelis-archive` 的现有内容和历史。
- 你希望 `Kaelis-archive` 继续作为 AI Agent OS 项目存在。
- 你愿意在套磁信中修改 GitHub 链接。

### 优势
- **零风险**：不覆盖任何已有仓库。
- **更清晰**：`research-portfolio` 名称直接表明用途。
- **更灵活**：未来可以持续维护学术作品，与工程项目分开。

### 步骤 1：在 GitHub 上创建新仓库

1. 登录 GitHub：https://github.com
2. 点击右上角 `+` → `New repository`
3. 填写信息：
   - **Repository name**: `research-portfolio`
   - **Description**: `Academic research portfolio: molecular docking, multi-omics, virtual cell modeling, and metabolomics`
   - **Visibility**: `Public`
   - **Initialize this repository with**: 全部**不勾选**（因为本地已有 README）
4. 点击 `Create repository`

### 步骤 2：本地推送到新仓库

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

### 步骤 3：验证新仓库

浏览器打开：https://github.com/Alex-conder/research-portfolio

确认：
- [ ] 首页显示 `Huang Teng — Academic Research Portfolio`
- [ ] 5 个项目目录清晰可见
- [ ] README 内容正确
- [ ] 没有敏感文件泄露

### 步骤 4：更新套磁信中的 GitHub 链接

将之前所有套磁信中的：

```
https://github.com/Alex-conder/Kaelis-archive
```

替换为：

```
https://github.com/Alex-conder/research-portfolio
```

需要修改的文件：
- `d:/vscode workplace/硕士套磁信_黄腾_四导师独立风格版_终稿.md`
- `d:/vscode workplace/硕士套磁信_黄腾_四导师终稿.docx`
- 已经发送过的邮件（如果已发送，无需追回，后续邮件用新链接）

### 步骤 5：处理原有 `Kaelis-archive`

有三种选择：

**选择 1：保留 `Kaelis-archive` 不变**
- 最简单。教授如果误点旧链接，看到的是 AI Agent OS 项目，但不影响新链接。

**选择 2：将 `Kaelis-archive` 的 README 修改为跳转页**
- 在 `Kaelis-archive` README 顶部添加：
  ```markdown
  > 📁 My academic research portfolio has moved to [Alex-conder/research-portfolio](https://github.com/Alex-conder/research-portfolio).
  ```

**选择 3：将 `Kaelis-archive` 设为私有或归档（Archive）**
- 如果 `Kaelis-archive` 不再维护，可以：
  - Settings → Danger Zone → Archive this repository
  - 或 Settings → Danger Zone → Change repository visibility → Private

---

## 方案对比

| 维度 | 方案 A：强制推送覆盖 | 方案 B：新建 research-portfolio |
|---|---|---|
| 风险 | 高（丢失原仓库历史） | 低（不破坏原仓库） |
| 操作复杂度 | 简单 | 中等（需改链接） |
| 套磁信链接 | 无需修改 | 需全部替换 |
| 仓库命名 | 保留 Kaelis-archive | 更专业的 research-portfolio |
| 推荐度 | ⭐⭐ | ⭐⭐⭐⭐⭐ |

---

## 最终建议

**强烈建议采用方案 B**：
1. 新建 `research` 仓库并推送作品集。
2. 将套磁信中的 GitHub 链接全部替换为新地址。
3. 将 `Kaelis-archive` 的 README 添加跳转提示，或将其设为私有/归档。
4. 将 `Kaelis-main` 设为私有。

这样既保护了原有项目，又为教授提供了一个清晰、专业的学术作品集入口。
