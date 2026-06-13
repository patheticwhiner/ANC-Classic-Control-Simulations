---
name: git-workflow
description: "Standard git workflow for this project — feature branches, never commit directly to master"
metadata: 
  node_type: memory
  type: project
  originSessionId: 47cabe08-1aa7-42be-be0c-1f5fd39c9e53
---

# Git 提交规范

## 原则

- **永远不在 master 上直接 commit** — master 只通过 PR/merge 接收变更
- **feature branch 上可以自由 amend、rebase** — 保持提交历史干净
- **提交前确认不包含无关文件** — `git status` 检查，只 stage 本次修改的文件

## 标准流程

```bash
# 1. 从最新 master 创建 feature 分支
git checkout master
git pull --rebase
git checkout -b feature/<简短描述>

# 2. 修改代码后，选择性 stage
git add <修改的文件>
# 确认: git status — 不应有无关联的文件

# 3. 提交 (conventional commits 格式)
git commit -m "fix: <简述>"
# 如需详细 body:
# git commit -m "fix: <简述>" -m "<详细说明>"

# 4. 推送到远程 feature 分支
git push -u origin feature/<简短描述>

# 5. 在 GitHub 上创建 PR，合并后清理本地
git checkout master
git pull --rebase
git branch -D feature/<简短描述>
```

## 与公司 Gerrit 流程的对应

| 公司流程 | GitHub 等价操作 |
|:---|:---|
| `git checkout -b feature` | 相同 |
| `git push HEAD:refs/for/master` | Push feature branch → 创建 GitHub PR → merge |
| `git branch -D feature` | 相同（PR 合并后执行） |

## Commit Message 格式

```
<type>: <简短描述>

<详细说明 (可选)>
```

常用 type：`fix`（修 bug）、`feat`（新功能）、`refactor`（重构）、`docs`（文档）。

## 反模式

- ❌ 在 master 上直接 `git commit`
- ❌ `git add .` 或 `git add -A`（容易误 stage 无关文件）
- ❌ `git push --force` 到 master 或已共享的分支
- ❌ 一个 commit 包含多个不相关的修改
