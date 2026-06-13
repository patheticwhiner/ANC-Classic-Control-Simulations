# 项目文档索引

> 生成: 2026-06-13 | 文档总数: 38 (不含 .claude/memory/)

## 使用指南

| 标签 | 含义 |
|:---|:---|
| 🎯 理论 | 定理-证明-引理链, 面向数学推导研究者 |
| 🔧 工程 | 代码/参数/仿真注意事项, 面向工程师 |
| 🎯🔧 兼有 | 同时包含严格推导和工程直觉 |
| ✅ 活跃 | 内容完整, 与当前项目状态一致 |
| ⚠️ 过时 | 内容不再适用或已被取代 |
| 📝 占位 | 仅有框架/标题, 内容待补充 |

---

## 项目顶层

| 文件 | 读者 | 状态 | 用途 |
|:---|:---|:---|:---|
| `README.md` | 🎯🔧 | ✅ | 项目主入口 — 概述、目录结构、四个 demo 简介 |
| `CONTEXT.md` | 🎯🔧 | ✅ | 术语表 — 规范术语定义 (抑制/收敛/稳态/测试基准等) |
| `DOCUMENT_INDEX.md` | 🎯🔧 | ✅ | 本文件 — 全项目文档索引 |

## 计划与架构

| 文件 | 读者 | 状态 | 用途 |
|:---|:---|:---|:---|
| `plans/Global_Work_Plan.md` | 🎯🔧 | ✅ | 四阶段全局计划 (Phase 1→4), 含当前状态和决策框架 |
| `plans/unified_eMOPSO_architecture.md` | 🔧 | 📝 | eMOPSO 统一入口架构设计 (TODO 未完成) |

## 测试套件

| 文件 | 读者 | 状态 | 用途 |
|:---|:---|:---|:---|
| `tests/README.md` | 🔧 | ✅ | tests/ 使用说明 — 信号生成 → 批量运行 → 结果解读 |

## 数据集

| 文件 | 读者 | 状态 | 用途 |
|:---|:---|:---|:---|
| `dataset/README.md` | 🔧 | ✅ | dataset/ 目录索引 — 可用模型、DataManager 用法 |
| `dataset/About_Plant_Performance.md` | 🎯🔧 | ✅ | Phase 1 被控对象性能上限分析 — Bode/Poisson/延迟三层约束 |

---

## demo1 — Youla 极点配置 + Q-RLS

| 文件 | 读者 | 状态 | 用途 |
|:---|:---|:---|:---|
| `demo1_RST/AboutRST.md` | 🎯 | ✅ | RST 控制器理论基础 — Diophantine 方程、灵敏度塑形 |
| `demo1_RST/About_RST_Bezout_Barrier.md` | 🎯🔧 | ✅ | **关键壁垒**: ARMAX(30,30,30,22) @48kHz Sylvester 矩阵病态 (条件数 ~1e16), bezoutd 不可行 |
| `demo1_RST/Review_demo1_Publishability.md` | 🎯 | ✅ | 发表可行性审稿 (Carmona2000 + Landau2005 方法) |

## demo2 — LQG + Youla-Q

| 文件 | 读者 | 状态 | 用途 |
|:---|:---|:---|:---|
| `demo2_LQG/AboutLQGpt1.md` | 🎯🔧 | ✅ | **主文档**: LQG-ANC 完整推导 — LQ/LQR/LQG/Kalman/ARMAX 辨识/实验 |
| `demo2_LQG/AboutLQGpt2.md` | 🔧 | 📝 | 占位: 辨识模型 ANC 问题记录 |
| `demo2_LQG/AboutLQR.md` | 🎯🔧 | ⚠️ | LQR 测试用例 (弹簧-质量-阻尼器) — 建议合并入 AboutLQGpt1 附录 |
| `demo2_LQG/Qian.md` | 🎯 | ✅ | **参考文献**: 钱梵梵 (2022) 博士论文完整推导 — 增广系统 + LQG + Youla-Q + RLS |
| `demo2_LQG/About_Q_Convergence.md` | 🎯🔧 | ✅ | **关键分析**: Q-RLS 收敛三条件 (干净残差/合格观测器/充足容量) + 裸 LQG 失败诊断 |
| `demo2_LQG/Review_demo2_Publishability.md` | 🎯 | ✅ | 发表可行性审稿 — 指出 Q 自适应未完成、无基线对比 |

## demo3 — Jafari F(z)+FIR+NLMS

| 文件 | 读者 | 状态 | 用途 |
|:---|:---|:---|:---|
| `demo3_Robust/Theory_Foundations.md` | 🎯 | ✅ | **理论基础**: 自适应控制 (Swapping Lemma/Laguerre/RLS/投影) + H∞ (DGKF/γ-迭代) |
| `demo3_Robust/Derivation_Problem1_Jafari_ContinuousCC.md` | 🎯🔧 | ✅ | 题1: 连续时间 CC — F(s) 频谱展平 + Laguerre 基 + RLS (28.9 dB) |
| `demo3_Robust/Derivation_Problem2_Jafari_DiscreteAVC.md` | 🎯🔧 | ✅ | 题2: 离散时间 AVC — FIR + H∞ 插值 + H∞优化 (25 dB) |
| `demo3_Robust/Derivation_Problem3_HinfSynthesis.md` | 🎯🔧 | ✅ | 题3: H∞ 混合灵敏度 — Bai (1997) 耳机模型 + W1/W2/W3 + hinfsyn |
| `demo3_Robust/Derivation_Problem4_RealAcousticPath.md` | 🎯🔧 | ✅ | **题4 (核心)**: ARMAX 真实声学路径 — F(z)展平 + FIR H∞插值 + 抽取RLS (23 dB) |
| `demo3_Robust/ExperimentReport_TimeVarying.md` | 🔧 | ✅ | 时变频实验报告: RLS 崩溃 → NLMS → 22 dB 全频段 |
| `demo3_Robust/Review_Problem4_Publishability.md` | 🎯 | ✅ | 题4 发表审稿: 建议 InterNoise → IEEE SPL, 核心贡献=高采样率 RLS 共线性诊断 |

## demo4 — ε-MOPSO RST

| 文件 | 读者 | 状态 | 用途 |
|:---|:---|:---|:---|
| `demo4_Robust/README.md` | 🔧 | ✅ | demo4 入口 — 目录结构、六种运行模式、快速开始 |
| `demo4_Robust/About_PSO_Foundations.md` | 🎯 | ✅ | **权威 PSO 理论**: 单目标→多目标→ε-MOPSO, 含 Jury 收敛性证明 (3h 阅读) |
| `demo4_Robust/About_RST_eMOPSO_spec.md` | 🎯🔧 | ✅ | RST+ε-MOPSO 规范 — 附录 B: ARMAX 适用性边界 |
| `demo4_Robust/About_ARMAX_Debugging.md` | 🎯🔧 | ✅ | ARMAX 调试全记录 — Toeplitz 递归放大 + ε-支配失效 |
| `demo4_Robust/About_MOEA_algorithms.md` | 🎯 | ⚠️ | 四种 MOEA 对比 (伪代码) — NSGA-II/MODE 未实现, PSO 内容与 PSO_Foundations 重复 |
| `demo4_Robust/About_Solution.md` | 🎯🔧 | ⚠️ | 柔性传动系统 RST 案例 — 引用他人工作未标注出处 |
| `demo4_Robust/Review_demo4_Publishability.md` | 🎯 | ✅ | 发表审稿: 当前不可发表, 建议教育/软件/arXiv 路线 |

## 信号激励

| 文件 | 读者 | 状态 | 用途 |
|:---|:---|:---|:---|
| `signal_excitation/AboutExcitation.md` | 🔧 | ✅ | 激励信号说明 — PRBS vs 带限白噪声, 实验方案 |

---

## 合并计划

### 已识别碎片

| 碎片组 | 涉及文件 | 处理 |
|:---|:---|:---|
| **demo2 LQR 理论** | `AboutLQR.md` (LQR 测试用例) + `AboutLQGpt1.md` (主文档) | `AboutLQR.md` 内容拉入 `AboutLQGpt1.md` 附录, 原文件加 "→ 见 AboutLQGpt1.md" 重定向 |
| **demo2 占位** | `AboutLQGpt2.md` (空) + `AboutLQGpt1.md` | `AboutLQGpt2.md` 的问题记录移入 `AboutLQGpt1.md` §待解决问题, 原文件删除 |
| **demo4 PSO 基础** | `About_PSO_Foundations.md` (权威) + `About_MOEA_algorithms.md` (重复) + `About_RST_eMOPSO_spec.md` (重复) | `About_MOEA_algorithms.md` 加已弃用声明; `About_RST_eMOPSO_spec.md` PSO 章节改为 "见 PSO_Foundations" 引用 |
| **demo4 ARMAX 失败** | `About_ARMAX_Debugging.md` + `About_RST_eMOPSO_spec.md` 附录 B | 附录 B 改为 "详细调试历程见 About_ARMAX_Debugging.md", 去重 |
| **demo4 案例** | `About_Solution.md` | 加标注 "引用来源待核实, 不可直接用于发表" |

---

## 阅读路线

### 🎯 理论研究者路线

1. `CONTEXT.md` → 术语对齐 (5min)
2. `plans/Global_Work_Plan.md` → 项目全景 (15min)
3. `dataset/About_Plant_Performance.md` → 被控对象约束 (30min)
4. **demo3**: `Theory_Foundations.md` → 公共理论 (2h) → `Derivation_Problem4` (核心, 1h) → 其余题解
5. **demo2**: `Qian.md` → LQG+Youla 框架 (1h) → `About_Q_Convergence.md` → Q 收敛条件 (30min)
6. **demo1**: `AboutRST.md` → `About_RST_Bezout_Barrier.md` (病态分析)
7. **demo4**: `About_PSO_Foundations.md` (2h) → `About_ARMAX_Debugging.md`

### 🔧 工程师路线

1. `README.md` → 项目入口 (5min)
2. `tests/README.md` → 如何跑测试 (10min)
3. `dataset/README.md` → 模型加载 (5min)
4. `plans/Global_Work_Plan.md` §Phase 2 → 测试基准设计 (10min)
5. **demo3**: `Derivation_Problem4` → 代码映射 + `ExperimentReport_TimeVarying` → 调参经验
6. **demo2**: `AboutLQGpt1.md` → 仿真设置 → `About_Q_Convergence.md` → 已知壁垒
7. `demo1_RST/About_RST_Bezout_Barrier.md` → 为什么极点配置失败
8. `demo4_Robust/About_ARMAX_Debugging.md` → 为什么 ε-MOPSO 失败
