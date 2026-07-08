# demo5_MarinoTomei — 自适应输出调节：从内模原理到非线性未知参数估计

## 概述

本目录实现了基于 **Marino & Tomei** 系列工作的自适应非线性输出调节控制器仿真，覆盖三组递进的理论主题：

| 主题                     | 理论来源                | 核心问题                                          |
| ------------------------ | ----------------------- | ------------------------------------------------- |
| **自适应输出调节** | Marino & Tomei (2011)   | 未知系统参数 + 已知频率扰动下的渐近跟踪与扰动抑制 |
| **频率自适应估计** | Marino & Tomei (2016)   | 未知频率多正弦扰动下的在线频率辨识与抵消          |
| **内模原理基础**   | Francis & Wonham (1976) | 并联 vs 级联内模结构的频域特性对比                |

核心控制问题：对二阶非线性系统

$$
\begin{aligned}
\dot{x}_1 &= x_2 + a_1 x_1 + \frac{1}{\beta}u + w_1(t) \\
\dot{x}_2 &= \frac{1}{\beta}b_2 u \\
y &= x_1
\end{aligned}
$$

在参数 $a_1, b_2, \beta$ 未知、扰动 $w_1$ 频率已知的条件下，设计自适应控制器 $u = -ke - \hat{\chi}_1 - \mu_1\hat{\theta}_1 - \mu_2\hat{\theta}_2$ 使输出 $y \to y_r = -w_2$。完整理论推导见 **[About_MarinoTomei.md](About_MarinoTomei.md)**。

## 目录结构

```
demo5_MarinoTomei/
├── README.md                              ← 本文件
├── About_MarinoTomei.md                   ← 完整理论推导（系统模型 → Lyapunov 稳定性）
├── Review_demo5_Publishability.md         ← 学术发表可行性评估
├── SIMULATION_REPORT.md                   ← ★ 仿真总结报告（本文档）
│
├── MarinoTomei_adaptive_regulation.m (integrator='ode45')                    ← ★ 自适应控制器 (ode45)
├── MarinoTomei_adaptive_regulation.m (integrator='euler')                    ← 同上 (Euler 定步长)
├── MarinoTomei_regulator_benchmark.m (integrator='rk4')                ← 自适应控制器 + 调节器方程基准 (RK4)
├── MarinoTomei_regulator_benchmark.m (integrator='ode45')              ← 同上 (ode45)
│
├── MarinoTomei_2016_adaptive_freq_est.m   ← ★ 2016 频率自适应估计 + 扰动抵消
├── MarinoTomei_2023_unstable.m            ← 2023 实验：不稳定传函系统
├── MarinoTomei_internal_model.m           ← 内模原理：并联 vs 级联传函频域分析
│
├── demos/                                 ← 辅助/教学演示
│   ├── demo_simple_adaptive.m               最简自适应控制器（教学入门）
│   ├── demo_regulator_solution.m            调节器方程求解 + 理论轨迹分析
│   ├── demo_ss_disturbance.m                扰动状态空间模型正确性验证
│   └── direct_regulator_solution.m          调节器方程数值求解函数
│
└── docs/                                  ← 文档与参考资料
    ├── MarinoTomei_2016_Derivation.md         ★ 2016 频率自适应估计完整数学推导
    ├── papers/                               ★ 所有 tex 文档 + 配套图片 (media/)
    ├── drawioNotebook/                       系统框图与调节器示意图
    └── simulinkModel/Tomei2023.mlx           Simulink 模型
```

## 仿真脚本差异对比

### 核心仿真：同一控制器，四种实现1






四个 `MarinoTomei_*.m` 脚本实现的是**同一个自适应输出调节控制器**（Marino & Tomei, 2011），差异仅在于**数值求解方式**和**是否附带调节器方程基准**：

|                                   | ode45 求解器                  | Euler/RK4 求解器            |
| --------------------------------- | ----------------------------- | --------------------------- |
| **仅控制器**                | `MarinoTomei_adaptive_regulation.m (integrator='ode45')`       | `MarinoTomei_adaptive_regulation.m (integrator='euler')`     |
| **控制器 + 调节器方程基准** | `MarinoTomei_regulator_benchmark.m (integrator='ode45')` | `MarinoTomei_regulator_benchmark.m (integrator='rk4')` |

**行对行对比：**

| 维度                 | `MarinoTomei_adaptive_regulation` (ode45)                                                  | `MarinoTomei_adaptive_regulation` (euler)      | `MarinoTomei_regulator_benchmark` (rk4)                            | `MarinoTomei_regulator_benchmark` (ode45) |
| -------------------- | ---------------------------------------------------------------------- | -------------------------- | ---------------------------------------------------- | --------------------------- |
| **控制律**     | $u = -ke - \hat{\chi}_1 - \mu_1\hat{\theta}_1 - \mu_2\hat{\theta}_2$ | ← 完全相同                | ← 完全相同                                          | ← 完全相同                 |
| **状态滤波器** | 赫尔维茨$D_{4\times4}$ + $\xi_1, \xi_2$                            | ← 完全相同                | ← 完全相同                                          | ← 完全相同                 |
| **参数投影**   | 平滑投影算子 ($\epsilon_r=0.1$)                                      | ← 完全相同                | ← 完全相同                                          | ← 完全相同                 |
| **参考控制器** | $R_c \in \mathbb{R}^{5\times5}$                                      | ← 完全相同                | 无（改用调节器方程）                                 | 无（改用调节器方程）        |
| **求解器**     | `ode45` 变步长                                                       | Euler 定步长 ($dt=0.01$) | RK4 (被控对象) + Euler (滤波器)                      | `ode45` 变步长            |
| **调节器方程** | 无                                                                     | 无                         | ✅ 求解$\Pi R = A\Pi + \frac{1}{\beta}b\Gamma + P$ | ✅ 同左                     |
| **输出图表数** | 2 张 (误差+估计 / 输出+扰动)                                           | 2 张 (同上)                | 4 张 (外加理论轨迹对比 + 切换点放大)                 | 4 张 (同左)                 |
| **运行速度**   | 快 (~2s)                                                               | 最快 (~1s)                 | 中等 (~5s，含调节器方程求解)                         | 慢 (~10s)                   |
| **适用场景**   | 日常验证                                                               | 单步调试、教学演示         | 需要与理论最优解对标                                 | 高精度对标                  |

> **选择指南**：快速验证用 `MarinoTomei_adaptive_regulation` (ode45)；需要理解每步计算用 `MarinoTomei_adaptive_regulation` (euler)（可在循环中设断点）；需要回答"控制器离理论最优还有多远"用 `MarinoTomei_regulator_benchmark` 系列。

### 其他仿真

| 脚本                                     | 与核心仿真的本质区别                                                                                                                                                  |
| ---------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `MarinoTomei_2023_unstable.m`          | **完全不同的系统模型**：被控对象是不稳定传函 $P(s)=\frac{10(s+2)}{(s+10)(s+1)(s-1)}$，控制策略为极点配置风格，与 Marino-Tomei 2011 输出调节框架无直接理论关联 |
| `MarinoTomei_internal_model.m`         | **不是控制器仿真**：它是频域分析工具，对比并联/级联谐振器组的传递函数、脉冲响应和频谱特性，服务于内模原理的直观理解                                             |
| `MarinoTomei_2016_adaptive_freq_est.m` | **不同理论问题**：2016 年论文处理的是**频率未知**的自适应估计与抵消（$\dot{\hat{\theta}} = \varepsilon \cdot w_2 \cdot y$），而核心仿真假设频率已知     |

### demos/ 辅助脚本

| 脚本                            | 用途                                  | 与主仿真的关系                                                           |
| ------------------------------- | ------------------------------------- | ------------------------------------------------------------------------ |
| `demo_simple_adaptive.m`      | 最简实现（不同 D 矩阵 + 硬边界投影）  | `MarinoTomei_adaptive_regulation.m (integrator='ode45')` 的**教学简化版**，去掉了参考控制器 $R_c$ |
| `demo_regulator_solution.m`   | 调节器方程求解 + 三区间理论轨迹可视化 | `MarinoTomei_regulator_benchmark.m` 中调节器方程求解的**独立验证版**         |
| `demo_ss_disturbance.m`       | 验证扰动状态空间模型与解析表达式一致  | 确保所有仿真中$w_1(t), w_2(t)$ 的生成是正确的                          |
| `direct_regulator_solution.m` | 调节器方程数值求解的独立函数          | 被`demo_regulator_solution.m` 调用                                     |

## 公共仿真参数与扰动场景

所有主仿真共用以下设置：

| 参数                             | 值                                           | 说明                               |
| -------------------------------- | -------------------------------------------- | ---------------------------------- |
| $a_1, b_2, \beta$              | 0.5, 0.5, 0.5                                | 未知系统参数（仿真中用于生成状态） |
| $\omega_1, \omega_2, \omega_3$ | 1, 0.5, 4 rad/s                              | 扰动频率                           |
| $k, g$                         | 7.1, 100                                     | 控制增益、参数自适应增益           |
| $D$                            | `[-2 1 0 0; -3 0 1 0; -2 0 0 1; -1 0 0 0]` | 赫尔维茨状态滤波器                 |
| 仿真时长                         | 150 s                                        | 覆盖三个扰动区间                   |

| 时间区间       | $w_1(t)$              | $w_2(t)$    | 测试目的                   |
| -------------- | ----------------------- | ------------- | -------------------------- |
| $[0, 50)$    | $\sin t + 0.2\sin 4t$ | $\sin 0.5t$ | 双频扰动 + 正弦跟踪        |
| $[50, 100)$  | $\sin t$              | $\sin 0.5t$ | 扰动频率切换瞬态           |
| $[100, 150]$ | $\sin t$              | $0$         | 纯调节问题（参考信号消失） |

## 阅读路径

1. **📖 [About_MarinoTomei.md](About_MarinoTomei.md)** — 2011 理论推导（系统模型 → 状态滤波器 → 投影算子 → Lyapunov 稳定性）
2. **📖 [docs/MarinoTomei_2016_Derivation.md](docs/MarinoTomei_2016_Derivation.md)** — 2016 理论推导（未知频率多正弦干扰的自适应输出反馈补偿，含平均化分析）
3. **🔬 [demos/demo_simple_adaptive.m](demos/demo_simple_adaptive.m)** — 最简实现，快速理解控制律结构
4. **🔬 [demos/demo_ss_disturbance.m](demos/demo_ss_disturbance.m)** — 理解扰动如何建模为状态空间
5. **★ [MarinoTomei_adaptive_regulation.m (integrator='ode45')](MarinoTomei_adaptive_regulation.m (integrator='ode45'))** — 完整自适应控制器（推荐首选）
6. **📊 [MarinoTomei_regulator_benchmark.m (integrator='rk4')](MarinoTomei_regulator_benchmark.m (integrator='rk4'))** — 实际轨迹 vs 理论最优轨迹对比
7. **📝 [Review_demo5_Publishability.md](Review_demo5_Publishability.md)** — 学术发表可行性评估

## 快速开始

```matlab
% 在 demo5_MarinoTomei/ 目录下运行:

% 完整自适应控制器（推荐首选）
MarinoTomei_adaptive_regulation

% 单步调试版（可在循环内设断点）
MarinoTomei_adaptive_regulation

% 理论最优基准对比
MarinoTomei_regulator_benchmark

% 不稳定系统实验
MarinoTomei_2023_unstable

% 内模原理频域分析
MarinoTomei_internal_model

% 频率估计器 (2016)
MarinoTomei_2016_adaptive_freq_est
```

### 预期输出

- **输出跟踪**: $y = x_1$ 渐近收敛到 $y_r = -w_2$
- **参数估计**: $\hat{\theta}_1 \to 1.25$, $\hat{\theta}_2 \to 0.25$
- **区间切换**: $t=50$ s 和 $t=100$ s 处出现瞬态，之后快速恢复
- **bench 系列额外输出**: 实际轨迹与 $\Pi w$ 理论最优轨迹的偏差曲线

## 依赖

- MATLAB R2020b+
- Control System Toolbox（`MarinoTomei_2023_unstable.m`）
- Signal Processing Toolbox（`MarinoTomei_internal_model.m`）

## 参考文献

1. Marino, R. & Tomei, P. (2011). "An adaptive learning control for nonlinear systems with unknown control directions." *Automatica*.
2. Marino, R. & Tomei, P. (2016). "Adaptive disturbance cancellation for nonlinear systems with unknown frequencies." *IEEE TAC*.
3. Marino, R. & Tomei, P. (2023). "Adaptive output regulation for nonlinear systems with unknown parameters." (见 `docs/papers/`)
4. Francis, B.A. & Wonham, W.M. (1976). "The internal model principle of control theory." *Automatica*, 12(5), 457–465.
