# PSO 基础：从单目标到 ε-MOPSO 的完整数学推导

> **本文档为下列读者定制**：自动化/控制工程背景，对优化理论有模糊概念，听说过 PSO 但不懂细节，完全未接触过多目标优化（Pareto 等概念），且需要严格的收敛性证明。如果你不是上述读者，绪论中的时间估算可能不适用于你。

---

## 绪论：为你定制的学习路径

### 0.1 你的起点

根据你告诉我的，你在打开这份文档之前的状态是：

| 领域 | 你的状态 | 这意味着 |
|:---|:---|:---|
| 优化问题的形式化表述 | 有模糊概念（知道 $\min f(x)$ 这种形式） | §1 需要为你把模糊印象**精确化**——定义、符号、术语的严格写法 |
| 粒子群优化（PSO） | 听过但不懂细节 | §2–§3 是全新的，从直觉 → 数学模型逐步搭建 |
| 多目标优化（Pareto 等） | 完全没学过 | §5 是你和这份文档之间**最重要的知识增量**——需要从零建立多目标优化的基本语言 |
| 收敛性证明 | 需要严格证明（Lemma/Theorem/Proof） | §4.3–4.5 和 §7 不是"可选阅读"——它们是为你写的核心内容 |

### 0.2 你为什么需要这份文档

你面前有一个具体的工程问题：

> 给定被控对象 $G(z^{-1}) = z^{-d} B(z^{-1}) / A(z^{-1})$。它有一个不稳定零点 $\beta = -1.79$（模长 1.79，单位圆外，非最小相位）。根据 Zames-Francis 积分等式，输出灵敏度函数 $S_{yp}(e^{j\omega})$ 在物理上**必须**满足：
>
> $$
> \int_{-\pi}^{\pi} \log|S_{yp}(e^{j\omega})| \cdot \frac{|\beta|^2 - 1}{|e^{j\omega} - \beta|^2} d\omega = \pi \log|B_a^{-1}(\beta)|
> $$
>
> 这不是设计指南——这是数学恒等式。其物理后果（水床效应）：你不能在所有频率上同时压低灵敏度。与此同时，时间延迟 $d$ 引入了**第二个独立的**积分约束。两个约束互为矛盾——改善一个必然恶化另一个。
>
> 再加上频域鲁棒性要求：$|S_{yp}(e^{j\omega})|$ 在 $\omega \in [0, \pi]$ 上必须钻进模值裕度和延迟裕度划定的"安全套筒"（200 个频率点 → 400 个不等式约束）。

现在逐一审视你手头的工具：

- **极点配置** → 让你把极点放在想要的位置，但不处理灵敏度函数在**所有频率上**的形状约束。
- **$H_\infty$ 优化** → 需要先选加权函数 $W(z)$，但 $W(z)$ 正是你要通过优化来**确定**的东西——鸡和蛋。
- **梯度下降** → 目标函数是 $\int_{-\pi}^{\pi} \log|S_d(e^{j\omega}, \boldsymbol{\theta})| \cdot K(\omega) d\omega$。频率响应对多项式根 $\boldsymbol{\theta}$ 的导数没有闭式解，搜索空间高度多峰。

**三个工具全部失效。** PSO 就是第四个工具。一群粒子在 $\boldsymbol{\theta}$ 空间中飞行，每颗粒子携带一组候选的 $(X, Y)$ 多项式根，$\varepsilon$-MOPSO 的 Archive 收敛到 Pareto 前沿——你从中选一个拐点解，通过 Bezout 方程反解出 $R, S, T$。

### 0.3 知识衔接：你已有的控制背景不是负担，是武器

你不是从零学 PSO。下表告诉你每一章的新概念和你已有知识的**精确对应关系**：

| 本文档章节 | 你已会的（控制背景） | 要新学的（PSO 概念） | 新旧知识的精确映射 |
|:---|:---|:---|:---|
| §1 优化问题表述 | LQR 的代价函数 $J = \int (x^T Q x + u^T R u)dt$；你在调 $Q$ 和 $R$ 时已经在做优化 | 一般优化问题的标准形式：决策变量 $\mathbf{x}$、目标函数 $f(\mathbf{x})$、约束 $g_j(\mathbf{x}) \le 0$ | 目标函数就是代价函数的推广。LQR 是优化，RST 参数整定也是优化——只是 $f$ 不再是二次型 |
| §2–§3 PSO 运动学 | 离散状态空间模型 $\mathbf{x}_{k+1} = \mathbf{A}\mathbf{x}_k + \mathbf{B}\mathbf{u}_k$ | 粒子速度和位置更新方程 $\mathbf{v}_{k+1} = w\mathbf{v}_k + c_1 r_1 (\mathbf{p} - \mathbf{x}_k) + c_2 r_2 (\mathbf{g} - \mathbf{x}_k)$ | PSO 本质上就是一个**受随机激励的二阶离散动力学系统**。把 $\mathbf{B}\mathbf{u}_k$ 替换为认知项 + 社会项，你就得到了 PSO |
| §4 收敛性分析 | Jury 稳定性判据：离散系统 $\mathbf{z}_{k+1} = \mathbf{M}\mathbf{z}_k$ 稳定 $\iff$ $\mathbf{M}$ 的所有特征值在单位圆内 $\iff$ 谱半径 $\rho(\mathbf{M}) < 1$ | 确定性 PSO 的收敛条件：$w < 1$，$\varphi > 0$，$2w - \varphi + 2 > 0$ | 剥离随机性后，PSO 简化为 $\mathbf{z}_{k+1} = \mathbf{M}\mathbf{z}_k$。**你用 Jury 判据就能证明它的收敛性**——不需要学任何新工具 |
| §5 MOPSO | 标量代价函数 → 唯一最优解（你一直是这样做的） | Pareto 支配 $\prec$、Pareto 前沿 $\mathcal{PF}$、外部 Archive $\mathcal{A}$ | Bode 图中的增益裕度和相位裕度本身就是一对矛盾目标——你已经理解"两个指标不能同时最优"，只是没有用"Pareto"这个名字叫它 |
| §6–§7 ε-MOPSO | ADC/DAC 的量化（quantization）：连续信号 → 离散电平；有限字长效应的精度权衡 | ε-支配 $\prec_\varepsilon$、Box 索引：连续 Pareto 前沿 → 离散网格 | ε 就是目标空间的"量化步长"。你完全理解"精度 vs. 存储开销"——这就是 $\varepsilon$ 参数的含义 |
| §8 衔接 RST | 极点配置法、Bezout 方程、灵敏度函数整形 | PSO 粒子到 RST 决策变量的映射 | 算法输出的 $\boldsymbol{\theta}$ 就是 $X$ 和 $Y$ 多项式的根——代回 Bezout 方程就得到 $R, S, T$ |

### 0.4 你的阅读路线

根据 §0.1 的基线，**你不需要跳过任何章节**。但阅读深度不同：

| 章节 | 你的策略 | 预计时间 |
|:---|:---|:---|
| **§1** 优化问题的标准表述 | 精读。你在"模糊概念"到"精确定义"之间有一段路。特别注意定义 1.1–1.3 中每个符号的**严格含义**——这是后续所有讨论的语言基础 | 20 min |
| **§2–§3** PSO 的直觉与数学模型 | 精读 §2（直觉）→ 通读 §3（方程），初读时暂不纠缠 Algorithm 1 的逐行伪代码。**重点**：式 (1)(2) 必须能闭眼写出来 | 25 min |
| **§4** 参数与收敛性 | §4.1–4.2 快速通读（参数直觉，10 min）→ **精读 §4.3–4.5**（这是你要求的部分：Lem 4.1 → Lem 4.2 → Thm 4.1 → Cor 4.1 → Thm 4.2，每一步都跟下来）。特别注意 **Jury 判据在 Thm 4.1 证明中的用法**——这是你已有知识和新知识的交汇点 | 45 min |
| **§5** MOPSO | 这是你**从未接触过的新领域**。精读全文。§5.1 的三个定义（支配、最优解、前沿）是多目标优化的基本语法——不理解它们，§6 根本读不下去 | 30 min |
| **§6** ε-支配与 ε-MOPSO | 精读。§6.2 的定义 (19) 和 §6.3 的 Box 索引 (20) 是核心。Algorithm 2 和 3 的伪代码是 `eMOPSO_core.m` 的精确镜像 | 30 min |
| **§7** ε-MOPSO 理论性质 | 精读。Thm 7.1（Archive 有界性）直接解释了为什么 $\varepsilon$ 是你需要在代码里调的参数。证明只用了基本的组合计数 | 15 min |
| **§8** 衔接 RST | 精读。这是所有理论知识落回到你控制问题的最后一公里。读完后立即打开 `run_RST_eMOPSO.m` 对照 | 15 min |

**总预计时间**：约 3 小时（建议分两次：第一次 §1–§5，第二次 §6–§8）。

### 0.5 文档标签说明

正文中用两种标签辅助导航：

| 标签 | 含义 | 对你的适用性 |
|:---|:---|:---|
| 🎯 **控制工程重点** | 与 RST 设计直接关联 | 全部必读（你的场景） |
| 📚 深入阅读 | 理论证明与形式化分析 | 对你**同样必读**（因为你要求了严格证明）——此标签仅用于提醒那些只关心"代码能不能跑"的读者可以跳过。你的阅读路线不受此标签限制 |

---

## 目录

1. [优化问题的标准表述](#1-优化问题的标准表述)
2. [PSO 的起源与直觉](#2-pso-的起源与直觉)
3. [标准 PSO 的数学模型](#3-标准-pso-的数学模型)
4. [参数分析与收敛性理论](#4-参数分析与收敛性理论)
5. [从单目标到多目标：MOPSO](#5-从单目标到多目标mopso)
6. [ε-支配与 ε-MOPSO](#6-ε-支配与-ε-mopso)
7. [ε-MOPSO 的理论性质](#7-ε-mopso-的理论性质)
8. [与 RST 控制器的衔接](#8-与-rst-控制器的衔接)
9. [参考文献](#9-参考文献)

---

## 1. 优化问题的标准表述

> 🎯 **控制工程重点**。你设计控制器时已经在做优化了——你选极点位置来让阶跃响应"好看"，你调 $Q$ 和 $R$ 矩阵来让 LQR 代价函数最小。本节把这种直觉形式化。

### 1.1 无约束连续优化

**定义 1.1（无约束单目标优化问题）**
给定目标函数 $f: \mathbb{R}^D \to \mathbb{R}$，寻找全局最优解 $\mathbf{x}^*$：

$$
\mathbf{x}^* = \underset{\mathbf{x} \in \mathbb{R}^D}{\arg\min}\; f(\mathbf{x})
$$

即 $\forall \mathbf{x} \in \mathbb{R}^D$，有 $f(\mathbf{x}^*) \le f(\mathbf{x})$。

> 🎯 **与你已知道的类比**：LQR 就是求解 $\mathbf{u}^* = \arg\min_{\mathbf{u}} \int (x^T Q x + u^T R u) dt$。这里的 $f(\mathbf{x})$ 可以是任何函数——你不必知道它的导数，甚至不必知道它的表达式（只要给定输入能算出输出）。

> **为什么不用梯度下降？** 在 RST 优化问题中，$f(\boldsymbol{\theta})$ 涉及频率积分 $\int_{-\pi}^{\pi} \log|S_d(e^{j\omega}, \boldsymbol{\theta})| \cdot K(\omega) d\omega$。这个函数对 $\boldsymbol{\theta}$ 的梯度没有闭式解，且存在大量局部极小值（因为 $X$ 和 $Y$ 多项式的根可以以不同方式组合出相似的频率响应）。

### 1.2 带约束优化

**定义 1.2（约束单目标优化问题）**

$$
\begin{aligned}
\underset{\mathbf{x} \in \mathbb{R}^D}{\min}\quad & f(\mathbf{x}) \\
\text{s.t.}\quad & g_j(\mathbf{x}) \le 0, \quad j = 1, \dots, J \\
& h_q(\mathbf{x}) = 0, \quad q = 1, \dots, Q \\
& \mathbf{x} \in \mathcal{S} = [\mathbf{x}_{\min}, \mathbf{x}_{\max}]^D
\end{aligned}
$$

其中 $\mathcal{S}$ 是搜索空间的超矩形边界（Box Constraint）。PSO 天然支持此类边界约束。

> 🎯 **RST 问题中的约束**：$g_j$ 正是 §0.2 中提到的"灵敏度函数必须钻进的安全套筒"——模值裕度下界和延迟裕度上界。每个频率点 $\omega_k$ 都贡献一对约束 $g_1(\omega_k) \le 0$ 和 $g_2(\omega_k) \le 0$。典型的频率离散化（200 点）意味着 **400 个约束**——用传统方法（如 SQP）处理非常繁琐，而 PSO 的惩罚函数法将其轻松吸收。

### 1.3 多目标优化（MOP）

**定义 1.3（多目标优化问题）**

$$
\begin{aligned}
\underset{\mathbf{x} \in \mathcal{S}}{\min}\quad & \mathbf{F}(\mathbf{x}) = \big(f_1(\mathbf{x}), f_2(\mathbf{x}), \dots, f_M(\mathbf{x})\big)^T \\
\text{s.t.}\quad & g_j(\mathbf{x}) \le 0, \quad j = 1, \dots, J
\end{aligned}
$$

多目标优化的核心矛盾：**不存在单一解同时最小化所有目标**。

> 🎯 **RST 问题中的多目标性**：每个不稳定零点 $\beta_i$ 对应一个 Zames-Francis 积分等式——这是一个独立的目标 $f_i$（你必须让灵敏度函数在 $\beta_i$ 对应的 Poisson 核加权积分上匹配理论值）。时间延迟 $d$ 对应另一个目标 $f_d$。这些目标**天然矛盾**：改善零点 $\beta_1$ 的匹配精度会恶化延迟积分的匹配精度，反之亦然。你不可能同时完美满足——只能找 Pareto 前沿上的折衷。

---

## 2. PSO 的起源与直觉

> 🎯 **控制工程重点**。本节零数学，纯直觉。理解"一群粒子为何能飞向最优解"是后续所有数学推导的地基。

### 2.1 仿生学动机

PSO 由 Kennedy 和 Eberhart 于 1995 年提出[1]，灵感来自鸟群和鱼群的集体觅食行为：

- 个体（鸟）在空间中飞行，每个个体知道自身到过的最佳位置
- 个体之间共享信息，知道整个群体到过的最佳位置
- 每个个体根据**自身经验**和**群体经验**来调整飞行方向与速度

### 2.2 从生物学到数学的类比

| 生物学概念 | 数学映射 |
|:---|:---|
| 鸟群 | 粒子群（Population of Particles） |
| 鸟的位置 | 决策变量向量 $\mathbf{x}_i \in \mathbb{R}^D$ |
| 鸟的速度 | 速度向量 $\mathbf{v}_i \in \mathbb{R}^D$ |
| 食物质量 | 目标函数值 $f(\mathbf{x}_i)$ |
| 个体记忆 | 个体最优 $\mathbf{pbest}_i$ |
| 群体知识 | 全局最优 $\mathbf{gbest}$ |

### 2.3 核心思想的形式化

每个粒子在每次迭代中整合三股"力量"：

1. **惯性**（Inertia）：保持原有运动趋势
2. **认知**（Cognitive）：被自己的最佳经验吸引
3. **社会**（Social）：被群体的最佳经验吸引

$$
\text{新速度} = \text{惯性项} + \text{认知项} + \text{社会项}
$$

---

## 3. 标准 PSO 的数学模型

> 🎯 **控制工程重点**。式 (1)(2) 就是你需要理解的全部——这两个方程定义了 PSO 的全部行为。如果你熟悉离散状态空间模型 $\mathbf{x}_{k+1} = \mathbf{A}\mathbf{x}_k + \mathbf{B}\mathbf{u}_k$，PSO 的更新方程本质上就是同一个东西：一个受随机激励的二阶离散动力学系统。

### 3.1 运动学方程

设 $N$ 个粒子在 $D$ 维空间中搜索。第 $i$ 个粒子在迭代 $k$ 时的状态由两个 $D$ 维向量描述：

- 位置 $\mathbf{x}_k^i \in \mathbb{R}^D$
- 速度 $\mathbf{v}_k^i \in \mathbb{R}^D$

**定义 3.1（标准 PSO 更新方程）**

$$
\boxed{\mathbf{v}_{k+1}^i = w \mathbf{v}_k^i + c_1 r_{1,k}^i \odot (\mathbf{pbest}_k^i - \mathbf{x}_k^i) + c_2 r_{2,k}^i \odot (\mathbf{gbest}_k - \mathbf{x}_k^i)} \tag{1}
$$

$$
\boxed{\mathbf{x}_{k+1}^i = \mathbf{x}_k^i + \mathbf{v}_{k+1}^i} \tag{2}
$$

其中：
- $\odot$：逐元素乘法（Hadamard product）
- $w \in [0, 1]$：**惯性权重**（Inertia Weight）
- $c_1, c_2 > 0$：**加速系数**（Acceleration Coefficients）
- $r_{1,k}^i, r_{2,k}^i \sim U(0,1)$：独立均匀随机变量（逐维独立采样）
- $\mathbf{pbest}_k^i$：粒子 $i$ 的历史最优位置
- $\mathbf{gbest}_k$：全局最优位置

### 3.2 pbest 与 gbest 的更新

**定义 3.2（pbest 更新规则）**

$$
\mathbf{pbest}_{k+1}^i = \begin{cases}
\mathbf{x}_{k+1}^i, & \text{if } f(\mathbf{x}_{k+1}^i) < f(\mathbf{pbest}_k^i) \\
\mathbf{pbest}_k^i, & \text{otherwise}
\end{cases} \tag{3}
$$

**定义 3.3（gbest 更新规则）**

$$
\mathbf{gbest}_{k+1} = \underset{i \in \{1,\dots,N\}}{\arg\min}\; f(\mathbf{pbest}_{k+1}^i) \tag{4}
$$

### 3.3 算法伪代码

```
Algorithm 1: Standard PSO
──────────────────────────────────────────
Input:  f: ℝᴰ → ℝ (目标函数)
        N (种群大小), Kₘₐₓ (最大迭代)
        w, c₁, c₂ (参数)
        lb, ub (搜索边界)
Output: gbest (找到的最优解)

1: for i = 1 to N do
2:     xᵢ ~ Uniform(lb, ub)       ▷ 随机初始化位置
3:     vᵢ ~ Uniform(-|ub-lb|, |ub-lb|)  ▷ 随机初始化速度
4:     pbestᵢ ← xᵢ
5: end for
6: gbest ← argmin_i f(pbestᵢ)
7:
8: for k = 1 to Kₘₐₓ do
9:     for i = 1 to N do
10:        r₁, r₂ ~ U(0,1)
11:        vᵢ ← w·vᵢ + c₁·r₁·(pbestᵢ - xᵢ) + c₂·r₂·(gbest - xᵢ)
12:        xᵢ ← xᵢ + vᵢ
13:        ▷ 边界处理
14:        xᵢ ← clamp(xᵢ, lb, ub)
15:        if f(xᵢ) < f(pbestᵢ) then
16:            pbestᵢ ← xᵢ
17:        end if
18:    end for
19:    gbest ← argmin_i f(pbestᵢ)
20: end for
```

### 3.4 速度钳制（Velocity Clamping）

为防止粒子"飞出"搜索空间，实践中常对速度幅值设限：

**定义 3.4（速度钳制）**

$$
v_{i,d} \leftarrow \text{clamp}(v_{i,d}, -V_{\max}, V_{\max}), \quad V_{\max} = \gamma \cdot (ub_d - lb_d)
$$

其中 $\gamma \in (0, 1]$，典型值为 $0.1 \sim 0.5$。

---

## 4. 参数分析与收敛性理论

> 🎯 **控制工程重点（前半）**，📚 深入阅读（后半）。§4.1–4.2 告诉你参数怎么选，**直接对应你运行 `run_RST_eMOPSO.m` 时需要调的 `w_max, w_min, c1, c2`**。§4.3–4.5 的收敛性证明使用你熟悉的离散系统分析工具（特征值、Jury 判据、谱半径）——对你来说这是必读内容。

### 4.1 各参数的直观角色

| 参数 | 太大 | 太小 | 推荐范围 |
|:---|:---|:---|:---|
| $w$ | 全局搜索过强，不收敛 | 快速陷入局部最优 | $[0.4, 0.9]$ |
| $c_1$ | 粒子"自私"，过于分散 | 忽略自身经验 | $[1.5, 2.5]$ |
| $c_2$ | 粒子"从众"，早熟收敛 | 忽略群体经验 | $[1.5, 2.5]$ |

### 4.2 自适应惯性权重

**定义 4.1（线性递减惯性权重）**

$$
w_k = w_{\max} - (w_{\max} - w_{\min}) \cdot \frac{k}{K_{\max}} \tag{5}
$$

直觉：早期 $w$ 较大，鼓励探索；后期 $w$ 较小，鼓励精化。典型设定：$w_{\max} = 0.9, w_{\min} = 0.4$。

### 4.3 确定性 PSO 的收敛性分析

> 📚 深入阅读（对你必读）。以下分析将 PSO 的随机项剥离后，得到一个二阶线性差分方程。收敛条件通过 **Jury 稳定性判据**导出——这恰好是你学离散控制系统时用过的工具。

为进行理论分析，我们首先**剥离随机性**，考察确定性版本的 PSO。

**定义 4.2（确定性 PSO）**

将 $r_{1,k}^i$ 和 $r_{2,k}^i$ 替换为其期望值 $\frac{1}{2}$，且假设 $\mathbf{pbest}$ 和 $\mathbf{gbest}$ 固定（即粒子已找到当前最优区域，仅考察局部收敛行为）：

$$
\mathbf{v}_{k+1} = w \mathbf{v}_k + \frac{c_1}{2}(\mathbf{p} - \mathbf{x}_k) + \frac{c_2}{2}(\mathbf{g} - \mathbf{x}_k) \tag{6}
$$

其中 $\mathbf{p} = \mathbf{pbest}$，$\mathbf{g} = \mathbf{gbest}$ 为常数向量。

**定义 4.3（一维简化）**

由于各维度解耦，分析一维情况即足矣。令 $p, g, x_k, v_k \in \mathbb{R}$，并定义凝聚最优：

$$
\phi = \frac{c_1 p + c_2 g}{c_1 + c_2} \tag{7}
$$

记 $\varphi = \frac{c_1 + c_2}{2}$。式 $(6)$ 重写为：

$$
v_{k+1} = w v_k + \varphi(\phi - x_k) \tag{8}
$$

连同 $x_{k+1} = x_k + v_{k+1}$，构成二阶线性差分方程组。

**Lemma 4.1（状态空间表示）**

定义状态向量 $\mathbf{z}_k = \begin{bmatrix} x_k - \phi \\ v_k \end{bmatrix}$。确定性 PSO 的动力学可写为：

$$
\mathbf{z}_{k+1} = \mathbf{M} \mathbf{z}_k, \quad \mathbf{M} = \begin{bmatrix}
1 - \varphi & w \\
-\varphi & w
\end{bmatrix} \tag{9}
$$

*证明*：

由 $x_{k+1} = x_k + v_{k+1} = x_k + w v_k + \varphi(\phi - x_k)$，

$
\begin{aligned}
x_{k+1} - \phi &= x_k - \phi + w v_k - \varphi(x_k - \phi) \\
&= (1 - \varphi)(x_k - \phi) + w v_k
\end{aligned}
$

$
v_{k+1} = w v_k - \varphi(x_k - \phi)
$

写成矩阵形式即得 $(9)$。$\square$

**Lemma 4.2（矩阵 $\mathbf{M}$ 的特征值）**

特征方程为：

$$
\lambda^2 - (1 - \varphi + w)\lambda + w = 0 \tag{10}
$$

两个特征值为：

$$
\lambda_{1,2} = \frac{1 - \varphi + w \pm \sqrt{(1 - \varphi + w)^2 - 4w}}{2} \tag{11}
$$

*证明*：$\det(\mathbf{M} - \lambda\mathbf{I}) = (1-\varphi-\lambda)(w-\lambda) + \varphi w = \lambda^2 - (1-\varphi+w)\lambda + w = 0$。$\square$

**Theorem 4.1（确定性 PSO 的收敛条件）**

确定性 PSO 的粒子位置 $x_k$ 收敛到凝聚点 $\phi$ 的**充要条件**是矩阵 $\mathbf{M}$ 的谱半径 $\rho(\mathbf{M}) < 1$。这等价于：

$$
\boxed{\begin{cases}
w < 1 \\
\varphi > 0 \\
2w - \varphi + 2 > 0
\end{cases}} \tag{12}
$$

即参数 $(w, \varphi)$ 必须落在以下三角形区域内：

$$
\varphi \in (0, 2w + 2), \quad w \in (-1, 1)
$$

当 $w \ge 0$（实践中总是如此），条件简化为：

$$
\varphi < 2(w + 1) \quad \text{且} \quad w < 1
$$

*证明*：

对二维离散线性系统 $\mathbf{z}_{k+1} = \mathbf{M}\mathbf{z}_k$，$\mathbf{z}_k \to \mathbf{0}$ 当且仅当 $\rho(\mathbf{M}) = \max\{|\lambda_1|, |\lambda_2|\} < 1$。

根据 Jury 稳定性判据（或直接分析特征值），特征多项式 $P(\lambda) = \lambda^2 + a_1\lambda + a_0$（其中 $a_1 = -(1-\varphi+w), a_0 = w$）的所有根在单位圆内的充要条件为：

1. $P(1) > 0 \implies 1 - (1-\varphi+w) + w > 0 \implies \varphi > 0$
2. $(-1)^2 P(-1) > 0 \implies 1 + (1-\varphi+w) + w > 0 \implies 2w - \varphi + 2 > 0$
3. $|a_0| < 1 \implies |w| < 1$

三条联立即得 $(12)$。$\square$

**Corollary 4.1（收敛到 pbest/gbest 加权平均）**

由式 $(7)$，粒子收敛的位置 $\phi$ 正是 $\mathbf{pbest}$ 和 $\mathbf{gbest}$ 的 $c_1 : c_2$ 加权平均。这意味着：

- 若 $c_1 \gg c_2$：粒子倾向于收敛到自身历史最优（强调探索）
- 若 $c_2 \gg c_1$：粒子倾向于收敛到全局最优（强调开发）

---

### 4.4 引入随机性后的收敛行为

**定义 4.4（随机 PSO 的一阶稳定条件）**[3]

考虑 $r_1, r_2 \sim U(0,1)$。对式 $(1)$ 取期望，得到与确定性分析相同的形式但 $\bar{\varphi} = \frac{c_1 + c_2}{2}$。因此，**期望意义下**的收敛条件同 Theorem 4.1。

**Theorem 4.2（广义收敛条件）**[2]

对于标准 PSO（含随机项），$\mathbb{E}[x_k] \to \phi$ 的充分条件为：

$$
w > \frac{1}{2}(c_1 + c_2) - 1 \quad \text{且} \quad w < 1 \tag{13}
$$

这被称为 **Trelea 收敛准则**[3]。

### 4.5 压缩因子 PSO（Clerc-Kennedy 变体）

**定义 4.5（压缩因子）**

Clerc 和 Kennedy[2] 提出了带压缩因子的 PSO，保证收敛性而无需速度钳制：

$$
\mathbf{v}_{k+1}^i = \chi \left[ \mathbf{v}_k^i + c_1 r_1 (\mathbf{pbest}_k^i - \mathbf{x}_k^i) + c_2 r_2 (\mathbf{gbest}_k - \mathbf{x}_k^i) \right] \tag{14}
$$

其中压缩因子：

$$
\chi = \frac{2}{\left|2 - \varphi - \sqrt{\varphi^2 - 4\varphi}\right|}, \quad \varphi = c_1 + c_2 > 4 \tag{15}
$$

典型参数：$c_1 = c_2 = 2.05$，$\varphi = 4.1$，得 $\chi \approx 0.7298$。此参数组合等价于标准 PSO 取 $w = 0.7298, c_1 = c_2 = 1.4962$。

---

## 5. 从单目标到多目标：MOPSO

> 🎯 **控制工程重点**。这是最关键的一步——你的 RST 问题天生就是多目标的（§0.2 中解释过）。从单目标 PSO 到 MOPSO，你只需要理解三个概念：Pareto 支配、外部 Archive、gbest 选择策略。注意：你告诉我你完全没学过这部分——所以请务必精读。

### 5.1 多目标优化的基本概念

**定义 5.1（Pareto 支配）**

解 $\mathbf{x}_1$ **Pareto-支配** $\mathbf{x}_2$（记作 $\mathbf{x}_1 \prec \mathbf{x}_2$），当且仅当：

1. $\forall m \in \{1, \dots, M\},\; f_m(\mathbf{x}_1) \le f_m(\mathbf{x}_2)$
2. $\exists k \in \{1, \dots, M\},\; f_k(\mathbf{x}_1) < f_k(\mathbf{x}_2)$

> 通俗地说：$\mathbf{x}_1$ 在所有目标上都不差于 $\mathbf{x}_2$，且在至少一个目标上严格更优。

**定义 5.2（Pareto 最优解）**

解 $\mathbf{x}^* \in \mathcal{S}$ 称为 **Pareto 最优解**，若不存在 $\mathbf{x} \in \mathcal{S}$ 使得 $\mathbf{x} \prec \mathbf{x}^*$。

> 通俗地说：你无法在不牺牲某个目标的情况下改善另一个目标——这就是"最优"的含义。

**定义 5.3（Pareto 前沿）**

所有 Pareto 最优解在目标空间的像构成 **Pareto 前沿（Pareto Front, PF）**：

$$
\mathcal{PF} = \{\mathbf{F}(\mathbf{x}^*) \mid \mathbf{x}^* \text{ 是 Pareto 最优解}\}
$$

### 5.2 单目标 PSO → 多目标 MOPSO 的关键改造

| 单目标 PSO 组件 | 多目标 MOPSO 改造 |
|:---|:---|
| **pbest 更新** | 若新解支配原 pbest，替换；若互不支配，随机选择 |
| **gbest 选择** | 不再是"唯一的全局最优"，而是从外部存档中按多样性策略挑选 |
| **外部存档（Archive）** | 存储搜索过程中找到的所有非支配解 |

### 5.3 MOPSO pbest 更新规则

**定义 5.4（MOPSO pbest 更新）**

$$
\mathbf{pbest}_{k+1}^i = \begin{cases}
\mathbf{x}_{k+1}^i, & \text{if } \mathbf{x}_{k+1}^i \prec \mathbf{pbest}_k^i \\
\mathbf{pbest}_k^i, & \text{if } \mathbf{pbest}_k^i \prec \mathbf{x}_{k+1}^i \\
\text{RandomChoice}(\mathbf{x}_{k+1}^i, \mathbf{pbest}_k^i), & \text{otherwise（互不支配）}
\end{cases} \tag{16}
$$

### 5.4 外部存档（External Archive）

**定义 5.5（Archive 的支配更新规则）**

设 $\mathcal{A}$ 为存档集合。当新解 $\mathbf{x}$ 到来时：

1. 若 $\exists \mathbf{a} \in \mathcal{A}$ 使 $\mathbf{a} \prec \mathbf{x}$：丢弃 $\mathbf{x}$（被存档支配）
2. 否则，将 $\mathbf{x}$ 加入 $\mathcal{A}$，并删除 $\mathcal{A}$ 中所有被 $\mathbf{x}$ 支配的解
3. 若 $|\mathcal{A}| > A_{\max}$：执行**截断**（删除最拥挤区域的解）

**截断机制的多样性维护**：MOPSO 采用**自适应网格（Adaptive Grid）**。将目标空间划分为 $K \times K \times \cdots \times K$（$M$ 维）的超立方体网格，优先删除粒子最密集网格中的解。

### 5.5 gbest 选择：Sigma 方法

**定义 5.6（Sigma 向量）**[4]

对于 $M$ 维目标向量 $\mathbf{f}$，其 Sigma 向量 $\boldsymbol{\sigma} \in \mathbb{R}^{M-1}$ 定义为：

$$
\sigma_m(\mathbf{f}) = \frac{f_1^2 - f_{m+1}^2}{\sum_{j=1}^M f_j^2}, \quad m = 1, \dots, M-1 \tag{17}
$$

或更简洁地（本文实现所用版本）：

$$
\sigma_m(\mathbf{f}) = f_1^2 - f_{m+1}^2, \quad m = 1, \dots, M-1 \tag{18}
$$

**gbest 选择**：为当前粒子计算 $\boldsymbol{\sigma}$，在 Archive 中寻找 $\boldsymbol{\sigma}$ 欧氏距离最近的解作为该粒子的 $\mathbf{gbest}$。这确保每个粒子被引导至 Pareto 前沿上离其"最近"的区域。

---

## 6. ε-支配与 ε-MOPSO

> 🎯 **控制工程重点**。ε 是你运行 `run_RST_eMOPSO.m` 时看到的最重要的参数之一（`options.epsilon`）。它控制 Pareto 前沿的"分辨率"——$\varepsilon$ 越小，解越密但计算越慢；$\varepsilon$ 越大，解越稀但跑得快。理解 ε-支配的定义和 Box 索引机制，你就理解了为什么这个算法能在有限时间内给出有限个解。

### 6.1 动机：为什么要引入 ε-支配？

标准 MOPSO 的 Archive 面临两个问题：

1. **无限增长**：理论上 Pareto 最优解集合可以是无穷大（连续前沿）
2. **支配压力衰减**：当 Archive 中解的数量趋于无穷时，新解几乎不可能支配任何已有解，Archive 更新退化

**ε-支配** 提供了一种"放松"的支配概念，允许算法在用户指定的精度 $\varepsilon$ 下工作，同时自然限制 Archive 大小。

### 6.2 ε-支配的严格定义

**定义 6.1（ε-支配）**[5]

给定 $\varepsilon > 0$，解 $\mathbf{x}_1$ **ε-支配** $\mathbf{x}_2$（记作 $\mathbf{x}_1 \prec_\varepsilon \mathbf{x}_2$），当且仅当：

$$
\boxed{\forall m \in \{1, \dots, M\},\; \frac{f_m(\mathbf{x}_1)}{1 + \varepsilon} \le f_m(\mathbf{x}_2) \quad \text{且} \quad \exists k,\; \frac{f_k(\mathbf{x}_1)}{1 + \varepsilon} < f_k(\mathbf{x}_2)} \tag{19}
$$

> **直觉**：$\mathbf{x}_1$ 不一定要在每个目标上都严格不差于 $\mathbf{x}_2$——允许最多 $(1+\varepsilon)$ 倍的松弛。这在数学上将连续的 Pareto 前沿"离散化"。

**Lemma 6.1（ε-支配蕴含 Pareto 支配的近似）**

若 $\mathbf{x}_1 \prec \mathbf{x}_2$（严格 Pareto 支配），则 $\mathbf{x}_1 \prec_\varepsilon \mathbf{x}_2$ 对所有 $\varepsilon \ge 0$ 均成立。但逆命题不成立：$\mathbf{x}_1 \prec_\varepsilon \mathbf{x}_2$ 不能推出 $\mathbf{x}_1 \prec \mathbf{x}_2$。

*证明*：正向显然（$\frac{f(\mathbf{x}_1)}{1+\varepsilon} \le f(\mathbf{x}_1) \le f(\mathbf{x}_2)$）。反例：取 $M=1, \varepsilon=0.1, f(\mathbf{x}_1)=10, f(\mathbf{x}_2)=10.5$，则 $10/1.1 \approx 9.09 \le 10.5$，故 $\mathbf{x}_1 \prec_\varepsilon \mathbf{x}_2$，但 $f(\mathbf{x}_1) < f(\mathbf{x}_2)$ 即 $\mathbf{x}_1$ 严格 Pareto 支配 $\mathbf{x}_2$（此时两者均为真）。更富信息量的反例：$M=2, \mathbf{f}(\mathbf{x}_1)=(10, 10), \mathbf{f}(\mathbf{x}_2)=(10.5, 9.5), \varepsilon=0.1$。$10/1.1 \approx 9.09 < 10.5$，$10/1.1 \approx 9.09 < 9.5$，满足 ε-支配；但 $f_2(\mathbf{x}_1)=10 > 9.5=f_2(\mathbf{x}_2)$，不满足 Pareto 支配。$\square$

### 6.3 Box 索引机制

ε-支配的核心数据结构是 **Box 索引**，它将目标空间离散化为一个网格。

**定义 6.2（Box 索引向量）**

对于目标向量 $\mathbf{f} \in \mathbb{R}_{>0}^M$，其在 ε-网格中的 Box 索引向量 $\mathbf{B} \in \mathbb{Z}^M$ 为：

$$
\boxed{B_m = \left\lfloor \frac{\log(f_m)}{\log(1 + \varepsilon)} \right\rfloor, \quad m = 1, \dots, M} \tag{20}
$$

- $\lfloor \cdot \rfloor$：向下取整（floor）
- 要求 $f_m > 0$（实践中设为 $\max(f_m, 10^{-12})$）

**Lemma 6.2（Box 索引的基本性质）**

1. **单调性**：若 $f_m^{(1)} \le f_m^{(2)}$，则 $B_m^{(1)} \le B_m^{(2)}$
2. **ε-区间**：Box $\mathbf{B}$ 内所有解在第 $m$ 目标上的值满足 $f_m \in [(1+\varepsilon)^{B_m}, (1+\varepsilon)^{B_m+1})$
3. **支配与 Box 的关系**：若 $\mathbf{f}^{(1)}$ ε-支配 $\mathbf{f}^{(2)}$，则 $\forall m, B_m^{(1)} \le B_m^{(2)}$ 且 $\exists k, B_k^{(1)} < B_k^{(2)}$（即 $\mathbf{B}^{(1)}$ "box-支配" $\mathbf{B}^{(2)}$）

*证明*：

(1) $f_m \mapsto \lfloor \frac{\log f_m}{\log(1+\varepsilon)} \rfloor$ 是单调不减的。

(2) 由 $B_m = \lfloor \frac{\log f_m}{\log(1+\varepsilon)} \rfloor$，有 $\frac{\log f_m}{\log(1+\varepsilon)} \in [B_m, B_m+1)$，即 $\log f_m \in [B_m\log(1+\varepsilon), (B_m+1)\log(1+\varepsilon))$，取指数得证。

(3) 由定义 (19) 和 (1) 直接导出。$\square$

### 6.4 ε-MOPSO 的 Archive 更新

**Algorithm 2: ε-Archive 更新**

```
Input:  A (当前 Archive), new_x, new_f
        ε (ε值), Aₘₐₓ (最大容量)
Output: 更新后的 A

1:  new_box ← box_index(new_f, ε)
2:  for each aⱼ ∈ A do
3:      arch_box ← box_index(aⱼ.f, ε)
4:      if new_box == arch_box then
5:          ▷ 同一 Box：保留目标值之和更小的
6:          if sum(new_f) < sum(aⱼ.f) then
7:              标记 aⱼ 待删除
8:          else
9:              return A  ▷ 新解不优，直接拒绝
10:         end if
11:     else if new_x ε-dominates aⱼ then
12:         标记 aⱼ 待删除
13:     else if aⱼ ε-dominates new_x then
14:         return A  ▷ 新解被支配，直接拒绝
15:     end if
16: end for
17: 删除所有标记的解
18: 将 new_x 加入 A
19: if |A| > Aₘₐₓ then
20:     从最拥挤的 Box 中随机删除一个解
21: end if
```

### 6.5 ε-MOPSO 完整算法

**Algorithm 3: ε-MOPSO**

```
Input:  F: ℝᴰ → ℝᴹ (多目标函数)
        N (种群大小), Kₘₐₓ (最大迭代)
        ε (ε值), wₘₐₓ, wₘᵢₙ, c₁, c₂
        lb, ub
Output: A (最终 Archive)

1:  ▷ 初始化
2:  for i = 1 to N do
3:      xᵢ ~ Uniform(lb, ub)
4:      vᵢ ← 0
5:      pbestᵢ ← xᵢ
6:  end for
7:  A ← InitializeArchive({xᵢ}, ε)   ▷ 用 ε-支配初始化
8:
9:  for k = 1 to Kₘₐₓ do
10:     w ← wₘₐₓ - (wₘₐₓ - wₘᵢₙ) × k/Kₘₐₓ
11:     for i = 1 to N do
12:         ▷ 选择 gbest（Sigma 方法）
13:         gbest ← SelectLeader(A, pbestᵢ)
14:
15:         ▷ 更新速度和位置
16:         r₁, r₂ ~ U(0,1)
17:         vᵢ ← w·vᵢ + c₁·r₁·(pbestᵢ - xᵢ) + c₂·r₂·(gbest - xᵢ)
18:         xᵢ ← xᵢ + vᵢ
19:         边界反射处理
20:
21:         ▷ 评估
22:         f ← F(xᵢ)
23:
24:         ▷ 更新 pbest（ε-支配判断）
25:         if xᵢ ε-dominates pbestᵢ then
26:             pbestᵢ ← xᵢ
27:         else if pbestᵢ ε-dominates xᵢ then
28:             ▷ 保持不变
29:         else
30:             ▷ 互不支配，随机选择
31:             if rand() < 0.5 then pbestᵢ ← xᵢ
32:         end if
33:
34:         ▷ 更新 Archive
35:         A ← ε-ArchiveUpdate(A, xᵢ, f, ε)
36:     end for
37: end for
```

---

## 7. ε-MOPSO 的理论性质

> 📚 深入阅读（对你必读）。以下定理保证了算法行为的可预测性和理论完备性。

### 7.1 Archive 容量的有限性

**Theorem 7.1（ε-Archive 的有界性）**

设目标空间的下界为 $\mathbf{f}_{\min} > \mathbf{0}$，上界为 $\mathbf{f}_{\max}$（由搜索空间的有界性和目标函数的连续性保证）。则 ε-Archive 中解的个数存在上界：

$$
|A_\varepsilon| \le \prod_{m=1}^M \left( \left\lfloor \frac{\log f_{\max,m} - \log f_{\min,m}}{\log(1+\varepsilon)} \right\rfloor + 1 \right) \tag{21}
$$

*证明*：根据 Archive 更新规则，每个 Box 内至多保留一个解。Box 的总数为各维 Box 索引范围的乘积。第 $m$ 维上 $B_m$ 的最小值为 $\lfloor \frac{\log f_{\min,m}}{\log(1+\varepsilon)} \rfloor$，最大值为 $\lfloor \frac{\log f_{\max,m}}{\log(1+\varepsilon)} \rfloor$，故独立 Box 的总数为 $\prod_{m=1}^M (B_m^{\max} - B_m^{\min} + 1)$，即式 $(21)$。$\square$

**Corollary 7.1（ε 对 Archive 大小的控制）**

当 $\varepsilon \to 0$，Archive 上界 $\to \infty$（退化为标准 Pareto Archive）；当 $\varepsilon$ 增大，上界迅速缩小。这是 $\varepsilon$ 作为"精度-效率权衡"参数的数学直觉。

### 7.2 ε-Pareto 前沿的近似精度

**定义 7.1（加性 ε-支配）**

某些文献中使用加性形式：$\mathbf{x}_1$ 加性 ε-支配 $\mathbf{x}_2$ 若 $\forall m, f_m(\mathbf{x}_1) - \varepsilon \le f_m(\mathbf{x}_2)$。本文采用的乘性形式 $(19)$ 等价于目标空间中的**相对误差**保证，对于不同量级的目标值更为公平。

**Theorem 7.2（乘性 ε-支配的近似保证）**

对于 Archive 中任意两个不同 Box 中的解 $\mathbf{a}_1, \mathbf{a}_2$，它们之间不存在 ε-支配关系。且 Archive 在目标空间中的"覆盖间隙"受 $\varepsilon$ 控制：

$$
\forall \mathbf{a} \in \mathcal{A}, \forall m, \quad \frac{f_m^{\text{true}}}{f_m(\mathbf{a})} \in \left[\frac{1}{1+\varepsilon}, 1+\varepsilon\right]
$$

即 Archive 中的每个解在真实 Pareto 前沿上，每个目标值的相对偏差不超过 $\varepsilon$（在乘性意义下）。

### 7.3 时间复杂度和空间复杂度

**Theorem 7.3（计算复杂度）**

单次迭代中 ε-MOPSO 的计算复杂度为：

- speed/position 更新：$O(ND)$
- pbest 更新：$O(N M)$
- Archive 更新：$O(N \cdot |A| \cdot M)$（每个粒子与每个存档解比较）

总复杂度：$O(K_{\max} \cdot (ND + N \cdot |A| \cdot M))$。由于 $|A|$ 受 $\varepsilon$ 限制（Theorem 7.1），算法是 **多项式时间** 的。

---

## 8. 与 RST 控制器的衔接

### 8.1 RST 优化问题的 PSO 编码

在本工程中，PSO 粒子的位置直接编码 RST 控制器设计中期望灵敏度函数的可调参数：

| PSO 概念 | RST 映射 |
|:---|:---|
| 粒子位置 $\mathbf{x}_i$ | 决策变量 $\boldsymbol{\theta} = (x_1, \dots, x_{n_X}, y_1, \dots, y_{n_Y})$ |
| 搜索空间 $\mathcal{S}$ | $x_j, y_j \in [-0.95, 0.95]$（根在单位圆内） |
| 目标函数 $\mathbf{F}$ | Zames-Francis 积分匹配误差 + Bezout 残差 |
| 约束 $g_j$ | 频域模值裕度 / 延迟裕度上下界 |

### 8.2 约束处理

约束通过**外部静态惩罚**融入目标函数。参见 [`penalty_function.m`](penalty_function.m)：

$$
\varphi_m(\boldsymbol{\theta}) = f_m(\boldsymbol{\theta}) + \Lambda \cdot \sum_{j} \overline{v_j^\kappa(\boldsymbol{\theta})}
$$

其中 $v_j(\boldsymbol{\theta}) = \max(0, g_j(\boldsymbol{\theta}))$ 只对违例部分惩罚。

### 8.3 完整求解链路

```
被控对象 G(z) → 提取不稳定零/极点 → 设 P_D, H_R, H_S
    → 构建优化命题 (RST_objective.m + RST_constraints.m)
    → ε-MOPSO 搜索 (eMOPSO_core.m)
    → 从 Pareto 前沿选拐点解
    → Bezout 反解 R,S (postprocess_RST.m)
    → 频域/时域验证
```

参见 [`RST_eMOPSO_spec.md`](RST_eMOPSO_spec.md) 和 [`Solution.md`](Solution.md) 获取 RST 端的完整数学推导。

---

## 9. 参考文献

[1] Kennedy, J. and Eberhart, R., "Particle Swarm Optimization", *Proc. IEEE Int. Conf. on Neural Networks*, pp. 1942–1948, 1995.

[2] Clerc, M. and Kennedy, J., "The Particle Swarm—Explosion, Stability, and Convergence in a Multidimensional Complex Space", *IEEE Trans. Evolutionary Computation*, 6(1), pp. 58–73, 2002.

[3] Trelea, I.C., "The Particle Swarm Optimization Algorithm: Convergence Analysis and Parameter Selection", *Information Processing Letters*, 85(6), pp. 317–325, 2003.

[4] Mostaghim, S. and Teich, J., "Strategies for Finding Good Local Guides in Multi-objective Particle Swarm Optimization (MOPSO)", *Proc. IEEE Swarm Intelligence Symposium*, pp. 26–33, 2003.

[5] Laumanns, M., Thiele, L., Deb, K., and Zitzler, E., "Combining Convergence and Diversity in Evolutionary Multiobjective Optimization", *Evolutionary Computation*, 10(3), pp. 263–282, 2002.

[6] Sierra, M.R. and Coello Coello, C.A., "Improving PSO-Based Multi-objective Optimization Using Crowding, Mutation and ε-Dominance", *EMO 2005*, LNCS 3410, pp. 505–519, 2005.

[7] Coello Coello, C.A. and Lechuga, M.S., "MOPSO: A Proposal for Multiple Objective Particle Swarm Optimization", *Proc. CEC 2002*, pp. 1051–1056.

[8] van den Bergh, F. and Engelbrecht, A.P., "A Study of Particle Swarm Optimization Particle Trajectories", *Information Sciences*, 176(8), pp. 937–971, 2006.

[9] Parsopoulos, K.E. and Vrahatis, M.N., *Particle Swarm Optimization and Intelligence: Advances and Applications*, IGI Global, 2010.

---

> 📖 **阅读路径建议**：本文档 → [`MOEA_algorithms.md`](MOEA_algorithms.md)（四种算法横向对比）→ [`RST_eMOPSO_spec.md`](RST_eMOPSO_spec.md)（RST 端数学建模）→ [`Solution.md`](Solution.md)（案例实战）→ [`README.md`](README.md)（代码使用手册）
