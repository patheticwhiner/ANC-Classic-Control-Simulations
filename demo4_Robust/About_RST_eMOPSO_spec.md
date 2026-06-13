# RST 控制器整定与 ε-MOPSO 算法实现规范

本文档为基于 ε-MOPSO 的 RST 控制系统参数整定与优化算法的数学建模与程序实现指南。请严格依据以下数学定义与伪代码编写仿真程序。

---

## 一、 控制系统建模与决策变量定义 (System & Variables)

### 1.1 受控对象与固定参数
受控对象 (Plant) 离散传递函数：
$$
G(z^{-1}) = z^{-d} \frac{B(z^{-1})}{A(z^{-1})}
$$
- 必须已知 $A(z^{-1})$ 和 $B(z^{-1})$ 的系数，以及纯延迟步数 $d$。
- 系统开环不稳定极点记为 $\alpha_i$。
- 系统开环不稳定零点记为 $\beta_i = r_i e^{j\phi_i}$ （其中模长 $r_i > 1$）。

用户需预先指定的基准多项式：
- 闭环主导极点多项式 $P_D(z^{-1})$。
- RST控制器的固定部分：$H_R(z^{-1}) = 1 + z^{-1}$， $H_S(z^{-1}) = 1 - z^{-1}$。

### 1.2 决策变量与搜索空间
优化算法的搜索空间粒子（决策变量）为向量 $\boldsymbol{\theta}$，它代表待求多项式 $X$ 和 $Y$ 的实数零点：
$$
\boldsymbol{\theta} = (x_1, x_2, \dots, x_{n_X}, y_1, y_2, \dots, y_{n_Y})^T
$$
程序中对应多项式构建方式：
- $$
  X(z^{-1}, \boldsymbol{\theta}) = \prod_{i=1}^{n_X} (1 - x_i z^{-1})
  $$
- $$
  Y(z^{-1}, \boldsymbol{\theta}) = \prod_{i=1}^{n_Y} (1 - y_i z^{-1})
  $$

**待优化的期望灵敏度函数**：
$$
S_d(z^{-1}, \boldsymbol{\theta}) = \frac{A(z^{-1})}{P_D(z^{-1})} \frac{X(z^{-1}, \boldsymbol{\theta})}{Y(z^{-1}, \boldsymbol{\theta})}
$$

---

## 二、 多目标优化问题建模 (MOP Formulation)

### 2.1 目标函数 (Objective Functions)
依据 Zames-Francis 积分等式构建目标函数。每个开环不稳定零点 $\beta_i$ 对应一个目标 $f_i$，纯延迟 $d$ 对应目标 $f_d$。积分域为 $\omega \in [-\pi, \pi]$，对应 $z = e^{j\omega}$。

**目标 1：针对不稳定零点 $\beta_i = r_i e^{j\phi_i}$**
$$
f_i(\boldsymbol{\theta}) = \left| \int_{-\pi}^{\pi} \log|S_d(e^{j\omega}, \boldsymbol{\theta})| \frac{r_i^2 - 1}{1 - 2r_i \cos(\omega - \phi_i) + r_i^2} d\omega - \pi \log|B_a^{-1}(\beta_i)| \right|
$$
*(其中
$$
B_a(z) = \prod \frac{z - \alpha_k}{1 - \bar{\alpha}_k z}
$$
为系统不稳定极点构成的 Blaschke 乘积)*

**目标 2：针对纯时间延迟 (相当于无穷远零点)**
$$
f_d(\boldsymbol{\theta}) = \left| \int_{-\pi}^{\pi} \log|S_d(e^{j\omega}, \boldsymbol{\theta})| d\omega - \pi \sum_k \log|\alpha_k| \right|
$$

### 2.2 不等式约束 (Constraints)
系统实际的输出灵敏度函数为 $S_{yp}(z^{-1}, \boldsymbol{\theta})$（在本文寻优阶段，可直接令 $S_{yp} \approx S_d$ 进行频率模板验证）。对所有频率 $\omega \in [0, \pi]$ 必须满足以下绝对值约束（$\leq 0$）：

- **下界约束 (Modulus Margin):**
  $$
  g_1(e^{j\omega}, \boldsymbol{\theta}) = 1 - |1 - e^{-j\omega}|^{-1} - |S_d(e^{j\omega}, \boldsymbol{\theta})| \leq 0
  $$
- **上界约束 (Delay Margin):**
  $$
  g_2(e^{j\omega}, \boldsymbol{\theta}) = |S_d(e^{j\omega}, \boldsymbol{\theta})| - 1 - |1 - e^{-j\omega}|^{-1} \leq 0
  $$

---

## 三、 PSO 基本算法实现 (Standard PSO)

> 📎 PSO 完整数学推导（含 Jury 收敛性证明）见 **`About_PSO_Foundations.md`**。本节仅给规范摘要。

设种群大小为 $n_{pop}$，最大迭代次数为 $k_{max}$。第 $i$ 个粒子的位置为 $\mathbf{x}^i_k$（即参数 $\boldsymbol{\theta}$），速度为 $\mathbf{v}^i_k$。

**运动学更新方程：**
$$
\mathbf{x}^i_{k+1} = \mathbf{x}^i_k + \mathbf{v}^i_{k+1}
$$
$$
\mathbf{v}^i_{k+1} = w_k \mathbf{v}^i_k + c_1 r_{1,k}^i (\mathbf{pbest}^i_k - \mathbf{x}^i_k) + c_2 r_{2,k}^i (\mathbf{gbest}_k - \mathbf{x}^i_k)
$$

- $r_{1,k}^i, r_{2,k}^i \sim U(0,1)$ 为随机数。
- $c_1 = 1.5, c_2 = 1.7$。
- **自适应惯性权重 $w_k$：**
  $$
  w_{k+1} = w_{max} - (w_{max} - w_{min}) \frac{k}{k_{max}}
  $$
  (其中 $w_{max} = 0.9, w_{min} = 0.4$)

---

## 四、 约束处理机制 (Static Penalty Function)

使用外部静态惩罚技术处理 $g_1, g_2$ 约束。将违反约束的惩罚项 $\delta(\mathbf{x})$ 加到目标函数上：
$$
\varphi_m(\mathbf{x}) = f_m(\mathbf{x}) + \delta(\mathbf{x}) \quad (m \text{ 为目标维度索引})
$$

**惩罚项计算：**
$$
\delta(\mathbf{x}) = \sum_{j=1}^{q} \Lambda_j v_j^\kappa(\mathbf{x})
$$
其中：
- $q$ 为约束总数。
- $\Lambda_j$ 为大常数惩罚系数（如 $10^4$）。
- $\kappa = 2$。
- 约束惩罚函数
  $$
  v_j(\mathbf{x}) = \max(0, g_j(\mathbf{x}))
  $$
  *(如果在某些频点违反，求所有频点违反量之和或最大值)*

---

## 五、 ε-MOPSO 核心机制 (ε-Domination & Archive)

> 📎 ε-支配的完整定义与性质证明见 **`About_PSO_Foundations.md`** §多目标优化。

### 5.1 ε-支配定义 (ε-Domination Rule)
对于两个解（粒子） $\boldsymbol{\theta}_1$ 和 $\boldsymbol{\theta}_2$，存在 $m_f$ 个目标。我们称 **$\boldsymbol{\theta}_1$ $\varepsilon$-支配 $\boldsymbol{\theta}_2$** (记作 $\boldsymbol{\theta}_1 \prec_\varepsilon \boldsymbol{\theta}_2$)，当且仅当：
1. $$
   \forall m \in \{1, 2, \dots, m_f\}, \quad \frac{f_m(\boldsymbol{\theta}_1)}{1 + \varepsilon} \leq f_m(\boldsymbol{\theta}_2)
   $$
2. $$
   \exists m \in \{1, 2, \dots, m_f\}, \quad \frac{f_m(\boldsymbol{\theta}_1)}{1 + \varepsilon} < f_m(\boldsymbol{\theta}_2)
   $$

### 5.2 档案盒划分 (Algorithm 2: Box Indexing)
将目标空间按 $\varepsilon$ 对数划分网格。对于任意输入的目标向量 $\mathbf{f} = (f_1, f_2, \dots, f_{m_f})$，计算其所在的 Box 索引向量 $\mathcal{B}$：

```text
For all m = 1, 2, ..., m_f do:
    B_m = floor( log(f_m) / log(1 + ε) )
End For
Output B = (B_1, ..., B_{m_f})  // B 是该解的网格坐标
```

# 附录 A：多目标优化测试函数 (Test functions for multiobjective optimization)

## 1. 测试函数 F1 (Test Function F1)

$$
\begin{cases}
f_1(\boldsymbol{x}) = x_1 \\
f_2(\boldsymbol{x}) = g(\boldsymbol{x})h(\boldsymbol{x}) \\
n = 2
\end{cases} \quad (A.1)
$$

**变量约束与辅助函数：**
$$
\boldsymbol{x} = (x_1, x_2, \dots, x_n)
$$
$$
0 \leq x_1 \leq 1, \quad -30 \leq x_2 \leq 30
$$
$$
g(\boldsymbol{x}) = 11 + x_2^2 - 10\cos(2\pi x_2)
$$
$$
h(\boldsymbol{x}) = 
\begin{cases} 
1 - \sqrt{\frac{f_1(\boldsymbol{x})}{g(\boldsymbol{x})}}, & \text{if } f_1(\boldsymbol{x}) \leq g(\boldsymbol{x}) \\ 
0, & \text{otherwise} 
\end{cases}
$$

---

## 2. 测试函数 F2 (Test Function F2)

$$
\begin{cases}
f_1(\boldsymbol{x}) = x_1 \\
f_2(\boldsymbol{x}) = g(x_2, x_3, \dots, x_n)h(f_1, g) + 1 \\
n = 30
\end{cases} \quad (A.2)
$$

**变量约束与辅助函数：**
$$
0 \leq \boldsymbol{x} = (x_1, x_2, \dots, x_n) \leq 1
$$
$$
g(\boldsymbol{x}) = 1 + \frac{9}{n - 1} \sum_{i=2}^{n} x_i
$$
$$
h(f_1, g) = 1 - \sqrt{\frac{f_1}{g}} - \left(\frac{f_1}{g}\right)\sin(10\pi f_1)
$$

---

## 3. 测试函数 F3 (Test Function F3)

$$
\begin{cases}
f_1(\boldsymbol{x}) = (1 + x_3)\cos\left(\frac{\pi x_1}{2}\right)\cos\left(\frac{\pi x_2}{2}\right) \\
f_2(\boldsymbol{x}) = (1 + x_3)\cos\left(\frac{\pi x_1}{2}\right)\sin\left(\frac{\pi x_2}{2}\right) \\
f_3(\boldsymbol{x}) = (1 + x_3)\sin\left(\frac{\pi x_1}{2}\right) \\
n = 3
\end{cases} \quad (A.3)
$$

**变量约束：**
$$
0 \leq \boldsymbol{x} = (x_1, x_2, x_3) \leq 1
$$

---

## 4. 测试函数 F4 (Test Function F4)

$$
\begin{cases}
f_1(\boldsymbol{x}) = x_1 \\
f_2(\boldsymbol{x}) = x_2 \\
f_3(\boldsymbol{x}) = 3.5 - \sum_{i=1}^{n} x_i \sin(n\pi x_i) \\
n = 3
\end{cases} \quad (A.4)
$$

**变量约束：**
$$
0 \leq \boldsymbol{x} = (x_1, x_2, x_3) \leq 1
$$
*(注：根据上下文与标准多目标测试函数文献还原，对原 PDF OCR 乱码处的连乘/累加表示进行了规范化纠正整理。)*


---

# 附录 B：Zames-Francis 框架的适用性边界——以 ARMAX 声学管道模型为例

> 🎯 **读者定位**：本附录写给未来接手此项目的同组研究者。你大概率已经实操过 ARMAX 系统辨识（理解 AIC/BIC 定阶、模型验证），但对 Zames-Francis 积分约束仅知其结论而不熟悉推导细节。因此：
> - **§B.1–B.2** 从你熟悉的辨识输出（零极点图）出发，论证问题根源——**这是你接手后最先需要理解的内容**。
> - **§B.3** 给出 Poisson 核退化的数学直觉——不要求你独立推导，但需要看懂"为什么 |β|→1 时目标函数会病态"。
> - **§B.4** 给出四种解决方案的量化对比——**这是你执行下一步实验的菜单**。
> - 如果你对 §B.3 的 Poisson 核推导感到陌生，回顾本文档 **§2.1** 的目标函数定义和 **§4.2** 的 Poisson 核公式 `(r²−1)/(1−2r cos(ω−φ)+r²)`——注意当 r→1 时分子趋近于 0 而分母趋近于 `2−2cos(ω−φ)`，核函数退化为 δ-脉冲。

---

## B.1 问题：ARMAX(30,30,30,22) 模型真的有 9 个"非最小相位零点"吗？

> 📎 完整调试历程（四次 Bezout 迭代、nX 提升轨迹、运行时间数据）见 **`About_ARMAX_Debugging.md`**。本附录聚焦根因分析和解决方案对比。

### B.1.1 实测数据

对 `armax_30303022` 模型（48 kHz 采样）的 B 多项式求根，得到 29 个零点。其中 9 个满足 `|z| > 1`（按定义属于"不稳定零点"）：

| 编号 | 零点 z | 模值 |z| | 分类 |
|:---|:---|:---|:---|
| β₁, β₂ | −0.40 ± 1.03j | **1.102** | 弱不稳定 |
| β₃, β₄ | −0.73 ± 0.70j | **1.009** | ⚠️ 几乎在单位圆上 |
| β₅, β₆ | 0.84 ± 0.87j | **1.208** | 弱不稳定 |
| β₇, β₈ | 1.22 ± 0.24j | **1.240** | 弱不稳定 |
| β₉ | 1.005 | **1.005** | ⚠️ 几乎在单位圆上 |

同时，A 多项式 30 个极点**全部在单位圆内**（开环稳定）。

### B.1.2 物理矛盾

声学管道是被动、因果的物理系统。其连续时间传递函数的所有零点均位于左半 s-平面（最小相位）。在 48 kHz 采样率下离散化后，真正的零点应当全部落在单位圆内。

**9 个 |z| > 1 的零点中，有 5 个模值 ≤ 1.01**——对于双精度浮点，|z| = 1.009 与 |z| = 0.991 之间的差异（~2%）远小于辨识方差。这意味着：**如果换一组激励信号重新辨识，这 5 个零点有一半概率会落在单位圆内**。

> 🎯 **对你而言**：这不难理解。你做辨识时见过类似现象——高阶 MA 部分的零点对训练数据的微小扰动极度敏感。下面的 §B.2 解释为什么这在 48 kHz + ARMAX(30,30,30) 配置下是**结构性必然**。

---

## B.2 根因：ARMAX(30,30,30,22) 伪非最小相位零点的三因素解释

### B.2.1 过参数化：MA 部分有远超物理需要的自由度

ARMAX 模型结构：

$$
A(z^{-1}) y(k) = z^{-d} B(z^{-1}) u(k) + C(z^{-1}) e(k)
$$

在 48 kHz 下，声学管道的物理动力学集中在 **0–500 Hz** 的窄带内（对应归一化频率 ω ∈ [0, 0.02π]）。这意味着 Nyquist 频率（24 kHz）以上的频段几乎全是噪声和抗混叠滤波器的残余。

ARMAX(30,30,30) 的 B 多项式有 30 个参数（实际 29 个零点）。物理上，拟合 0–500 Hz 的动力学最多需要 4–8 个零点。**剩余的 21–25 个 MA 自由度**被迫去"解释"：
- 高频测量噪声的频谱特征
- 抗混叠滤波器的相位畸变
- 扬声器/传声器的非线性残余

这些自由度缺乏物理约束，零点的模值在辨识过程中"漂移"——部分恰好落在单位圆外。

> 🎯 **对你而言的类比**：这和你用 AIC 选阶时看到的陡降点原理一致——AIC 在"真实动力学阶数"处出现拐点。如果你对 48 kHz 的声学管道数据跑一次 AIC，你会发现最优阶数远小于 30。这意味着 ARMAX(30,30,30) 是显著**过参数化**的——而 MA 过参数化的必然产物就是近单位圆的伪零点。

### B.2.2 高采样率放大效应

离散化零点的模值 |z| 与连续时间零点 s 的关系为：

$$
|z| = |e^{s T_s}|
$$

在 48 kHz（T_s = 20.8 μs）下：
- 一个连续时间零点 s = −100 rad/s → |z| = e^{−0.0021} ≈ **0.998**
- 所有真实零点都被压缩到 |z| ∈ (0.995, 1.0) 的极窄环形区域内

辨识算法需要在 29 个零点中区分模值为 0.995（真实）和 1.005（artifact）的零点——**信噪比不允许这种精度**。

### B.2.3 Poisson 核的退化：当 |β| → 1 时目标函数奇异

Zames-Francis 框架中，每个不稳定零点 β_i = r_i e^{jφ_i} 对应的 Poisson 核为：

$$
K_i(\omega) = \frac{r_i^2 - 1}{1 - 2 r_i \cos(\omega - \phi_i) + r_i^2}
$$

**关键观察**：当 r → 1⁺ 时，分子 r² − 1 → 0⁺，而分母（在 ω ≈ φ_i 处）→ (1 − r)² → 0⁺。

Poisson 核的极限行为：

$$
\lim_{r \to 1^+} K(\omega) = 2\pi \cdot \delta(\omega - \phi_i)
$$

即：

| r 值 | Poisson 核行为 | 对目标函数的影响 |
|:---|:---|:---|
| r = 1.5 | 光滑钟形，半宽 ~0.5 rad | 积分鲁棒，数值稳定 |
| r = 1.1 | 尖锐峰，半宽 ~0.05 rad | 需要密集频率采样 |
| r = 1.01 | 近 δ-脉冲，半宽 ~0.005 rad | **500 点频率离散严重欠采样** |
| r = 1.005 | 几乎就是 δ(ω−φ) | 积分值完全取决于恰好落在峰值处的那 1–2 个频点 |

对于 β₉（r = 1.005），即使用 500 个频点，Poisson 核的有效宽度只覆盖 2–3 个频点——**梯形积分的误差可以轻松达到 100%**。

> 🎯 **对你而言**：你不必独立推出这个极限。你只需要知道：当你看到 `RST_objective.m` 里某个 f_β 的目标值是 `5.23` 而另一个是 `0.87` 时，差异可能**不是优化意义上的 Pareto trade-off，而是纯数值噪声**——两个目标都在被 Poisson 核的欠采样误差主导。

---

## B.3 对 ε-MOPSO 优化框架的量化影响

### B.3.1 目标函数空间的退化

在 11 个目标中（9 个 f_β + f_delay + f_bezout）：

| 目标 | 对应零点 | Poisson 核状态 | 数值可信度 |
|:---|:---|:---|:---|
| f_β₁–f_β₄ | |z| = 1.10–1.24 | 正常 | ✅ 可信 |
| f_β₅–f_β₉ | |z| = 1.005–1.01 | **近 δ-脉冲** | ❌ 被数值噪声主导 |
| f_delay | — | 无核 | 正常 | ✅ 可信 |
| f_bezout | — | — | 正常 | ✅ 可信 |

**结论**：11 维目标空间中，至少 5 个维度（45%）的有效信息含量接近于零。

### B.3.2 Pareto 支配在高噪声维度上的失效

ε-MOPSO 的 Box 索引基于 `log(f_m)`。当 f_m 的值被数值噪声主导时：
- 两个本应等价的解，因噪声差异被分配到**不同的 Box**
- ε-支配判断在噪声维度上随机地接受/拒绝候选解
- Archive 的多样性被"假多样性"填充——保留了大量在噪声维度上不同但在真正重要的维度上等价的解

**这解释了 `About_ARMAX_Debugging.md` 中观察到的现象**：nX=50 时，Pareto 前沿上找到了"系数小但鲁棒性弱"的解——这个解在噪声维度上表现好（Poisson 积分碰巧落在低值），但在真实的鲁棒性指标（GM）上退化。

### B.3.3 搜索空间效率损失

62 维决策变量空间中，ε-MOPSO 需要花费大量迭代预算在 5 个无信息目标维度上"区分"粒子——这解释了为什么 400 次迭代后 Archive 仍未饱和（收敛慢）。

---

## B.4 解决方案：四种路径的量化对比

### 方案 A：硬阈值筛选 β

直接丢弃 |β| < 1.05 的不稳定零点，只保留物理上有意义的。

| 参数 | 变更前 | 变更后 |
|:---|:---|:---|
| 不稳定零点数 | 9 | **4**（仅 |z| > 1.05） |
| 目标数 | 11 | **6** |
| 搜索难度 | 62 维 / 11 目标 | 62 维 / **6 目标** |
| Poisson 核鲁棒性 | 5/9 退化 | 全部正常 |

**代价**：丢弃了模型的 5 个零点信息——如果这些零点恰好是真实非最小相位零点，会导致灵敏度整形约束不足。

**实现**：在 `run_eMOPSO.m` 第 89 行过滤 `abs(beta) > 1.05`。

### 方案 B：模型降阶后重新提取 β

用 `balred` 或 `reduce` 将 ARMAX(30,30,30) 降阶到物理上有意义的阶数（建议 6–8 阶），再从降阶模型中提取 β。

| 参数 | 变更前 | 变更后 |
|:---|:---|:---|
| 模型阶数 | 30 | **6–8** |
| 零点总数 | 29 | **5–7** |
| 期望不稳定零点 | ？ | **≤2** |
| Bezout 矩阵维度 | 83×83 | **~20×20** |
| R/S 系数爆炸 | 严重 | **基本消除** |

**代价**：模型保真度下降——降阶模型的高频响应与原始 ARMAX 模型存在偏差。

**实现**：
```matlab
sys_full = tf(B_delayed, A, Ts, 'Variable', 'z^-1');
sys_red = balred(sys_full, 6);
[B_red, A_red] = tfdata(sys_red, 'v');
```

### 方案 C：频域加权替代 Poisson 积分

放弃 Zames-Francis 逐点 β 匹配，改为在频域上直接约束灵敏度函数的**包络形状**：

$$
f_{\text{env}}(\boldsymbol{\theta}) = \int_{0}^{\pi} \left| |S_d(e^{j\omega})| - W^{-1}(\omega) \right|_+ d\omega
$$

其中 W⁻¹(ω) 是用户指定的频域权重上界（而非从 β 导出的 Poisson 积分约束）。

**代价**：失去了 Zames-Francis 的严密理论支撑，频域权重的选择变为主观工程判断。

### 方案 D：H∞ 加权灵敏度整形（框架级变更）

延续 `About_Solution.md` 中 Q-参数化路线，直接将问题建模为：

$$
\min_{Q} \| W_S S_{yp} \|_\infty + \| W_T T_{yp} \|_\infty
$$

完全脱离 Zames-Francis 积分框架和 Bezout 反解。

**代价**：实现复杂度最高，需要 Control System Toolbox 的 H∞ 求解器。

### 方案对比总结

| 维度 | 方案 A<br>硬阈值 | 方案 B<br>模型降阶 | 方案 C<br>频域包络 | 方案 D<br>H∞ 整形 |
|:---|:---|:---|:---|:---|
| 改动量 | 1 行代码 | ~30 行 | RST_objective 重写 | 框架重写 |
| ZF 理论完整性 | 部分保留 | 部分保留 | 放弃 | 放弃 |
| 运行时间改善 | 微小 | **显著**（8×） | 中等 | 未知 |
| R/S 系数爆炸 | 不解决 | **解决** | 不解决 | 解决 |
| 发表论证压力 | 低 | 低 | 中 | 高 |
| **推荐优先级** | ⭐⭐ | ⭐⭐⭐ | ⭐ | ⭐⭐ |

> 🎯 **对你而言的下一步**：方案 B（模型降阶）是最优起点——它同时解决零点 artifact、运行时间、和 R/S 系数爆炸三个问题，且代码改动最小。降阶后如果 β 数量仍然偏高（>3），再叠加方案 A 的硬阈值。
>
> 如果要执行方案 B，你需要关注的代码入口是 `run_eMOPSO.m:75` 的 `load_armax_full()` 函数——在加载后插入 `balred` 即可。降阶后的系统信息（B, A, d, Ts）无缝替换现有 `sys_info`，不需要修改任何优化层代码。

---

## B.5 终端衔接：回到你的代码

| 本文档概念 | 你代码中的映射 | 受影响的文件 |
|:---|:---|:---|
| 不稳定零点 β | `sys_info.beta` — 在 `run_eMOPSO.m:89` 由 `abs(sys_zeros) > 1` 筛选 | `run_eMOPSO.m:88-99` |
| Poisson 核积分 | `RST_objective.m:113-124` 中 `trapz(omega, log_S_d_mag .* K)` | `RST_objective.m:103-125` |
| Blaschke 乘积 | `RST_objective.m:107` 中 `blaschke_product(beta_i, alpha)` — 对于稳定系统（无 α）退化为 1 | `RST_objective.m:221-237` |
| 目标函数缩放 ×1000 | `wrapped_objective()` 中的 `f = RST_objective(theta, sys_info) * 1000` | `run_eMOPSO.m:439` |
| 近单位圆 β 的 Poisson 核欠采样 | `sys_info.nFreq = 500` — 对于 r=1.005 的 β，500 点严重不足 | `run_eMOPSO.m:188` |
| 模型降阶入口 | `load_armax_full()` — 方案 B 的修改点 | `run_eMOPSO.m:275-309` |

**读完本附录后的建议行动顺序**：
1. 先跑一次 `run_eMOPSO('RST_toy')` 确认框架在低阶教学模型上正常
2. 实施方案 B：在 `load_armax_full()` 末尾加 `balred`，降阶到 8 阶
3. 观察降阶后模型的不稳定零点数量（预期 ≤ 2 个）
4. 如果仍有个别 |β| ≈ 1.01 的零点，叠加方案 A 的硬阈值
5. 重新运行 `run_eMOPSO('RST_armax')`，对比运行时间和控制器性能
