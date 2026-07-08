# Marino & Tomei (2016) — 完整数学推导：未知频率多正弦干扰的自适应输出反馈补偿

> **论文出处**: Marino, R. & Tomei, P. (2016). "Adaptive disturbance rejection for unknown stable linear systems." *Transactions of the Institute of Measurement and Control*, 38(6): 640–647.
>
> **本文档定位**: 应试级完整推导，覆盖问题制定→补偿器设计→稳定性求解→结论的全过程。
> 对应仿真脚本: `MarinoTomei_2016_adaptive_freq_est.m`

---

## 第一部分：问题制定 (Problem Formulation)

### 1. 系统模型

考虑单输入单输出 (SISO) 线性时不变系统：

$$y(s) = P(s) [u(s) + d(s)] \tag{1}$$

其中：
- $y, u, d \in \mathbb{R}$ 分别是输出、输入和干扰
- **假设 1 (稳定性)**：$P(s)$ 是适当的 (proper) 且稳定的（所有极点均在左半平面）
- **假设 2 (符号已知)**：静态增益 $P(0)$ 的符号已知；在干扰频率处，$\text{Re}[P(j\omega_i)] \neq 0$ 时其符号已知（Case A），或 $\text{Re}[P(j\omega_i)] = 0$ 时 $\text{Im}[P(j\omega_i)]$ 的符号已知（Case B）

### 2. 干扰模型（外系统 Exosystem）

干扰 $d(t)$ 是带偏置的多正弦信号：

$$d(t) = d_0 + \sum_{i=1}^q A_i \sin(\omega_i t + \phi_i) \tag{2}$$

其中幅度 $A_i > 0$、相位 $\phi_i$ 和频率 $\omega_i > 0$ **均未知**。频率互不相等 ($\omega_i \neq \omega_j$ for $i \neq j$)。

其状态空间描述（外系统）为：

$$\begin{cases} \dot{w}_0 = 0 \\ \dot{w}_{2i-1} = w_{2i} \\ \dot{w}_{2i} = -\omega_i^2 w_{2i-1} \end{cases} \quad,\quad d = w_0 + \sum_{i=1}^q w_{2i-1} \tag{3}$$

外系统矩阵为分块对角形式：

$$\dot{w} = S w, \quad S = \text{diag}\left(0, \begin{bmatrix} 0 & 1 \\ -\omega_1^2 & 0 \end{bmatrix}, \ldots, \begin{bmatrix} 0 & 1 \\ -\omega_q^2 & 0 \end{bmatrix}\right) \tag{4}$$

---

## 第二部分：补偿器设计 (Compensator Design)

基于**内模原理 (Internal Model Principle)**，补偿器必须包含干扰发生器的结构。由于频率 $\omega_i$ 未知，引入自适应律更新频率平方的估计值 $\hat{\theta}_i \approx \omega_i^2$。

### 1. 控制律与观测器（Case A：$\text{Re}[P(j\omega_i)] \neq 0$）

构造如下自适应补偿器：

$$\begin{aligned} \dot{\hat{w}}_0 &= k \cdot \text{sign}[P(0)] \cdot y \\ \dot{\hat{w}}_{2i-1} &= \hat{w}_{2i} + k \cdot \text{sign}\{\text{Re}[P(j\omega_i)]\} \cdot y \\ \dot{\hat{w}}_{2i} &= -\hat{\theta}_i \hat{w}_{2i-1} \\ u &= -\hat{w}_0 - \sum_{i=1}^q \hat{w}_{2i-1} \end{aligned} \tag{5}$$

其中 $k > 0$ 是观测器增益，$i = 1, 2, \ldots, q$。

### 2. 参数自适应律

利用输出 $y$ 与估计状态 $\hat{w}_{2i}$ 的相关性更新频率估计：

$$\dot{\hat{\theta}}_i = \varepsilon \cdot \text{sign}\{\text{Re}[P(j\omega_i)]\} \cdot \hat{w}_{2i} \cdot y \tag{6}$$

其中 $\varepsilon > 0$ 是足够小的自适应步长。

**Case B 备选**（当 $\text{Re}[P(j\omega_i)] = 0$ 时）：将 $\hat{w}_{2i-1}$ 和 $\hat{w}_{2i}$ 的角色互换，使用 $\text{sign}\{\text{Im}[P(j\omega_i)]\}$ 代替。本文后续分析仅针对 Case A。

### 3. 仿真实现注解

在 `MarinoTomei_2016_adaptive_freq_est.m` 中，被控对象取为 $P(s) = 1/(s+1)^2$。此时：

$$\begin{aligned} \text{sign}[P(0)] &= \text{sign}[1] = +1 \\ \text{sign}\{\text{Re}[P(j0.5)]\} &= \text{sign}[+0.48] = +1 \\ \text{sign}\{\text{Re}[P(j2.0)]\} &= \text{sign}[-0.12] = -1 \end{aligned}$$

代码中的交替符号（`dw1 = w2 + k*y`, `dw3 = w4 - k*y`, `dtheta1 = +ε*w2*y`, `dtheta2 = -ε*w4*y`）恰好是 sign 函数的实例化结果。

---

## 第三部分：稳定性求解与数学推导 (Mathematical Derivation)

求解的核心在于证明误差系统在平衡点处的**局部指数稳定性**。

### 1. 定义误差变量

定义状态误差和参数误差：

$$\tilde{w}_i = w_i - \hat{w}_i, \quad \tilde{\theta}_i = \omega_i^2 - \hat{\theta}_i \tag{7}$$

由式 (3) 和 (5) 推导误差动力学方程。对第 $i$ 个频率分量：

$$\begin{aligned} \dot{\tilde{w}}_{2i} &= \dot{w}_{2i} - \dot{\hat{w}}_{2i} \\ &= -\omega_i^2 w_{2i-1} + \hat{\theta}_i \hat{w}_{2i-1} \\ &= -\omega_i^2 w_{2i-1} + \hat{\theta}_i (w_{2i-1} - \tilde{w}_{2i-1}) \\ &= -\hat{\theta}_i \tilde{w}_{2i-1} - \tilde{\theta}_i w_{2i-1} \end{aligned} \tag{8}$$

式 (8) 揭示了误差动力学的双线性结构：第一项 $-\hat{\theta}_i \tilde{w}_{2i-1}$ 是标称内模动力学，第二项 $-\tilde{\theta}_i w_{2i-1}$ 是参数误差与扰动状态的乘积——这是自适应系统中典型的**非线性参数化**特征。

### 2. 闭环传递函数分析

根据控制律 (5)，闭环系统的输出 $y(s)$ 可表示为：

$$y(s) = P(s) \left[u(s) + d(s)\right] = P(s) \left[-\hat{w}_0(s) - \sum_{i=1}^q \hat{w}_{2i-1}(s) + d(s)\right] \tag{9}$$

将观测器动力学 (5) 变换到频域。对 $\hat{w}_{2i-1}, \hat{w}_{2i}$ 子系统：

$$\begin{aligned} s\hat{w}_{2i-1} &= \hat{w}_{2i} + k \cdot s_{\omega i} \cdot y \\ s\hat{w}_{2i} &= -\hat{\theta}_i \hat{w}_{2i-1} \end{aligned} \tag{10}$$

其中 $s_{\omega i} \triangleq \text{sign}\{\text{Re}[P(j\omega_i)]\}$。消去 $\hat{w}_{2i}$：

$$\hat{w}_{2i-1}(s) = \frac{k \cdot s_{\omega i} \cdot s}{s^2 + \hat{\theta}_i} y(s), \quad \hat{w}_{2i}(s) = -\frac{k \cdot s_{\omega i} \cdot \hat{\theta}_i}{s^2 + \hat{\theta}_i} y(s) \tag{11}$$

同理，$\hat{w}_0(s) = \frac{k \cdot \text{sign}[P(0)]}{s} y(s)$。

将以上代入 (9)，得到 $y$ 满足的闭环特征方程：

$$y(s) = P(s) \left[-\frac{k \cdot s_0}{s} - \sum_{i=1}^q \frac{k \cdot s_{\omega i} \cdot s}{s^2 + \hat{\theta}_i} + 1\right]^{-1} d(s) \tag{12}$$

定义灵敏度函数类 $G(s)$，使得 $y = G(s) \prod_{i=1}^q (s^2 + \hat{\theta}_i) \cdot (\text{扰动项})$。当 $\hat{\theta}_i \to \omega_i^2$ 时，闭环传递函数在 $s = \pm j\omega_i$ 处具有零点——这正是内模原理产生**扰动抵消**的频域解释。

### 3. 平均化理论 (Averaging Theory) 处理

当自适应增益 $\varepsilon$ 足够小 ($0 < \varepsilon \ll 1$) 时，$\hat{\theta}$ 的演化速度远慢于状态变量。这是**双时间尺度**系统：状态以 $O(1)$ 时间常数演化，而参数以 $O(1/\varepsilon)$ 时间常数演化。

考虑"冻结"参数 $\hat{\theta}$ 下的稳态响应。设 $y_{ss}(t)$ 和 $\hat{w}_{2i,ss}(t)$ 为固定 $\hat{\theta}$ 时系统 (5) 的稳态解。自适应律 (6) 的平均系统为：

$$\dot{\tilde{\theta}}_{av} = \varepsilon \cdot f_{av}(\tilde{\theta}_{av}) \tag{13}$$

其中平均向量场定义为：

$$f_{av,i}(\tilde{\theta}) = -\frac{1}{T_i} \int_0^{T_i} s_{\omega i} \cdot y_{ss}(t) \cdot \hat{w}_{2i,ss}(t) \, dt \tag{14}$$

$T_i = 2\pi/\omega_i$ 为第 $i$ 个扰动分量的周期。

### 4. 关键推导：稳态响应的解析计算

当 $\hat{\theta}$ 固定时，将 $y_{ss}$ 和 $\hat{w}_{2i,ss}$ 用扰动模态展开。对每个频率 $\omega_i$，$y_{ss}$ 包含 $\sin(\omega_i t)$ 和 $\cos(\omega_i t)$ 分量。由频域关系 (11)：

$$\hat{w}_{2i,ss}(s) = -\frac{k \cdot s_{\omega i} \cdot \hat{\theta}_i}{s^2 + \hat{\theta}_i} y_{ss}(s) \tag{15}$$

在稳态正弦激励下，$y_{ss}(t) = |G_y(j\omega_i)| A_i \sin(\omega_i t + \phi_i + \angle G_y(j\omega_i))$。代入上式并利用 $\frac{1}{s^2 + \hat{\theta}_i}$ 在频率 $\omega_i$ 处的频率响应，得到：

$$y_{ss}(t) = |P(j\omega_i)| A_i \sin(\omega_i t + \psi_i) + \text{其他频率分量} \tag{16}$$

$$\hat{w}_{2i,ss}(t) = -k \cdot s_{\omega i} \cdot \hat{\theta}_i \cdot \frac{|P(j\omega_i)| A_i}{\omega_i^2 - \hat{\theta}_i} \cos(\omega_i t + \psi_i) + \cdots \tag{17}$$

其中 $\psi_i$ 是相位偏移。将 (16)-(17) 代入平均积分 (14)，利用三角函数的正交性（$\int_0^{2\pi/\omega_i} \sin(\omega_i t + \alpha) \cos(\omega_i t + \alpha) dt = 0$，且 $\int_0^{2\pi/\omega_i} \sin^2(\omega_i t) dt = \pi/\omega_i$）：

$$\begin{aligned} f_{av,i}(\tilde{\theta}) &= -\frac{1}{T_i} \int_0^{T_i} s_{\omega i} \cdot \left[|P| A_i \sin(\omega_i t + \psi_i)\right] \cdot \left[-k \cdot s_{\omega i} \cdot \hat{\theta}_i \cdot \frac{|P| A_i}{\omega_i^2 - \hat{\theta}_i} \cos(\omega_i t + \psi_i)\right] dt \\ &= -\frac{k \cdot \hat{\theta}_i \cdot |P(j\omega_i)|^2 A_i^2}{2(\hat{\theta}_i - \omega_i^2)} + O(\varepsilon) \end{aligned} \tag{18}$$

注意 $(s_{\omega i})^2 = 1$（符号的平方恒为 1）。

当 $\hat{\theta}_i \approx \omega_i^2$ 时，令 $\tilde{\theta}_i = \omega_i^2 - \hat{\theta}_i$，则 $\hat{\theta}_i - \omega_i^2 = -\tilde{\theta}_i$。代入得到平均动力学的线性化形式：

$$\boxed{f_{av,i} = -\frac{k \cdot \omega_i^2 \cdot |P(j\omega_i)|^2 A_i^2}{2} \cdot \tilde{\theta}_i \triangleq -\kappa_i \tilde{\theta}_i} \tag{19}$$

由于 $k > 0$、$\omega_i^2 > 0$、$|P(j\omega_i)|^2 > 0$、$A_i^2 > 0$，有 $\kappa_i > 0$。

### 5. 指数收敛性

将 (19) 代入平均系统 (13)：

$$\dot{\tilde{\theta}}_{av,i} = -\varepsilon \kappa_i \tilde{\theta}_{av,i} \tag{20}$$

这是一阶线性衰减动力学，解为 $\tilde{\theta}_{av,i}(t) = \tilde{\theta}_{av,i}(0) \cdot e^{-\varepsilon \kappa_i t}$。

由平均化理论的**一阶近似定理**（Khalil, 2002, Theorem 10.4）：对足够小的 $\varepsilon$，实际误差 $\tilde{\theta}_i(t)$ 与平均系统解 $\tilde{\theta}_{av,i}(t)$ 之差在有限时间区间上为 $O(\varepsilon)$。结合 (20) 的指数稳定性，存在 $\varepsilon^* > 0$ 使得对所有 $\varepsilon \in (0, \varepsilon^*)$，$\tilde{\theta}_i(t) \to 0$ 指数收敛。

---

## 第四部分：结论 (Convergence Results)

### 1. 参数收敛

通过平均化分析证明了：

$$\lim_{t \to \infty} \hat{\theta}_i(t) = \omega_i^2, \quad i = 1, 2, \ldots, q \tag{21}$$

收敛速率约为 $O(e^{-\varepsilon \kappa_i t})$，其中 $\kappa_i = \frac{1}{2} k \omega_i^2 |P(j\omega_i)|^2 A_i^2$。

### 2. 输出收敛

当频率估计准确 ($\hat{\theta}_i \to \omega_i^2$) 时，补偿器变为精确内模。由于系统 $P(s)$ 稳定，且闭环矩阵在平衡点 $\tilde{\theta} = 0$ 处是 Hurwitz 的——其特征值为 $P(s)$ 的极点（左半平面）加上 $s = 0, \pm j\omega_1, \ldots, \pm j\omega_q$（虚轴上的一对内模极点）——根据内模原理，闭环系统满足**输出调节**性质：

$$\lim_{t \to \infty} y(t) = 0 \tag{22}$$

### 3. 总体结论

存在足够小的 $k > 0$ 和 $\varepsilon > 0$，使得闭环系统 (1)-(5)-(6) 在平衡点

$$\{\hat{w} = w, \hat{\theta}_i = \omega_i^2, y = 0\}$$

处是**局部指数稳定**的。

---

## 应试得分要点

| 要点 | 说明 |
|:---|:---|
| **内模原理 (IMP)** | 补偿器结构 $[s(s^2+\hat{\theta}_i)]^{-1}$ 是干扰发生器 $[s(s^2+\omega_i^2)]^{-1}$ 的估计拷贝 |
| **符号函数的作用** | $\text{sign}\{\text{Re}[P(j\omega_i)]\}$ 决定反馈极性，源自 MRAC 的严格正实条件——符号错误导致梯度方向反转、参数发散 |
| **双时间尺度** | 状态快速 ($\sim 1/k$)、参数慢速 ($\sim 1/\varepsilon$)，$\varepsilon \ll k$ 是分离条件 |
| **平均化方法** | 处理非线性参数化（$\hat{\theta}_i$ 出现在 $\sin/\cos$ 的自变量中实为频率）的标准严谨手段 |
| **正交性积分** | 平均化计算中 $\int \sin(\omega_i t)\cos(\omega_i t) = 0$ 和 $\int \sin^2 = 1/2$ 是关键技巧 |
| **持续性激励** | $A_i > 0$ 保证扰动持续存在，为参数收敛提供必要条件 |
| **Case A vs Case B** | 本文仅推导 Case A ($\text{Re}[P(j\omega_i)] \neq 0$)；当实部为零时需切换至 Im 符号和正交补偿器结构 |

---

## 与仿真代码的对应关系

| 数学符号 | 代码变量 (plant 模式) | 所在行 |
|:---|:---|:---|
| $\hat{w}_0$ | `y(2)` = w0 | 228 |
| $\hat{w}_1$ | `y(3)` = w1 | 229 |
| $\hat{w}_2$ | `y(4)` = w2 | 230 |
| $\hat{w}_3$ | `y(5)` = w3 | 231 |
| $\hat{w}_4$ | `y(6)` = w4 | 232 |
| $\hat{\theta}_1$ | `y(7)` = theta1 | 233 |
| $\hat{\theta}_2$ | `y(8)` = theta2 | 234 |
| $y$ | `y(1)` = x1 | 222 |
| $k$ | `k` = 0.5 | 21 |
| $\varepsilon$ | `epsilon` = 0.05 | 22 |
| $u = -\hat{w}_0 - \hat{w}_1 - \hat{w}_3$ | `u = -y(2) - y(3) - y(5)` | 220 |

---

## 参考文献

1. Marino, R. & Tomei, P. (2016). "Adaptive disturbance rejection for unknown stable linear systems." *Transactions of the Institute of Measurement and Control*, 38(6): 640–647.
2. Khalil, H.K. (2002). *Nonlinear Systems*, 3rd ed. Prentice Hall. (平均化理论: Theorem 10.4)
3. Francis, B.A. & Wonham, W.M. (1976). "The internal model principle of control theory." *Automatica*, 12(5): 457–465.
4. Marino, R. & Tomei, P. (2011). "An adaptive learning regulator for uncertain minimum phase systems with undermodeled unknown exosystems." *Automatica*, 47(4): 739–747.
