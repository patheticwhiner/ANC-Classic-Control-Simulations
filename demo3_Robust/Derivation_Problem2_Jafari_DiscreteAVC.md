# 题2：未知频率周期扰动的离散时间鲁棒自适应抑制

> **Jafari & Ioannou (2015 TAC / 2015 JVC) 离散时间 AVC 控制器 · 考试答题级题解**
>
> **前置阅读**：[Theory_Foundations.md](Theory_Foundations.md) Part A（自适应控制）和 Part B（H∞ 鲁棒控制）。本题同时用到两部分理论。

---

## 绪论

### 你的知识基线

| 领域 | 你的状态 | 这意味着 |
|:---|:---|:---|
| 离散时间系统 | 能跟上推导 | $z^{-1}$ 延迟算子、$z = e^{j\omega T_s}$ 频率映射不陌生 |
| FIR 滤波器 | 熟悉 | $\sum \theta_k z^{-k}$ 对你不是新概念 |
| H∞ 优化 | 有模糊概念 | 已通过 [Theory_Foundations.md Part B](Theory_Foundations.md) 建立基础 |
| RLS 自适应 | 有模糊概念 | 已通过 [Theory_Foundations.md §A.4](Theory_Foundations.md) 理解 |

**预计阅读时间**：80 分钟。

---

## 题目

> **已知**：
>
> **(1)** 离散时间被控对象（$T_s = 1/480$ s）：
>
> $$\boxed{G_0(z) = \frac{-0.00146\,(z - 0.1438)(z - 1)}{(z - 0.7096)(z^2 - 0.04369z + 0.01392)}}$$
>
> **(2)** 扰动信号（单频正弦 + 有界噪声）：
>
> $$d(k) = \sin(\omega_0 k) + v(k), \qquad \omega_0 = 0.0521 \text{ rad/sample } (\approx 25 \text{ rad/s } \approx 3.98 \text{ Hz})$$
>
> $|v(k)| \le 0.1$，标准差 $\sigma_v = 0.02$。
>
> **(3)** 控制架构：
>
> $$u(k) = -K(z,\theta)[y(k)]$$
>
> **求**：
> (a) $G_0(z)$ 的关键特征分析（极点、零点、在扰动频率处的增益）
> (b) 固定 FIR 控制器的 H∞ 最优设计（含插值约束）
> (c) 自适应 FIR 控制器的 RLS 在线更新律
> (d) 鲁棒稳定性条件的数值验证
> (e) $N$ 对性能的影响和频率偏移鲁棒性

---

## 解

### 第1步 · 审题

本题与题1的关键差异：

| 特征 | 题1（连续 CC） | 题2（离散 AVC） |
|:---|:---|:---|
| 时间域 | 连续 $s$ | 离散 $z$，$T_s = 1/480$ |
| 基函数 | Laguerre $\frac{\lambda^{k+1}}{(s+\lambda)^{k+1}}$ | FIR $z^{-k}$（延迟线） |
| 扰动 | 多频（2 个正弦） | 单频（1 个正弦） |
| 固定控制器 | $F(s)$ 频谱展平 | H∞ 优化 + 插值约束 |
| 理论依赖 | Part A 全部 | Part A + Part B |

**核心困难**：$G_0(z)$ 在 $\omega_0 = 0.0521$ rad/sample 处的增益极小（约 $2.3 \times 10^{-4}$，即 −73 dB），意味着控制器需要极大的增益（≈4300）才能在该频率产生有效抵消。这带来了参数尺度问题。

---

### 第2步 · 被控对象特征分析

#### 2.1 零极点

$$G_0(z) = \frac{-0.00146(z - 0.1438)(z - 1)}{(z - 0.7096)(z^2 - 0.04369z + 0.01392)}$$

- **极点**：$p_1 = 0.7096$（实极点，低频主导），$p_{2,3} = 0.02185 \pm j0.1157$（复共轭，$|p_{2,3}| \approx 0.118$，对应约 8.8 Hz 的阻尼模态）
- **零点**：$z_1 = 1$（积分器——$G_0(e^{j0}) = 0$，DC 处完全阻塞），$z_2 = 0.1438$（实零点，稳定）
- **DC 增益**：$G_0(1) = 0$（因分子有因子 $(z-1)$）

#### 2.2 扰动频率处的增益

在 $\omega_0 = 0.0521$ rad/sample：

$$G_0(e^{j\omega_0}) \approx 2.30 \times 10^{-4} \; \angle \; -91.7^\circ$$

即幅值约 −72.8 dB。这意味着要使 $|G_0(e^{j\omega_0})K(e^{j\omega_0})| = 1$（插值约束），控制器在扰动频率处需要约 $|K(e^{j\omega_0})| \approx 4350$（+72.8 dB）的增益。这解释了为何代码中需要 $\theta_{\max}$ 很大（$10^7$）且 $P(0)$ 很大（$10^6 I$）。

---

### 第3步 · 固定 FIR 控制器的 H∞ 最优设计

#### 3.1 控制器参数化

取 $N$ 阶 FIR 滤波器：

$$
\boxed{K(z, \theta) = \sum_{k=0}^{N-1} \theta_k z^{-k} = \theta^T \alpha(z)} \tag{2.1}
$$

其中 $\theta = [\theta_0, \dots, \theta_{N-1}]^T$，$\alpha(z) = [1, z^{-1}, \dots, z^{-(N-1)}]^T$。

本题取 $N = 30$（提供充足自由度以同时满足插值约束和最小化 H∞ 范数）。$\theta_k$ 均为**实数**。

#### 3.2 优化问题形式化

**目标**：在所有满足插值约束的 FIR 系数中，最小化闭环传递函数 $G_0(z)K(z,\theta)$ 的 H∞ 范数。

（注：IMC 架构下闭环输出 $y = (1 - G_0K)[d]$，$G_0K$ 在扰动频率处应为 1 以实现完全抵消，同时 $\|G_0K\|_\infty$ 应尽可能小以限制噪声放大。）

**插值约束**（强制在扰动频率处完全抵消）：

$$
\boxed{G_0(e^{j\omega_0}) \cdot K(e^{j\omega_0}, \theta) = 1} \tag{2.2}
$$

这是一个复数等式约束。因 $\theta$ 为实数，将其拆分为实部和虚部两个实约束：

$$
\begin{bmatrix} \operatorname{Re}\{G_0(e^{j\omega_0})\alpha^T(e^{j\omega_0})\} \\ \operatorname{Im}\{G_0(e^{j\omega_0})\alpha^T(e^{j\omega_0})\} \end{bmatrix} \theta = \begin{bmatrix} 1 \\ 0 \end{bmatrix} \tag{2.3}
$$

记 $A_{eq} \theta = b_{eq}$（$A_{eq} \in \mathbb{R}^{2 \times N}$）。

**H∞ 范数约束**（对频率网格 $\{\omega_k\}_{k=1}^{N_f}$，$N_f = 2000$）：

$$|G_0(e^{j\omega_k}) \cdot K(e^{j\omega_k}, \theta)| \le \gamma, \qquad k = 1, \dots, N_f$$

这给出了 $N_f$ 个二阶锥约束 [Theory_Foundations.md 引理 B.1](Theory_Foundations.md)：

$$\sqrt{\big(\operatorname{Re}\{G_0(e^{j\omega_k})K(e^{j\omega_k})\}\big)^2 + \big(\operatorname{Im}\{G_0(e^{j\omega_k})K(e^{j\omega_k})\}\big)^2} \le \gamma$$

#### 3.3 标准形式与求解

综上，优化问题为：

$$
\boxed{\begin{aligned}
\underset{\theta \in \mathbb{R}^N, \; \gamma \ge 0}{\min} \quad & \gamma \\
\text{s.t.} \quad & A_{eq} \theta = b_{eq} \quad \text{(插值约束)} \\
& |G_0(e^{j\omega_k}) \cdot \alpha^T(e^{j\omega_k})\theta| \le \gamma, \quad k = 1, \dots, N_f \\
& \|\theta\|_2 \le R_0 \quad \text{(参数范数上界，可选)}
\end{aligned}} \tag{2.4}
$$

这是一个**二阶锥规划（SOCP）**——凸优化问题，全局最优解唯一。

三种求解策略：

| 方法 | 工具 | 用途 |
|:---|:---|:---|
| **最小二乘（LS）** | `Aeq \ beq` | 仅满足插值约束，取最小范数解。$\gamma_{LS}$ 通常较大但 $\|\theta\|$ 最小——在实践中降噪效果最好！ |
| **CVX** | `cvx_begin ... cvx_end` | 直接建模求解 SOCP，得理论最优 $\gamma$ |
| **YALMIP + MOSEK/SeDuMi** | `sdpvar` + `optimize` | 工业级 SOCP 求解器，$\gamma$ 与 CVX 一致 |

**关键发现**：对此 $G_0$（极低的扰动频率增益），不加 $\|\theta\|$ 约束的 H∞ 最优解会产生 $\|\theta\| \approx 7 \times 10^5$ 的极端参数——在非 $\omega_0$ 频率放大了宽带噪声，**反而增噪**。LS 最小范数解在实践中给出最好的实际降噪效果。因此在代码中 LS 解被选为最终固定控制器。

#### 3.4 结果

$$N = 30: \quad \gamma_{LS} \approx 0.0137, \quad \|\theta_{LS}\| \approx 7.1 \times 10^4$$

小 $\gamma$ 意味着 $|G_0K| \le 0.0137$ 在所有频率成立——在扰动频率处增益为 1（插值），在其他所有频率处被严格抑制到 < 0.014。

---

### 第4步 · 自适应 FIR 控制器

#### 4.1 参数化与控制律

与固定控制器相同的 FIR 结构，但参数在线更新：

$$K(z, \theta(k)) = \sum_{i=0}^{N-1} \theta_i(k) z^{-i} \tag{2.5}$$

控制律（IMC 结构）：

$$\boxed{u(k) = -K(z, \theta(k))[y(k)] = -\sum_{i=0}^{N-1} \theta_i(k) \cdot y(k-i)} \tag{2.6}$$

#### 4.2 回归向量的构造

在 IMC 架构中 $y = d + G_0[u]$，辅佐信号 $z = y - G_0[u]$。控制目标是使 $K(z,\theta^*) = G_0^{-1}(z)$（在扰动频率处）。

依 [Theory_Foundations.md §A.3](Theory_Foundations.md) 构造线性参数化。离散时间 FIR 情况下，回归向量直接取滤波后的延迟信号：

$$
\boxed{\phi(k) = G_0(z)F(z) \cdot \begin{bmatrix} y(k) \\ y(k-1) \\ \vdots \\ y(k-N+1) \end{bmatrix}} \tag{2.7}
$$

或者等价地（代码中的实现）：

$$\phi_i(k) = G_0(z)F(z)[y(k-i)], \qquad i = 0, \dots, N-1$$

其中 $F(z)$ 在此问题中为常数增益 $F = 100$（Paper 的早期版本未使用频谱展平滤波器 $F(s)$，改用常数增益补偿 $G_0$ 的低通衰减）。

#### 4.3 RLS 自适应律

直接引用 [Theory_Foundations.md §A.4 定理 A.3](Theory_Foundations.md)（离散时间 RLS）：

$$
\begin{aligned}
m_s^2(k) &= 1 + \gamma_0 \phi^T(k)\phi(k), \qquad \gamma_0 = 1 \\[4pt]
\varepsilon(k) &= \frac{y(k) - \theta^T(k-1)\phi(k)}{m_s^2(k)} \\[4pt]
K(k) &= \frac{P(k-1)\phi(k)}{m_s^2(k) + \phi^T(k)P(k-1)\phi(k)} \\[4pt]
P(k) &= P(k-1) - K(k)\phi^T(k)P(k-1) \\[4pt]
\theta(k) &= \operatorname{Proj}_\Theta\Big(\theta(k-1) + K(k)\varepsilon(k)\Big)
\end{aligned} \tag{2.8}
$$

初始化：$\theta(0) = 0$，$P(0) = 10^6 I_N$（匹配参数尺度），$\theta_{\max} = 10^7$（参考固定解 $\|\theta_{LS}\| \approx 7 \times 10^4$，取宽松上界）。

投影 $\operatorname{Proj}_\Theta$ 的定义见 [Theory_Foundations.md §A.5](Theory_Foundations.md)，球形约束 $\Theta = \{\theta : \|\theta\| \le \theta_{\max}\}$。

---

### 第5步 · 鲁棒稳定性

#### 5.1 闭环传递函数

IMC 结构中：$y = d + G_0[u] = d - G_0K[y]$，故：

$$y = \frac{1}{1 + G_0K}[d], \qquad u = -\frac{K}{1 + G_0K}[d]$$

$K(s,\theta)$ 参数变化时闭环的鲁棒稳定性由 [Theory_Foundations.md 引理 A.3](Theory_Foundations.md) 保证。

#### 5.2 参数尺度与稳定性

固定解 $\|\theta_{LS}\| \approx 7.1 \times 10^4$ 看似大，但因 $G_0$ 在扰动频率处的增益仅 $2.3 \times 10^{-4}$，$|G_0K| \le 0.0137 < 1$ 在所有频率成立——Nyquist 轨迹远在单位圆内，稳定裕度充裕。

---

### 第6步 · $N$ 的影响与频率鲁棒性

#### 6.1 FIR 阶数 $N$ 的扫描

代码中测试了 $N = 5, 10, 15, 20, 30, 50$：

| $N$ | $\gamma$ | 抑制 (dB) | 备注 |
|:---|:---|:---|:---|
| 5 | 较大 | 差 | 自由度不足以同时满足插值和抑制带外增益 |
| 10 | 中等 | ~10 dB | 基本可用 |
| 15 | 小 | ~20 dB | 接近饱和 |
| 30 | 0.0137 | ~25 dB | **推荐**——$N > 30$ 改善不显著 |

**理论解释**：插值约束消耗 2 个实自由度（复数等式 → 2 个实等式），剩余 $N - 2$ 个自由度用于在全频带最小化 $\|G_0K\|_\infty$。$N$ 越大，带外抑制越充分——但边际收益递减。

#### 6.2 频率偏移鲁棒性

当实际扰动频率 $\omega$ 偏离设计频率 $\omega_0$ 时，理论抑制（由 $|1 - G_0(e^{j\omega})K(e^{j\omega})|$ 决定）随频偏增大而下降。代码扫描了 ±30% 频偏范围：

- 在 $\omega_0 \pm 10\%$ 内，抑制仍然 > 15 dB
- 超过 ±20% 偏移，抑制显著退化

这验证了适应性的必要性：若扰动频率未知或漂移，固定控制器不再最优，必须自适应在线跟踪。

---

### 第7步 · 数值结果

运行 `JafariJVC_DiscreteAVC.m`（50 s 仿真，控制器在 10 s 后开启）：

| 方案 | $\gamma$ | $\|\theta\|$ | 抑制 (dB) | 注 |
|:---|:---|:---|:---|:---|
| 固定 LS | 0.0137 | $7.1 \times 10^4$ | ~25 dB | 最小范数，实际最优 |
| 固定 H∞ | $< 0.01$ | $\approx 7 \times 10^5$ | 差 | 极端参数放大宽带噪声 |
| 自适应 RLS | — | 收敛至 ~$10^4$ | ~25 dB | 在线跟踪，无需已知 $\omega_0$ |

**自适应 vs. 固定的关键优势**：当 $t \ge 30$ s 时额外加入二次频率分量（85, 125 rad/s），自适应控制器自动适应并压制，而固定控制器对此完全失效（因其设计仅针对 $\omega_0$）。

---

### 第8步 · 代码映射

#### 固定控制器部分

| 数学量 | 代码行 (`JafariJVC_DiscreteAVC.m`) | MATLAB |
|:---|:---|:---|
| $G_0(z)$ | L11 | `G0 = (-0.00146*(z-0.1438)*(z-1)) / (...)` |
| $N = 30$ | L43 | `N_fir = 30` |
| $\Phi$ 矩阵 | L53–56 | `Phi_grid(:,k) = z_grid(:) .^ (-k)` |
| $A_{eq}, b_{eq}$ | L61–62 | `Aeq = [real(G_interp*Phi_interp); imag(...)]` |
| LS 解 | L65 | `theta_ls = Aeq \ beq` |
| $\gamma_{LS}$ | L66–67 | `mag_ls = abs(G_grid(:) .* (Phi_grid * theta_ls))` |
| H∞ 优化 | L77–102 | `yalmip(...)` + `optimize(Constraints, gamma_opt)` |
| 闭环仿真 | L149–157 | `u_fixed(k) = -sum(theta_fixed'.*x_fixed(...))` |

#### 自适应控制器部分

| 数学量 | 代码行 (`JafariJVC_DiscreteAVC.m`) | MATLAB |
|:---|:---|:---|
| $N = 30$ | L263 | `N_adapt = 30` |
| $P(0)$ | L269 | `P = eye(N_adapt) * 1e6` |
| $\theta_{\max}$ | L268 | `theta_max = 1e7` |
| 回归向量 $\phi$ | L299+ | `phi_reg(i) = ... G0(z)F(z)[z(k-i)]` |
| RLS 更新 | L299+ | `K = P*phi / (m2 + phi'*P*phi)` |

---

### 终 · 回到你的问题

| 子问题 | 答案 |
|:---|:---|
| (a) $G_0(z)$ 特征 | 三阶离散系统，有 $(z-1)$ 积分因子（DC 零增益），扰动频率增益 ≈ −73 dB，控制需要极大增益 |
| (b) 固定 H∞ 设计 | $N = 30$ FIR + 插值约束 $G_0(e^{j\omega_0})K(e^{j\omega_0}) = 1$ → SOCP 凸优化。LS 最小范数解 $\gamma \approx 0.014$，实践中优于理论 H∞ 解 |
| (c) 自适应律 | 离散 RLS [式 (2.8)]，$P(0) = 10^6 I$，$\theta_{\max} = 10^7$，投影保证有界 |
| (d) 鲁棒稳定性 | $\|G_0K\|_\infty \le 0.014 < 1$，Nyquist 裕度充裕 |
| (e) $N$ 的影响 | $N = 30$ 饱和；频率偏移 ±10% 内保留 >15 dB 抑制，超出后显著退化——验证了自适应的必要性 |

---

> 📖 **相关文档**：
> - 公共理论：[Theory_Foundations.md](Theory_Foundations.md)（Part A + Part B）
> - 连续时间版本：[题1 题解](Derivation_Problem1_Jafari_ContinuousCC.md)
> - H∞ 标准方法：[题3 题解](Derivation_Problem3_HinfSynthesis.md)
> - 代码 bug 详情：[DebugReport_JafariJVC_Discrete.md](Derivation_Problem2_Appendix_DebugReport.md)

---

# 附录

## 附录A · 离散域代码四重缺陷诊断

> 来源：`Derivation_Problem2_Appendix_DebugReport.md`

对照 2015 TAC/JVC 论文，原版 `JafariJVC_Discrete.m` 存在以下四重缺陷：

**缺陷 1：回归向量与控制律的因果顺序错误**

正确的因果顺序（论文 §4）：$y(k) \to z(k) \to \phi(k) \to \text{RLS} \to u(k) \to G_0$。原版代码将 RLS 更新置于控制律计算之前，导致第 $k$ 步使用了尚未计算的控制信息——形成代数环。

**缺陷 2：差分方程维度不匹配**

$G_0(z)$ 的递推式 $y(k) = -\sum a_i y(k-i) + \sum b_i u(k-i)$ 要求分子分母延迟对齐。原版代码中 $b$ 向量的延迟前导零未正确处理，导致 $u$ 和 $y$ 的时基错位。

**缺陷 3：RLS 归一化信号路径错误**

归一化信号 $m_s^2(k) = 1 + \gamma_0\phi^T(k)\phi(k)$ 应基于**更新前**的 $\theta(k-1)$ 计算先验误差 $e_0(k) = z(k) - \theta^T(k-1)\phi(k)$。原版混用了先验和后验误差。

**缺陷 4：符号翻转**

IMC 结构中 $u = -K(z)[z]$ 的负号在 FIR 实现中被遗漏——控制器变成了正反馈。

**修复验证**：修正后四重缺陷的修复版与独立验证的 `JafariTAC_DiscreteAVC.m`（不同论文的独立实现）输出一致。

---

## 附录B · MIMO 推广

> 来源：`Derivation_Problem2_Appendix_MIMO.md`（对应 Jafari & Ioannou, IEEE TAC, 2016）

### B.1 MIMO 问题制定

$$y(k) = G(z)[u(k)] + d(k) = (I + \Delta_m(z))G_0(z)[u(k)] + d(k)$$

其中 $y(k) \in \mathbb{R}^{n_y}$，$u(k) \in \mathbb{R}^{n_u}$。$G_0(z)$ 为 $n_y \times n_u$ 标称传递函数矩阵。

### B.2 MIMO 与 SISO 的关键差异

| | SISO（题2） | MIMO（TAC 2016） |
|:---|:---|:---|
| 插值约束 | $G_0(e^{j\omega_0})K(e^{j\omega_0}) = 1$（标量） | $G_0(e^{j\omega_0})K(e^{j\omega_0}) = I_{n_y}$（矩阵恒等式） |
| 控制器 | FIR 滤波器 | FIR 矩阵 $\Theta(z) \in \mathbb{R}^{n_u \times n_y}[z^{-1}]$ |
| RLS 回归 | 向量 $\phi \in \mathbb{R}^N$ | 矩阵 $\Phi \in \mathbb{R}^{N \times n_y}$ |
| 内模 | $D_s(z) = z^2 - 2\cos(\omega_0 T_s)z + 1$ | 对角阵 $\mathbf{D}_s(z)$，每通道独立内模 |

### B.3 MIMO 插值约束

复数矩阵等式拆为 $2n_y^2$ 个实约束，RLS 的协方差矩阵 $P \in \mathbb{R}^{N n_u \times N n_u}$——维数爆炸。论文中通过将 MIMO 解耦为 $n_y$ 个 MISO 子问题来规避。
