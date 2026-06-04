针对该论文中提出的**未知周期干扰的鲁棒自适应抑制问题**，其核心思想是结合“内模原理/多项式方程”与“鲁棒自适应滤波”。由于直接提取干扰的频率可能存在鲁棒性问题，作者将问题巧妙地转化为**FIR滤波器的过参数化设计**和**带有鲁棒修正的递归最小二乘(RLS)参数估计**。

以下是完整的系统问题制定与数学推导过程，整理为Markdown文档：

---

# 鲁棒自适应未知周期干扰抑制：问题制定与数学推导

## 1. 系统与干扰建模

### 1.1 被控对象模型
考虑一个带有未建模动态的稳定离散时间单输入单输出 (SISO) 系统：
$$
\begin{aligned}
y(k) &= G(z)[u(k)] + d(k) \\
     &= G_0(z)\big(1 + \Delta_m(z)\big)[u(k)] + d(k)
\end{aligned}
$$
其中：
*   $G_0(z) = \frac{Z_0(z)}{R_0(z)}$ 为已知的标称模型（$R_0(z)$ 是 $n$ 阶首一 Hurwitz 多项式）。
*   $\Delta_m(z)$ 为未知的乘性未建模动态。
*   $d(k)$ 为外界干扰。

### 1.2 干扰模型
外界干扰 $d(k)$ 假设由未知的周期信号和有界随机噪声组成：
$$
d(k) = d_s(k) + v(k) = \sum_{i=1}^{n_\omega} \kappa_i \sin(\omega_i T_s k + \varphi_i) + v(k)
$$
由于正弦信号 $d_s(k)$ 可被线性差分方程生成，其对应的内部生成多项式（零化多项式）已知为：
$$
D_s(z) = \prod_{i=1}^{n_\omega} (z^2 - 2\cos(\omega_i T_s)z + 1)
$$
**控制目标**：设计控制输入 $u(k)$，在保证闭环系统鲁棒稳定的前提下，完全抑制正弦干扰 $d_s(k)$，并最小化噪声 $v(k)$ 对输出 $y(k)$ 的放大效应。

---

## 2. 闭环控制结构与问题转化

文章采用了一种类似于内模控制(IMC) / Q-参数化的结构，引入辅助信号 $\zeta(k)$ 和前馈FIR滤波器 $F(z)$：
$$
\zeta(k) = y(k) - G_0(z)[u(k)]
$$
控制器结构定义为：
$$
u(k) = -F(z)[\zeta(k)]
$$
其中 $F(z)$ 被参数化为 $N$ 阶FIR滤波器：
$$
F(z, \theta) = \frac{\Theta(z)}{z^N} = \sum_{i=0}^{N-1} \theta_i z^{i-N} = \theta^T \alpha(z)
$$
其中，参数向量 $\theta = [\theta_0, \dots, \theta_{N-1}]^T$，延迟向量 $\alpha(z) = [z^{-N}, z^{1-N}, \dots, z^{-1}]^T$。

**闭环传递函数推导**：
将 $u = -F \zeta$ 代入 $\zeta$ 的定义式中：
$$
\zeta = y - G_0 u \implies \zeta = \big(G u + d\big) - G_0 u = d + (G - G_0)u = d + G_0 \Delta_m u
$$
进一步代入 $u = -F \zeta$：
$$
\zeta = d - G_0 \Delta_m F \zeta \implies \zeta = \frac{1}{1 + G_0 \Delta_m F} d
$$
此时系统输出 $y$ 为：
$$
y = \zeta + G_0 u = \zeta - G_0 F \zeta = (1 - G_0 F) \zeta
$$
因此，从干扰 $d$ 到输出 $y$ 的闭环敏感度函数 $H(z)$ 为：
$$
y(k) = H(z)[d(k)] = \frac{1 - G_0(z) F(z, \theta)}{1 + \Delta_m(z) G_0(z) F(z, \theta)} \big[d_s(k) + v(k)\big]
$$
**稳定性约束**：由小增益定理，为保证系统在存在未建模动态 $\Delta_m(z)$ 时依然稳定，需满足：
$$
\|\Delta_m(z)G_0(z)F(z, \theta)\|_\infty < 1
$$

---

## 3. 理想情况(已知频率)：问题制定与最优求解方程

在理想情况下，如果干扰频率 $\omega_i$（即 $D_s(z)$）已知，我们希望 $y(k)$ 完全抑制 $d_s(k)$。
从上述 $H(z)$ 公式可知，当且仅当敏感度函数的分子 $(1 - G_0 F)$ 的零点包含 $D_s(z)$ 的所有根时，正弦干扰被完全抑制。

展开 $1 - G_0 F$：
$$
1 - G_0 F = 1 - \left(\frac{Z_0}{R_0}\right) \left(\frac{\Theta}{z^N}\right) = \frac{z^N R_0 - Z_0 \Theta}{z^N R_0}
$$
要让该分式的分子含有 $D_s(z)$ 作为因子，必然存在一个多项式 $L(z)$ 满足丢番图方程 (Diophantine Equation)：
$$
z^N R_0(z) - Z_0(z) \Theta(z) = D_s(z) L(z)
$$
即：
$$
Z_0(z)\Theta(z) + D_s(z)L(z) = z^N R_0(z)
$$
**最优性能问题制定**：
当滤波器阶数 $N > n_d$（过参数化）时，上述方程存在无穷多解。为了减少噪声 $v(k)$ 的放大，应最小化 $G_0 F$ 的峰值。因此构成约束优化问题：
$$
\begin{aligned}
\min_{\theta} \quad & J(\theta) = \|G_0 F\|_\infty = \left\| \frac{Z_0 \Theta}{z^N R_0} \right\|_\infty \\
\text{s.t.} \quad & Z_0 \Theta + D_s L = z^N R_0 \quad (\text{完全抑制条件})
\end{aligned}
$$
令这个优化的最优解对应的参数向量为 $\theta^*$，在理想情况下，$\theta^*$ 即为我们寻找的“最优抗扰滤波参数”。

---

## 4. 实际情况(未知频率)：自适应问题制定

在实际中，干扰频率未知，因此 $\theta^*$ 未知，需要设计自适应律 $\hat{\theta}(k)$ 实时估计。
此时控制器修改为带有增益 $\tau_0$ 的时变自适应形式：
$$
u(k) = -\tau_0 \hat{\theta}^T(k-1)\alpha(z)[\zeta(k)] = -\tau_0 \hat{\theta}^T(k-1) w(k)
$$
其中 $w(k) = \alpha(z)[\zeta(k)]$ 为滤波器的状态向量。

### 4.1 提取线性参数模型 (Linear Parametric Model)
为了使用自适应算法，必须将不可测的 $\theta^*$ 提取到可测的回归模型中。
令系统存在的最优滤波参数向量为 $\theta^*$（使得增益 $1 - \tau_0 G_0 F(z, \theta^*)$ 在干扰频率处为零）。
将 $\zeta(k)$ 作用于该理想函数：
$$
\big(1 - \tau_0 G_0(z)F(z, \theta^*)\big)[\zeta(k)] = \big(1 - \tau_0 G_0(z)F(z, \theta^*)\big)[d_s(k) + v(k) + u_\Delta(k)]
$$
因为 $\theta^*$ 的定义即为了抵消 $d_s(k)$，所以 $(1 - \tau_0 G_0 F(z, \theta^*))[d_s(k)] = 0$。剩余项定义为不可测残差 $\eta(k)$：
$$
\eta(k) = \big(1 - \tau_0 G_0(z)F(z, \theta^*)\big)[v(k) + u_\Delta(k)]
$$
同时展开等式左侧：
$$
\begin{aligned}
\big(1 - \tau_0 G_0(z)F(z, \theta^*)\big)[\zeta(k)] &= \zeta(k) - \tau_0 G_0(z)\big[\theta^{*T}\alpha(z)\zeta(k)\big] \\
&= \zeta(k) - \theta^{*T} \big(\tau_0 G_0(z)[w(k)]\big)
\end{aligned}
$$
定义回归向量 $\phi(k) = \tau_0 G_0(z)[w(k)]$（它是已知或可测的），则有极其优雅的参数化等式：
$$
\zeta(k) - \theta^{*T}\phi(k) = \eta(k) \implies \zeta(k) = \theta^{*T}\phi(k) + \eta(k)
$$
这就是用于参数更新的**标准线性参数模型**（等式 11）。

### 4.2 Swapping Lemma处理可测输出
真实输出 $y(k)$ 与回归向量的关系需要利用交换引理(Swapping Lemma)处理时变的 $\hat{\theta}$：
$$
G_0(z)[\hat{\theta}^T w] = \hat{\theta}^T G_0(z)[w] + \eta_{swap}
$$
其中 $\eta_0 = \tau_0 \eta_{swap}$ 包含参数变化量 $\Delta\hat{\theta}(k)$ 导致的暂态误差。代入 $y = \zeta + G_0 u$：
$$
y = \zeta - \tau_0 G_0(z)[\hat{\theta}^T w] = \zeta - \tau_0\hat{\theta}^T G_0(z)[w] - \eta_0 = \zeta - \hat{\theta}^T\phi - \eta_0
$$

---

## 5. 鲁棒自适应控制律推导与求解

由于未知干扰存在多个频率分量，为了保证充分的自由度同时满足 $H_\infty$ 最小化约束，论文设定滤波器阶数 $N > 2\bar{n}_\omega$（即过参数化）。
过参数化导致**持续激励条件 (Persistent Excitation, PE) 丢失**。在失去PE条件且存在建模误差 $\eta(k)$ 的情况下，普通的自适应律会导致参数漂移(Parameter Drift)。

因此，问题求解最终落脚于**具有归一化和参数投影的鲁棒递归最小二乘算法 (Robust RLS with Projection)**。

### 5.1 误差与归一化设计
构建预测信号：
$$
\hat{\zeta}(k) = \hat{\theta}^T(k-1)\phi(k)
$$
根据 4.2 节推导可知：
$$
\zeta(k) - \hat{\zeta}(k) = \zeta(k) - \hat{\theta}^T(k-1)\phi(k) = y(k) + \eta_0(k)
$$
忽略暂态误差 $\eta_0$，可以直接利用被控对象输出 $y(k)$ 作为误差信号，归一化误差定义为：
$$
\varepsilon(k) = \frac{y(k) + \eta_0(k)}{m^2(k)} \approx \frac{y(k)}{m^2(k)}
$$
归一化信号：
$$
m^2(k) = 1 + c_0 \phi^T(k)\phi(k) \quad (c_0 > 0)
$$

### 5.2 鲁棒RLS参数更新求解公式
应用鲁棒RLS更新滤波器参数 $\hat{\theta}(k)$：
**1. 协方差矩阵更新：**
$$
P(k) = P(k-1) - \frac{P(k-1)\phi(k)\phi^T(k)P(k-1)}{m^2(k) + \phi^T(k)P(k-1)\phi(k)}
$$
*(实际应用中需附带协方差重置 (Covariance Resetting) 防止协方差矩阵奇异或衰减过小)*

**2. 参数向量更新与投影：**
$$
\hat{\theta}(k) = \text{proj} \Big( \hat{\theta}(k-1) + P(k)\varepsilon(k)\phi(k) \Big)
$$
其中，投影算子 $\text{proj}(\cdot)$ 确保参数向量 $\hat{\theta}(k)$ 始终限制在紧集 $\mathcal{S}(\theta_{\max}) = \{x \in \mathbb{R}^N \mid x^T x \le \theta_{\max}^2\}$ 内，以此硬约束切断了参数漂移的可能性。

### 5.3 闭环求解总结 (Theorem 2保证)
通过上述制定与求解机制：
1. **控制律合成**：每个时刻解算出 $\hat{\theta}(k)$ 后，合成控制信号 $u(k) = -\tau_0 \hat{\theta}^T(k-1)\alpha(z)[\zeta(k)]$。
2. **渐近性能**：证明结果指出，闭环系统全局稳定，且输出能量界由建模误差规模和噪声方差决定：
$$
\frac{1}{T} \sum_{i=k}^{k+T-1} |y(i)|^2 \le c(\mu_\Delta^2 + v_0^2) + \frac{c}{T}
$$
无噪声和未建模动态时，$y(k) \to 0$。这表明系统能够在无需预先显式计算干扰频率的条件下，自适应求解出逼近最优的 $\theta^*$ 进行干扰抵消。