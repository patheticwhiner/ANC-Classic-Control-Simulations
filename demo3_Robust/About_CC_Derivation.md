# Jafari 2017 JVC 自适应 CC 控制器：算法推导、加速与实时性分析

## 读者定位

**论文复现验证者 & 实时实现工程师**——你已阅读 Jafari & Ioannou 2017 JVC 论文，关心：(1) 修正版代码与论文公式的逐项对应关系，(2) 从 O(N²) 到 O(N) 的加速推导，(3) 实时嵌入式部署的可行性。

---

## 1. 系统与控制结构

### 1.1 问题制定

被控对象（Eq. 1）：

$$
y(t) = G(s)[u(t)] + d(t), \quad G(s) = G_0(s)(1 + \Delta_m(s))
$$

其中 $G_0(s) = \frac{0.5(s-0.2)}{s^2 + s + 1.25}$，$\Delta_m(s) = -\mu s$（$\mu=0.001$）。

扰动 $d(t)$ 由未知频率的正弦分量加宽带噪声组成。

### 1.2 控制器架构（Eq. 6–8）

$$
\begin{aligned}
z(t) &= y(t) - G_0(s)[u(t)] \\[4pt]
u(t) &= -F(s)\big[K(s, \theta(t))[z(t)]\big]
\end{aligned}
$$

辅佐信号 $z(t)$ 剥离了标称模型 $G_0$ 的贡献。当 $G = G_0$（无建模误差）时，$z(t) = d(t)$。

### 1.3 自适应滤波器参数化（Eq. 3）

$K(s,\theta)$ 用 Laguerre 基函数展开：

$$
K(s, \theta) = \sum_{k=0}^{N-1} \theta_k \frac{\lambda^{N-k}}{(s+\lambda)^{N-k}}
= \theta^T \Lambda(s)
$$

其中 $\Lambda(s) = \left[\frac{\lambda^N}{(s+\lambda)^N}, \frac{\lambda^{N-1}}{(s+\lambda)^{N-1}}, \dots, \frac{\lambda}{s+\lambda}\right]^T$。

记 $\Lambda_i(s) = \frac{\lambda^{N-i+1}}{(s+\lambda)^{N-i+1}}$（$i=1,\dots,N$），即 $\Lambda_i$ 是 $N-i+1$ 个一阶节 $\frac{\lambda}{s+\lambda}$ 的级联。

### 1.4 F(s) 滤波器设计（Eq. 13）

$F(s)$ 的设计目标是将 $G_0(s)F(s)$ 整形为近似全通的频谱平坦化滤波器。对于含 RHP 零点 $s=0.2$ 的 $G_0(s)$，将该零点关于虚轴镜像反射至 $s=-0.2$，得到：

$$
F(s) = \kappa_0 \frac{2\alpha^2(s^2 + s + 1.25)}{(s+\alpha)^2(s+0.2)}
$$

其中 $\kappa_0 = 0.5$，$\alpha = 500$。

---

## 2. 自适应律推导

### 2.1 线性参数模型（Eq. 12）

利用 Swapping Lemma，定义回归向量：

$$
\phi(t) = G_0(s)F(s)\Lambda(s)[z(t)]
$$

则有线性参数化形式：

$$
z(t) = \theta^{*T}\phi(t) + \eta(t)
$$

其中 $\theta^*$ 为最优参数，$\eta(t)$ 为噪声/建模误差项。

### 2.2 连续时间 RLS 自适应律（Eq. 15–17）

归一化信号：$m_s^2(t) = 1 + \gamma_0 \phi^T(t)\phi(t)$

估计误差：$\varepsilon(t) = \frac{z(t) - \theta^T(t)\phi(t)}{m_s^2(t)}$

协方差更新（微分方程）：

$$
\dot{P}(t) = -P(t)\frac{\phi(t)\phi^T(t)}{m_s^2(t)}P(t)
$$

参数更新（带投影）：

$$
\dot{\theta}(t) = \text{proj}\Big(P(t)\varepsilon(t)\phi(t)\Big)
$$

### 2.3 Euler 离散化（数字实现）

采样周期 $T_s = 10^{-4}$ s，前向 Euler 近似：

$$
\begin{aligned}
P(k+1) &= P(k) - T_s \cdot P(k)\frac{\phi(k)\phi^T(k)}{m_s^2(k)}P(k) \\[4pt]
\bar{\theta}(k+1) &= \theta(k) + T_s \cdot P(k)\varepsilon(k)\phi(k) \\[4pt]
\theta(k+1) &= \begin{cases}
\bar{\theta}(k+1), & \|\bar{\theta}\| \le \theta_{\max} \\
\theta_{\max} \cdot \frac{\bar{\theta}(k+1)}{\|\bar{\theta}(k+1)\|}, & \text{otherwise}
\end{cases}
\end{aligned}
$$

控制律：

$$
u(k) = -F(z)\big[\theta^T(k)\Lambda(z)[z(k)]\big]
$$

**关键**：控制律中使用 $\Lambda(z)[z(k)]$（仅经过基函数滤波），而非 $\phi(k)$（还经过了 $G_0F$）。这是原版代码的核心 bug。

---

## 3. 原版代码四个缺陷的数学分析

### 缺陷 1：频率单位（Hz → rad/s）

原版将论文的 $\omega = 70, 187$ rad/s 误作 $f = 70, 187$ Hz，导致实际扰动角频率变为 $\omega' = 2\pi f \approx 440, 1175$ rad/s。$F(s)$ 的截止频率 $\alpha = 500$ rad/s，在 1175 rad/s 处增益已衰减至不足 −10 dB，控制器在扰动频率处无有效增益。

### 缺陷 2：参数个数 $N=6 \to 20$

论文要求 $N \ge 2n_{\max}$（$n_{\max}=5$），取 $N=20$。$N=6$ 意味着 FIR 滤波器仅有 6 个自由度，不足以在多个频率处构造零点。

### 缺陷 3：控制律信号混淆

原版：$u = -F[\theta^T\phi] = -F[\theta^T G_0 F \Lambda[z]]$

修正：$u = -F[\theta^T \Lambda[z]]$

差异是 $\phi$ 中额外的 $G_0F$。当 $\theta$ 收敛至 $\theta^T G_0 F \Lambda \approx I$ 时，原版退化为 $u \approx -F[z]$（固定增益反馈），丧失频率选择性。

### 缺陷 4：高阶多项式数值崩溃

$\Lambda_1(s) = \frac{\lambda^{20}}{(s+\lambda)^{20}}$ 的系数跨 50+ 个数量级，在 IEEE 754 double 精度下 Tustin 离散化后 `den(1)` 被舍入为零，差分方程除以零产生 NaN。

---

## 4. 从 $O(N^2)$ 到 $O(N)$ 的加速推导

### 4.1 原始算法复杂度

设 $N$ 个基函数，每个 $\Lambda_i$ 的级联深度为 $d_i = N-i+1$。总级联节数：

$$
\sum_{i=1}^{N} d_i = \sum_{i=1}^{N} (N - i + 1) = \frac{N(N+1)}{2}
$$

$N=20$ 时，每样本需 $\frac{20 \times 21}{2} = 210$ 次一阶 IIR 滤波。

加上每 $\Lambda_i$ 输出经过 $F$ 和 $G_0$（各 1 次），总计：

$$
\text{滤波次数/样本} = \frac{N(N+1)}{2} + 2N = \frac{N(N+5)}{2}
$$

$N=20$ 时为 $\frac{20 \times 25}{2} = 250$ 次。

### 4.2 共享级联链优化

**关键洞察**：对所有 $i=1,\dots,N$，$\Lambda_i$ 是级联链 $G_1 \to G_1^2 \to \dots \to G_1^N$（其中 $G_1(s)=\frac{\lambda}{s+\lambda}$）的不同深度输出。

具体地：

$$
\begin{aligned}
\Lambda_1[z] &= G_1^N[z] \quad\text{(经过全部 N 个一阶节)} \\
\Lambda_2[z] &= G_1^{N-1}[z] \quad\text{(经过前 N-1 个)} \\
&\vdots \\
\Lambda_N[z] &= G_1^1[z] \quad\text{(仅经过第 1 个)}
\end{aligned}
$$

因此，只需将 $z(k)$ 依次通过 $N$ 个一阶节 $G_1$，记录每个中间输出即可得到全部 $\Lambda_i[z]$：

```
y ← z(k)
for s = 1..N:
    y ← G₁[y]    // 第 s 个一阶节
    Λ_{N-s+1}[z] ← y   // 第 s 级输出 = Λ_{N-s+1}[z]
```

这一步将 $\Lambda$ 计算从 $\frac{N(N+1)}{2}$ 次降至 $N$ 次。

优化后总滤波次数：

$$
\text{滤波次数/样本} = N + 2N = 3N
$$

$N=20$ 时为 60 次。加速比：

$$
\frac{N(N+5)/2}{3N} = \frac{N+5}{6} \xrightarrow{N=20} \frac{25}{6} \approx 4.17\times
$$

### 4.3 实测加速数据

| 版本 | Λ 滤波/样本 | 总滤波/样本 | 10s 仿真耗时 | 样本/s | 加速比 |
|------|-----------|-----------|-------------|--------|--------|
| IIR_filter (原) | 210 | 250 | ~150 s | ~670 | 1× |
| filter() 独立级联 | 210 | 250 | 56.7 s | 1,763 | 2.6× |
| filter() 共享级联 | **20** | **60** | **15.8 s** | **6,347** | **9.5×** |

15.8 s 仿真 10 s → 实时因子 **1.58×**（仅比实时慢 58%）。

所有版本输出 −28.9 dB 抑制，$\|\theta\| = 2.51$，结果完全一致。

### 4.4 进一步优化路径

| 方案 | 方法 | 预期加速 | 代价 |
|------|------|---------|------|
| 降 $N$ 到 12 | $N \ge 2n_{\max} = 10$ 仍满足 | 1.7× | 频率选择性略微降低 |
| 降采样 $f_s$ | $f_s = 2000$ Hz 仍满足 Nyquist | 5× | 不能处理 >1000 Hz 扰动 |
| MEX/C 编译 | 核心循环用 C 重写 | 10–50× | 开发成本 |
| 组合：$N$=12 + $f_s$=2k + MEX | 全部叠加 | **>100×** | 最高开发成本 |

---

## 5. 实时嵌入式部署可行性分析

### 5.1 每样本计算量

当前算法（$N=20$, $f_s=10^4$）每样本主要操作：

| 操作 | 次数 | 每次运算量 | 总浮点运算 |
|------|------|----------|-----------|
| Λ 级联 (filter 1阶) | 20 | 2M + 1A | ~60 FLOP |
| F, G0 phi路径 (filter 2-3阶) | 40 | 5M + 4A | ~360 FLOP |
| P 更新 ($N \times N$) | 1 | $O(N^2)$ | ~800 FLOP |
| θ 更新 ($N \times 1$) | 1 | $O(N)$ | ~40 FLOP |
| **总计** | | | **~1260 FLOP** |

在 10 kHz 采样率下：$1260 \times 10^4 = 12.6$ MFLOP/s。这对嵌入式处理器（ARM Cortex-M7 @ 480 MHz ≈ 960 MMACS）而言**完全可行**。

### 5.2 实时性条件

设处理器能力 $C$ (FLOP/s)，采样周期 $T_s$：

$$
\frac{\text{FLOP/样本}}{C} < T_s \quad\Rightarrow\quad C > \frac{1260}{10^{-4}} = 12.6 \text{ MFLOP/s}
$$

现代单片机（如 STM32H7）单精度浮点 >200 MFLOP/s，留有 **16× 裕度**。

### 5.3 降 $N$ 和降采样后的最小配置

对于 $n_{\max}=2$（两个扰动频率），取 $N=6$，$f_s=2000$ Hz：

- 滤波次数/样本 = $3N = 18$ 次
- FLOP/样本 ≈ 300
- 所需：$300 \times 2000 = 0.6$ MFLOP/s
- 最低可在 **Cortex-M0+ @ 48 MHz** 上运行

---

## 6. 参数调优方法论

### 6.1 两个自由参数

| 参数 | 作用 | 增大效果 | 减小效果 |
|------|------|---------|---------|
| $\kappa_0$ | $F(s)$ 增益标量 | 抑制↑，控制能量↑，鲁棒裕度↓ | 抑制↓，控制能量↓，鲁棒裕度↑ |
| $\theta_{\max}$ | 参数投影界 | 允许更强的自适应 | 强制限制闭环增益，保证鲁棒性 |

### 6.2 鲁棒稳定性条件（定理 3）

$$
\kappa_0 \cdot \bar{\theta}_M \cdot \|\bar{F}G_0\Delta_m\|_1 < 1
$$

其中 $\bar{F} = F/\kappa_0$，$\bar{\theta}_M \ge \|\theta\|$ 为投影界（取 $\theta_{\max}$）。

对于当前系统：$\|\bar{F}G_0\Delta_m\|_\infty \approx 0.25$。

实际 $\|\theta\| \approx 2.5$，取 $\theta_{\max}=5$ 满足条件：
$$
0.5 \times 5 \times 0.25 = 0.625 < 1 \;\checkmark
$$

### 6.3 推荐配置

| 场景 | $\kappa_0$ | $\theta_{\max}$ | 抑制 | 控制 RMS | 鲁棒 LHS |
|------|-----------|----------------|------|---------|----------|
| 🥇 性能优先 | 0.50 | 5 | −28.9 dB | 193 | 0.62 ✅ |
| 🥈 均衡 | 0.20 | 10 | −16.9 dB | 168 | 0.50 ✅ |
| 🥉 低控制/实时 | 0.10 | 10 | −8.8 dB | 126 | 0.25 ✅ |

---

## 7. 各版本实现对比

### 7.1 算法级优化

**O(N²) → O(N) 的关键洞察**：$\Lambda_i$ 构成一条级联链的连续快照。

```
原实现（独立级联）：           优化实现（共享链）：
z → G₁ → G₁ → G₁ → Λ₁[z]     z → G₁ → y₁ = Λ₂₀[z]
z → G₁ → G₁ → Λ₂[z]               → G₁ → y₂ = Λ₁₉[z]
z → G₁ → Λ₂₀[z]                    → G₁ → ...
                                  → G₁ → y₂₀ = Λ₁[z]
(210 次 filter)                (20 次 filter)
```

### 7.2 数值稳定性级联实现

**原因**：$\Lambda_i(s) = \lambda^k/(s+\lambda)^k$ 中 $k$ 最大为 20，展开后多项式系数跨 $500^{20} \approx 10^{54}$ 量级，超过 double 精度。级联一阶节避免系数展开，每节系数维持在 $O(1)$。

### 7.3 版本总结

| 文件 | Λ 算法 | 特点 |
|------|--------|------|
| `Jafari_AdaptiveCC.m` | 单高阶 TF | 原版（含全部 4 个缺陷） |
| `Jafari_AdaptiveCC_Fixed.m` | 独立 IIR_filter 级联 | 修正版，正确但慢 |
| `Jafari_AdaptiveCC_VFast.m` | **共享级联链** | **推荐生产版本** |

---

## 附录：符号表

| 符号 | 含义 | 值 |
|------|------|-----|
| $G_0(s)$ | 标称被控对象 | $\frac{0.5(s-0.2)}{s^2+s+1.25}$ |
| $F(s)$ | 频谱平坦化滤波器 | $\kappa_0\frac{2\alpha^2(s^2+s+1.25)}{(s+\alpha)^2(s+0.2)}$ |
| $\Lambda_i(s)$ | Laguerre 基函数 | $\lambda^{N-i+1}/(s+\lambda)^{N-i+1}$ |
| $N$ | 基函数个数 | 20 |
| $\lambda$ | 基函数极点 | 500 rad/s |
| $\kappa_0$ | F 增益 | 0.5 |
| $\alpha$ | F 带宽参数 | 500 rad/s |
| $T_s$ | 采样周期 | $10^{-4}$ s |
| $\gamma_0$ | 归一化常数 | 1 |
| $P(0)$ | 初始协方差 | $500 \cdot I_N$ |
| $\theta_{\max}$ | 投影界 | 5 |
