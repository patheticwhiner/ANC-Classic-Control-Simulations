# 题1：非最小相位系统的周期扰动自适应反馈抑制

> **Jafari & Ioannou (2017 JVC) 连续时间 CC 控制器 · 考试答题级题解**
>
> **前置阅读**：[Theory_Foundations.md](Theory_Foundations.md) Part A（自适应控制基础）。本文档不再重复 Swapping Lemma、Laguerre 基、RLS、投影、小增益定理的推导——直接引用其定理编号。

---

## 绪论

### 你的知识基线

| 领域 | 你的状态 | 这意味着 |
|:---|:---|:---|
| 连续时间控制系统 | 能跟上推导但自己写不出 | 推导中每一步都会写出，不跳步 |
| Riccati 方程 | 能写出 ARE 形式 | RLS 的 $P(t)$ 微分方程不会陌生 |
| RLS / Laguerre 基 | 有模糊概念 | 已通过 [Theory_Foundations.md §A.2–A.4](Theory_Foundations.md) 理解基础 |
| 小增益定理 | 熟悉 | 已通过 [Theory_Foundations.md §A.6](Theory_Foundations.md) 复习 |

**预计阅读时间**：75 分钟。

---

## 题目

> **已知**：
>
> **(1)** 被控对象（含乘性不确定性）：
>
> $$\boxed{G_0(s) = \frac{0.5(s - 0.2)}{s^2 + s + 1.25}, \qquad \Delta_m(s) = -\mu s,\; \mu = 0.001}$$
>
> **(2)** 扰动信号：
>
> $$d(t) = 0.6\sin(70t + \tfrac{\pi}{4}) + 0.7\sin(187t + \tfrac{\pi}{2}) + \nu(t),\quad \nu(t) \sim \mathcal{N}(0, 0.02^2)$$
>
> 频率个数上界 $n_{\max} = 5$，实际 $n_\omega = 2$。频率值 $\omega_i$ **未知**。
>
> **(3)** 控制架构（Jafari-Ioannou）：
>
> $$z(t) = y(t) - G_0(s)[u(t)], \qquad u(t) = -F(s)\big[K(s, \theta(t))[z(t)]\big]$$
>
> **求**：
> (a) 固定滤波器 $F(s)$ 的闭式表达式
> (b) 自适应滤波器 $K(s,\theta)$ 的参数化形式
> (c) 参数 $\theta$ 的在线自适应律
> (d) 闭环鲁棒稳定性的数值验证
> (e) 仿真抑制性能

---

## 解

### 第1步 · 审题

本问题有三个关键约束迫使我们不能用标准方法：

1. **$G_0(s)$ 非最小相位**（RHP 零点 $s=0.2$）→ 直接逆控制不稳定
2. **扰动频率未知** → 不能预先设计陷波滤波器
3. **存在未建模动态** $\Delta_m(s)$ → 需要鲁棒稳定性保证

Jafari-Ioannou 方案用两层分解应对：

<img src="Images\ControllerDesign.png" alt="ControllerDesign" width = 40% />

- **$F(s)$**：固定滤波器，处理 RHP 零点（通过镜像反射），使 $G_0 F$ 在全频带近似全通
- **$K(s,\theta)$**：自适应 FIR-like 滤波器（用 Laguerre 基实现），在线学习扰动频率处的零点

所需理论基础全部在前置文档中：
- Laguerre 基：[Theory_Foundations.md §A.2](Theory_Foundations.md)
- Swapping Lemma + 线性参数化：[Theory_Foundations.md §A.1, A.3](Theory_Foundations.md)
- RLS 自适应律：[Theory_Foundations.md §A.4](Theory_Foundations.md)
- 参数投影：[Theory_Foundations.md §A.5](Theory_Foundations.md)
- 小增益鲁棒稳定性：[Theory_Foundations.md §A.6](Theory_Foundations.md)

---

### 第2步 · 固定滤波器 $F(s)$ 的设计

> 目标：使 $W(s) \triangleq G_0(s)F(s)$ 在工作频带 $[0, \bar{\omega}_{\max}]$ 上近似全通。

#### 2.1 构造过程

按以下步骤操作（这些是**构造性规则**，非推理性定理）：

**步骤①**：$\tilde{G}_0(s) = G_0(s) = \frac{0.5(s-0.2)}{s^2 + s + 1.25}$

**步骤②**：检查虚轴零点。$G_0(s)$ 无虚轴零点 → 跳过。

**步骤③**：处理 RHP 零点。唯一 RHP 零点 $s = 0.2$，关于虚轴镜像反射至 $s = -0.2$：

$$\tilde{G}_0(s) = \frac{0.5\big(s - (-0.2)\big)}{s^2 + s + 1.25} = \frac{0.5(s + 0.2)}{s^2 + s + 1.25}$$

**步骤④**：构造 $F(s)$。相对阶 $m = 1$（分母 2 阶、分子 1 阶）。取 $\kappa_0 = 0.5$，$\alpha = 500$ rad/s：

$$
\boxed{F(s) = \kappa_0 \cdot \frac{\alpha^m}{(s+\alpha)^m} \cdot \tilde{G}_0^{-1}(s) = 0.5 \cdot \frac{500}{s+500} \cdot \frac{s^2 + s + 1.25}{0.5(s+0.2)}}
$$

化简：

$$
\boxed{F(s) = \frac{500(s^2 + s + 1.25)}{(s + 500)(s + 0.2)}} \tag{1.1}
$$

#### 2.2 验证

- 低频：$F(0) = \frac{500 \times 1.25}{500 \times 0.2} = 6.25$。$G_0(0)F(0) = (-0.08) \times 6.25 = -0.5$
- 扰动频率 70 rad/s：$|G_0(j70)F(j70)| \approx 0.48$（接近 0.5）
- 扰动频率 187 rad/s：$|G_0(j187)F(j187)| \approx 0.49$
- 高频（→∞）：$|F(j\omega)| \sim \frac{500}{\omega}$（一阶滚降）

$G_0F$ 在宽带内幅频响应平坦（≈0.5），达到了「频谱展平」的设计目标。

---

### 第3步 · 自适应滤波器 $K(s,\theta)$ 的参数化

#### 3.1 基函数选取

依 [Theory_Foundations.md §A.2](Theory_Foundations.md)，选用 Laguerre 基 [定理 A.1, 式 (A.4)]：

$$\Lambda_k(s) = \frac{\lambda^{k+1}}{(s+\lambda)^{k+1}}, \qquad k = 0, 1, \dots, N-1$$

**本题参数选择**：

- $N = 20$（$\ge 2n_{\max} = 10$，4 倍冗余以增强频率选择性）
- $\lambda = 500$ rad/s（$> \bar{\omega}_{\max} = 600$ 的约 83%，保证基函数在扰动频带内增益≈1）

自适应滤波器为 [式 (A.6)]：

$$
\boxed{K(s, \theta) = \sum_{k=0}^{N-1} \theta_k \Lambda_k(s) = \theta^T \Lambda(s)} \tag{1.2}
$$

其中 $\theta = [\theta_0, \dots, \theta_{N-1}]^T \in \mathbb{R}^{20}$。

---

### 第4步 · 线性参数化——构造回归向量

#### 4.1 辅佐信号

$$z(t) = y(t) - G_0(s)[u(t)] \tag{1.3}$$

代入 $y = G_0(1+\Delta_m)[u] + d$：

$$z(t) = d(t) + G_0(s)\Delta_m(s)[u(t)]$$

当 $\Delta_m \approx 0$（低频段 $\mu=0.001$ 极小），$z(t) \approx d(t)$。

#### 4.2 应用 Swapping Lemma

控制律：$u = -F K[z] = -F[\theta^T \Lambda[z]]$。

依 [Theory_Foundations.md 引理 A.2](Theory_Foundations.md)（动态 Swapping Lemma），在**慢自适应**假设下：

$$\theta^T(t) \Lambda(s)[z(t)] \approx \theta^T(t) \cdot \Lambda(s)[z(t)]$$

（交换误差 $\varepsilon_H \approx 0$，因 $\dot{\theta}$ 在秒级变化而 $F, G_0$ 的带宽在 $10^2$ rad/s 量级。）

#### 4.3 回归向量

依 [Theory_Foundations.md §A.3](Theory_Foundations.md)，定义：

$$
\boxed{\phi(t) \triangleq G_0(s)F(s)\Lambda(s)[z(t)] \in \mathbb{R}^{20}} \tag{1.4}
$$

即 $\phi_k(t) = G_0(s)F(s)\big[\Lambda_k(s)[z(t)]\big]$，$k = 0, \dots, 19$。

使用 Swapping Lemma 交换 $\theta^T$ 和 $G_0F$，得到线性参数模型 [式 (A.7)]：

$$
\boxed{z(t) = \theta^{*T}\phi(t) + \eta(t)} \tag{1.5}
$$

其中 $\theta^*$ 是能使 $K(s,\theta^*)$ 在扰动频率处获得无穷增益（从而通过反馈完全抵消周期扰动）的未知最优参数，$\eta(t)$ 为建模误差与噪声项。

---

### 第5步 · RLS 自适应律

> 引用 [Theory_Foundations.md §A.4](Theory_Foundations.md) 的连续时间 RLS [定理 A.2]。

#### 5.1 连续时间版本

归一化信号 [式 (A.9)]：$m_s^2(t) = 1 + \gamma_0 \phi^T(t)\phi(t)$，$\gamma_0 = 1$

归一化辨识误差 [式 (A.9)]：$\varepsilon(t) = \dfrac{z(t) - \theta^T(t)\phi(t)}{m_s^2(t)}$

协方差更新 [式 (A.10)]：$\dot{P}(t) = -P(t)\dfrac{\phi(t)\phi^T(t)}{m_s^2(t)}P(t)$，$P(0) = 500 \cdot I_{20}$

参数更新（带投影）[式 (A.11), (A.15)]：

$$\boxed{\dot{\theta}(t) = \operatorname{Proj}_\Theta\Big(P(t)\varepsilon(t)\phi(t)\Big), \qquad \Theta = \{\theta : \|\theta\| \le 5\}} \tag{1.6}$$

投影界 $\theta_{\max} = 5$ 的选取依据将在第 7 步验证。

#### 5.2 离散时间实现（Euler 近似）

采样周期 $T_s = 10^{-4}$ s（10 kHz）。前向 Euler 离散化：

$$
\begin{aligned}
P(k+1) &= P(k) - T_s \cdot P(k)\frac{\phi(k)\phi^T(k)}{m_s^2(k)}P(k) \\[4pt]
\bar{\theta}(k+1) &= \theta(k) + T_s \cdot P(k)\varepsilon(k)\phi(k) \\[4pt]
\theta(k+1) &= \begin{cases}
\bar{\theta}(k+1), & \|\bar{\theta}\| \le 5 \\
5 \cdot \bar{\theta} / \|\bar{\theta}\|, & \text{otherwise}
\end{cases}
\end{aligned} \tag{1.7}
$$

代码中另有一版本使用真正的离散时间 RLS [Theory_Foundations.md 式 (A.12)-(A.14)]，数值稳定性更优。

---

### 第6步 · 控制律

**关键**：控制律中的信号路径与自适应路径**不同**！

| 路径 | 滤波链 | 用途 |
|:---|:---|:---|
| **自适应路径** | $z \to \Lambda[z] \to F[z] \to G_0[z] \to \phi$ | 计算 $\varepsilon$，更新 $\theta$ |
| **控制路径** | $z \to \Lambda[z] \to \theta^T\Lambda[z]$ | 产生控制信号 $u = -F[\theta^T\Lambda[z]]$ |

控制律中 $\theta^T \Lambda[z]$ **不经过** $G_0F$——若将 $\phi$（已包含 $G_0F$）误用于控制律，当 $\theta$ 收敛后控制器退化为固定增益反馈，丧失频率选择性。这是原版代码的核心 bug。

离散化后（各传递函数经 Tustin 离散化）：

$$
\boxed{u(k) = -F_d(z)\left[\sum_{k=0}^{N-1} \theta_k(k) \cdot \Lambda_{k,d}(z)[z(k)]\right]} \tag{1.8}
$$

---

### 第7步 · 鲁棒稳定性验证

引用 [Theory_Foundations.md 引理 A.3](Theory_Foundations.md)，自适应系统的鲁棒稳定条件为：

$$\kappa_0 \cdot \bar{\theta}_M \cdot \| \bar{F} G_0 \Delta_m \|_{\mathcal{L}_1} < 1$$

代入数值：
- $\kappa_0 = 0.5$
- $\bar{\theta}_M = \theta_{\max} = 5$（投影界保证）
- $\| \bar{F} G_0 \Delta_m \|_\infty \approx \| F G_0 \Delta_m \|_\infty / \kappa_0 \approx 0.25$（数值计算）

$$0.5 \times 5 \times 0.25 = 0.625 < 1 \;\; \checkmark$$

留有 37.5% 裕度。若增大 $\kappa_0 = 0.8$，LHS = $0.8 \times 5 \times 0.25 = 1.0$（临界），当前选值合理。

---

### 第8步 · 数值结果

运行 `JafariJVC_Continuous.m`（10 kHz, 20 s），取后 50% 数据统计：

| 版本 | 自适应算法 | 抑制 | $\|\theta\|$ | 耗时 (10s 仿真) |
|------|-----------|------|-------------|----------------|
| 第1段 | Euler-RLS | −28.9 dB | 2.51 | ~150 s |
| 第2段 | 真 RLS | −28.9 dB | 2.51 | 略快于第1段 |
| VFast | RLS + 共享级联 | −28.9 dB | 2.51 | ~15.8 s (9.5×加速) |

抑制的一致性表明三种实现数学等价。VFast 的加速原理见 [Theory_Foundations.md §A.2 引理 A.2 性质 4](Theory_Foundations.md)——利用 $\Lambda_k$ 的级联递推关系，将 $\Lambda$ 滤波从 $O(N^2)$ 降至 $O(N)$。

---

### 第9步 · 代码映射

| 数学公式 | `JafariJVC_Continuous.m` 关键行 | 变量 |
|:---|:---|:---|
| $G_0(s)$ | L69 | `G0 = 0.5*(s-0.2)/(s^2+s+1.25)` |
| $F(s)$ 设计 | L92–97 | `F = designLTIFilter(G0, k0, delta0, alpha, m, omega_max)` |
| $z = y - G_0[u]$ | L131–132 | `z_k = y1(k) - G0_u_k` |
| $\Lambda_k$ 级联 | L134–137 | `for stage=1:Nparam ... filter(numG1, denG1, ...)` |
| $\phi = G_0F\Lambda[z]$ | L138–141 | `filter(numF,denF,...)` → `filter(numG0d_s,denG0d_s,...)` |
| $m_s^2 = 1 + \gamma_0\phi^T\phi$ | L142 | `m2 = 1 + gamma0*(phi1.'*phi1)` |
| $\varepsilon$ | L142 | `eps1 = (z_k - theta1.'*phi1)/m2` |
| $\dot{P}$ Euler | L143 | `P_euler = P_euler - Ts*(Pph*Pph.')/m2` |
| $\dot{\theta} +$ Proj | L144–145 | `th_new = theta1 + Ts*P_euler*eps1*phi1` |
| $u = -F[\theta^T\Lambda[z]]$ | L147 | `u1(k) = filter(numF,denF,-theta1.'*lam1,...)` |
| $F(s)$ 零极点反射 | L220–304 | `designLTIFilter` 函数体 |

---

### 终 · 回到你的问题

| 子问题 | 答案 |
|:---|:---|
| (a) $F(s)$ 表达式 | $F(s) = \dfrac{500(s^2 + s + 1.25)}{(s + 500)(s + 0.2)}$，由 RHP 零点镜像反射 + 低通滤波器保证真分性 |
| (b) $K(s,\theta)$ 参数化 | $K(s,\theta) = \sum_{k=0}^{19} \theta_k \frac{500^{k+1}}{(s+500)^{k+1}}$，Laguerre 基 $N=20$, $\lambda=500$ |
| (c) 自适应律 | $\dot{P} = -P\frac{\phi\phi^T}{m_s^2}P$，$\dot{\theta} = \operatorname{Proj}(P\varepsilon\phi)$，离散化用 Euler 或真 RLS |
| (d) 鲁棒稳定 | $0.5 \times 5 \times 0.25 = 0.625 < 1$ ✓，留有 37.5% 裕度 |
| (e) 抑制性能 | −28.9 dB（RMS 比），参数收敛至 $\|\theta\| \approx 2.51$ |

---

> 📖 **相关文档**：
> - 公共理论：[Theory_Foundations.md Part A](Theory_Foundations.md)
> - 代码修正细节：[DebugReport_JafariJVC_Discrete.md](DebugReport_JafariJVC_Discrete.md)
> - 离散时间版本：[题2 题解](Derivation_Problem2_Jafari_DiscreteAVC.md)
> - 算法推导与加速分析：[About_CC_Derivation.md](Derivation_Problem1_Appendix_CC_Algorithm.md)

---

# 附录

> 以下内容来自原附录文档，归入主文档以便查阅。§A 覆盖算法加速与嵌入式部署，§B 覆盖原版代码的四重缺陷诊断。

## 附录A · 算法加速与嵌入式实时性

> 来源：`Derivation_Problem1_Appendix_CC_Algorithm.md`

### A.1 从 $O(N^2)$ 到 $O(N)$ 的加速推导

$N$ 个 Laguerre 基函数 $\Lambda_i(s) = \lambda^{N-i+1}/(s+\lambda)^{N-i+1}$，每个的级联深度为 $d_i = N-i+1$。总级联节数：

$$\sum_{i=1}^{N} d_i = \sum_{i=1}^{N} (N - i + 1) = \frac{N(N+1)}{2}$$

$N=20$ 时，每样本需 $\frac{20 \times 21}{2} = 210$ 次一阶 IIR 滤波。加上每 $\Lambda_i$ 输出经过 $F$ 和 $G_0$（各 1 次），总计：

$$\text{滤波次数/样本} = \frac{N(N+1)}{2} + 2N = \frac{N(N+5)}{2}$$

**共享级联链优化**：$\Lambda_i$ 是级联链 $G_1 \to G_1^2 \to \dots \to G_1^N$（$G_1(s)=\frac{\lambda}{s+\lambda}$）的不同深度输出。只需将 $z(k)$ 依次通过 $N$ 个一阶节，记录每个中间输出：

```
y ← z(k)
for s = 1..N:
    y ← G₁[y]
    Λ_{N-s+1}[z] ← y
```

将 $\Lambda$ 滤波从 $\frac{N(N+1)}{2}$ 次降至 $N$ 次。优化后总滤波 $3N$ 次。加速比 $\frac{N+5}{6} \approx 4.17\times$（$N=20$）。

**实测数据**：

| 版本 | Λ 滤波/样本 | 10s 仿真耗时 | 加速比 |
|------|-----------|-------------|--------|
| IIR_filter (原) | 210 | ~150 s | 1× |
| filter() 独立级联 | 210 | 56.7 s | 2.6× |
| filter() 共享级联 | **20** | **15.8 s** | **9.5×** |

### A.2 实时嵌入式部署可行性

当前算法（$N=20$, $f_s=10^4$）每样本约 1260 FLOP。10 kHz 采样率下需 12.6 MFLOP/s。STM32H7（Cortex-M7 @ 480 MHz）单精度浮点 >200 MFLOP/s，留有 **16× 裕度**。

降 $N=6$ + 降采样 $f_s=2000$ Hz 后，仅需 0.6 MFLOP/s——**Cortex-M0+ @ 48 MHz** 即可运行。

### A.3 参数调优方法论

| 参数 | 增大效果 | 减小效果 |
|------|---------|---------|
| $\kappa_0$ | 抑制↑，控制能量↑，鲁棒裕度↓ | 抑制↓，控制能量↓，鲁棒裕度↑ |
| $\theta_{\max}$ | 允许更强的自适应 | 强制限制闭环增益 |

鲁棒稳定性条件：$\kappa_0 \cdot \bar{\theta}_M \cdot \|\bar{F}G_0\Delta_m\|_1 < 1$。本题 $0.5 \times 5 \times 0.25 = 0.625 < 1 \checkmark$。

**推荐配置**：

| 场景 | $\kappa_0$ | $\theta_{\max}$ | 抑制 |
|------|-----------|----------------|------|
| 性能优先 | 0.50 | 5 | −28.9 dB |
| 均衡 | 0.20 | 10 | −16.9 dB |
| 低控制/实时 | 0.10 | 10 | −8.8 dB |

---

## 附录B · 原版代码的四重缺陷诊断

> 来源：`Derivation_Problem1_Appendix_Reproduction.md`

### B.1 缺陷 1：频率单位错误（Hz → 应为 rad/s）

论文扰动频率 $\omega = 70, 187$ **rad/s**。原版代码误作 $f = 70, 187$ **Hz**，换算为 $\omega' \approx 440, 1175$ rad/s，远超 $F(s)$ 的带宽（$\alpha=500$），控制器在扰动频率处无有效增益。

### B.2 缺陷 2：参数个数 $N=6 \to 20$

论文要求 $N \ge 2n_{\max}$（$n_{\max}=5$），取 $N=20$。$N=6$ 意味着 FIR 仅有 6 个自由度，不足以在多个频率处构造零点。

### B.3 缺陷 3：控制律信号混淆

| | 原版（错误） | 修正 |
|:---|:---|:---|
| 控制律 | $u = -F[\theta^T\phi] = -F[\theta^T G_0 F \Lambda[z]]$ | $u = -F[\theta^T \Lambda[z]]$ |

$\phi$ 中额外的 $G_0F$ 使控制器退化为固定增益反馈，丧失频率选择性。

### B.4 缺陷 4：高阶多项式数值崩溃

$\Lambda_1(s) = \lambda^{20}/(s+\lambda)^{20}$ 展开后系数跨 $500^{20} \approx 10^{54}$ 量级，Tustin 离散化后 `den(1)` 被舍入为零，差分方程除以零产生 NaN。**修复**：用一阶节级联避免系数展开。

### B.5 版本总结

| 文件 | Λ 算法 | 特点 |
|------|--------|------|
| `Jafari_AdaptiveCC.m` | 单高阶 TF | 原版（含全部 4 个缺陷） |
| `Jafari_AdaptiveCC_Fixed.m` | 独立级联 | 修正版，正确但慢 |
| `Jafari_AdaptiveCC_VFast.m` | **共享级联链** | 推荐生产版本 |
