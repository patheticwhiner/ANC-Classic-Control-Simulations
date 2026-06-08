# 理论基础：自适应鲁棒控制与 H∞ 综合

> **本文档是 `demo3_Robust` 下所有题解的公共理论地基。各题解直接引用本文档中的定理编号，不再重复推导。**

---

## 绪论

### 本文档的定位

`demo3_Robust` 目录下的三组仿真对应三个具体控制问题，但共享共同的理论工具箱：

| 工具箱 | 内容 | 被哪些题引用 |
|:---|:---|:---|
| **Part A** 自适应控制基础 | Swapping Lemma、Laguerre 基、线性参数化、RLS、投影、小增益定理 | 题1（连续 CC）、题2（离散 AVC） |
| **Part B** H∞ 鲁棒控制基础 | H∞ 范数、广义被控对象、混合灵敏度、DGKF 方法、γ-迭代 | 题2（H∞ 插值约束）、题3（H∞ 混合灵敏度） |

### 依赖关系

```
Part A（自适应控制）
  ├── A.1 Swapping Lemma          ← 题1§3, 题2§4 引用
  ├── A.2 Laguerre 基函数          ← 题1§2 引用
  ├── A.3 线性参数化               ← 题1§3, 题2§4 引用
  ├── A.4 RLS 自适应律             ← 题1§4, 题2§5 引用
  ├── A.5 参数投影算法             ← 题1§4, 题2§5 引用
  └── A.6 小增益定理与鲁棒稳定性    ← 题1§5, 题2§6 引用

Part B（H∞ 鲁棒控制）
  ├── B.1 H∞ 范数定义              ← 题2§3, 题3§2 引用
  ├── B.2 广义被控对象与混合灵敏度  ← 题3§3 引用
  ├── B.3 DGKF 方法                ← 题3§4 引用
  └── B.4 γ-迭代与中心控制器        ← 题3§4 引用
```

---

# Part A：自适应控制理论基础

---

## A.1 Swapping Lemma（交换引理）

> 🎯 **自查**：你是否熟悉「线性时不变算子的叠加原理」？若是，本节只需 5 分钟（核心公式 (A.3) 即可）。

### A.1.1 问题的提出

在自适应控制中，常见的困难是：参数向量 $\theta(t)$（时变）与信号 $x(t)$ 通过一个**动态系统**（传递函数）耦合在一起，你需要将 $\theta$ 从传递函数内部「提出」到外面，才能得到标准的线性参数化形式 $\theta^T \phi$。

### A.1.2 常参数情况

**引理 A.1（静态 Swapping Lemma）** 设 $H(s)$ 是任意稳定的线性时不变传递函数，$x(t) \in \mathbb{R}^N$ 为信号向量，$\theta \in \mathbb{R}^N$ 为**常向量**。则：

$$
\boxed{H(s)[\theta^T x(t)] = \theta^T H(s)[x(t)]} \tag{A.1}
$$

其中 $H(s)[x(t)]$ 表示对 $x(t)$ 的每个分量分别施加 $H(s)$。

*证明*：传递函数 $H(s)$ 是线性算子。令 $h(t) = \mathcal{L}^{-1}\{H(s)\}$ 为其冲激响应，则：

$$
\begin{aligned}
H(s)[\theta^T x](t) &= \int_0^t h(t-\tau) \left(\sum_{i=1}^N \theta_i x_i(\tau)\right) d\tau \\
&= \sum_{i=1}^N \theta_i \int_0^t h(t-\tau) x_i(\tau) d\tau \\
&= \sum_{i=1}^N \theta_i \cdot H(s)[x_i](t) = \theta^T H(s)[x](t)
\end{aligned}
$$

第二步将常数 $\theta_i$ 提出积分号，第三步利用 $H(s)$ 对每个分量的线性性。$\square$

**注释 A.1**：式 (A.1) 的本质是：**线性算子与常系数的内积可以交换顺序**。这是自适应控制中几乎所有线性参数化技巧的共同根源。

### A.1.3 时变参数情况

**引理 A.2（动态 Swapping Lemma）** 当 $\theta(t)$ 是时变信号时，存在交换误差项 $\varepsilon_H(t)$：

$$
H(s)[\theta^T(t) x(t)] = \theta^T(t) H(s)[x(t)] + \varepsilon_H(t) \tag{A.2}
$$

其中 $\varepsilon_H(t)$ 的显式表达式为：

$$
\varepsilon_H(t) = \int_0^t h(t-\tau)\big[\theta^T(\tau) - \theta^T(t)\big] x(\tau) d\tau \tag{A.3}
$$

*证明*：

$$
\begin{aligned}
H(s)[\theta^T x](t) &= \int_0^t h(t-\tau) \theta^T(\tau) x(\tau) d\tau \\
&= \int_0^t h(t-\tau) \theta^T(t) x(\tau) d\tau + \int_0^t h(t-\tau)[\theta^T(\tau) - \theta^T(t)] x(\tau) d\tau \\
&= \theta^T(t) H(s)[x](t) + \varepsilon_H(t)
\end{aligned}
$$

第一步到第二步：在积分内加减 $\theta^T(t)$；第二步到第三步：利用常参数情形的引理 A.1。$\square$

**推论 A.1（慢自适应近似）** 若 $\|\dot{\theta}(t)\|$ 在所有时刻均足够小（即参数变化比 $H(s)$ 的动态慢得多），则 $\varepsilon_H(t) \approx 0$，常参数情形的交换成立。

*论证*：由中值定理，$\theta(\tau) - \theta(t) = -\dot{\theta}(\xi)(t-\tau)$ 对某 $\xi \in (\tau, t)$。代入 (A.3)：

$$|\varepsilon_H(t)| \le \|\dot{\theta}\|_\infty \cdot \int_0^t |h(t-\tau)|(t-\tau) \cdot \|x(\tau)\| d\tau$$

当 $\|\dot{\theta}\|_\infty$ 很小且 $H(s)$ 稳定（$h(t)$ 指数衰减）时，此积分有界且很小。$\square$

---

## A.2 Laguerre 基函数

> 🎯 **自查**：你是否理解「函数空间的基」这一概念（例如 Fourier 级数是 $\{e^{jn\omega t}\}$ 的线性组合）？若是，Laguerre 基不过是换了一组基函数，本节约 15 分钟。

### A.2.1 定义与直觉

**定义 A.1（Laguerre 基函数）** 给定极点参数 $\lambda > 0$，第 $k$ 个 Laguerre 基函数（$k = 0, 1, 2, \dots$）定义为：

$$
\boxed{\Lambda_k(s) = \frac{\lambda^{k+1}}{(s + \lambda)^{k+1}} = \underbrace{\frac{\lambda}{s + \lambda} \cdot \frac{\lambda}{s + \lambda} \cdot \cdots \cdot \frac{\lambda}{s + \lambda}}_{k+1 \text{ 个相同的一阶低通滤波器级联}}} \tag{A.4}
$$

**直觉 1（级联低通的物理画面）**：$\Lambda_k(s)$ 是一串完全相同的一阶低通滤波器 $\frac{\lambda}{s+\lambda}$ 的级联。这就像 FIR 滤波器用一串延迟 $z^{-1}$ 来近似任意频率响应一样——Laguerre 基用一串一阶低通滤波器来近似任意传递函数。

| 对比维度 | FIR 基 | Laguerre 基 |
|:---|:---|:---|
| 基函数 | $z^{-k}$（延迟线） | $\frac{\lambda^{k+1}}{(s+\lambda)^{k+1}}$（低通级联） |
| 参数域 | 离散时间 | 连续时间（可离散化） |
| 物理直觉 | 时间平移 | 时间+频率平滑 |
| 近似能力 | 任意 FIR 频率响应 | 任意稳定严格真传递函数 |

**直觉 2（为什么用 Laguerre 而非其他基？——一个摄影类比）**：想象你要拍一张照片来描述一个传递函数的频率响应。

- **FIR 基（延迟线）** 就像用一架**固定焦距的相机**在时间轴上每隔 $T_s$ 拍一张——你在时间上获得了精确的分辨率，但在频率上只能用更多的点数（更高的 N）来逐渐逼近。
- **Laguerre 基** 就像用一架**可以同时调整焦距和曝光时间的相机**——$\lambda$ 是焦距（决定了你关注哪个频段），$k$ 是多重曝光次数。低 $k$ 的基函数拍"广角"照片（宽频带、低分辨率），高 $k$ 的基函数拍"长焦"照片（窄频带、高分辨率）。把不同焦距的照片叠在一起，你就得到了全貌。

这种"多分辨率"特性是 Laguerre 基比 FIR 基在连续时间域更高效的根本原因：你不需要 $N = 100$ 个延迟线来近似一个二阶系统——通常 $N = 6 \sim 20$ 个 Laguerre 基函数就能提供足够的频率选择性。

**直觉 3（级联递推——$O(N^2) \to O(N)$ 加速的物理图景）**：性质 4 的级联递推 $\Lambda_k = \frac{\lambda}{s+\lambda}\Lambda_{k-1}$ 意味着你可以像**流水线**一样计算所有基函数输出：

```
z(t) ─→ [G₁] ─→ [G₁] ─→ [G₁] ─→ ... ─→ [G₁]
         ↓       ↓       ↓               ↓
        Λₙ₋₁    Λₙ₋₂    Λₙ₋₃           Λ₀
```

你把信号 $z(t)$ 推进一条 $N$ 个相同一阶低通滤波器串联的流水线，每经过一级就"抽出"一个样本——这就是那个级的 $\Lambda_k[z]$。$N$ 个基函数的输出只需 $N$ 次滤波操作，而非 $N(N+1)/2$ 次。这就是代码中 `Jafari_AdaptiveCC_VFast.m` 比原始版本快 9.5 倍的根本原因。

**引理 A.2（基本性质）** Laguerre 基函数满足：

1. **DC 增益**：对所有 $k$，$\Lambda_k(0) = 1$（每个基函数在稳态下通过所有信号，不衰减）
2. **−3 dB 频率**：$\Lambda_k(j\omega)$ 的 −3 dB 截止频率为 $\omega_c = \lambda\sqrt{2^{1/(k+1)} - 1}$。$k$ 越大，截止频率越低
3. **严格真性**：$\Lambda_k(s)$ 的相对阶为 $k+1$
4. **级联递推**：$\Lambda_{k}(s) = \frac{\lambda}{s+\lambda} \Lambda_{k-1}(s)$（其中 $\Lambda_{-1}(s) \triangleq 1$），这一性质是实现 $O(N)$ 加速的关键

### A.2.2 完备性

**定理 A.1（Laguerre 基在 $\mathcal{H}_2$ 中的完备性）** 对于任意稳定严格真（strictly proper）的传递函数 $H(s) \in \mathcal{H}_2$，存在唯一的系数序列 $\{\theta_k\}_{k=0}^\infty$ 使得：

$$
H(s) = \sum_{k=0}^{\infty} \theta_k \Lambda_k(s) \tag{A.5}
$$

在 $\mathcal{H}_2$ 范数 $\|H\|_2^2 = \frac{1}{2\pi}\int_{-\infty}^{\infty} |H(j\omega)|^2 d\omega$ 的意义下收敛。并且 $\{\Lambda_k(s)\}_{k=0}^\infty$ 构成 $\mathcal{H}_2$ 的一组**正交基**：

$$
\langle \Lambda_k, \Lambda_\ell \rangle_{\mathcal{H}_2} = \frac{1}{2\pi}\int_{-\infty}^{\infty} \Lambda_k(j\omega) \overline{\Lambda_\ell(j\omega)} d\omega = 0 \quad (k \ne \ell)
$$

*证明梗概*：通过变量代换 $w = \frac{s-\lambda}{s+\lambda}$（从 $s$ 平面到 $w$ 平面的 Möbius 变换），将 $\Lambda_k(s)$ 映射为 $w^{k+1}$，后者恰为 $|w|<1$ 单位圆盘上 Hardy 空间的标准正交基 $\{1, w, w^2, \dots\}$。由 Hardy 空间理论知此组基是完备的。详细证明参见标准泛函分析教材（如 Sz.-Nagy & Foias, *Harmonic Analysis of Operators on Hilbert Space*）。$\square$

### A.2.3 截断近似

在工程实践中，级数 (A.5) 需要在有限项数 $N$ 处截断：

$$
\boxed{K(s, \theta) = \sum_{k=0}^{N-1} \theta_k \Lambda_k(s) = \theta^T \Lambda(s)} \tag{A.6}
$$

其中 $\theta = [\theta_0, \dots, \theta_{N-1}]^T \in \mathbb{R}^N$，$\Lambda(s) = [\Lambda_0(s), \dots, \Lambda_{N-1}(s)]^T$。

**参数 $N$ 的选择原则**：若要精确抑制 $n_{\max}$ 个不同频率的正弦扰动，每个频率需要一对共轭零点（2 个自由度），因此 $N \ge 2n_{\max}$。实践中 $N = 2n_{\max} \sim 4n_{\max}$ 提供冗余。

**$\lambda$ 的选择原则**：$\lambda$ 应大于最大可能的扰动频率 $\bar{\omega}_{\max}$，以确保所有 Laguerre 基函数在所需频带内保持单位增益附近。典型取值 $\lambda \in [1.5\bar{\omega}_{\max}, 5\bar{\omega}_{\max}]$。

---

## A.3 线性参数化

> 🎯 **自查**：你知道「线性回归」吗？$y = \beta_0 + \beta_1 x_1 + \cdots + \beta_p x_p$ 就是线性参数化。自适应控制要做的就是把非线性控制律重写成这种形式。本节约 10 分钟。

### A.3.1 什么是线性参数化

**定义 A.2（线性参数模型，Linear Parametric Model, LPM）** 考虑未知参数向量 $\theta^* \in \mathbb{R}^N$。若测量信号 $z(t)$ 可以写为：

$$
\boxed{z(t) = \theta^{*T} \phi(t) + \eta(t)} \tag{A.7}
$$

其中 $\phi(t) \in \mathbb{R}^N$ 是**已知的**回归向量（由可测信号经过滤波得到），$\eta(t)$ 是建模误差/噪声项，则称 $z(t)$ 被**线性参数化**。

**为什么线性参数化如此重要？** 因为 (A.7) 将未知参数 $\theta^*$ 与已知信号 $\phi$ 的关系写成**内积**形式。一旦有了它，所有标准参数估计算法（梯度下降、RLS、LMS）都可以直接套用。

### A.3.2 构造线性参数化的一般策略

从形如 $z = H(s, \theta)[u] + d$ 的非线性模型出发：

1. **剥离确定性部分**：通过类似 IMC 的结构，构造辅佐信号 $z' = z - \text{(已知模型贡献)}$
2. **应用 Swapping Lemma**（引理 A.1/A.2）：将 $\theta$ 从传递函数内部交换到外部
3. **定义回归向量** $\phi$：将 $\theta$ 前面的所有滤波操作打包为一个向量信号

### A.3.3 辨识误差与归一化

**定义 A.3（辨识误差）** 给定参数估计 $\theta(t)$，定义**先验辨识误差**（a priori identification error）：

$$
\boxed{\hat{z}(t) = \theta^T(t)\phi(t), \qquad e(t) = z(t) - \hat{z}(t)} \tag{A.8}
$$

**定义 A.4（归一化辨识误差）** 为增强对大幅值信号的鲁棒性，引入归一化：

$$
\boxed{m_s^2(t) = 1 + \gamma_0 \phi^T(t)\phi(t), \qquad \varepsilon(t) = \frac{e(t)}{m_s^2(t)}} \tag{A.9}
$$

其中 $\gamma_0 > 0$ 是归一化增益。式 (A.9) 中的 "1" 防止 $\phi = 0$ 时除以零。

**注释 A.2**：归一化是自适应控制鲁棒化的关键步骤。它保证了即使 $\phi(t)$ 无界（如系统不稳定），$\varepsilon(t)$ 仍然有界。在无归一化的情况下，自适应律可能因 $\phi$ 的瞬态大值而发散（称为**参数漂移** / parameter drift）。

---

## A.4 递归最小二乘（RLS）自适应律

> 🎯 **自查**：你是否理解最小二乘的基本思想——$\hat{\theta} = \arg\min \sum (z_i - \theta^T \phi_i)^2$？RLS 是它的递归（在线）版本。本节约 25 分钟。

### A.4.1 从批量最小二乘到递归最小二乘

**问题**：给定直到时刻 $t$ 的数据 $\{z(\tau), \phi(\tau)\}_{\tau=0}^t$，求 $\theta$ 使加权累积平方误差最小：

$$
J_t(\theta) = \int_0^t e^{-\beta(t-\tau)} \frac{(z(\tau) - \theta^T\phi(\tau))^2}{m_s^2(\tau)} d\tau + e^{-\beta t} (\theta - \theta_0)^T P_0^{-1} (\theta - \theta_0)
$$

其中 $\beta \ge 0$ 是指数遗忘因子（$\beta = 0$ 为无遗忘），$P_0 > 0$ 为初始参数协方差矩阵。

**定理 A.2（连续时间 RLS）** 令 $\beta \to 0$（无遗忘），最小化 $J_t(\theta)$ 的参数估计 $\theta(t)$ 满足以下微分方程：

$$
\boxed{\dot{P}(t) = -P(t) \frac{\phi(t)\phi^T(t)}{m_s^2(t)} P(t), \qquad P(0) = P_0 = p_0 I_N} \tag{A.10}
$$

$$
\boxed{\dot{\theta}(t) = P(t)\varepsilon(t)\phi(t)} \tag{A.11}
$$

其中 $\varepsilon(t) = (z(t) - \theta^T(t)\phi(t))/m_s^2(t)$ 为归一化辨识误差。

*推导*：令 $P^{-1}(t) = \int_0^t \frac{\phi(\tau)\phi^T(\tau)}{m_s^2(\tau)} d\tau + P_0^{-1}$（即信息矩阵，information matrix）。则最优解为 $\theta(t) = P(t) \int_0^t \frac{\phi(\tau)z(\tau)}{m_s^2(\tau)} d\tau + P(t)P_0^{-1}\theta_0$。

对 $P^{-1}(t)$ 求导：$\frac{d}{dt}P^{-1}(t) = \frac{\phi(t)\phi^T(t)}{m_s^2(t)}$。

利用恒等式 $\frac{d}{dt}(P^{-1}) = -P^{-1}\dot{P}P^{-1}$，得：

$$-P^{-1}\dot{P}P^{-1} = \frac{\phi\phi^T}{m_s^2}$$

左乘 $P$、右乘 $P$：$\dot{P} = -P\frac{\phi\phi^T}{m_s^2}P$。此即 (A.10)。

对 $\theta(t)$ 的表达式求导（利用 $\dot{P}$ 的表达式和 Leibniz 积分法则）即得 (A.11)。$\square$

**注释 A.3（$P(t)$ 的单调性）**：由 (A.10)，$\dot{P} \le 0$（半负定），故 $P(t)$ 随时间单调递减。这与 Kalman 滤波器中协方差矩阵的行为一致：随着更多数据被纳入，参数估计的不确定性持续下降。

**注释 A.4（RLS vs. 梯度下降）**：比较：

| 方法 | 更新律 | 步长 |
|:---|:---|:---|
| 梯度下降 | $\dot{\theta} = -\eta \nabla_\theta(\frac{1}{2}\varepsilon^2) = \eta \varepsilon \phi$ | 标量 $\eta$ |
| RLS | $\dot{\theta} = P \varepsilon \phi$ | 矩阵 $P$（**曲率自适应**） |

RLS 的步长 $P$ 自动沿「信息量大」的方向放大更新、沿「信息量小」的方向缩小更新。这就像牛顿法之于梯度下降：前者使用了二阶导数（曲率）信息。

**直觉 4（$P$ 矩阵的几何图景——「各向异性学习率」）**：标量学习率 $\eta$ 就像一个**球**——你在参数空间的每个方向上迈出相同的步长。但真实的问题从来不是各向同性的。

想象你在山脊上行走寻找山谷的最低点。山脊的方向（$\phi$ 的方向）你已经看得很清楚——这就是「信息量大」的方向，你应该大步走。垂直于山脊的方向你很模糊——这就是「信息量小」的方向，你应该小步走以免掉进沟里。

$P$ 矩阵就是一个**椭球**——它在信息量大的方向拉长（大步长），在信息量小的方向压扁（小步长）。并且这个椭球会随着你探索而不断变形：每当你得到一个新的 $\phi(k)$，你就「削掉」$P$ 在 $\phi$ 方向上的一部分不确定度（$\dot{P} = -P\frac{\phi\phi^T}{m_s^2}P$），那个方向的步长就缩小。

**直觉 5（$P$ 与 Kalman 滤波器的关系——你不是第一次见到它）**：如果你学过 Kalman 滤波，$P$ 就是**参数估计误差的协方差**。Kalman 滤波器中有：

$$\hat{x}_{k|k} = \hat{x}_{k|k-1} + K_k(y_k - C\hat{x}_{k|k-1})$$

RLS 中有：

$$\theta(k) = \theta(k-1) + K(k)(z(k) - \phi^T(k)\theta(k-1))$$

这两条公式**完全相同**！唯一的区别是 Kalman 中的「状态预测」在 RLS 中变成了「参数不变」（$\theta$ 是常数，所以预测就是上一时刻的值）。Woodbury 恒等式 (A.12)-(A.13) 不过是 Kalman 协方差更新在没有过程噪声时的退化形式。这意味着你对 Kalman 滤波器的所有直觉——「$P$ 大=不确定性大=学得快」「随着数据增多 $P$ 单调减小」「$\phi$ 的方向贡献信息最多」——在 RLS 中全部适用。

### A.4.2 离散时间 RLS

**定理 A.3（离散时间 RLS）** 对于离散时间线性参数模型 $z(k) = \theta^{*T}\phi(k) + \eta(k)$，最小化 $\sum_{i=1}^k \lambda^{k-i} (z(i) - \theta^T\phi(i))^2$（$\lambda \in (0,1]$ 为遗忘因子）的递归解为：

$$
\boxed{K(k) = \frac{P(k-1)\phi(k)}{\lambda + \phi^T(k)P(k-1)\phi(k)}} \tag{A.12}
$$

$$
\boxed{P(k) = \frac{1}{\lambda}\left[I - K(k)\phi^T(k)\right] P(k-1)} \tag{A.13}
$$

$$
\boxed{\theta(k) = \theta(k-1) + K(k)\left(z(k) - \theta^T(k-1)\phi(k)\right)} \tag{A.14}
$$

其中 $P(k)$ 为协方差矩阵，$K(k)$ 为 Kalman 增益向量。

*推导*：这是标准 Kalman 滤波的确定性等价形式。从矩阵求逆引理（Woodbury identity）出发：

$$(A + BC)^{-1} = A^{-1} - A^{-1}B(I + CA^{-1}B)^{-1}CA^{-1}$$

应用到 $P^{-1}(k) = \lambda P^{-1}(k-1) + \phi(k)\phi^T(k)$，即得 (A.12)-(A.13)。参数更新 (A.14) 是加权最小二乘的正规方程。$\square$

**注释 A.5（归一化在离散 RLS 中的作用）**：当 $m_s^2(k)$ 出现在分母中时，(A.12) 的分母变为 $\lambda m_s^2(k) + \phi^T P \phi$。这确保了即使 $\phi(k)$ 很大，$K(k)$ 也不会失控。

---

## A.5 参数投影算法

> 🎯 **自查**：你知道约束优化的基本概念吗（可行域、梯度投影）？投影算法就是在线参数估计的约束处理。

### A.5.1 为什么需要投影

无约束 RLS 的 $\theta(t)$ 可能漂移到远离真实值的区域，原因包括：
- 持续激励（PE）条件不满足时，协方差矩阵 $P(t)$ 在某些方向不衰减
- 噪声和未建模动态可能将 $\theta$ 推离有效区域

投影将 $\theta$ 始终约束在一个先验已知的紧凸集内，保证参数有界。

### A.5.2 投影算子的定义

**定义 A.5（投影算子）** 令 $\Theta \subset \mathbb{R}^N$ 为紧凸集，$\delta\Theta$ 为其边界，$\Theta^o$ 为其内部。对于从 $\theta \in \Theta$ 出发的搜索方向 $\tau \in \mathbb{R}^N$，投影算子定义为：

$$
\boxed{\operatorname{Proj}_\Theta(\theta, \tau) = \begin{cases}
\tau, & \text{if } \theta \in \Theta^o \\[4pt]
\tau, & \text{if } \theta \in \delta\Theta \text{ and } \tau^T \nabla g(\theta) \le 0 \\[6pt]
\tau - \Gamma \dfrac{\nabla g \nabla g^T}{\nabla g^T \Gamma \nabla g} \tau, & \text{otherwise}
\end{cases}} \tag{A.15}
$$

其中 $g(\theta) \le 0$ 是定义凸集 $\Theta = \{\theta : g(\theta) \le 0\}$ 的约束函数，$\Gamma$ 是一个正定矩阵（通常取 $\Gamma = P$ 或 $\Gamma = I$）。

**直觉 6（投影的几何图景——「弹珠碰到碗壁」）**：用二维平面（$N = 2$）来想。$\Theta$ 是一个圆盘（$\theta_1^2 + \theta_2^2 \le \theta_{\max}^2$）。你在圆盘内放一颗弹珠（当前参数 $\theta$），你要沿着方向 $\tau$ 推动它。

- **情况 1**：弹珠在圆盘**内部**（$\theta \in \Theta^o$）。你怎么推都行——$\operatorname{Proj} = \tau$（直接沿 $\tau$ 走）。

- **情况 2**：弹珠在圆盘边界上，但你想把它**往圆盘里面推**（$\tau^T \nabla g \le 0$，即 $\tau^T\theta \le 0$，意味着 $\tau$ 与半径方向 $\theta$ 的夹角 $\ge 90^\circ$）。没问题——$\operatorname{Proj} = \tau$，弹珠沿着边界滑回内部。

- **情况 3**：弹珠在圆盘边界上，你想把它**往圆盘外面推**（$\tau^T\theta > 0$，$\tau$ 与 $\theta$ 同向）。你不能让它出去——所以你把 $\tau$ 中指向外面的分量（即 $\tau$ 在 $\theta$ 方向上的投影 $\frac{\tau^T\theta}{\|\theta\|^2}\theta$）砍掉，只保留垂直于 $\theta$ 的切向分量。弹珠沿着边界**滑动**而不离开圆盘。

这就是为什么投影后的更新方向是 $\tau - \frac{\tau^T\theta}{\|\theta\|^2}\theta$——它等价于「沿着圆盘边界的切线方向走」。物理上这就是离心力被碗壁的法向反力抵消后，物体只能沿碗壁滑动的道理。

**特例：球形约束**。当 $\Theta = \{\theta : \|\theta\| \le \theta_{\max}\}$ 时（最常用的形式），$g(\theta) = \frac{1}{2}(\theta^T\theta - \theta_{\max}^2)$，$\nabla g(\theta) = \theta$。投影简化为：

$$
\operatorname{Proj}_\Theta(\theta, \tau) = \begin{cases}
\tau, & \|\theta\| < \theta_{\max} \text{ or } (\|\theta\| = \theta_{\max} \text{ and } \tau^T\theta \le 0) \\[6pt]
\tau - \dfrac{\theta\theta^T}{\theta^T\theta} \tau = \tau - \dfrac{\tau^T\theta}{\|\theta\|^2}\theta, & \|\theta\| = \theta_{\max} \text{ and } \tau^T\theta > 0
\end{cases}
$$

**定理 A.4（投影的性质）** 若 $\theta^* \in \Theta$（真实参数在可行域内），则对任意 $\tau$：

$$
(\theta - \theta^*)^T \Gamma^{-1}(\operatorname{Proj}_\Theta(\theta, \tau) - \tau) \le 0 \tag{A.16}
$$

*证明思路*：此不等式保证了投影「不会伤害」Lyapunov 稳定性分析——在 Lyapunov 函数的导数中，投影项总是产生非正贡献。详见 Ioannou & Sun (2012), §4.4。$\square$

---

## A.6 小增益定理与鲁棒稳定性

> 🎯 **自查**：你理解 Nyquist 判据吗？小增益定理是 Nyquist 判据在范数框架下的推广——只要开环增益的某种范数小于 1，闭环就稳定。

### A.6.1 $\mathcal{L}_2$ 小增益定理

**定理 A.5（小增益定理，Small Gain Theorem）** 设 $H_1(s)$ 和 $H_2(s)$ 均为稳定传递函数，且 $\|H_1 H_2\|_\infty < 1$。则反馈互联系统 $(I - H_1 H_2)^{-1}$ 是 $\mathcal{L}_2$-稳定的。

*证明*（简要）：由 Nyquist 判据，闭环稳定的条件是 $H_1(j\omega)H_2(j\omega)$ 的 Nyquist 轨迹不环绕 $-1$ 点。若 $\sup_\omega |H_1(j\omega)H_2(j\omega)| < 1$，则轨迹始终在单位圆内，不可能与负实轴上的 $[-1, 0]$ 段相交，故不环绕 $-1$。$\square$

### A.6.2 自适应系统的鲁棒稳定条件

**引理 A.3（自适应控制中的小增益条件）** 对于形如 $u = -F K(\theta)[z]$（其中 $K(\theta)$ 是自适应滤波器）的控制器结构，若：

$$
\boxed{\kappa_0 \cdot \bar{\theta}_M \cdot \|F G_0 \Delta_m\|_{\mathcal{L}_1} < 1} \tag{A.17}
$$

成立，则闭环系统是 $\mathcal{L}_2$-稳定的。其中 $\kappa_0$ 是 $F(s)$ 的增益标量，$\bar{\theta}_M$ 是 $\|\theta(t)\|$ 的上界（由投影保证），$\|\cdot\|_{\mathcal{L}_1} = \int_0^\infty |h(t)| dt$ 是冲激响应的绝对积分。

此条件将**自适应参数界**（$\bar{\theta}_M$）和**未建模动态界**（$\|\cdot\|_{\mathcal{L}_1}$）耦合为一个可直接检查的数值不等式。投影参数 $\theta_{\max}$ 的选择必须使 (A.17) 成立。

---

# Part B：H∞ 鲁棒控制理论基础

---

## B.1 H∞ 范数

> 🎯 **自查**：你理解传递函数的 Bode 图吗？H∞ 范数精确刻画了「Bode 幅频曲线的最高峰」——系统在所有频率上对正弦输入的最大放大倍数。

### B.1.1 定义

**定义 B.1（H∞ 范数）** 对于稳定的传递函数矩阵 $G(s) \in \mathcal{RH}_\infty^{p \times q}$，其 H∞ 范数定义为：

$$
\boxed{\|G\|_\infty \triangleq \sup_{\omega \in \mathbb{R}} \bar{\sigma}\big(G(j\omega)\big)} \tag{B.1}
$$

其中 $\bar{\sigma}(\cdot)$ 表示矩阵的最大奇异值。对 SISO 系统，简化为 $\|G\|_\infty = \sup_{\omega \in \mathbb{R}} |G(j\omega)|$——即 Bode 幅频曲线的峰值。

### B.1.2 物理意义：能量增益

**引理 B.1（H∞ 范数的时域解释）** H∞ 范数等于系统在 $\mathcal{L}_2$ 信号空间上的**诱导范数**（induced norm）：

$$
\|G\|_\infty = \sup_{u \in \mathcal{L}_2, u \ne 0} \frac{\|y\|_2}{\|u\|_2}, \qquad \|x\|_2^2 \triangleq \int_0^\infty x^T(t)x(t) dt \tag{B.2}
$$

即：$\|G\|_\infty$ 是系统对所有可能的能量有限输入信号，输出能量相对输入能量的**最大放大倍数**。

*证明梗概*：通过 Parseval 定理将时域 $\mathcal{L}_2$ 范数转换为频域积分，再利用奇异值分解取上确界。$\square$

**推论 B.1（性能要求的 H∞ 表述）**：

- **扰动抑制**：$\|S\|_\infty$ 小 $\iff$ 在所有频率上对扰动的最大放大倍数小
- **鲁棒稳定性**：$\|T\|_\infty$ 小 $\iff$ 对高频未建模动态的峰值增益受限

### B.1.3 范数的基本性质

**引理 B.2（H∞ 范数的次乘性）** $\|G_1 G_2\|_\infty \le \|G_1\|_\infty \cdot \|G_2\|_\infty$。

**引理 B.3（小增益定理的 H∞ 表述）** 若 $H_1, H_2 \in \mathcal{RH}_\infty$ 且 $\|H_1 H_2\|_\infty < 1$，则 $(I - H_1 H_2)^{-1} \in \mathcal{RH}_\infty$。

---

## B.2 广义被控对象与混合灵敏度

> 🎯 **自查**：你理解标准反馈控制回路的框图吗？广义被控对象不过是把它重新画成「一切输入到一切输出」的矩阵传递函数。

### B.2.1 标准 H∞ 控制问题

**定义 B.2（广义被控对象）** 将你要控制的系统、权重函数、以及它们之间的连接关系打包成广义被控对象 $P(s)$：

$$
\begin{bmatrix} z \\ y \end{bmatrix} = P(s) \begin{bmatrix} w \\ u \end{bmatrix} =
\begin{bmatrix} P_{11}(s) & P_{12}(s) \\ P_{21}(s) & P_{22}(s) \end{bmatrix}
\begin{bmatrix} w \\ u \end{bmatrix}
$$

其状态空间实现为：

$$
\boxed{P(s) = \begin{bmatrix} A & B_1 & B_2 \\ C_1 & D_{11} & D_{12} \\ C_2 & D_{21} & D_{22} \end{bmatrix} =
\left[\begin{array}{c|cc} A & B_1 & B_2 \\ \hline C_1 & D_{11} & D_{12} \\ C_2 & D_{21} & D_{22} \end{array}\right]} \tag{B.3}
$$

其中 $w$ 是外部输入（扰动、噪声、参考），$u$ 是控制输入，$z$ 是性能输出（希望抑制的信号），$y$ 是测量输出。

**定义 B.3（标准 H∞ 最优控制问题）** 求控制器 $K(s)$（$u = K(s) y$）使得：

1. 闭环系统内稳定（internal stability）
2. 从 $w$ 到 $z$ 的闭环传递函数 $T_{zw}(s) = \mathcal{F}_\ell(P, K)$ 的 H∞ 范数最小化：

$$
\boxed{\gamma_{\text{opt}} = \inf_{K \text{ stabilizing}} \| \mathcal{F}_\ell(P, K) \|_\infty} \tag{B.4}
$$

其中 $\mathcal{F}_\ell(P, K) = P_{11} + P_{12}K(I - P_{22}K)^{-1}P_{21}$ 称为下线性分式变换（lower LFT）。

**定义 B.4（次优 H∞ 控制问题）** 给定 $\gamma > 0$，求稳定化控制器 $K(s)$ 使得 $\| \mathcal{F}_\ell(P, K) \|_\infty < \gamma$。（$\gamma = 1$ 是最常见的目标——若成功，则所有加权性能目标满足。）

### B.2.2 混合灵敏度问题

**定义 B.5（混合灵敏度）** 对标准反馈回路（被控对象 $G(s)$，控制器 $K(s)$），定义：

- **灵敏度函数**：$S(s) = (I + G(s)K(s))^{-1}$（从扰动到输出的传递函数）
- **互补灵敏度函数**：$T(s) = G(s)K(s)(I + G(s)K(s))^{-1} = I - S(s)$
- **控制灵敏度函数**：$K(s)S(s)$（从扰动到控制输入）

混合灵敏度问题是通过三个频域权重函数 $W_1(s), W_2(s), W_3(s)$ 同时约束三者：

$$
\boxed{\left\| \begin{bmatrix} W_1 S \\ W_2 K S \\ W_3 T \end{bmatrix} \right\|_\infty < 1} \tag{B.5}
$$

**权重函数的设计直觉**：

| 权重 | 形状 | 作用 |
|:---|:---|:---|
| $W_1(s)$ | **低通**（低频高增益） | 迫使 $S$ 在低频段很小 → 好的扰动抑制 |
| $W_2(s)$ | **高通或常数** | 限制控制能量，防止执行器饱和 |
| $W_3(s)$ | **高通**（高频高增益） | 迫使 $T$ 在高频段很小 → 对未建模动态鲁棒 |

在 MATLAB 中，`augw(G, W1, W2, W3)` 自动将 $G$ 和三个权重连接为广义被控对象 $P$。

**直觉 7（广义被控对象的物理意义——「一张写了所有关系的接线图」）**：这是理解 H∞ 时最容易卡住的概念。用一句话说：**广义被控对象就是把「你的被控对象 + 所有的权重滤波器 + 它们的连接关系」打包成一个矩阵传递函数。**

不要在脑中默念「$P_{11}, P_{12}, P_{21}, P_{22}$」——用物理信号来理解。回到这张图：

```
              ┌───────────────────────────────┐
              │        广义被控对象 P(s)        │
              │                               │
    w ───────┤                               ├──────→ z
  （扰动）    │    ┌─────┐    ┌─────┐    ┌───┐ │   （性能指标）
              │    │ W₁  │←───│  S  │←───┤   │ │
              │    └─────┘    └─────┘    │   │ │
              │                          │ G │ │
              │    ┌─────┐    ┌─────┐    │   │ │
    u ───────┤    │ W₂  │←───│ KS  │←───┤   ├─┼──────→ y
  （控制输入）│    └─────┘    └─────┘    │   │ │   （测量输出）
              │    ┌─────┐    ┌─────┐    └───┘ │
              │    │ W₃  │←───│  T  │←─────────┤
              │    └─────┘    └─────┘           │
              └───────────────────────────────┘
```

每一根信号线的物理含义：

| 信号 | 物理含义 | 来自哪里 | 去哪里 |
|:---|:---|:---|:---|
| $w$ | **外部世界进入系统的所有扰动**。在 ANC 中就是噪声 $d(t)$ | 外部 | 进入 $G$（作为扰动）和 $P_{11}$ |
| $u$ | **你设计的控制器的输出**。在 ANC 中就是扬声器的驱动信号 | 控制器 $K$ | 进入 $G$ 的输入端 |
| $z$ | **所有你想"压制"的信号的加权组合**。$z = [W_1 S w, \; W_2 K S w, \; W_3 T w]^T$ | 被控对象 + 权重 | 这是你要最小化 $\|z\|_2 / \|w\|_2$ 的对象 |
| $y$ | **控制器能测量到的信号**。在 ANC 中就是误差传声器信号 | 被控对象输出 + $w$ | 输入到控制器 $K$ |

**关键洞察**：$P(s)$ 的四个子块各有物理含义：

- **$P_{11}$**：$w \to z$ 的**直接**传递（不经过控制器）——代表即使 $u = 0$，扰动也会通过权重通路产生性能损失
- **$P_{12}$**：$u \to z$——控制输入通过被控对象 + 权重对性能的影响。你想让 $K$ 通过这条路**抵消** $P_{11} w$
- **$P_{21}$**：$w \to y$——扰动到测量输出的传递。这就是「噪声是怎么被传感器看到的」
- **$P_{22}$**：$u \to y$——这就是你的**被控对象本身** $G(s)$！控制输入如何影响测量输出

所以 LFT $\mathcal{F}_\ell(P, K) = P_{11} + P_{12}K(I - P_{22}K)^{-1}P_{21}$ 的物理含义是：

$$\text{闭环性能} = \underbrace{P_{11}}_{\text{直接馈通}} + \underbrace{P_{12}K(I - P_{22}K)^{-1}P_{21}}_{\text{控制器通过被控对象施加的抵消作用}}$$

H∞ 的目标就是选 $K$ 使得这两项的 H∞ 范数最小。

**数值例子（题 3 的 $P_{aug}$）**：在题 3 中调用 `augw(G, W1, 1e-3, W3)`，返回的 $P_{aug}$ 是你原来的二阶 $G(s)$ 加上 $W_1$（二阶低通）和 $W_3$（一阶高通）的状态串联——总共约 5-7 个状态。`hinfsyn` 求解后返回的 $K(s)$ 也是相似阶数的动态控制器。你不需要手动写出 $P_{aug}$ 的 4 个子块——`augw` 替你做了所有连线。

---

## B.3 DGKF 方法（双 Riccati 方程解法）

> 🎯 **自查**：你已掌握代数 Riccati 方程 $A^T P + PA - PBR^{-1}B^T P + Q = 0$（来自 LQR）。DGKF 是它的 H∞ 推广。

### B.3.1 简化假设

为简化表述，对 $P(s)$ 施加以下简化假设：

| 条件 | 含义 |
|:---|:---|
| $(A, B_1)$ 可稳，$(C_1, A)$ 可检测 | 标准能稳/能观条件 |
| $(A, B_2)$ 可稳，$(C_2, A)$ 可检测 | 保证存在稳定化控制器 |
| $D_{12}^T [C_1 \; D_{12}] = [0 \; I]$ | 正交归一化（无交叉项） |
| $\begin{bmatrix} B_1 \\ D_{21} \end{bmatrix} D_{21}^T = \begin{bmatrix} 0 \\ I \end{bmatrix}$ | 正交归一化（无交叉项） |
| $D_{11} = 0, D_{22} = 0$ | 无前馈（常见于混合灵敏度） |

这些假设在实际问题中可通过**回路变换**（loop-shifting）实现，MATLAB 的 `hinfsyn` 自动处理。

### B.3.2 DGKF 定理

**定理 B.1（DGKF 定理，Doyle-Glover-Khargonekar-Francis, 1989）** 在上述简化假设下，存在一个稳定化控制器 $K(s)$ 使 $\| \mathcal{F}_\ell(P, K) \|_\infty < \gamma$ 的**充要条件**是以下三个条件同时成立：

**条件 1（控制 Riccati 方程）**：存在半正定解 $X_\infty \ge 0$ 满足

$$
\boxed{A^T X_\infty + X_\infty A + C_1^T C_1 + X_\infty(\gamma^{-2} B_1 B_1^T - B_2 B_2^T) X_\infty = 0} \tag{B.6}
$$

且 $A + (\gamma^{-2}B_1 B_1^T - B_2 B_2^T) X_\infty$ 是稳定矩阵（所有特征值具有负实部）。

**条件 2（滤波 Riccati 方程）**：存在半正定解 $Y_\infty \ge 0$ 满足

$$
\boxed{A Y_\infty + Y_\infty A^T + B_1 B_1^T + Y_\infty(\gamma^{-2} C_1^T C_1 - C_2^T C_2) Y_\infty = 0} \tag{B.7}
$$

且 $A + Y_\infty(\gamma^{-2}C_1^T C_1 - C_2^T C_2)$ 是稳定矩阵。

**条件 3（耦合条件）**：

$$
\boxed{\rho(X_\infty Y_\infty) < \gamma^2} \tag{B.8}
$$

其中 $\rho(\cdot)$ 表示矩阵的谱半径（最大特征值的模）。

**注释 B.1（与 LQR/LQE 的类比）**：
- (B.6) 是 LQR Riccati 方程的推广——当 $\gamma \to \infty$（不加 H∞ 约束），$\gamma^{-2}B_1B_1^T \to 0$，退化为标准 LQR 的 ARE：$A^T X + XA - XB_2 B_2^T X + C_1^T C_1 = 0$
- (B.7) 是 Kalman 滤波 Riccati 方程的推广
- (B.8) 是 LQG 中分离原理在 H∞ 框架下的类比——$X_\infty Y_\infty$ 的谱半径度量了控制和估计之间的耦合程度

**直觉 8（为什么需要两个 Riccati 方程？——「控制」和「观测」不再分离）**：这是 H∞ 与 LQG 最根本的区别。

在 **LQG** 中，你可以**先设计 LQR**（不管噪声），**再设计 Kalman 滤波器**（不管控制），然后把两者串联——分离原理保证这依然是最优的。两个 Riccati 方程是**独立求解**的。

在 **H∞** 中，分离原理**不成立**。因为 H∞ 优化的不是期望值，而是**最坏情况**能量增益。控制器的状态反馈增益 $F_\infty = -B_2^T X_\infty$ 和观测器注入增益 $L_\infty = -Y_\infty C_2^T$ 通过 $\gamma$ **相互耦合**——你设计控制器时必须同时考虑观测器，设计观测器时必须同时考虑控制器。这就是为什么需要**两个** Riccati 方程，且它们被耦合条件 $\rho(X_\infty Y_\infty) < \gamma^2$ 绑定在一起。

**直觉 9（谱半径条件 $\rho(X_\infty Y_\infty) < \gamma^2$ 的物理含义——「控制与估计不能同时太强」）**：

$X_\infty$ 度量了「状态反馈控制有多强」——$X_\infty$ 大意味着 LQR 代价中的状态惩罚项 $C_1^T C_1$ 很大，你需要很强的控制来压制它。

$Y_\infty$ 度量了「状态估计有多精确」——$Y_\infty$ 大意味着过程噪声 $B_1 B_1^T$ 很大（或传感器 $C_2$ 很弱），你需要很强的观测器来重建状态。

耦合条件 $\rho(X_\infty Y_\infty) < \gamma^2$ 说：**你不能同时拥有极强的控制和极强的估计。**

这是一个非常深刻的物理约束。它本质上说的是：

> 如果控制器严重依赖精确的状态估计（$Y_\infty$ 大）来产生大幅度的控制动作（$X_\infty$ 大），那么任何估计误差都会被放大——在最坏情况下，这会导致闭环不稳定。$\gamma$ 越大（你对性能越宽松），这对 $(X_\infty, Y_\infty)$ 的绑定就越松；$\gamma$ 越小（你要求性能越高），它们就越被挤压。在 $\gamma \to \gamma_{\text{opt}}$ 的极限下，$\rho(X_\infty Y_\infty) \to \gamma_{\text{opt}}^2$——饱和。

这就是为什么 `hinfsyn` 找不到 $\gamma < 1$ 的控制器时你需要「放宽权重」：放宽 $W_1$ 等价于减小 $C_1$（减小了 $X_\infty$），放宽 $W_3$ 等价于允许更大的 $T$（减小了高频约束），两者都会松动 $\rho(X_\infty Y_\infty)$ 的上限。

### B.3.3 中心控制器的构造

**定理 B.2（中心控制器的状态空间实现）** 当 DGKF 三个条件满足时，一个可行的次优控制器 $K(s)$ 具有状态空间实现：

$$
\boxed{K(s) = \begin{bmatrix} A_K & B_K \\ C_K & D_K \end{bmatrix}}
$$

其中：

$$
\begin{aligned}
A_K &= A + \gamma^{-2} B_1 B_1^T X_\infty + B_2 F_\infty + Z_\infty L_\infty C_2 \\
B_K &= -Z_\infty L_\infty \\
C_K &= F_\infty \\
D_K &= 0
\end{aligned} \tag{B.9}
$$

且 $F_\infty = -B_2^T X_\infty$（状态反馈增益），$L_\infty = -Y_\infty C_2^T$（输出注入增益），$Z_\infty = (I - \gamma^{-2} Y_\infty X_\infty)^{-1}$。

*证明*：参见 Doyle, Glover, Khargonekar & Francis (1989)。关键是验证通过此 $K(s)$，闭环满足 $\|T_{zw}\|_\infty < \gamma$。$\square$

---

## B.4 γ-迭代与二分搜索

> 🎯 你不需要手动求解 Riccati 方程——MATLAB 的 `hinfsyn` 自动完成一切。但理解其内部 γ-迭代逻辑有助于你诊断「gamma ≥ 1」的失败。

### B.4.1 从次优到最优：γ-迭代

DGKF 定理给出的是**次优**条件（对给定 $\gamma$）。最优 H∞ 范数 $\gamma_{\text{opt}}$ 通过以下二分搜索确定：

```
Algorithm B.1: γ-Iteration (hinfsyn 内部逻辑)
────────────────────────────────────────
Input:  P(s)（广义被控对象）, ε（容差，~1e-6）
Output: K_opt(s), γ_opt

1: γ_low  ← 某个下界（足够大使得条件成立）
2: γ_high ← 某个上界（足够小使得条件不成立）
3: while (γ_high − γ_low) / γ_low > ε:
4:     γ ← (γ_low + γ_high) / 2
5:     尝试求解 DGKF 条件 (B.6)-(B.8) 对当前 γ
6:     if 三个条件都满足:
7:         γ_high ← γ    ▷ 可以更小
8:     else:
9:         γ_low ← γ     ▷ 需要更大
10: end while
11: γ_opt ← γ_high
12: 用 γ_opt 代入 (B.9) 构造 K_opt(s)
```

### B.4.2 γ 值的判读

| γ 值 | 含义 | 应对策略 |
|:---|:---|:---|
| $\gamma < 1$ | **设计成功**。所有加权性能目标满足 | 可以尝试增大 $W_1$ 的增益（更严格的性能要求） |
| $\gamma \approx 1$ | 刚好满足，接近性能极限 | 可以接受，或略微放宽某个权重 |
| $\gamma > 1$ | **设计失败**。在给定权重下不存在满足所有要求的稳定化控制器 | 放宽 $W_1$ 的带宽/增益，或降低 $W_3$ 的高频要求 |
| 无解 | ARE 求解失败或谱半径条件不满足 | 检查 $P(s)$ 的能稳/能观性；放宽权重 |

### B.4.3 与 LQG 的终极对比

| 维度 | LQG | H∞ (DGKF) |
|:---|:---|:---|
| 不确定性模型 | 白噪声（已知统计特性） | 范数有界（未知统计特性） |
| 性能指标 | $\mathbb{E}[\int (x^T Q x + u^T R u)dt]$ | $\|T_{zw}\|_\infty$（最坏情况能量增益） |
| Riccati 方程 | 1 个（LQR ARE）| 2 个（控制 + 滤波）+ 耦合条件 |
| 分离原理 | 严格成立 | 不成立——控制器和观测器耦合（通过 $\gamma^2$） |
| 鲁棒性 | 极差（LQG 没有鲁棒性保证！） | 显式保证（通过 $W_3 T$ 约束） |

---

## 参考文献

1. Ioannou, P.A. & Sun, J. (2012). *Robust Adaptive Control*. Dover Publications. （Part A 的主要参考）
2. Doyle, J.C., Glover, K., Khargonekar, P.P., & Francis, B.A. (1989). "State-space solutions to standard $H_2$ and $H_\infty$ control problems." *IEEE Trans. Automatic Control*, 34(8), 831–847. （Part B 的原始论文）
3. Zhou, K., Doyle, J.C., & Glover, K. (1996). *Robust and Optimal Control*. Prentice Hall. （H∞ 综合的权威教材）
4. Jafari, S. & Ioannou, P.A. (2017). "Adaptive feedback control of periodic disturbances in the presence of time-varying feedback path." *Journal of Vibration and Control*, 23(5), 826–841.

---

> 📖 **阅读下一步**：
> - Part A 的直接应用 → [题1：连续时间 CC 自适应控制](Derivation_Problem1_Jafari_ContinuousCC.md)
> - Part A + Part B 的结合应用 → [题2：离散时间 AVC 鲁棒自适应](Derivation_Problem2_Jafari_DiscreteAVC.md)
> - Part B 的标准应用 → [题3：H∞ 混合灵敏度控制](Derivation_Problem3_HinfSynthesis.md)
