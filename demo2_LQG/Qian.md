这是一个关于“管道内有界随机噪声与多频率音调噪声的主动反馈控制”的完整学术推导文档。本推导严格遵循文中提及的 **LQG + Youla (Q) 参数化自适应控制** 框架。

---

# 管道噪声主动反馈衰减控制系统：问题制定与数学求解推导

## 一、 问题制定 (Problem Formulation)

### 1.1 被控对象（次级通路）模型

在管道主动噪声控制（ANC）系统中，声波通过扬声器（次级声源）产生的补偿声场抵消原始噪声。设定次级通路的离散时间状态空间模型为 $\Sigma_p$：

$$
\Sigma_p : \begin{cases} x(k+1) = Ax(k) + Bu(k) \\ e(k) = y(k) = Cx(k) + \omega(k) + v(k) \end{cases} \quad \text{--- (1)}
$$

其中：

* $x(k) \in \mathbb{R}^n$ 为系统状态。
* $u(k)$ 为控制器输出信号（输入到补偿扬声器）。
* $e(k)$ 为误差信号（传感器/麦克风采集的残余噪声）。
* $\omega(k)$ 为等效噪声信号，包含音调噪声 $\omega_1(k)$ 和有界随机噪声 $\omega_2(k)$。
* $v(k)$ 为测量白噪声。

### 1.2 噪声特性描述

* **音调噪声 ($\omega_1$)**：具有未知且时变的频率 $B_\ell$。
  $$
  \omega_1(k) = \sum_{\ell=1}^{s} A_\ell(k) \sin[B_\ell(k)k + \phi_\ell(k)]
  $$
* **有界随机噪声 ($\omega_2$)**：由白噪声 $n(k)$ 通过成形滤波器 $\Sigma_{\omega2}$ 产生：
  $$
  \Sigma_{\omega2} : \begin{cases} x_\omega(k+1) = A_\omega x_\omega(k) + B_\omega n(k) \\ \omega_2(k) = C_\omega x_\omega(k) \end{cases} \quad \text{--- (2)}
  $$

### 1.3 增广系统构建

将被控对象与随机噪声模型合并，得到增广状态空间模型 $\Sigma$：

$$
\Sigma : \begin{cases} \tilde{x}(k+1) = \tilde{A}\tilde{x}(k) + \tilde{B}_1 u(k) + \tilde{B}_2 n(k) \\ e(k) = \tilde{C}\tilde{x}(k) + \omega_1(k) + v(k) \end{cases} \quad \text{--- (3)}
$$

其中 $\tilde{x}(k) = [x(k)^T, x_\omega(k)^T]^T$。

---

## 二、 第一步：基础 LQG 控制器设计 (Base LQG Design)

### 2.1 目标函数

设计最优控制器 $u(k)$ 以最小化二次型性能指标 $J$，主要用于衰减有界随机噪声 $\omega_2$：

$$
J = \lim_{N \to \infty} \frac{1}{N} E \sum_{k=0}^{N-1} (\tilde{x}_k^T Q_e \tilde{x}_k + u_k^T R_e u_k) \quad \text{--- (4)}
$$

其中 $Q_e = \tilde{C}^T \tilde{C}$，$R_e = 10^{-4}$（标准实验设定）。

### 2.2 代数黎卡提方程 (ARE) 求解

根据**分离原理**，LQG 控制器由状态反馈增益 $K$ 和观测器增益 $L$ 组成：

1. **求解反馈增益 $K$**：

   $$
   K = -(R_e + \tilde{B}_1^T P \tilde{B}_1)^{-1} \tilde{B}_1^T P \tilde{A}
   $$

   其中 $P$ 是以下 ARE 的解：

   $$
   P = \tilde{A}^T P \tilde{A} - \tilde{A}^T P \tilde{B}_1 (R_e + \tilde{B}_1^T P \tilde{B}_1)^{-1} \tilde{B}_1^T P \tilde{A} + Q_e
   $$
2. **求解观测器增益 $L$**：

   $$
   L = X \tilde{C}^T (\tilde{C} X \tilde{C}^T + V_e)^{-1}
   $$

   其中 $X$ 是估计误差协方差矩阵，由对应的 ARE 求解。

### 2.3 结论：基础控制器 $K_{LQG}$

$$
K_{LQG} = \left[ \begin{array}{c:c} \tilde{A} + \tilde{B}_1 K + L\tilde{C} & -L \\ \hdashline K & 0 \end{array} \right] \quad \text{--- (5)}
$$

---

## 三、 第二步：Youla (Q) 参数化与自适应调节

为了在不破坏系统稳定性的前提下抑制频率未知的音调噪声 $\omega_1$，引入 Youla 参数化。

### 3.1 控制器结构变换

根据 Youla 定理，所有稳定控制器可表示为 $J$ 块与稳定参数 $Q$ 的连接。
**$J$ 块（系统内部逻辑）**：

$$
\begin{cases} \hat{x}(k+1) = (\tilde{A} + L\tilde{C} + \tilde{B}_1 K)\hat{x}(k) - Ly(k) + \tilde{B}_1 y_Q(k) \\ u(k) = K\hat{x}(k) + y_Q(k) \\ r(k) = y(k) - \tilde{C}\hat{x}(k) \end{cases} \quad \text{--- (6)}
$$

其中 $r(k)$ 是用于自适应调节的残差信号。

### 3.2 自适应 Youla 参数 $Q(z)$

设定 $Q(z)$ 为一个线性参数化形式：

$$
Q(z) = \left( \sum_{i=1}^{n_q} \theta_i z^{1-i} \right) H(z) \quad \text{--- (7)}
$$

* $\theta = [\theta_1, \dots, \theta_{n_q}]^T$ 为待调参数矢量。
* $H(z)$ 为带通滤波器，用于限制自适应动作在噪声频率范围内。

---

## 四、 第三步：参数估计算法 (RLS求解)

利用递归最小二乘法 (RLS) 在线修正 $\theta(k)$。

### 4.1 误差转换

定义修改后的误差信号 $\tilde{e}(k)$，使其与参数偏差呈线性关系：

$$
\tilde{e}(k) = \phi^T(k) \tilde{\theta}(k) + e^0(k)
$$

其中 $\phi(k)$ 是经 $T_{12}$ 和 $H(z)$ 滤波后的回归矢量。

### 4.2 RLS 迭代公式

1. **计算先验误差** $\tilde{e}(k)$。
2. **计算增益矩阵 $P(k)$**（引入遗忘因子 $\lambda$）：
   $$
   P(k+1) = \frac{1}{\lambda} \left[ P(k) - \frac{P(k)\phi(k)\phi^T(k)P(k)}{1 + \phi^T(k)P(k)\phi(k)} \right] \quad \text{--- (8)}
   $$
3. **更新参数矢量**：
   $$
   \hat{\theta}(k+1) = \hat{\theta}(k) + P(k)\phi(k) \frac{\tilde{e}(k)}{1 + \phi^T(k)P(k)\phi(k)} \quad \text{--- (9)}
   $$

---

## 五、 标准实验模型数据 (Experimental Data)

根据文中的实验验证部分，以下为实现该控制器的关键参数：

| 参数类别                       | 具体数值 / 描述                                                                   |
| :----------------------------- | :-------------------------------------------------------------------------------- |
| **系统阶数**             | 识别得到的次级通路模型$\Sigma_p$ 为 **40阶** ARX模型                      |
| **采样频率**             | $f_s = 5000 \text{ Hz}$                                                         |
| **音调噪声频率范围**     | $200 \sim 650 \text{ Hz}$                                                       |
| **随机噪声带宽**         | $600 \sim 650 \text{ Hz}$ (Band-limited)                                        |
| **Youla 参数阶数**       | $n_q = 6$ (Case 1/2) 或 $n_q = 8$ (Case 3)                                    |
| **遗忘因子 $\lambda$** | 固定$\lambda = 0.98$ 或 变遗忘因子 $\lambda(k) = 0.9996\lambda(k-1) + 0.0004$ |
| **权重函数 $H(z)$**    | 截止频率为$200$ 和 $650 \text{ Hz}$ 的带通滤波器 (BPF)                        |
| **LQG 惩罚参数**         | $Q_e = \tilde{C}^T\tilde{C}$, $R_e = 10^{-4}$                                 |
| **噪声方差设定**         | $N_e = 200$, $V_e = 10^{-2}$                                                  |

---

## 六、 控制器结论 (Conclusion)

通过上述推导，得出的控制器结论如下：

1. **稳定性保证**：由于采用 Youla 参数化框架，只要 $Q(z)$ 是稳定的且满足受限增益条件 ($|Q\Delta|_\infty < 1$)，无论自适应参数如何变化，闭环系统始终保持全局稳定。
2. **双重抑制能力**：
   * **基础 LQG 部分** 通过对随机噪声状态的估计和最优反馈，实现了对 $600-650\text{Hz}$ 宽带噪声的预压制。
   * **自适应 $Q$ 部分** 基于内模原理 (IMP)，通过 RLS 算法使 $Q$ 动态收敛于音调噪声频率的逆模型，实现对 $350\text{Hz}$、$430\text{Hz}$ 等离散频率噪声的完全消除（实验显示衰减达 **-30dB 至 -35dB**）。
3. **动态鲁棒性**：引入带通滤波器 $H(z)$ 和变遗忘因子，有效解决了自适应算法在非目标频段的溢出问题，并提升了对时变频率（线性扫频噪声）的跟踪速度。

---

**End of Document**
