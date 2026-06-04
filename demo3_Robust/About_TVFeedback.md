基于提供的论文《Adaptive feedback suppression of unknown periodic components of acoustic noises with time-varying characteristics》，针对该被控系统，以下是问题制定与求解的完整数学推导过程。本文档直接对应您所要求的各个核心部分。

# 具有时变特征的未知声学周期噪声自适应反馈抑制推导

## 1. 系统模型 (被控对象与问题制定)
根据论文 Section 3 (Problem formulation) 和 Section 6 (Numerical simulation)，系统的控制目标是通过生成控制信号 $u(t)$ 来使得麦克风位置的误差信号 $y(t)$ 最小化，以抑制未知的声学扰动 $d(t)$。

带建模误差的系统实际动态方程（参见 Eq. 1 及 Fig. 2）为：
$$y(t) = G(s)[u(t)] + d(t)$$
$$G(s) = G_0(s)(1 + \Delta_m(s))$$

其中各个参数定义如下：
*   **$G_0(s)$**：已知标称被控对象（Nominal plant model）的传递函数，包含扬声器、声学路径和麦克风的动态特性。假定具有稳定的极点，但可能存在非最小相位零点（右半平面零点）。
*   **$\Delta_m(s)$**：未建模动态（乘性不确定性），代表高频范围的误差或输入延迟。
*   **$d(t)$**：未知的周期性窄带声学干扰信号，由未知时变频率、相位和幅值的正弦波以及宽带噪声 $v(t)$ 组成（参见 Eq. 2）：
    $$d(t) = \sum_{i=1}^{n_f(t)} a_i(t) \sin(\omega_i(t)t + \varphi_i(t)) + v(t)$$
*   **具体仿真中的 $G_0(s)$ 参数 (Section 6)**：
    $$G_0(s) = \frac{0.5(s - 0.2)}{s^2 + s + 1.25}$$
    其具有一个右半平面不稳定零点 $s = 0.2$。未建模动态设为 $\Delta_m(s) = -\mu s$ （$\mu = 0.001$）。

---

## 2. 回归向量 $\phi$ 的构造
为了在线自适应地估计未知扰动的频率并进行抑制，系统设计了一个参数化的自适应控制器（Eq. 6 - 8）。

**控制架构：**
控制输入由下式生成：
$$u(t) = -F(s)[K(s, \theta(t))[z(t)]]$$
其中辅助信号 $z(t)$ 通过剥离标称控制输入的影响得到：
$$z(t) = y(t) - G_0(s)[u(t)]$$

**自适应滤波器 $K(s, \theta(t))$ 与基向量 $\Lambda(s)$：**
滤波器采用以下结构（Eq. 3）：
$$K(s, \theta(t)) = \sum_{k=0}^{N-1} \theta_k(t) \frac{\lambda^{N-k}}{(s+\lambda)^{N-k}}$$
我们可以将其写为参数与信号内积的形式：$K(s, \theta(t))[z(t)] = \theta^T(t)\Lambda(s)[z(t)]$，其中：
$$\Lambda(s) = \left[ \frac{\lambda^N}{(s+\lambda)^N}, \dots, \frac{\lambda}{s+\lambda} \right]^T$$

**回归向量 $\phi(t)$ 的推导 (Eq. 12)：**
利用 Swapping Lemma，并引入 LTI 滤波器 $F(s)$，令 $G_0'(s) = F(s)G_0(s)$，我们可以得到参数化的线性模型：
$$z(t) = \theta^{*T}(t)\phi(t) + \eta(t)$$
在此模型中，**回归向量 $\phi(t)$** 是完全可测的，其构造方式为将 $z(t)$ 依次通过基向量 $\Lambda(s)$、标称模型 $G_0(s)$ 和滤波器 $F(s)$：
$$\phi(t) = G_0'(s)\Lambda(s)[z(t)] = F(s)G_0(s)\Lambda(s)[z(t)]$$

---

## 3. LTI 滤波器 $F(s)$ 的设计 
**注：** 本文不要求严格正实(SPR)条件，而是利用鲁棒最小二乘(RLS)直接对参数化模型辨识。$F(s)$ 的主要**设计目标是“频谱平坦化 (flatten the spectrum)”**：保证回归向量 $\phi(t)$ 在扰动频率范围内 $[0, \bar{\omega}_{max}]$ 具有足够大的增益，避免由于 $G_0(s)$ 在特定频率衰减导致自适应律停滞。

**设计步骤与约束 (Section 5.1)：**
1. 初始令 $\tilde{G}_0(s) = G_0(s)$。
2. 若 $\tilde{G}_0(s)$ 在虚轴有零点 $s = j\omega_0$，将其微调为 $s = -\delta_0 + j\omega_0$ ($\delta_0 > 0$)。
3. 若 $\tilde{G}_0(s)$ 有右半平面不稳定零点 $s = \sigma_0 \pm j\omega_0$ ($\sigma_0 > 0$)，将其沿虚轴镜像反射至左半平面： $s = -\sigma_0 \pm j\omega_0$。
4. **$F(s)$ 的最终定义 (Eq. 13)：**
   $$F(s) = \kappa_0 \frac{\alpha^m}{(s+\alpha)^m} \tilde{G}_0^{-1}(s) \triangleq \kappa_0 \bar{F}(s)$$
   其中：$m>0$ 是大于 $G_0(s)$ 相对阶数的整数，$\alpha > 0$ 是保证截止频率足够大的低通滤波器常数，$\kappa_0$ 用于调节整体增益。

**鲁棒稳定性设计约束 (Eq. 18)：**
设计 $F(s)$ 时，需要满足未建模动态存在下的鲁棒小增益条件定理：
$$\kappa_0 \bar{\theta}_M ||\bar{F}(s)G_0(s)\Delta_m(s)||_1 < 1$$
*这意味着：提升 $\kappa_0$ 能加快自适应收敛速度，但过大会放大高频处未建模动态 $\Delta_m(s)$ 的影响，导致闭环失稳。*

**仿真中的特定 $F(s)$ (Section 6)：**
对于 $G_0(s) = \frac{0.5(s-0.2)}{s^2+s+1.25}$，其零点 $0.2$ 被镜像为 $-0.2$，取 $m=2$，故 $\tilde{G}_0^{-1}(s) = \frac{2(s^2+s+1.25)}{s+0.2}$，最终得到：
$$F(s) = \kappa_0 \frac{2\alpha^2(s^2 + s + 1.25)}{(s+\alpha)^2(s+0.2)}$$

---

## 4. 自适应律设计 (RLS / Gradient)
采用带死区和参数投影 (Projection) 的鲁棒递归最小二乘算法 (Robust RLS) 在线估计未知参数 $\theta(t)$ (Section 5.2)。

**估计模型与误差定义：**
* 预测器： $\hat{z}(t) = \theta^T(t)\phi(t)$
* 归一化信号： $m_s^2(t) = 1 + \gamma_0 \phi^T(t)\phi(t)$
* 估计误差： $\epsilon(t) = \frac{z(t) - \hat{z}(t)}{m_s^2(t)}$

**参数更新方程 (Eq. 16 - 17)：**
协方差矩阵 $P(t)$ 的微分方程（带有协方差重置修改，当 $\lambda_{min}(P) \le \rho_1$ 时重置 $P(t_r^+) = \rho_0 I$ 防止奇异）：
$$\dot{P}(t) = -P(t) \frac{\phi(t)\phi^T(t)}{m_s^2(t)} P(t)$$

参数向量 $\theta(t)$ 的更新方程 (带有投影算子保证参数有界性 $||\theta|| \le \theta_{max}$)：
$$\dot{\theta}(t) = \text{proj} \left( P(t) \epsilon(t) \phi(t) \right)$$
具体投影算子公式如下（定义凸集边界面 $g(\theta) = \theta^T \theta - \theta_{max}^2 \le 0$）：
$$
\dot{\theta}(t) = 
\begin{cases} 
P(t)\epsilon(t)\phi(t), & \text{if } \theta \in \mathbb{S}_0 \text{ or if } \theta \in \delta\mathbb{S} \text{ and } (P\epsilon\phi)^T \nabla g \le 0 \\ 
P(t)\epsilon(t)\phi(t) - P(t) \frac{\nabla g \nabla g^T}{\nabla g^T P(t) \nabla g} P(t)\epsilon(t)\phi(t), & \text{otherwise}
\end{cases}
$$

---

## 5. 仿真参数表
论文 Section 6 在验证所提理论时，使用了以下具体仿真参数构建被控对象与控制器：

| 参数类别          | 符号                   | 具体数值 (Value)                | 物理含义/备注                                              |
| :---------------- | :--------------------- | :------------------------------ | :--------------------------------------------------------- |
| **被控对象**      | $G_0(s)$               | $\frac{0.5(s-0.2)}{s^2+s+1.25}$ | 标称传递函数，含 $s=0.2$ 右半平面不稳定零点                |
| **未建模动态**    | $\Delta_m(s)$          | $-0.001 s$                      | 乘性高频建模误差 (参数 $\mu=0.001$)                        |
| **扰动信号**      | $\omega_1$, $\omega_2$ | $70$ rad/s, $187$ rad/s         | 周期性扰动频率，后续变频场景中 $\omega_1$ 变为 $100$ rad/s |
| **宽带噪声**      | $v(t)$                 | 均值 $0$, 标准差 $0.02$         | 截断的高斯白噪声，界限为 $|v(t)| \le 0.1$                  |
| **自适应滤波器**  | $N$                    | $20$                            | 阶数设计参数 ($N \ge 2n_{max}$，其中最大频数 $n_{max}=5$)  |
|                   | $\lambda$              | $500$                           | 基函数极点位置参数                                         |
| **$F(s)$ 滤波器** | $\kappa_0$             | $0.5$ ($-6$ dB)                 | 低频段增益控制设计常数                                     |
|                   | $\alpha$               | $500$                           | 极点常数，使得在频率 $[0, 600]$ rad/s 内有足够带宽         |
| **自适应律/RLS**  | $P(0)$                 | $500 I_{20 \times 20}$          | 初始协方差矩阵                                             |
|                   | $\theta(0)$            | $0_{20 \times 1}$               | 参数向量初始值                                             |
|                   | $T_s$                  | $10^{-4}$ s                     | 仿真使用的系统采样周期步长                                 |
|                   | $\gamma_0$             | 常数 $> 0$                      | 归一化设计常数（论文未显式列出具体值，理论需为正数）       |