# 离散时间 SISO 系统鲁棒自适应干扰衰减：问题制定、求解与复现诊断

> 对应原书第 6.6.1 节。基于 Jafari, Ioannou, Fitzpatrick, Wang (IEEE TAC, 2015)。

---

## 1. 问题制定

### 1.1 被控对象

考虑带有未建模动态的离散时间 SISO 系统：

$$
y(k) = G(z)[u(k)] + d(k), \qquad G(z) = G_0(z)\big(1 + \Delta_m(z)\big)
$$

**标称模型**（AVC 主动悬架）：

$$
G_0(z) = \frac{-0.00146(z - 0.1438)(z - 1)}{(z - 0.7096)(z^2 - 0.04369z + 0.01392)}, \qquad T_s = \frac{1}{480}\text{ s}
$$

**未建模动态**（高频主导）：

$$
\Delta_m(z) = \frac{-0.0001}{(z + 0.99)^2}
$$

### 1.2 扰动模型

$$
d(k) = d_s(k) + v(k) = \sum_{i=1}^{n_\omega} A_i \sin(\omega_i T_s k + \varphi_i) + v(k)
$$

其中 $v(k) \sim \mathcal{N}(0, \sigma_v^2)$，$\sigma_v = 0.02$。

**仿真扰动配置**（$0 \le t < 30$ s）：

| $i$ | $\omega_i$ (rad/s) | $\omega_i T_s$ (rad) | $A_i$ | $\varphi_i$ |
|-----|---------------------|----------------------|-------|-------------|
| 1 | 25 | 0.0521 | 0.7 | $\pi/3$ |
| 2 | 225 | 0.4688 | 0.5 | $\pi/4$ |

$t \ge 30$ s 后新增 $\omega_3 = 85$、$\omega_4 = 125$ rad/s。

**控制目标**：设计 $u(k)$ 使闭环系统鲁棒稳定，完全抑制 $d_s(k)$ 并最小化 $v(k)$ 的放大。

---

## 2. 控制器结构

### 2.1 Youla 参数化形式

引入辅佐信号剥离标称模型的贡献：

$$
\zeta(k) = y(k) - G_0(z)[u(k)] \tag{1}
$$

控制律采用 FIR 滤波器 $K(z,\theta)$ 经预补偿器 $F(z)$ 的级联形式：

$$
u(k) = -F(z)\Big[K\big(z, \hat{\theta}(k-1)\big)[\zeta(k)]\Big] \tag{2}
$$

### 2.2 FIR 滤波器参数化

$$
K(z, \hat{\theta}) = \sum_{i=0}^{N-1} \hat{\theta}_i z^{i-N} = \hat{\theta}^\top \alpha(z) \tag{3}
$$

其中基向量：

$$
\alpha(z) = \begin{bmatrix} z^{-N} & z^{-(N-1)} & \cdots & z^{-1} \end{bmatrix}^\top
$$

时域实现：$\alpha(z)[\zeta(k)] = \begin{bmatrix} \zeta(k-N) & \zeta(k-N+1) & \cdots & \zeta(k-1) \end{bmatrix}^\top \triangleq w(k)$。

### 2.3 预补偿器 $F(z)$

在此例中 $F(z) = 100$（静态增益），作用是将回归信号的幅度提升到合理范围。

---

## 3. 参数辨识模型

### 3.1 回归向量的构造

将标称模型 $G_0(z)$ 和预补偿器 $F(z)$ 纳入回归向量：

$$
\phi(k) = G_0(z)F(z)\alpha(z)[\zeta(k)] = F \cdot G_0(z)[w(k)] \tag{4}
$$

分量形式（$i = 1, \dots, N$）：

$$
\phi_i(k) = F \cdot G_0(z)\big[\zeta(k - (N-i+1))\big]
$$

### 3.2 线性参数模型

利用 Swapping Lemma，系统可表示为关于 $\theta^*$ 的线性回归：

$$
\zeta(k) = \theta^{*\top}\phi(k) + \eta(k) \tag{5}
$$

其中 $\eta(k)$ 包含噪声、未建模动态以及参数暂态误差。

---

## 4. 鲁棒自适应律

### 4.1 归一化 RLS

归一化信号与估计误差：

$$
m^2(k) = 1 + \gamma_0 \phi^\top(k)\phi(k), \qquad \gamma_0 = 1 \tag{6}
$$

$$
\varepsilon(k) = \frac{\zeta(k) - \hat{\theta}^\top(k-1)\phi(k)}{m^2(k)} \tag{7}
$$

### 4.2 协方差与参数更新

$$
P(k) = P(k-1) - \frac{P(k-1)\phi(k)\phi^\top(k)P(k-1)}{m^2(k) + \phi^\top(k)P(k-1)\phi(k)} \tag{8}
$$

$$
\hat{\theta}(k) = \text{proj}\Big(\hat{\theta}(k-1) + P(k)\varepsilon(k)\phi(k)\Big) \tag{9}
$$

投影算子 $\text{proj}(\cdot)$ 保证 $\|\hat{\theta}\| \le \theta_{\max}$。

### 4.3 协方差重置（Covariance Resetting）

**这是复现成功的关键机制。** 当干扰被阶段性抑制后，持续激励（PE）条件丧失，$P(k)$ 迅速收缩至近零。若不处理，$t=30$ s 新扰动出现时算法失去自适应能力。

触发条件（Eq. 6.61）：

$$
\lambda_{\min}\big(P(k)\big) \le \beta_1 \;\Longrightarrow\; P(t_r) = \beta_0 I, \qquad \beta_0 > \beta_1 > 0
$$

工程实现可用 $\text{trace}(P)/N \le \beta_1$ 近似触发，$\beta_0 = P(0)$ 即可。

### 4.4 控制律合成

每时刻用更新后的参数生成控制：

$$
u(k) = -F \cdot \hat{\theta}^\top(k-1) w(k) = -F \sum_{i=1}^{N} \hat{\theta}_i(k-1) \zeta(k-N+i-1) \tag{10}
$$

**注意**：书中**不对 $u(k)$ 施加限幅**。

---

## 5. 闭环分析

### 5.1 理想情况（$G = G_0$, $\Delta_m = 0$）

此时 $\zeta(k) = d(k)$。闭环输出：

$$
y(k) = G_0(z)[u(k)] + d(k) = -F \cdot G_0(z) K(z)[d(k)] + d(k)
$$

频域：$Y(e^{j\omega}) = \big(1 - F \cdot G_0(e^{j\omega}) \cdot K(e^{j\omega})\big) D(e^{j\omega})$

完全抑制条件：

$$
F \cdot G_0(e^{j\omega}) \cdot K(e^{j\omega}) = 1 \quad \text{（在所有扰动频率处）} \tag{11}
$$

### 5.2 鲁棒情况（$G = G_0(1+\Delta_m)$）

$$
\zeta(k) = d(k) + G_0(z)\Delta_m(z)[u(k)]
$$

闭环敏感度：

$$
H(z) = \frac{1 - F G_0 K}{1 + F G_0 \Delta_m K}
$$

鲁棒稳定性（小增益定理）：$\|F G_0 \Delta_m K\|_\infty < 1$。

---

## 6. 完整仿真参数

| 参数 | 符号 | 值 | 来源 |
|------|------|-----|------|
| 采样周期 | $T_s$ | $1/480$ s | book §6.6.1 |
| FIR 阶数 | $N$ | $50$ | Table 4.1 |
| 预补偿器 | $F(z)$ | $100$（静态） | Table 4.1 |
| 归一化增益 | $\gamma_0$ | $1$ | Eq. (6) |
| 初始协方差 | $P(0)$ | $20 I_{50}$ | Table 4.1 |
| 参数初值 | $\hat{\theta}(0)$ | $0$ | Table 4.1 |
| 投影界 | $\theta_{\max}$ | 建议 $\ge 100$（避免约束收敛方向） | — |
| 协方差重置 $\beta_0$ | — | $20$ | Eq. (6.61) |
| 协方差重置 $\beta_1$ | — | $0.5$ | Eq. (6.61) |
| 仿真时长 | $T_{\text{end}}$ | $50$ s | Fig. 6.8 |

---

## 7. 复现诊断：单频 vs 多频的性能断裂

### 7.1 实验结果

| 场景 | 方法 | 抑制 | 备注 |
|------|------|------|------|
| 单频 25 rad/s | 离线 LS | **−28.6 dB** ✅ | — |
| 单频 25 rad/s | 在线 RLS | **−27.5 dB** ✅ | — |
| 双频 25+225 | 离线 LS | −3.7 dB ❌ | RCND = $2.6\times 10^{-19}$ |
| 双频 25+225 | 在线 RLS | −3.7 dB ❌ | — |
| 四频 (book 全场景) | 在线 RLS | −2.5 dB ❌ | — |
| 四频 | 凸优化 MOSEK | **−12.8 dB** | 频域约束优化 |

### 7.2 根因分析

**问题不在自适应算法，而在参数化模型与被控对象的匹配。**

$G_0(z)$ 在扰动频率处增益极低：

| $\omega$ (rad/s) | $|G_0(e^{j\omega T_s})|$ | 所需 $|K|$ (F=100) |
|-------------------|--------------------------|---------------------|
| 25 | 0.0088 (−41 dB) | 1.14 |
| 225 | 0.0013 (−58 dB) | 8.0 |

FIR 基函数 $\alpha(z) = [z^{-N}, \dots, z^{-1}]$ 在 $\omega T_s = 0.052$ rad（25 rad/s）处，$N=50$ 个列向量跨越的相位范围仅 $[-2.6, -0.05]$ rad。如此窄的相位跨度使得 $\Phi$ 矩阵的列近线性相关。

**单频**：仅需匹配一个复数增益，近奇异仍可行。**多频**：需要同时匹配多个复数增益，$\Phi^\top\Phi$ 的条件数爆炸（$10^{-19}$），最小二乘解退化为数值噪声。

### 7.3 为什么 MOSEK 可行

MOSEK 在**频域**直接求解约束优化：

$$
\min_{\theta} \max_{\omega \in [0,\pi]} |G_0(e^{j\omega})F(e^{j\omega},\theta)|
$$

绕过了时域 $\Phi$ 矩阵的条件数问题。频域约束天然分离了不同频率分量。

### 7.4 复现建议

1. **单频场景**：RLS 自适应方案完全可行（−28 dB），适合演示基本收敛行为
2. **多频场景**：需采用频域约束优化（MOSEK/CVX）或切换至 2017 JVC 的 $F(s)$ 频谱平坦化结构
3. **协方差重置**必不可少——即使单频，无此机制则 30s 新扰动无法适应
4. **$u(k)$ 不应限幅**——书中未施加限幅

---

## 8. 与 2017 JVC 方案的关系

| | 2015 TAC（本文） | 2017 JVC |
|---|---|---|
| 系统域 | 纯离散 $z$ | 连续 $s$ → 离散化 |
| 控制器 | $u = -F \cdot \theta^\top \alpha(z)[\zeta]$ | $u = -F(s)[\theta^\top \Lambda(s)[z]]$ |
| 基函数 | $\alpha(z) = [z^{-N}, \dots, z^{-1}]$（延迟） | $\Lambda(s) = [\lambda^N/(s+\lambda)^N, \dots]$（Laguerre） |
| 预补偿器 | $F=100$（静态增益） | $F(s) = \kappa_0 \frac{2\alpha^2(s^2+s+1.25)}{(s+\alpha)^2(s+0.2)}$（频谱平坦化） |
| 条件数 | 极差（多频） | 良好（$F(s)$ 整形后 $\Phi$ 列不相关） |
| 适用 | 单频/离线优化 | 多频自适应 |

---

## 9. 代码实现要点

```matlab
% 关键参数
N = 50;  F_gain = 100;  gamma0 = 1;  P0 = 20;
beta0 = 20;  beta1 = 0.5;  theta_max = 100;

% 仿真循环每步
zeta_k = y(k) - G0_filter(u_prev);                % (1)
w_vec  = flipud(buffer_zeta);                      % α(z)[ζ]
for i = 1:N
    phi_vec(i) = F_gain * G0_filter(zeta_delayed); % (4)
end
pe = zeta_k - theta'*phi_vec;                       % 预测误差
m2 = 1 + gamma0*(phi_vec'*phi_vec);                 % (6)
eps = pe / m2;                                      % (7)
P = P - (P*phi)*(P*phi)'/(m2 + phi'*P*phi);        % (8)
if trace(P)/N < beta1, P = beta0*eye(N); end        % (6.61)
theta = proj(theta + P*eps*phi, theta_max);         % (9)
u(k) = -F_gain * (theta' * w_vec);                  % (10) — 不限幅!
```

---

## 参考资料

- Jafari, Ioannou, Fitzpatrick, Wang. "Robustness and Performance of Adaptive Suppression of Unknown Periodic Disturbances." *IEEE TAC*, 60(8):2166–2171, 2015.
- 本书第 6.6.1 节：Discrete-Time SISO Robust Adaptive Disturbance Attenuation
- 本书 Eq. (6.13)–(6.21), Eq. (6.61), Figure 6.8
