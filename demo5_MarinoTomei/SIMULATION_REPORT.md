# demo5_MarinoTomei 仿真总结报告

> 基于 Marino & Tomei 系列论文的自适应输出调节仿真实验全览与结论。
> 报告日期: 2026-07-08

---

## 一、仿真全景图

本目录共 **11 个 .m 仿真脚本**，覆盖 Marino & Tomei 三个研究阶段：

| 阶段            | 论文年份 | 核心问题                                     | 对应脚本                                                                  |
| :-------------- | :------- | :------------------------------------------- | :------------------------------------------------------------------------ |
| **阶段2** | 2011     | 被控对象不完全已知，外系统频率已知但参数未知 | `MarinoTomei_adaptive_regulation`, `MarinoTomei_regulator_benchmark`, `demo_simple_adaptive` |
| **阶段3** | 2016     | 被控对象未知，外系统频率和参数均未知         | `MarinoTomei_2016_adaptive_freq_est`                                    |
| **阶段3** | 2023     | 不稳定最小相位系统的高增益自适应调节         | `MarinoTomei_2023_unstable`                                             |

辅助脚本:

- **理论验证**: `demo_regulator_solution`, `demo_ss_disturbance`, `direct_regulator_solution`
- **结构分析**: `MarinoTomei_internal_model`1111

---

## 二、各仿真详细结论

### 2.1 Marino & Tomei (2011) — 自适应输出调节

**四份实现, 同一控制器, 不同求解器/基准:**

| 脚本                          | 求解器          | 额外功能          | 运行速度 |
| :---------------------------- | :-------------- | :---------------- | :------- |
| `MarinoTomei_adaptive_regulation.m (integrator='ode45')`       | ode45 (变步长)  | 参考控制器$R_c$ | ~2s      |
| `MarinoTomei_adaptive_regulation.m (integrator='euler')`       | Euler (dt=0.01) | 参考控制器$R_c$ | ~1s      |
| `MarinoTomei_regulator_benchmark.m (integrator='ode45')` | ode45           | + 调节器方程基准  | ~10s     |
| `MarinoTomei_regulator_benchmark.m (integrator='rk4')`   | RK4 (dt=0.01)   | + 调节器方程基准  | ~5s      |

#### 控制律

$$
u = -k e - \hat{\chi}_1 - \mu_1 \hat{\theta}_1 - \mu_2 \hat{\theta}_2
$$

其中控制器由三部分构成:

1. **比例误差反馈** $-k e$: 提供基本镇定
2. **自适应内模观测** $-\hat{\chi}_1$: 观测器估计外系统状态
3. **参数自适应补偿** $-\mu_1\hat{\theta}_1 - \mu_2\hat{\theta}_2$: regressor × 频率参数

#### 扰动场景 (三区间测试)

| 区间 | 时间       | $w_1(t)$              | $w_2(t)$    | 测试目的              |
| :--- | :--------- | :---------------------- | :------------ | :-------------------- |
| 1    | [0, 50)    | $\sin t + 0.2\sin 4t$ | $\sin 0.5t$ | 双频扰动 + 正弦跟踪   |
| 2    | [50, 100)  | $\sin t$              | $\sin 0.5t$ | 扰动频率成分突减      |
| 3    | [100, 150] | $\sin t$              | 0             | 纯调节 (参考信号消失) |

#### 仿真结论

1. **输出跟踪收敛性**: 系统输出 $y = x_1$ 在所有三个区间均渐近收敛到 $y_r = -w_2$, 区间切换处 ($t=50$s, $t=100$s) 出现瞬态后快速恢复。稳态误差趋近于零。
2. **参数估计收敛**:

   - $\hat{\theta}_1 \to \omega_1^2 + \omega_2^2 = 1.25$ (真值: 1.25)
   - $\hat{\theta}_2 \to \omega_1^2 \cdot \omega_2^2 = 0.25$ (真值: 0.25)
   - 参数收敛速度取决于自适应增益 $g=100$ 和持续激励条件
3. **投影算子的必要性**: 使用平滑投影算子 ($\epsilon_r=0.1$, $r_\Omega=2$) 确保参数估计保持在有界区域内, 防止参数漂移导致控制器失效。
4. **求解器对比**:

   - ode45 变步长与 Euler 定步长 ($dt=0.01$) 结果几乎一致
   - Euler 法适合嵌入式实现 (固定计算量), 但需注意数值稳定性
   - RK4 精度的 benchmark 系列额外提供与理论最优轨迹 $\Pi w$ 的偏差分析
5. **控制器与理论最优的差距**: benchmark 系列的调节器方程求解显示, 实际自适应控制器轨迹与理论最优 $\Pi w$ 间的偏差随参数收敛而减小。$\Pi$ 和 $\Gamma$ 矩阵的值依赖于当前活跃的扰动频率成分——不同区间的解存在显著差异 (Frobenius 范数 $\sim 10^{-1}$ 量级)。

---

### 2.2 Marino & Tomei (2016) — 自适应频率估计与扰动抵消

**脚本**: `MarinoTomei_2016_adaptive_freq_est.m`
**合并自**: `MarinoTomei2016/freq_estimator_core.m` + `MarinoTomei2016/run_freq_estimator.m`

提供两种运行模式:

| 模式          | 说明                                           | 状态数 |
| :------------ | :--------------------------------------------- | :----- |
| `plant`     | 估计器 + 二阶被控对象$P(s) = 1/(s+1)^2$ 闭环 | 9      |
| `estimator` | 纯频率估计器, 直接估计扰动$d(t)$             | 7      |

#### 控制问题

扰动信号 $d(t) = A_0 + \sin(\omega_1 t) + \sin(\omega_2 t)$ 的频率 $\omega_1=0.5$, $\omega_2=2.0$ **未知**, 需在线估计。

控制器: $u = -\hat{w}_0 - \hat{w}_1 - \hat{w}_3$  (抵消估计的扰动各分量)

频率参数自适应律:

$$
\dot{\hat{\theta}}_1 = \varepsilon \cdot w_2 \cdot y, \quad \dot{\hat{\theta}}_2 = -\varepsilon \cdot w_4 \cdot y
$$

#### 仿真结论

1. **频率估计精度**: 参数估计收敛到真值:

   - $\hat{\theta}_1 \to \omega_1^2 = 0.25$
   - $\hat{\theta}_2 \to \omega_2^2 = 4.00$
   - 初始估计 $\hat{\theta}_1(0)=0.3$, $\hat{\theta}_2(0)=4.2$ 均偏离真值, 验证了自适应机制的收敛能力
2. **Plant 模式 vs Estimator 模式**:

   - Plant 模式: 经过被控对象 $P(s)$ 的低通滤波效应, 输出 $y(t)$ 稳态收敛到零 (扰动被抵消)
   - Estimator 模式: 直接估计扰动 $\hat{d}(t)$, 估计误差迅速收敛, 验证了估计器本身的正确性
3. **自适应增益的权衡**:

   - $\varepsilon=0.05$ 较小时, 参数收敛慢但瞬态平滑
   - 增大 $\varepsilon$ 加速收敛但可能引起振荡
4. **与 2011 方法的本质区别**:

   - 2011: 频率已知, 自适应参数是外系统特征多项式的系数
   - 2016: 频率未知, 自适应参数直接是频率的平方, 对先验知识要求更少

---

### 2.3 Marino & Tomei (2023) — 不稳定系统的高增益自适应调节

**脚本**: `MarinoTomei_2023_unstable.m`

#### 系统模型

四个可选示例, 均为**不稳定最小相位系统**:

| 示例 | 传递函数$P(s)$                 | 相对阶 |
| :--- | :------------------------------- | :----- |
| 1    | $10/[(s+10)(s+1)(s-1)]$        | 3      |
| 2    | $10(s+2)/[(s+10)(s+1)(s-1)]$   | 2      |
| 3    | $10(s+2)^2/[(s+10)(s+1)(s-1)]$ | 1      |
| 4    | 示例 1 + 未建模扰动              | 3      |

系统特征值: $\{-10, -1, +1\}$ (右半平面极点 $+1$ 导致开环不稳定)

#### 控制器结构

采用**高增益线性补偿** + **自适应内模**结构:

- 高增益补偿: $v = -K y$, $K = \bar{\rho} - \underline{\rho} + 1$
- 频率参数自适应: $\dot{\hat{\theta}} = -\varepsilon \cdot \text{error\_signal} \cdot v$

#### 仿真结论

1. **不稳定系统的镇定**: 即使被控对象含 RHP 极点 $(s=+1)$, 高增益自适应控制器仍能镇定系统, 验证了 2023 年论文的理论结论。
2. **相对阶的影响**:

   - 相对阶越低 (示例 3, $\rho=1$), 高增益反馈更直接, 响应更快
   - 相对阶越高 (示例 1, $\rho=3$), 需更多滤波级, 响应更平缓
3. **未建模扰动的鲁棒性** (示例 4): 增加 $0.1\sin(t+0.2)$ 未建模扰动后, 控制器仍维持稳定, 但参数估计出现小幅波动。
4. **实现注意事项**: 当前实现使用 $A \backslash (B u + d)$ 计算状态, 这在数值上等价于求解 $Ax = Bu+d$ 而非真实的 ODE 积分——对不稳定系统这可能引入数值误差。**建议改用 ode45 求解完整状态方程。**

---

### 2.4 理论验证仿真

#### 2.4.1 调节器方程求解 (`demo_regulator_solution.m` + `direct_regulator_solution.m`)

**问题**: 求解调节器方程 $\Pi R = A\Pi + \frac{1}{\beta}b\Gamma + P$, $C\Pi = q$

**方法**: 将 Sylvester 方程 + 输出约束构造为联合线性方程组, 直接数值求解。

**结论**:

- 三个扰动区间对应三组不同的 $(\Pi, \Gamma)$ 解
- 解之间存在显著差异 (Frobenius 范数 $\sim 10^{-1}$), 说明理想控制输入 $\Gamma w$ 取决于当前活跃的频率成分
- 残差验证精度达到 $10^{-15}$ 量级 (机器精度), 验证了解法的正确性
- $\Pi$ 矩阵将外系统状态 $w$ 映射为理想系统状态 $x_r = \Pi w$

#### 2.4.2 扰动状态空间验证 (`demo_ss_disturbance.m`)

**问题**: 验证多正弦扰动的状态空间模型 $\dot{w} = R w$ 是否与解析表达式 $w_1(t) = \sin(\omega_1 t) + 0.2\sin(\omega_3 t)$ 一致。

**方法**: 构建 $6\times6$ 矩阵 $R$ 表示三个正弦模态, 用 RK4 求解并与直接 $\sin()$ 计算对比。

**结论**:

- 状态空间模型与解析表达式**完全一致**, 最大误差 $\sim 10^{-15}$ (数值精度)
- 区间切换处 ($t=50$s) 无额外误差, 验证了分段输出矩阵 $P_i$, $q_i$ 的正确性
- 这一验证为 benchmark 系列中"理论最优轨迹 $\Pi w$"的计算提供了正确性保证

---

### 2.5 内模结构分析 (`MarinoTomei_internal_model.m`)

**问题**: 并联 vs 级联内模结构的频域特性差异。

**方法**: 构建 $q$ 个并联/级联的谐振器组:

$$
\text{并联: } \sum_{i=1}^{q} \frac{s^2 + 2\zeta_i\omega_i s + \theta_i}{s^2 + \theta_i}, \quad \text{级联: } \prod_{i=1}^{q} \frac{s^2 + 2\zeta_i\omega_i s + \theta_i}{s^2 + \theta_i}
$$

**结论**:

- 并联结构: 各谐振峰独立、可单独调节, 适合频率间隔较大的场景
- 级联结构: 谐振峰相互作用, 整体增益更高, 带宽更宽
- Bode 图显示: 并联结构的相位特性在谐振峰之间更平缓
- FFT 频谱清晰显示谐振峰位置 $\omega_0 = \sqrt{\theta_i}$, 与设计参数一致
- **对 ANC 应用的启示**: 并联结构更适合多频噪声控制 (各频率独立调节)

---

## 三、跨仿真综合结论

### 3.1 自适应输出调节的有效性

Marino & Tomei 系列方法的核心优势在于**最小化先验假设**:

- **2011**: 只需频率数量 $m$ 和传递函数在扰动频率处的符号
- **2016**: 只需频率数量 $m$, 频率本身可以未知
- **2023**: 只需相对阶上界 $\bar{\rho}$ 和高频增益符号

### 3.2 关键的符号假设

所有 2011 年仿真依赖 **$\text{sgn}(P(j\omega_i))$ 已知**这一假设。这一假设源于模型参考自适应控制 (MRAC) 理论——反馈增益的符号决定了参数自适应律中梯度方向的正确性。符号错误会导致参数向错误方向更新, 使系统发散。

### 3.3 自适应增益的调参策略

| 参数                    | 作用           | 增大效果                         | 典型值     |
| :---------------------- | :------------- | :------------------------------- | :--------- |
| $k$                   | 误差反馈增益   | 加速输出收敛, 但过大导致振荡     | 7.1        |
| $g$ / $\varepsilon$ | 参数自适应增益 | 加速参数收敛, 但过大导致瞬态恶化 | 100 / 0.05 |
| $\epsilon_r$          | 投影平滑参数   | 越大边界过渡越平滑               | 0.1        |
| $r_\Omega$            | 参数边界       | 约束参数估计范围                 | 2.0        |

### 3.4 数值实现建议

1. **首选 ode45**: 变步长自适应, 精度和速度兼顾, `MarinoTomei_adaptive_regulation.m (integrator='ode45')` 是最佳起点
2. **定步长用于嵌入式**: `MarinoTomei_adaptive_regulation.m (integrator='euler')` 和 `bench_rk4.m` 展示固定步长可行性
3. **投影算子不可省略**: 无投影约束时参数估计可能漂移出稳定区域

### 3.5 当前仿真的局限性

1. **未覆盖 2013 年方法**: Marino & Tomei (2013) 的鲁棒自适应观测器和完整滤波变换控制器尚未实现
2. **仅 SISO 系统**: 所有仿真限于单输入单输出, 未涉及 MIMO
3. **连续域设计**: 未涉及离散化实现, 与 ANC 实际应用的数字控制器仍有差距
4. **2023 仿真的数值问题**: `MarinoTomei_2023_unstable.m` 中 $A \backslash (Bu+d)$ 不是标准 ODE 求解

---

## 四、文件索引

| 文件                                     | 行数 | 类型                   | 关联论文 |
| :--------------------------------------- | :--- | :--------------------- | :------- |
| `MarinoTomei_2016_adaptive_freq_est.m` | 280  | 自适应频率估计         | 2016     |
| `MarinoTomei_adaptive_regulation.m (integrator='ode45')`                  | 265  | 自适应输出调节         | 2011     |
| `MarinoTomei_adaptive_regulation.m (integrator='euler')`                  | 272  | 自适应输出调节 (Euler) | 2011     |
| `MarinoTomei_regulator_benchmark.m (integrator='ode45')`            | ~350 | 调节 + 基准 (ode45)    | 2011     |
| `MarinoTomei_regulator_benchmark.m (integrator='rk4')`              | ~350 | 调节 + 基准 (RK4)      | 2011     |
| `MarinoTomei_2023_unstable.m`          | 173  | 不稳定系统自适应       | 2023     |
| `MarinoTomei_internal_model.m`         | 97   | 内模结构分析           | —       |
| `demos/demo_simple_adaptive.m`         | 130  | 最简自适应 (教学)      | 2011     |
| `demos/demo_regulator_solution.m`      | 363  | 调节器方程求解         | —       |
| `demos/demo_ss_disturbance.m`          | 227  | 扰动 SS 模型验证       | —       |
| `demos/direct_regulator_solution.m`    | 123  | 调节器方程求解器       | —       |

---

## 五、参考文献

1. Marino, R. & Tomei, P. (2011). "An adaptive learning regulator for uncertain minimum phase systems with undermodeled unknown exosystems." *Automatica*, 47(4): 739–747.
2. Marino, R. & Tomei, P. (2013). "Disturbance cancellation for linear systems by adaptive internal models." *Automatica*, 49(5): 1494–1500.
3. Marino, R. & Tomei, P. (2015). "Output Regulation for Unknown Stable Linear Systems." *IEEE TAC*, 60(8): 2213–2218.
4. Marino, R. & Tomei, P. (2016). "Adaptive disturbance rejection for unknown stable linear systems." *Transactions of the Institute of Measurement and Control*, 38(6): 640–647.
5. Marino, R. & Tomei, P. (2021). "Adaptive output regulation for minimum-phase systems with unknown relative degree." *Automatica*, 130: 109670.
6. Tomei, P. & Marino, R. (2023). "Adaptive regulation for minimum phase systems with unknown relative degree and uncertain exosystems." *Automatica*, 147: 110678.
