# JafariJVC_Discrete.m 调试报告：四重缺陷的根因分析与修复

> **读者**：能独立推导 RLS/LMS、能直接从 $G(z)=b(z)/a(z)$ 写出时域递推公式、已运行过 `JafariTAC_DiscreteAVC.m` 和 `JafariJVC_Continuous.m`。
>
> **目标**：(1) 理解每个 bug 的数学根因——为什么"看起来能跑"的代码结果不对；(2) 学习"从现象到根因"的系统调试方法论。

---

## §0 知识基线

| 领域 | 你的状态 | 本文档策略 |
|:---|:---|:---|
| RLS/LMS 推导 | 能独立推导 | 不展开 RLS 公式，聚焦适配问题 |
| 离散差分方程 | 能直接写出递推式 | 展示代码与公式的维度/结构偏差 |
| 项目脚本 | 已跑过 Continuous、TAC 等多个版本 | 用对比表展示离散域独有的问题 |
| 调试方法论 | 需要系统化 | §6 独立成章 |

**动机锚定**：你有能力独立写出 ANC 自适应仿真，但四个 bug 在代码层面的表现与其数学层面的错误之间隔着一层"看起来能跑通"的假象。Bug 3（控制律信号混淆）和 Bug 4（参数尺度）尤其隐蔽——代码不报错，收敛曲线看起来"在运行"，但 `‖θ‖` 始终被夹在投影边界上，抑制 ~0 dB。

**关联文件状态**：`JafariJVC_Continuous.m` 的代码修改（频率单位修正、Λ 基函数 + Euler/RLS 双段仿真）为本次会话前已存在的修改，在修复了文本编码问题后已恢复。本文档聚焦 `JafariJVC_Discrete.m` 的四重缺陷。

---

## §1 系统模型与控制器结构

### 1.1 被控对象与扰动

$$
G_0(z) = \frac{-0.00146(z-0.1438)(z-1)}{(z-0.7096)(z^2 - 0.04369z + 0.01392)}, \quad T_s = \frac{1}{480}\text{s}
$$

扰动：$d(k) = \sin(\omega_0 k) + v(k)$，$\omega_0 = 0.0521$ rad/sample（≈ 4 Hz）。

### 1.2 控制器与自适应律

$$
\begin{aligned}
z(k) &= y(k) - G_0(z)[u(k-1)] \\[4pt]
u(k) &= -F(z)\Big[K(z,\theta)[z(k)]\Big], \quad K(z,\theta) = \sum_{i=1}^{N} \theta_i z^{-i} \\[4pt]
\phi_i(k) &= G_0(z)F(z)[z(k-i)], \quad i = 1,\dots,N
\end{aligned}
$$

本示例中 $F(z) = 1$。归一化 RLS 递推：

$$
\begin{aligned}
m_s^2(k) &= 1 + \gamma_0\,\phi^T\phi \\[4pt]
\varepsilon(k) &= \frac{z(k) - \theta^T(k-1)\phi(k)}{m_s^2(k)} \\[4pt]
K(k) &= \frac{P(k-1)\phi(k)}{m_s^2(k) + \phi^T(k)P(k-1)\phi(k)} \\[4pt]
\theta(k) &= \text{Proj}_{\|\theta\|\le\theta_{\max}}\Big(\theta(k-1) + K(k)\,\varepsilon(k)\Big) \\[4pt]
P(k) &= P(k-1) - K(k)\big(P(k-1)\phi(k)\big)^T
\end{aligned}
$$

**关键区分**——控制律用原始延迟 $z(k-i)$，非回归向量 $\phi(k)$：

$$
u(k) = -F(z)\Big[\theta^T(k)\,z(k-\{1:N\})\Big]
$$

---

## §2 Bug 1：差分方程维度不匹配（第 180 行）

### 数学定义

$G_0(z)=\frac{b(z)}{a(z)}$ 的时域递推：

$$
y(k) = \frac{-\sum_{j=2}^{n_a} a_j\,y(k-j+1) + \sum_{j=1}^{n_b} b_j\,u(k-j+1)}{a_1}
$$

右侧是标量。**每个内积的左右操作数维度必须匹配。**

### 错误

```matlab
% 原代码
anti_t(k) = (-a(2:end)*anti_t(k-1:-1:k-length(a)+1)' + b*u_t)/a(1);
%                                                             ^^^^^
%                            b: 1×3 行向量, u_t: 标量 → b*u_t: 1×3 向量, 无法赋值给标量
```

同一文件中固定控制器部分（第 133–134 行）的**正确写法**就在不远处：

```matlab
% 正确（固定控制器仿真）
antinoise(k) = (-a(2:end)*antinoise(k-1:-1:k-length(a)+1)' + ...
                 b*u_fixed(k-1:-1:k-length(b))')/a(1);
```

### 修复

```matlab
u_hist = zeros(1, Nsim);      % 新增
u_hist(k) = u_t;
anti_t(k) = (-a(2:end)*anti_t(k-1:-1:k-length(a)+1)' + ...
              b*u_hist(k-1:-1:k-length(b))')/a(1);
```

**教训**：MATLAB 中 `row_vector * scalar` 合法且不报警，但结果是向量。差分方程实现时，MA 部分必须与等长的过去输入向量做内积。

---

## §3 Bug 2：`c2d` 用于已离散模型（第 240 行）

### 错误

```matlab
Fd = c2d(F, Ts, 'tustin');   % ❌ F = tf(1,1,Ts) 已是离散模型
Gd = c2d(G0, Ts, 'tustin');  % ❌ G0 = tf('z',Ts) 已是离散模型
```

`c2d`（continuous-to-discrete）拒绝离散模型输入。

### 修复

```matlab
Fd = F;    % 离散模型直接赋值
Gd = G0;
```

**教训**：`tf('z', Ts)` → 离散模型；`tf('s')` → 连续模型。从连续版（如 `JafariJVC_Continuous.m`）迁移代码到离散版时，`c2d`/`d2c` 的适用条件反转。

---

## §4 Bug 3：控制律信号混淆（第 313 行）

四个 bug 中**最隐蔽**的一个——代码不报错、RLS 收敛曲线看起来正常，但自适应完全失效。

### 数学上正确的关系

| 用途 | 表达式 | 信号路径 |
|:---|:---|:---|
| **回归向量**（RLS 辨识） | $\phi_i(k) = G_0(z)F(z)[z(k-i)]$ | 经过 G₀·F |
| **控制律**（生成 u） | $u(k) = -F(z)[\theta^T z(k-\{1:N\})]$ | $z$ 原始延迟，不经 G₀ |

辨识用 $\phi$ 是因为 Swapping Lemma：$z = \theta^{*T}\phi + \eta$。但控制器 $K(z,\theta) = \sum \theta_i z^{-i}$ 是 FIR 滤波器，必须直接作用于 $z$。

### 错误

```matlab
% 原代码
Kz_k = theta_adapt' * phi_reg;          % φ = G₀·F[z_delays]
[u_adapt(k), ~] = filter(numFd, denFd, -Kz_k, F_ctrl_zf);
```

当 $F(z)=1$：原代码实现 $u = -\theta^T \cdot G_0[z_{\text{delays}}]$，正确应实现 $u = -\theta^T \cdot z_{\text{delays}}$。

### 修复

```matlab
% 在回归向量循环中同步收集原始 z 延迟
for i = 1:N_adapt
    z_delay = 0;
    if k > i, z_delay = z_hist(k - i); end
    z_delays(i) = z_delay;                               % 保存原始延迟
    [f_out, ~] = filter(numFd, denFd, z_delay, ...);     % 继续计算 φ
    [phi_reg(i), ~] = filter(numGd, denGd, f_out, ...);
end
% ...
Kz_k = theta_adapt' * z_delays;   % ← 用 z_delays，非 phi_reg
```

### 为什么这是一个跨文件的结构性错误

`About_CC_Derivation.md`（第 116 行）记录了连续域 Laguerre 基版本（`Jafari_AdaptiveCC.m`）的完全相同错误：

> **关键**：控制律中使用 Λ(z)[z(k)]（仅经过基函数滤波），而非 φ(k)（还经过了 G₀F）。这是原版代码的核心 bug。

同一个结构错误以几乎相同的形式出现在离散域 FIR 版本中——两个文件由不同的控制器参数化方式写就，但犯了同一类信号路径混淆。

---

## §5 Bug 4：参数投影与协方差尺度（第 237–239 行）

### 现象

修复 Bug 3 后抑制仍 ~0 dB，`‖θ‖` 精确卡在投影边界 `θ_max=10` 上。

### 根因：|G₀| 在 ω₀ 处极低

$$
|G_0(e^{j\omega_0})| \approx 2.3 \times 10^{-4} \quad (\approx -73\text{ dB})
$$

G₀ 的分子含 $(z-1)$ 因子——在 DC 处有一个零点，严重衰减低频。使 $G_0(e^{j\omega_0})\cdot K(e^{j\omega_0})=1$ 需要：

$$
|K(e^{j\omega_0})| = \frac{1}{|G_0(e^{j\omega_0})|} \approx 4380
$$

对于 FIR 滤波器 $K(z,\theta)=\sum_{i=1}^{30}\theta_i z^{-i}$，$|K| \approx \frac{\|\theta\|\sqrt{N}}{2} \approx 2.74\|\theta\|$，因此需要 $\|\theta\| \gtrsim 1600$。

固定 H∞ 控制器实测 $\|\theta\| = 695,957$——而原代码 `θ_max=10` 比所需值小了 **~70,000 倍**。

$P(0)=500$ 使 RLS 每步更新量 $\Delta\theta_i \approx 0.11$，即使无投影也需 ~14,500 步（~30s）才能收敛到 $\|\theta\|=1600$。

### 与连续域版本的对比——为什么同一个参数在两个版本中表现截然不同

| 参数 | Continuous (JafariJVC_Continuous) | Discrete (JafariJVC_Discrete) |
|:---|:---|:---|
| $G_0$ | $\frac{0.5(s-0.2)}{s^2+s+1.25}$ | $\frac{-0.00146(z-0.1438)(z-1)}{(z-0.7096)(z^2-0.04369z+0.01392)}$ |
| $\vert G_0\vert$ @ disturbance | ~0.4（−8 dB） | ~2.3×10⁻⁴（−73 dB） |
| $F$ | 频谱平坦化滤波器（精心设计） | $F(z)=1$ |
| 所需 $\|\theta\|$ | ~2.5 | ~1600–700k |
| `θ_max=10` 足够？ | ✅ | ❌ 差 160× 以上 |

**核心洞察**：直接将连续域的 `θ_max=10` 搬到离散域而不重新计算 $|G_0(e^{j\omega_0})|$，是根本性的工程失误。固定 H∞ 控制器的 $\|\theta\| = 695,957$ 早已暗示了所需参数尺度——只是这个信息没有被传递到自适应部分的设计中。

### 修复

```matlab
% 原代码
theta_max = 10;
P = eye(N_adapt) * 500;

% 修复后
theta_max = 1e7;              % 匹配固定控制器 ||θ|| ≈ 7×10⁵
P = eye(N_adapt) * 1e6;       % 匹配参数尺度，确保收敛速度
```

---

## §6 调试方法论

### 6.1 诊断顺序

当自适应算法「能跑但无效」时：

| 步骤 | 检查项 | 本案例发现 |
|:---|:---|:---|
| **① 算法结构** | 三条信号路径是否与公式一致 | 控制律用了 φ 而非 z_delays |
| **② 参数范围** | `‖θ‖` 是否被投影截断 | θ_max=10 vs. 需要 1600+ |
| **③ 数值实现** | MATLAB API 语义、向量维度 | b\*u_t 维度不匹配、c2d 误用 |
| **④ 收敛速度** | P(0) 是否允许在仿真时长内收敛 | P(0)=500 过小 |

### 6.2 关键诊断信号

```
‖θ‖ = 10.000    ← 精确卡在投影边界 → θ_max 太小
tr(P) = 1.4e4    ← 从初始值几乎不下降 → 收敛未启动
控制RMS = 0.0    ← 控制信号完全为零 → 结构错误或参数截断
```

**法则**：如果 `‖θ‖` 在运行几秒后就精确等于 `θ_max` 且不再变化，投影界几乎肯定太小。

### 6.3 对比参考法

```
Continuous  ──→  抑制 −28.9 dB, ||θ|| ≈ 2.5, θ_max=5
                         │
          G₀ 增益差 ~1700×│F 滤波器差 ~50 dB
                         │
Discrete    ──→  抑制  −0.1 dB, ||θ|| = 10 (被夹), θ_max=10
```

逐项比对 G₀、F、θ_max、P(0)，差异自然指向根因。

### 6.4 参考脚本

| 脚本 | 用途 |
|:---|:---|
| `JafariTAC_DiscreteAVC.m` | 固定控制器仿真（CVX/MOSEK），‖θ‖≈7e5 — 提供了参数尺度的上界 |
| `demos/Jafari2015_DiscreteAVC.m` | 原始参考实现，含正确的自适应结构 |
| `JafariJVC_Continuous.m` | 不同 G₀/F 的连续版本，θ_max=5 能工作 — 提供对比基准 |
| `demos/debug_loop2_rootcause.m` | 三类 RLS bug 的最小诊断脚本 |
| `About_CC_Derivation.md` | 连续域的算法推导 + 第 116 行的控制律信号混淆记录 |

---

## §7 修改汇总与工程经验

### 修改清单

| # | 行号 | 类别 | 原代码 | 修复 |
|:---|:---|:---|:---|:---|
| 1 | 174, 180–183 | 差分方程 | `b*u_t`（向量×标量） | `b*u_hist(k-1:-1:k-length(b))'` |
| 2 | 240–243 | 模型转换 | `c2d(F, Ts, 'tustin')` | `F`（已是离散模型） |
| 3 | 284–288, 318 | 控制律结构 | `Kz_k = θ' * φ_reg` | `Kz_k = θ' * z_delays` |
| 4 | 237–239 | 参数尺度 | `θ_max=10, P(0)=500I` | `θ_max=1e7, P(0)=1e6I` |

### 修复后预期

| 指标 | 修复前 | 修复后（预期） |
|:---|:---|:---|
| 自适应抑制 | ~0 dB | 接近 −53.5 dB（固定 FIR 水平） |
| `‖θ‖` 稳态 | 10.0（被夹） | ~1e4–1e6 |
| `tr(P)` 终值 | 1.4e4（未收敛） | ≪ 1e4 |
| 控制 RMS | ~0 | 与固定控制器量级相当 |

### 工程经验

1. **让固定控制器的 $\|\theta\|$ 成为自适应控制器 $\theta_{\max}$ 的基准**——固定 H∞ 解已经给出了所需参数的量级上界，忽略它就是浪费已有信息。

2. **控制律和回归量用不同的信号路径**——这个错误在连续域 Laguerre 基版本和离散域 FIR 版本中出现了一模一样的模式。结构审查时，逐条画出 $z \to \phi \to \varepsilon \to \theta \to u \to y \to z$ 的完整信号流图，每个箭头标注其经过的滤波器链。

3. **MATLAB 的"看起来能跑"是调试自适应控制的最大敌人**——`row_vector * scalar` 不报错、`c2d` 对离散模型报错但错误信息可能被忽略、控制律用错信号不报错。唯一可靠的验证手段是对固定控制器和自适应控制器跑相同的闭环仿真，比较 RMS 抑制。

4. **不同域（连续 vs 离散）之间移植参数不会平移生效**——$|G_0|$ 在扰动频率处的差异可能跨越 3 个数量级（本案 −8 dB → −73 dB），直接影响 $\theta_{\max}$ 和 $P(0)$ 的合理量级。

---

*生成时间：2026-06-06* · *关联文件：`JafariJVC_Discrete.m`, `About_CC_Derivation.md`, `demos/debug_loop2_rootcause.m`*
