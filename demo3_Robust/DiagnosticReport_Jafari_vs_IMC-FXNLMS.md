# 诊断报告：Jafari NLMS 为何在此 plant 上劣于 IMC-FXNLMS

> **日期**：2026-07-07
> **模型**：ARMAX(16,60,4,43), $f_s=12000$ Hz, sim fit=89.9%
> **涉及脚本**：`JafariANC_SecPath_12k.m`

---

## 1. 问题陈述

在同一 IMC 架构、相同 NLMS 参数 ($\mu_0=0.05$, $N=128$, $\varepsilon=0.01$) 的条件下：

| 方法 | regressor 滤波 | 额外展平 | Part 1 固定音 | Part 2 扫频 | Part 3 FM |
|:---|:---|:---|:---|:---|:---|
| **IMC-FXNLMS** | $\phi = G_0(z)[z]$ | 无 | 胜出 | 胜出 | 胜出 |
| **Jafari NLMS** | $\phi = H(z)[z] = G_0F(z)[z]$ | $F(z)$ NMP反射+LPF | 劣于 | 劣于 | 劣于 |

**核心发现**：Jafari 方法在所有工况上劣于更简单的 IMC-FXNLMS 基线。

---

## 2. Jafari 方法的理论基础与前提假设

### 2.1 方法本质

Jafari 2015 TAC 的核心创新是为 IMC 架构引入**频谱展平滤波器** $F(z)$：

$$F(z) = \kappa_0 \cdot L(z) \cdot \frac{A(z^{-1})}{\tilde{B}(z^{-1})}$$

其中：
- $L(z) = \frac{1-\rho}{1-\rho z^{-1}}$：一阶低通，确保 $F(z)$ 因果可实现
- $A(z^{-1})$：对消 plant 极点，展平谐振峰
- $\tilde{B}(z^{-1})$：将 NMP 零点 $z_i$ ($|z_i|>1$) 反射为 $1/\bar{z}_i$ ($|1/\bar{z}_i|<1$) 的最小相位重构

滤波后的有效 plant：

$$H(z) = G_0(z) \cdot F(z) = z^{-d} \cdot \frac{B(z)}{L(z) \cdot \tilde{B}(z^{-1})}$$

由于 $|B(e^{j\omega})| = |\tilde{B}(e^{j\omega})|$（全通因子保持幅频），理论上 $|H(e^{j\omega})| \approx 1/|L(e^{j\omega})|$——在 LPF 通带内接近平坦。

### 2.2 隐含假设

Jafari 方法隐含以下前提：

| # | 假设 | 出处 |
|:--|:---|:---|
| H1 | NMP 零点数量少（1-2 个） | 全通相位累积可控 |
| H2 | Plant 阶数低（$\le$ 5 阶） | $F(z)$ 和 $H(z)$ 阶数不高，瞬态短 |
| H3 | 极点远离单位圆 | $A(z^{-1})$ 在 $F(z)$ 分子中不引入近不稳定零极点 |
| H4 | $F(z)$ 引入的额外群延迟不影响自适应跟踪 | $L(z)$ 一阶足够 |

**Jafari 原论文使用的 plant 是 3 阶、0-1 个 NMP 零点。**

---

## 3. 本 plant 违反的假设

### 3.1 22 个 NMP 零点 → 全通相位灾难

$B(z)$ 共 60 个零点，其中 22 个在单位圆外。NMP 反射 $z \to 1/\bar{z}$ 在幅度上无损，但**相位差不可消除**：

$$\angle H(e^{j\omega}) = -d\omega + \underbrace{\angle B(e^{j\omega}) - \angle \tilde{B}(e^{j\omega})}_{\text{22 个全通因子的净相位}}$$

每个 NMP 零点贡献的相位滞后随频率累积。22 个零点在大频率范围内产生远超 $180^\circ$ 的额外相位旋转，使 $H(z)$ **远非 SPR**（Strictly Positive Real）。Jafari 方法的收敛保证直接失效。

### 3.2 $\tilde{B}(z^{-1})$ 的 60 个极点 → regressor 污染

$F(z)$ 的分母含 $\tilde{B}(z^{-1})$，是 $B(z)$ 的全通等价——60 阶 FIR 多项式。$\tilde{B}(z^{-1})$ 的 59 个零点成为 $F(z)$ 的**极点**。当 plant 的 NMP 零点靠近单位圆时（$\max|z_i| = 16.78$），其倒数 $1/|z_i| \approx 0.06$ 深入单位圆内——这本身不是问题。

但反过来：plant 的 MP 零点可能**非常靠近单位圆**。$B(z)$ 有 38 个 MP 零点分布在单位圆内外。B̃ 将 NMP 零点替换为其倒数，倒数的模 $<1$，但如果原始 MP 零点中有模 $\approx 1$ 的，B̃ 的极点（即 B̃ 的零点模的倒数关系？不，B̃ 本身是多项式，其"极点"体现在 F(z) 分母中 B̃(z⁻¹) 的系数）……

实际上真正的问题是：**$\tilde{B}(z^{-1})$ 作为 60 阶 FIR 多项式，在 F(z) 的分母中产生 60 个极点**。对于单频扰动，这些极点被持续激励，产生长时 ringing，污染 regressor $\phi = H(z)[z]$。

### 3.3 Plant 极点在 $|p|=0.97$ → $F(z)$ 中的对消不干净

$A(z)$ 的模最大极点在 0.9696。$F(z)$ 将其放入分子（$A(z^{-1})$），与 $G_0(z)$ 分母中的 $A(z)$ 试图对消形成 $H(z)$。但在有限精度和时变信号下，对消不完全。$A(z)$ 仍出现在代码实现的 $H(z)$ 中（见 §4），引入额外 16 个极点。

### 3.4 阶数爆炸

| 滤波器 | Jafari 原论文（典型） | 本 plant |
|:---|:---|:---|
| $G_0(z)$ 阶数 | 3 | 16 (A) + 60 (B) = **76** |
| $F(z)$ 阶数 | ~5 | 16 (num) + **61** (den) |
| $H(z)$ 阶数（代码实现） | ~8 | **119** (num) + **77** (den) |
| $H(z)$ 阶数（理论化简后） | ~5 | **60** + delay（B̃ 的阶数） |

即使按理论化简，$H(z)$ 仍有 60+ 阶——远超 NLMS 的跟踪能力。

---

## 4. 代码实现问题：$A(z)$ 未消去

当前 `JafariANC_SecPath_12k.m` 中 H(z) 的计算：

```matlab
H_num = conv([zeros(1,d), B], F_num);   % F_num = A → z^{-d}·B·A
H_den = conv(A, F_den);                  % F_den = conv(L_den, B_tilde) → A·L·B̃
```

这保留了 $A(z)$ 在分子和分母中，名义上应彼此对消，但实际上：
1. 增加 16 个冗余极点，$H(z)$ 阶数从 61 → 77
2. 每个冗余极点都在 NLMS 的 regressor 滤波中产生瞬态
3. 扫频/FM 工况下频率连续变化 → 极点持续被重激励 → regressor 永久不干净

**正确做法**：直接消去 A(z)：
```matlab
H_reduced_num = [zeros(1,d), B];          % z^{-d}·B(z), 103 阶
H_reduced_den = conv(L_den, B_tilde);     % L(z)·B̃(z⁻¹), 61 阶
```

这将 $H(z)$ 从 77 阶降到 61 阶，但**不能解决根本问题**——61 阶仍然是 regressor 的严重污染源。

---

## 5. 为什么 IMC-FXNLMS 胜出

IMC-FXNLMS 仅用 $G_0(z)$ 滤波 $z$：

```
φ = G₀(z)[z_delays]
θ 学习: θ'·G₀(z)[z] ≈ z(k)
u = -θ'·z_delays
```

三条优势：

1. **阶数低**：$G_0(z)$ 的 regressor 滤波仅 16 (A) + 43 (delay) + 60 (B) = 103 阶 FIR 等效——虽然不低，但对比 H(z) 的 77 阶 IIR，收敛所需的瞬态要短。

2. **不引入人工极点**：$G_0(z)$ 的极点就是物理声学管道的真实极点。$F(z)$ 引入的 B̃ 极点是**人工构造**的，没有物理意义，在时变工况下持续产生 ringing。

3. **128 抽头 FIR 学到的是 plant 在该频率的增益/相位补偿**：对于单频，NLMS 只需 2 个有效抽头。额外的 126 个抽头用于拟合 plant 在邻频的响应。因为 $G_0(z)$ 的 regressor 没有 $F(z)$ 的额外相位旋转，FIR 更容易找到有效解。

**结论**：当 plant 足够复杂时，**不做展平比做坏的展平更好**。

---

## 6. 实验证据

### 6.1 三工况抑制量汇总

| 工况 | 固定 FIR | Jafari NLMS | IMC-FXNLMS | $\Delta$ (J−FX) |
|:---|:---|:---|:---|:---|
| Part 1: 固定 265 Hz | ~4 dB | 劣 | **胜** | 负值 |
| Part 2: 扫频 230→300 Hz | 随频率衰减 | 劣 | **胜** | 负值 |
| Part 3: FM $\pm$30 Hz | 偏离时衰减 | 劣 | **胜** | 负值 |

IMC-FXNLMS 在所有工况、所有时间段均优于或等于 Jafari NLMS。

### 6.2 定性观察

- **Jafari NLMS**：收敛慢，regressor 受 $H(z)$ 76 阶 IIR 瞬态污染，$\theta$ 权值持续振荡
- **IMC-FXNLMS**：收敛快，稳态抑制深，$\theta$ 权值平滑

---

## 7. 建议

### 7.1 短期——在此 plant 上放弃 Jafari 方法

对于 ARMAX(16,60,4,43)，IMC-FXNLMS 是更合适的自适应基线。当前 `JafariANC_SecPath_12k.m` 应保留 IMC-FXNLMS 作为主要对比对象，对 Jafari NLMS 添加注释说明其不适用性。

### 7.2 中期——换用符合 Jafari 假设的 plant

若需展示 Jafari 方法的理论优势，应使用：
- 低阶 plant（3-5 阶）
- 0-2 个 NMP 零点
- 极点 $|p| < 0.9$
- 原论文的 demo plant（TAC 2015 的 G₀(z)）或类似模型

### 7.3 长期——探索中间方案

对于复杂 plant，可考虑：
- **部分展平**：仅对消最显著的 NMP 零点（而非全部 22 个），减少 $F(z)$ 阶数
- **频段选择展平**：仅在设计频率邻域做展平，而非全频段
- **FIR 近似 F(z)**：用有限冲激响应近似 $F(z)$，避免 IIR 瞬态

---

## 8. 与方法论相关的关键教训

| 教训 | 说明 |
|:---|:---|
| **验证假设** | Jafari 方法隐含的「低阶少 NMP」假设在真实声学管道上不成立 |
| **退化基线比复杂方法更可靠** | IMC-FXNLMS 只有 $G_0$ 滤波，但因为没有引入人工极点，反而更好 |
| **「展平」是一把双刃剑** | 展平幅度响应可能以恶化相位响应为代价；对 22 NMP 的 plant，相位代价 > 幅度收益 |
| **阶数即代价** | $F(z)$ 的每一次阶数增加都在自适应 regressor 中引入延迟和瞬态——这两者抵消了幅度展平带来的理论收益 |

---

> **参考文献**
> - Jafari, S., Ioannou, P. A., & Rudd, B. (2015). Adaptive suppression of periodic disturbances in ANC systems. *IEEE TAC*.
> - Kuo, S. M., & Morgan, D. R. (1996). *Active Noise Control Systems*. Wiley.
