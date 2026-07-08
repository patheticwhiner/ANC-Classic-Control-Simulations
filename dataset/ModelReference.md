# 被控对象模型参考手册

> 本文档完整记录 `dataset/` 中所有合成理论模型的数学定义、来源文献、关键特性和使用方式。实测辨识模型（`armax_30303022` 等）参见 `sysid_models/`。

---

## 目录

1. [syn_TAC2015_3rd](#1-syn_tac2015_3rd--jafari-tac-2015-离散-avc)
2. [syn_JVC2017_3rd](#2-syn_jvc2017_3rd--jafari-jvc-2017-连续-cc)
3. [syn_JVC2017_6th](#3-syn_jvc2017_6th--jafari-jvc-2017-高阶示例)
4. [syn_Bai1997_4th](#4-syn_bai1997_4th--bai-1997-耳机-h)
5. [syn_Carmona2000_7th](#5-syn_carmona2000_7th--carmona-2000-管道-anc)
6. [syn_MassSpringDamper_2nd](#6-syn_massspringdamper_2nd--质量-弹簧-阻尼器)
7. [syn_RSTtoy_2nd](#7-syn_rsttoy_2nd--rst-教学模型)

---

## 1. syn_TAC2015_3rd — Jafari TAC 2015 离散 AVC

### 来源

Jafari, S., Ioannou, P. A., & Rudd, B. (2015). Adaptive suppression of periodic disturbances in ANC systems. *IEEE Transactions on Automatic Control*.

### 数学模型

$$G_0(z) = \frac{-0.00146(z - 0.1438)(z - 1)}{(z - 0.7096)(z^2 - 0.04369z + 0.01392)}$$

| 属性 | 值 |
|:---|---|
| 域 | 离散时间 |
| 阶数 | 3 (2 极点 + 1 零点 + 1 积分零点) |
| 采样周期 | $T_s = 1/480$ s ($f_s = 480$ Hz) |
| NMP 零点 | 1 个 ($z=1$, 积分型) |
| 极点 | $0.7096$, $0.0218 \pm 0.1157i$ |
| DC 增益 | ~$0.0017$ ($-55.4$ dB) |

### 关键特性

- **极低增益**：在典型扰动频率 (25 rad/s ≈ 4 Hz) 处 $\vert G_0\vert \approx -73$ dB
- 控制器需极高增益补偿 → H∞ 优化不加 $\Vert\theta\Vert$ 约束时解发散
- **NMP 零点 $z=1$**（积分型）在 $\omega=0$ 处引入相位反转 → 固定控制器带宽受限

### 使用方

| 脚本 | 用途 |
|:---|:---|
| `JafariTAC_DiscreteAVC.m` | TAC 2015 完整复现：H∞ 插值 + 自适应 RLS |
| `JafariJVC_DiscreteAVC.m` | JVC Example 6.1：固定 FIR + H∞ 优化 |
| `demos/Jafari2015_DiscreteAVC_Unified.m` | 统一对比：LS/H∞/RLS/论文复现 |

### 加载方式

```matlab
modelFile = fullfile('..', 'dataset', 'syn_TAC2015_3rd.mat');
load(modelFile, 'model');
G0 = model.G0;  Ts = model.Ts;  fs = model.fs;
```

---

## 2. syn_JVC2017_3rd — Jafari JVC 2017 连续 CC

### 来源

Jafari, S., & Ioannou, P. A. (2016). Robust adaptive attenuation of unknown periodic disturbances. *Journal of Vibration and Control*.

### 数学模型

$$G_0(s) = \frac{0.5(s - 0.2)}{s^2 + s + 1.25}$$

| 属性 | 值 |
|:---|---|
| 域 | 连续时间 |
| 阶数 | 3 (2 极点 + 1 零点) |
| 标称离散化 | $f_s = 10$ kHz (Tustin) |
| NMP 零点 | 1 个 ($s=0.2$, RHP) |
| 极点 | $-0.5 \pm 1.0i$ ($\omega_n \approx 1.12$, $\zeta \approx 0.45$) |
| DC 增益 | $-0.08$ ($-21.9$ dB) |

### 关键特性

- **轻阻尼谐振** ($\zeta=0.45$)：在 $\omega \approx 1$ rad/s 处有 ~$+6$ dB 峰
- **单个 NMP 零点** $s=0.2$：满足 Jafari F(s) 展平的前提假设
- **标称离散化 10 kHz**：自适应控制器在离散域设计（CC → CD via Tustin）

### 使用方

| 脚本 | 用途 |
|:---|:---|
| `JafariJVC_Continuous.m` | 连续域 F(s) + 自适应 RLS (Euler/NLMS) |
| `demos/Jafari_AdaptiveCC_Unified.m` | CC 统一对比：原始/修正/VFast |
| `demos/Jafari_AdaptiveCD_Unified.m` | CD 统一对比：α(z) 基 / Λ(z) 基 |

### 加载方式

```matlab
modelFile = fullfile('..', 'dataset', 'syn_JVC2017_3rd.mat');
load(modelFile, 'model');
G0 = model.G0;  % s-domain tf, 需 s = tf('s') 构建其他传函
```

---

## 3. syn_JVC2017_6th — Jafari JVC 2017 高阶示例

### 来源

同上 (Jafari & Ioannou, 2016 JVC)，用于展示 F(s) 展平对高阶系统的效果。

### 数学模型

$$G_0(s) = \frac{s(s^2 + 4)(s - 0.8)(s + 1.4)}{(s + 0.5)^3(s + 2)^2(s + 3)}$$

| 属性 | 值 |
|:---|---|
| 域 | 连续时间 |
| 阶数 | 6 (6 极点 + 5 零点) |
| NMP 零点 | 1 个 ($s=0.8$, RHP) |
| 极点 | $-0.5$ (三重), $-2$ (二重), $-3$ |
| 稳态 | 0 型系统（原点零点，零 DC 增益） |

### 使用方

| 脚本 | 用途 |
|:---|:---|
| `JafariJVC_Continuous.m` (§0 示例) | F(s) 设计教学演示 |

---

## 4. syn_Bai1997_4th — Bai 1997 耳机 H∞

### 来源

Bai, M. R., & Lee, D. (1997). Implementation of an active headset by using H∞ robust control. *IEEE Transactions on Speech and Audio Processing*.

### 数学模型

$$P(s) = 0.3921 \cdot \frac{(s + 3.0841)(s - 1.0320)(s + 0.4387)(s - 0.0034)}{(s - 0.6612 \pm 0.3483i)(s + 0.4426 \pm 0.3324i)}$$

| 属性 | 值 |
|:---|---|
| 域 | 连续时间 |
| 阶数 | 4 (4 极点 + 4 零点) |
| 标称离散化 | $f_s = 4$ kHz (Tustin) |
| RHP 零点 | 2 个 ($s=1.032$, $s=0.0034$) |
| RHP 极点 | 2 个 ($s=0.6612 \pm 0.3483i$) — **开环不稳定** |
| LHP 极点 | 2 个 ($s=-0.4426 \pm 0.3324i$) |

### 关键特性

- **开环不稳定** + **NMP 零点** → 经典 H∞ 混合灵敏度 benchmark
- H∞ 权重设计：$W_1$ 低频性能 (DC~400 Hz), $W_3$ 高频鲁棒性 (>1 kHz)
- 论文 Table II 提供 5 阶参考控制器

### 使用方

| 脚本 | 用途 |
|:---|:---|
| `Bai1997_Hinf.m` | H∞ 综合 + 论文控制器验证 + 带限白噪声测试 |
| `demos/Bai1997_MultiSine.m` | H∞ 多正弦扰动测试 |

### 加载方式

```matlab
modelFile = fullfile('..', 'dataset', 'syn_Bai1997_4th.mat');
load(modelFile, 'model');
P_zp = model.G0_zpk;   % zpk 形式 (推荐, 保留零极点语义)
P_tf = model.G0_tf;    % tf 形式
```

---

## 5. syn_Carmona2000_7th — Carmona 2000 管道 ANC

### 来源

Carmona, R., & Alvarado, V. M. (2000). Active noise control of a duct using robust control theory. *ASME J. of Vibration and Acoustics*.

### 数学模型

$$G_0(z) = z^{-6} \cdot \frac{B(z^{-1})}{A(z^{-1})}$$

| 属性 | 值 |
|:---|---|
| 域 | 离散时间 |
| 阶数 | 7 极点 + 13 零点 (含延迟) |
| 采样频率 | $F_s = 2000$ Hz |
| 延迟 | $d = 6$ 采样点 |
| $A(z^{-1})$ | $1 -1.3941z^{-1} -0.0389z^{-2} +1.2131z^{-3} -1.1895z^{-4} +0.0430z^{-5} +1.0517z^{-6} -0.6267z^{-7}$ |

$$B(z^{-1}) = \begin{aligned} &0.0304 +0.0709z^{-1} -0.0947z^{-2} -0.0170z^{-3} -0.0104z^{-4} -0.0787z^{-5} \\ &+0.0414z^{-6} +0.0380z^{-7} +0.0250z^{-8} +0.0366z^{-9} \\ &-0.0584z^{-10} +0.0540z^{-11} -0.0862z^{-12} -0.6267z^{-13} \end{aligned}$$

### 关键特性

- 声学管道典型模型：长延迟 + 高阶 FIR 分子 → 非最小相位特性
- 极点接近单位圆（$\max|p|=0.97$）
- 设计频率：100 Hz 和 200 Hz

### 使用方

| 脚本 | 用途 |
|:---|:---|
| `demo1_RST/Carmona2000.m` | RST 极点配置 + 灵敏度函数塑形 |
| `demo1_RST/Landau2005.m` | 自适应极点配置 (继承 Carmona 模型) |

### 加载方式

```matlab
modelFile = fullfile('..', 'dataset', 'syn_Carmona2000_7th.mat');
load(modelFile, 'model');
A_coeffs = model.A_coeffs;  B_coeffs = model.B_coeffs;
d = model.d_delay;  Fs = model.fs;  Ts = model.Ts;
```

---

## 6. syn_MassSpringDamper_2nd — 质量-弹簧-阻尼器

### 来源

标准控制理论教材模型。

### 数学模型

$$\dot{x} = \begin{bmatrix} 0 & 1 \\ -k/m & -b/m \end{bmatrix} x + \begin{bmatrix} 0 \\ 1/m \end{bmatrix} u, \quad y = \begin{bmatrix} 1 & 0 \\ 0 & 1 \end{bmatrix} x$$

| 属性 | 值 |
|:---|---|
| 域 | 连续时间 (状态空间) |
| 阶数 | 2 |
| 参数 | $m=1$ kg, $k=1$ N/m, $b=0.5$ N·s/m |
| 固有频率 | $\omega_n = \sqrt{k/m} = 1$ rad/s |
| 阻尼比 | $\zeta = b/(2\sqrt{mk}) = 0.25$ |

### 使用方

| 脚本 | 用途 |
|:---|:---|
| `demo2_LQG/LQRdemo.m` | LQR 状态反馈 + Kalman 滤波器演示 |

### 加载方式

```matlab
modelFile = fullfile('..', 'dataset', 'syn_MassSpringDamper_2nd.mat');
load(modelFile, 'model');
A = model.A; B = model.B; C = model.C; D = model.D;
sys = model.sys;  % ss 对象
```

---

## 7. syn_RSTtoy_2nd — RST 教学模型

### 来源

demo4_Robust MOEA benchmark，用于验证多目标进化算法在 RST 控制器设计中的效果。

### 数学模型

$$G_0(z) = z^{-1} \cdot \frac{0.2 + 0.15z^{-1}}{1 - 1.2z^{-1} + 0.45z^{-2}}$$

| 属性 | 值 |
|:---|---|
| 域 | 离散时间 |
| 阶数 | 2 (2 极点 + 1 零点 + 1 延迟) |
| 采样周期 | $T_s = 1$ s |
| 延迟 | $d = 1$ |
| 极点 | $0.6 \pm 0.3464i$ (稳定) |
| 零点 | $-0.75$ (最小相位) |

### 关键特性

- 故意简化的教学模型：低阶、稳定、无 NMP 零点
- 延迟 $d=1$ 保证 RST 多项式方程有解
- 用于 eMOPSO/NSGA-II 等算法的基准测试

### 使用方

| 脚本 | 用途 |
|:---|:---|
| `demo4_Robust/benchmark_MOEAs.m` | MOEA RST 设计 benchmark |
| `demo4_Robust/run_eMOPSO.m` | eMOPSO RST 优化 (case `rst_toy`) |

### 加载方式

```matlab
modelFile = fullfile('..', 'dataset', 'syn_RSTtoy_2nd.mat');
load(modelFile, 'model');
B = model.B_poly;  A = model.A_poly;
d = model.d_delay;  Ts = model.Ts;
```

---

## 8. syn_Ho2020_ALE — Ho 2020 窄带 ALE 前馈 ANC

### 来源

Ho, C. Y., Shyu, K. K., Chang, C. Y., & Kuo, S. M. (2020). Efficient Narrowband Noise Cancellation System Using Adaptive Line Enhancer. *IEEE/ACM Transactions on Audio, Speech, and Language Processing*, 28.

### 数学模型

**主通路** (primary path, 扰动从噪声源传播到误差传声器):

$$P(z) = z^{-5} - 0.3z^{-6} + 0.2z^{-7}$$

**次级通路** (secondary path, 抗噪声从次级扬声器到误差传声器):

$$S(z) = z^{-2} + 1.5z^{-3} - z^{-4}$$

| 属性 | 值 |
|:---|---|
| 域 | 离散时间 (FIR) |
| P(z) 阶数 | 7 (延迟 d=5 + 2 阶 FIR) |
| S(z) 阶数 | 4 (延迟 d=2 + 2 阶 FIR) |
| 采样频率 | $f_s = 4000$ Hz |
| S(z) NMP 零点 | 是 ($z \approx -1.5 \pm 0.5i$, 模 > 1) |
| 架构 | 前馈窄带 ANC + ALE |

### 关键特性

- **前馈结构**：P(z) 和 S(z) 均建模为纯 FIR（无极点），适合 FX-LMS 及其变体
- **S(z) 有 NMP 零点**：标准 FXLMS 收敛慢，需 ALE 预处理参考信号
- **论文贡献**：用 ALE 从参考麦克风提取窄带分量，减少对次级通路建模精度的依赖
- **适合验证**：FXLMS vs ALE-FXLMS vs Jafari IMC 在 NMP 次级通路下的对比

### 模型特殊性

与仓库中其他模型不同，此模型包含 **两条通路 P(z) + S(z)**，对应前馈 ANC 架构：

```
噪声源 → P(z) → [+] → 误差Mic
                  ↑
参考Mic → ALE → W(z) → S(z) ┘
```

### 使用方

| 脚本 | 用途 |
|:---|:---|
| (待编写) | Ho 2020 ALE-FXLMS 复现与对比 |

### 加载方式

```matlab
modelFile = fullfile('..', 'dataset', 'syn_Ho2020_ALE.mat');
load(modelFile, 'model');
P = model.P;   % 主通路 tf 对象
S = model.S;   % 次级通路 tf 对象
fs = model.fs;  Ts = model.Ts;
```

---

## 附录 A：统一加载接口

所有模型均可通过 `DataManager` 统一加载：

```matlab
info = DataManager('syn_TAC2015_3rd');
G0 = info.G0;  Ts = info.Ts;  fs = info.fs;
```

列出所有可用模型：

```matlab
DataManager()
```

## 附录 B：添加新模型的步骤

1. 在 `dataset/export_all_synthetic.m` 中添加模型定义
2. 运行脚本生成 `.mat` 文件
3. 在 `DataManager.m` 的 `sources` 注册表中添加条目
4. 更新 `README.md` 和本文档的模型清单
5. 更新所有引用该模型的脚本为 `load()` 方式
