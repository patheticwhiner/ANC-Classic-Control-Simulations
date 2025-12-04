## LQGdemo2 仿真实验详解

文件：`LQGdemo2.m`

本说明文档基于仓库中的 `AboutLQG.md` 并对 `LQGdemo2.m` 的理论基础、数学推导、实现细节与仿真步骤做详尽说明，便于读者理解、复现与改进该实验。

## 一、概述

本程序将被控对象（次级通路）与带限噪声生成器（干扰模型）耦合为一个增广离散时间状态空间系统，通过 LQG 方法（LQR + Kalman）设计主动噪声控制（ANC）。程序使用 MATLAB 内置的 `lqg` 函数一步到位设计控制器，也实现了逐步仿真验证闭环性能。

主要步骤：

- 构造原系统（plant）状态空间模型 $(A_f,B_f,C_f)$；
- 构造干扰/外系统（exosystem 或 bandlimited noise）模型 $(A_w,B_w,C_w)$；
- 将两者增广为整体系统并建立控制与干扰输入矩阵；
- 选择 LQ 性能权重矩阵 $Q,R$ 以及过程/测量噪声协方差 $Q_n,R_n$；
- 使用 `lqg` 设计 LQG 控制器，提取状态反馈 $K$ 和卡尔曼增益 $L$；
- 运行逐步仿真，记录并绘图展示控制信号、反噪声与残余噪声。

## 二、数学背景与关键公式

1) 系统与噪声模型（离散时间）：

$$
\begin{aligned}
 x_{sys}(k+1) & = A_f x_{sys}(k) + B_f u(k) \\
 x_w(k+1) & = A_w x_w(k) + B_w e(k) \\
 y(k) & = C_f x_{sys}(k) + C_w x_w(k) + v(k)
\end{aligned}
$$

其中 $e(k)$ 为白噪声驱动干扰过程（外系统过程噪声），$v(k)$ 为测量噪声。

将二者增广，记增广状态 $\tilde{x} = [x_{sys}; x_w]$，则增广系统为：

$$
\tilde{x}(k+1) = \underbrace{\begin{bmatrix}A_f & 0 \\
 0 & A_w\end{bmatrix}}_{\tilde{A}} \tilde{x}(k) + \underbrace{\begin{bmatrix}B_f \\
 0\end{bmatrix}}_{\tilde{B}} u(k) + \underbrace{\begin{bmatrix}0 \\
 B_w\end{bmatrix}}_{\tilde{G}} e(k)
$$

输出：

$$
y(k)=\underbrace{\begin{bmatrix}C_f & C_w\end{bmatrix}}_{\tilde{C}}\tilde{x}(k)+v(k)
$$

2) LQR（离散时间）

我们用离散 LQR 来设计最优状态反馈，使性能指标

$$J=\sum_{k=0}^{\infty} \left( x_k^\top Q x_k + u_k^\top R u_k \right)$$

最小化。常见选择 $Q=C^\top C$ 用以直接惩罚输出能量，$R$ 为控制能量权重（正定）。离散 LQR 可由 `dlqr` 求解，得到反馈律 $u(k)=-Kx(k)$。

3) Kalman 滤波器（离散时间 LQE）

假设过程噪声 $w(k)$（这里用于增广模型的 $e(k)$）和测量噪声 $v(k)$ 为零均值白噪声：

$$\mathbb{E}[w(k)w(k)^\top]=Q_n,\quad \mathbb{E}[v(k)v(k)^\top]=R_n,\quad \mathbb{E}[w(k)v(k)^\top]=0$$

离散稳态卡尔曼滤波器增益 $L$ 可由 `dlqe` 或 `kalman` 计算，状态估计递推：

$$\hat{x}(k+1)=\tilde{A}\hat{x}(k)+\tilde{B}u(k)+L\big(y(k)-\tilde{C}\hat{x}(k)\big)$$

4) LQG（分离性）

分离原理表明，当系统可控可观时，最优控制器可由 LQR 与 Kalman 滤波器串联组成：控制 $u=-K\hat{x}$，其中 $K$ 来自 LQR，$L$ 来自 Kalman。MATLAB 的 `lqg` 可在满足一定假设下一步完成这两个设计，返回的结构包含 `K` 与 `L`。

注意：LQG 的鲁棒性一般不如 H_\infty，当模型误差较大时需谨慎使用，可结合 LTR 或鲁棒控制方法改进。

## 三、`LQGdemo2.m` 代码实现详解（逐段说明）

下面按 `LQGdemo2.m` 的执行顺序逐段解释代码（行号可能因编辑器而异）：

1) 初始化与加载数据

```matlab
clear; close all; clc;
load('dataset\bandlimitedNoise.mat');  % 包含 Aw, Bw, Cw
load('dataset\systemIdentification.mat');  % 包含 Af, Bf, Cf
```

说明：程序依赖两个预先生成的 .mat 文件：

- `bandlimitedNoise.mat`：外系统（干扰）状态空间矩阵 `Aw, Bw, Cw`；
- `systemIdentification.mat`：原系统（次级通路或被控对象）状态空间矩阵 `Af, Bf, Cf`。

2) 构造增广系统

```matlab
n = size(Af, 1);
p = size(Aw, 1);
A_sys = Af; B_sys = Bf; C_sys = Cf;
A_dist = Aw; B_dist = Bw; C_dist = Cw;
A = blkdiag(A_sys, A_dist);
B = [B_sys; zeros(p, size(B_sys, 2))];
G = [zeros(n, size(B_dist, 2)); B_dist];
C = [C_sys, C_dist];
```

说明：构建增广矩阵 $\tilde{A}=\mathrm{blkdiag}(A_f,A_w)$，控制输入矩阵在增广模型上只作用于原系统状态，干扰噪声由 $G$（或 \tilde{G}）作用于外系统子空间；输出矩阵 $C$ 拼接两个子系统的测量输出。

3) 设计 LQG 控制器

```matlab
Q = C' * C;  % 状态权重
R = 10e-4;   % 控制权重
Qn = 2 * eye(size(G, 1));
Rn = 1e-4 * eye(size(C, 1));
sys = ss(A, [B G], C, 0, 1);
lqg_controller = lqg(sys, blkdiag(Q,R), blkdiag(Qn,Rn));
K = lqg_controller.K;
L = lqg_controller.L;
```

要点解释：

- 这里使用 `Q=C' * C` 的目的是直接以输出能量为被惩罚项（常用技巧，尤其用于 ANC 问题）；
- 标量 `R = 10e-4` 实际等价于 $1\times 10^{-3}$，取值决定控制器对输入能量的惩罚；
- `Qn` 是“过程噪声协方差”，维度为增广输入维度（这里是 `size(G,1)`），代表我们对系统过程不确定性的假设；
- `Rn` 是测量噪声协方差，通常取较小以反映传感器精度较高的假设；
- `ss(A,[B G],C,0,1)` 构建了离散时间状态空间模型（采样时间 1）；`lqg` 的第二个与第三个参数分别是 LQ 权重以及噪声协方差矩阵（放在块对角中，以匹配系统的输入/噪声结构）。

注意：`lqg` 要求传入的系统 `sys` 将控制输入与过程噪声放在同一个输入矩阵中，即 `u` 与 `w` 被视作并列的输入。

4) 仿真初始化

```matlab
N = 2000;
x_sys = zeros(n,1); x_dist = zeros(p,1); x_hat = zeros(n+p,1);
...  % 历史记录变量初始化
```

5) 闭环仿真主循环（逐步计算）

主要流程：

- 生成干扰过程噪声：`e_dist = sqrtm(Qn)*randn(...)`；
- 生成测量噪声：`v = sqrtm(Rn)*randn(...)`；
- 更新干扰模型状态：`x_dist = A_dist * x_dist + B_dist * e_dist`；
- 计算干扰信号：`d = C_dist * x_dist`；
- 计算控制律：`u = -K * x_hat`（仅使用估计状态）；
- 更新原系统状态：`x_sys = A_sys * x_sys + B_sys * u`；
- 计算测量输出：`y = C_sys * x_sys + d + v`；
- 卡尔曼预测与校正：

```matlab
x_hat = A * x_hat + [B G] * [u; e_dist];
x_hat = x_hat + L * (y - C * x_hat);
```

注意：这里把实际的过程噪声样本 `e_dist` 传入预测方程以做仿真（真实系统中这不可观测，但仿真时可用于生成真实系统轨迹）。卡尔曼滤波器在校正时使用测量残差 `y - C*x_hat` 与增益 `L` 更新估计。

6) 结果计算与绘图

- `d_history = C_dist * x_dist_history;` 干扰信号轨迹；
- `anti_history = C_sys * x_sys_history;` 反噪声信号（由控制器产生，目标与干扰相抵消）；
- 绘制控制信号、反噪声与残余输出。

## 四、实现细节与设计选择理由

1) 为什么选择增广模型？

增广将外部窄带干扰建模为低维线性系统（exosystem）。通过估计外系统的状态，控制器可以产生相位/幅度匹配的反噪声，从而在输出处实现消除或减小干扰波形。

2) 关于权重矩阵的选择：

- $Q=C^\top C$：直接将输出能量作为被惩罚项，符合 ANC 的目标（最小化残余输出）；
- $R$ 较小：允许控制器使用较大控制能量以增强消噪效果；实际工程中应结合执行器最大力矩/能量限制调整；
- $Q_n$ 与 $R_n$：分别体现过程扰动与测量噪声的统计假设。若测量噪声被低估（取值过小），卡尔曼滤波器会过分信任测量，导致估计噪声放大或不稳定；若对过程噪声估计过小，则滤波器对模型误差反应迟钝。

3) 使用 `lqg` 与分离设计的差异：

`lqg` 为一步到位的设计函数，内部组合了 LQR 与 Kalman 的求解，但其适用场景有一定限制（参见 MATLAB 文档与 `AboutLQG.md`）。如果需要更灵活的权重、积分器、或非控制的确定性输入，应手动分别调用 `dlqr`/`lqr` 与 `kalman`/`dlqe`，然后将它们串联（或用 `lqgreg`）。

## 五、工程注意点与潜在问题（Edge cases）

列出常见的 5 个边界/异常情况及建议：

1. 模型维度错误：加载的 `Af/Aw` 维度不匹配会导致 `blkdiag` 或矩阵拼接失败，确保 `Af`、`Aw`、`Bf`、`Bw`、`Cf`、`Cw` 维度正确。
2. 噪声协方差选择失当：`Rn` 取值过小会导致滤波器数值不稳定；`Qn` 取值过大可能导致滤波器过度依赖模型噪声，估计波动。
3. 控制能量受限：现实中执行器有饱和限制，应在仿真中加入饱和处理（clipping）或将该限制作为设计约束（增大 R）。
4. 离散化/采样问题：文件中 `ss(...,1)` 指定采样时间 1（单位与数据一致性很重要）；若原模型为连续时间，需要先准确离散化。
5. LQG 鲁棒性：LQG 在存在建模误差时可能不稳定，需要考虑 LTR/H_infty 或引入增益/相位裕度分析。

## 六、如何运行（在 MATLAB 中）

在 MATLAB 工作目录切换到 `demo2_LQG` 文件夹并运行脚本：

```matlab
cd('e:/Code/MATLAB_ANC/LQGbasedANC/demo2_LQG');
LQGdemo2;
```

如果 `dataset` 文件夹不在当前目录，请确保 MATLAB 路径包含 `dataset`，或者使用 `load` 时给出完整路径。

## 七、质量门（quick checks）

在完成设计/修改后请做以下快速检查：

1. 代码能否运行且无维度错误（运行 `LQGdemo2` 并观察是否报错）；
2. 绘图是否产生合理的波形：反噪声与干扰在频率/相位上应有抵消趋势；
3. 检查 `K` 与 `L` 的维度是否与增广系统一致；
4. 调整 `Q,R,Qn,Rn` 的数值并观察收敛性；
5. 如果闭环不稳定，先检查滤波器增益 `L` 与控制增益 `K` 是否来源于同一离散化采样时间。

## 八、改进建议与后续工作

1. 将 `lqg` 替换为手动组合 `dlqr` + `dlqe`/`kalman`，以便更细粒度地设置噪声矩阵与积分权重；
2. 在仿真中加入执行器饱和、延迟与量化模型，评估实际可实现性；
3. 尝试 LTR 或 H_\infty 设计以提高鲁棒性；
4. 添加参数扫描脚本，自动化搜索最佳 $Q,R,Q_n,R_n$；
5. 对生成的带限噪声与辨识数据做更多 Monte-Carlo 验证，统计控制效果的均值与方差。

## 九、参考（来自仓库与 MATLAB 文档）

- `AboutLQG.md`（仓库）
- MATLAB 文档：`lqg`, `dlqr`, `dlqe`, `kalman`, `ss`

---

文件创建者：自动生成文档（基于 `AboutLQG.md` 与 `LQGdemo2.m` 的内容）

要求覆盖检查（映射）：

- 理论背景/公式推导：Done（见第二节、第四节）；
- 代码实现介绍：Done（见第三节、五节）；
- 运行指导和质量门：Done（见第六、第七节）；
- 改进建议：Done（见第八节）。
