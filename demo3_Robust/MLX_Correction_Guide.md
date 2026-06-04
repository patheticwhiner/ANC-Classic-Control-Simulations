# MLX 实时脚本修正指南

对照 `demos/Jafari_adaptive_CC_fixed.m` 和 `Jafari_adaptive_CD_fixed.m`，修改以下 MLX 脚本中的代码块。

---

## 一、JafariJVC_Continuous.mlx（7 个代码块）

### Block 1: G₀(s) 示例定义 ✅ 无需修改
教学性代码块，展示 `designLTIFilter` 的使用。保留。

### Block 2: Λ(s) 滤波器组 Bode 图示例 ⚠️ 建议更新

**问题**：N=5 示例不够充分。建议增加说明文字指出论文实际用 N=20。

**修正**：只需修改 text 段落说明（非代码块），建议在文本中补充：
> "论文仿真中 N=20，此处以 N=5 便于可视化。"

### Block 3: 问题设置 🔴 必须修正

**问题 1**：`f = [70 187]` — 频率单位错误，论文使用 rad/s 非 Hz

**问题 2**：`Delta = epsi*s` 的参数名与论文不一致

**修正代码**（替换整个 Block 3）：
```matlab
clear; close all; clc;
%% 参数设置 (2017 JVC Section 6)
fs = 10000;              % 采样频率 [Hz], T_s = 1e-4 s
Ts = 1/fs;
Tsim = 20;               % 仿真总时间 [s]
Nsim = Tsim*fs;
t = (0:Nsim-1)*Ts;

%% 1. 标称模型 G0(s) 和不确定性
s = tf('s');
G0 = 0.5*(s-0.2)/(s^2+s+1.25);           % 标称模型（含 RHP 零点 s=0.2）
Delta_m = -0.001*s;                        % 乘性未建模动态 (论文 μ=0.001)
G = G0*(1+Delta_m);                        % 实际模型

%% 2. 扰动信号 d(t)
% ★ 论文 ω=[70, 187] rad/s (≈ 11.1, 29.8 Hz)，不是 Hz!
omega = [70, 187];        % rad/s
A = [0.6, 0.7];
phi0 = [0, pi/3];
d = zeros(1,Nsim);
for k = 1:length(omega)
    d = d + A(k)*sin(omega(k)*t + phi0(k));
end
rng(42);  d = d + 0.02*randn(1,Nsim);

%% 对象特性
figure; bodemag(G0); grid on; title('G_0(s) 频响');
```

### Block 4: F(s) 滤波器设计 ⚠️ 参数修正

**问题**：`delta0 = 0.` 导致 `designLTIFilter` 对 RHP 零点处理不正确（零点刚好在虚轴上时才有用）。代码应显式调用或直接用解析式。

**修正代码**（替换整个 Block 4）：
```matlab
%% 3. 固定滤波器 F(s) 设计 (Eq. 13)
% F(s) = κ₀ · 2α²(s²+s+1.25) / ((s+α)²(s+0.2))
% 设计原理: G₀ 的 RHP 零点 s=0.2 镜像反射至 s=-0.2
kappa0 = 0.5;    % -6 dB 增益
alpha_val = 500; % 带宽参数
m_val = 2;       % 相对度

F = kappa0 * 2*alpha_val^2*(s^2+s+1.25) / ((s+alpha_val)^2*(s+0.2));

% 也可使用通用设计函数
% F = designLTIFilter(G0, kappa0, 0.01, alpha_val, m_val, 600);

%% 验证 G₀F 频谱平坦化效果
figure;
bodemag(G0, F, F*G0); grid on;
legend('G_0(s)', 'F(s)', 'G_0(s)F(s)');
title('频谱平坦化: |G_0(jω)F(jω)| ≈ 1 在 [0, α] 范围内');
```

### Block 5: 简化自适应控制器 🔴 必须重写

**问题**：这段代码结构完全错误——
1. `phi_reg = -y(k-1:-1:k-Nparam)` — 用延迟的 y 而非 G₀FΛ[z]
2. `eps_k = y(k)` — 用 y 替代归一化预测误差
3. `theta = theta - K_gain*eps_k` — 符号错误且公式不对
4. 缺少 z=G₀[u] 和 u=−F[θᵀΛ[z]] 的结构

**修正代码**（替换整个 Block 5，作为教学简化版）：
```matlab
%% 自适应控制器结构示意 (简化版)
% 完整实现见 Block 6

% K(s,θ) = Σ θ_k · λ^{N-k}/(s+λ)^{N-k}
Nparam = 20;            % 论文 N=20
lambda = 500;           % 基函数极点 (rad/s)

% Λ(s) 基函数: 级联一阶节实现 (避免高阶系数数值崩溃)
G1 = tf(lambda, [1 lambda]);  % 一阶节 λ/(s+λ)

% 离散化
G1d = c2d(G1, Ts, 'tustin');
[numG1, denG1] = tfdata(G1d, 'v');

% 绘制几个 Λ_i(s) 的频响
figure; hold on;
for i = [1, 5, 10, 15, 20]
    order = Nparam - i + 1;
    [mag, ~, w] = bode(tf(lambda^order, poly(-lambda*ones(1,order))), logspace(0,4,500));
    semilogx(w/(2*pi), 20*log10(squeeze(mag)), 'DisplayName', sprintf('Λ_{%d}', i));
end
grid on; legend; xlabel('Hz'); ylabel('dB');
title('Λ_i(s) 基函数频响 (N=20, λ=500)');
```

### Block 6: 完整仿真 🔴 必须重写

**问题**（与原版 `Jafari_adaptive_CC.m` 相同的全部 bug）：
1. 频率单位 (Hz → rad/s)
2. N=20 但基函数实现为 FIR 延迟（应为 Laguerre 级联）
3. 控制律 `u = F[-θᵀφ]` 应为 `u = F[-θᵀΛ[z]]`
4. 回归向量使用手动差分方程（应统一用 `filter()`）
5. 自适应律使用离散 RLS 应改为 Euler 离散化连续 RLS
6. 缺少参数投影

**修正代码**（替换整个 Block 6，与 `Jafari_adaptive_CC_fixed.m` 对齐）：
```matlab
%% ====== 完整仿真 (2017 JVC CC 修正版) ======
% 离散化滤波器
Fd = c2d(F, Ts, 'tustin');
[numF, denF] = tfdata(Fd, 'v');
Gd = c2d(G, Ts, 'tustin');
[numG, denG] = tfdata(Gd, 'v');
G0d = c2d(G0, Ts, 'tustin');
[numG0, denG0] = tfdata(G0d, 'v');

% 自适应律参数
P = eye(Nparam) * 500;   % P(0)
gamma0 = 1.0;
theta_max = 5;           % 投影界
theta = zeros(Nparam, 1);

% 滤波器状态初始化
nF = max(length(numF), length(denF))-1;
nG0 = max(length(numG0), length(denG0))-1;
% Λ 级联状态 (共享链优化: N个一阶节)
zL = zeros(Nparam, 1);
% phi 路径状态
phi_F_zf = zeros(Nparam, nF);
phi_G0_zf = zeros(Nparam, nG0);
% 被控对象状态
G_zf = zeros(nG0, 1);
G0_zf = zeros(nG0, 1);
F_ctrl_zf = zeros(nF, 1);

% 信号存储
y = zeros(1,Nsim);  u = zeros(1,Nsim);
plant_out = zeros(1,Nsim);  z_hist = zeros(1,Nsim);
theta_hist = zeros(Nparam,Nsim);  P_tr = zeros(1,Nsim);

%% 仿真循环
lambda_vals = zeros(Nparam, 1);
phi_reg = zeros(Nparam, 1);

fprintf('===== 2017 JVC CC 仿真 =====\n');
for k = 2:Nsim
    % 1. y(k) = G[u(k-1)] + d(k)
    [yG, G_zf] = filter(numG, denG, u(k-1), G_zf);
    plant_out(k) = yG;  y(k) = yG + d(k);

    % 2. z(k) = y(k) - G₀[u(k-1)]
    [g0u, G0_zf] = filter(numG0, denG0, u(k-1), G0_zf);
    z_k = y(k) - g0u;  z_hist(k) = z_k;

    % 3. Λ[z] + F[Λ[z]] + G₀[F[Λ[z]]] → φ
    % ★ 共享级联链: 一遍扫描得全部 Λ_i[z]
    y_stage = z_k;
    for stage = 1:Nparam
        [y_stage, zL(stage)] = filter(numG1, denG1, y_stage, zL(stage));
        lambda_vals(Nparam - stage + 1) = y_stage;
    end
    for i = 1:Nparam
        [f_out, phi_F_zf(i,:)] = filter(numF, denF, lambda_vals(i), phi_F_zf(i,:)');
        [phi_reg(i), phi_G0_zf(i,:)] = filter(numG0, denG0, f_out, phi_G0_zf(i,:)');
    end

    % 4. 自适应律 (Euler 离散化连续 RLS)
    m2 = 1 + gamma0*(phi_reg'*phi_reg);
    pred_err = z_k - theta'*phi_reg;
    epsilon = pred_err / m2;
    P = P - Ts*(P*(phi_reg*phi_reg')*P)/m2;
    theta_new = theta + Ts*P*epsilon*phi_reg;
    if norm(theta_new) > theta_max
        theta = theta_new*(theta_max/norm(theta_new));
    else
        theta = theta_new;
    end
    theta_hist(:,k) = theta;  P_tr(k) = trace(P);

    % 5. ★ 控制律 u = -F[θᵀΛ[z]] (用 Λ[z], 非 φ!)
    Kz = theta'*lambda_vals;
    [u(k), F_ctrl_zf] = filter(numF, denF, -Kz, F_ctrl_zf);

    if mod(k, round(Nsim/8))==0
        fprintf('  t=%.1fs: ||θ||=%.2f, tr(P)=%.1e\n', t(k), norm(theta), trace(P));
    end
end

%% 性能评估
fi = round(0.5*Nsim):Nsim;
rms_d = sqrt(mean(d(fi).^2));
rms_y = sqrt(mean(y(fi).^2));
fprintf('\n抑制效果: %.1f dB\n', 20*log10(rms_y/rms_d));

%% 绘图
figure('Position', [100 100 1000 600]);
subplot(2,2,1); plot(t, d, 'Color', [.7 .7 .7]); hold on; plot(t, y, 'b-');
title(sprintf('输出 vs 扰动 (%.1f dB)', 20*log10(rms_y/rms_d))); legend('d','y'); grid on;
subplot(2,2,2); plot(t, theta_hist(1:5:end,:)'); title('θ(t) 演化'); grid on;
subplot(2,2,3); yyaxis left; plot(t, vecnorm(theta_hist)); ylabel('||θ||');
yyaxis right; semilogy(t, max(P_tr,1e-10)); ylabel('tr(P)');
title('收敛指标'); grid on; xlabel('t (s)');
subplot(2,2,4);
[Pyy, fp] = pwelch(y(fi).*hann(length(fi))', hann(512), 256, 512, fs);
[Pdd, ~] = pwelch(d(fi).*hann(length(fi))', hann(512), 256, 512, fs);
semilogy(fp(fp<=100), Pdd(fp<=100), 'Color', [.7 .7 .7]); hold on;
semilogy(fp(fp<=100), Pyy(fp<=100), 'b-');
xline(70/(2*pi), 'r--'); xline(187/(2*pi), 'r--');
title('稳态 PSD'); legend('d','y'); grid on; xlim([0 50]);
```

### Block 7: `designLTIFilter` 函数 ✅ 保留但增加注释
此函数是实现 F(s) 通用设计的工具，结构正确。建议在函数注释中增加：
> "对于论文中的 G₀(s)=0.5(s-0.2)/(s²+s+1.25)，等效于调用 `F=designLTIFilter(G0, 0.5, 0.01, 500, 2, 600)`"

---

## 二、JafariJVC_Discrete.mlx（1 个代码块，多个空块待填充）

### Block 1: 系统设置 ⚠️ 需小幅修正

**问题 1**：`omega0 = 0.0521` 是离散域频率 (rad/sample)，正确。
**问题 2**：`F = tf(1,1,Ts)` 为 F(z)=1，但论文的 2015 TAC 用的是 F=100（静态增益），且这里标注的是"Example 6.1"应使用 2015 TAC 的参数。

### 空代码块：需填充内容

| 块位置 | 标题 | 建议填充 |
|--------|------|---------|
| Block 2 | Example 1 固定滤波器 | MOSEK/CVX 凸优化求解最优 θ |
| Block 3 | (1) 阶次影响分析 | N=10,15,20,30,50 的对比 |
| Block 4 | (2) 鲁棒性分析 | 不同 Δₘ 幅度下的抑制效果 |
| Block 5 | Example 2 自适应滤波器 | 与 `Jafari_adaptive_CD_fixed.m` 对齐的完整自适应仿真 |

### Block 5 建议填充代码（自适应 CD 修正版）：

```matlab
%% 自适应 CD 控制器 (2017 JVC 离散实现)
% 与 Continuous 版结构完全一致，但明确标注为离散时间实现
% 参数和算法与 Jafari_adaptive_CD_fixed.m 对齐

% 离散化
G1d = c2d(tf(500, [1 500]), Ts, 'tustin');
[numG1, denG1] = tfdata(G1d, 'v');
Fd = c2d(F, Ts, 'tustin');  % F(s) 来自 Block 4
[numF, denF] = tfdata(Fd, 'v');
G0d = c2d(G0, Ts, 'tustin');
[numG0, denG0] = tfdata(G0d, 'v');

Nparam = 20;  P = eye(Nparam)*500;
theta = zeros(Nparam, 1);
zL = zeros(Nparam, 1);  % Λ 级联状态

% [仿真循环 — 与 Continuous Block 6 相同]
% [此处省略循环代码，详见 demos/Jafari_adaptive_CD_fixed.m]
```

---

## 三、JafariTAC_DiscreteAVC.mlx

此文件对应 2015 TAC 论文的 AVC 悬架模型，与 2017 JVC 是不同的仿真场景。代码块需要与 `Jafari2015_discrete_AVC_fixed.m` 对齐。

**主要差异**：
- 被控对象为离散 AVC 模型 (Ts=1/480)，含 (z-1) 因子
- F(z)=100 为静态增益（非动态滤波器）
- 基函数为纯延迟 α(z)=[z⁻ᴺ,...,z⁻¹]（非 Laguerre）
- 自适应律为标准离散 RLS（非 Euler 离散化连续 RLS）
- 需包含协方差重置机制

---

## 四、应用修正的方法

MATLAB 2025a 中打开 .mlx 文件 → 找到对应代码块 → 复制上方修正代码 → 替换原代码 → 保存。

每个被替换的代码块前已标注 `★` 以示修改项。
