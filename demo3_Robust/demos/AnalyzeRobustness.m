%% 控制能量与鲁棒性分析
% 分析 G₀ 的频率特性、控制信号频谱、以及鲁棒稳定性条件
clear; close all; clc;
addpath('../../functions');

%% ====== 系统定义 ======
s = tf('s');
fs = 10000; Ts = 1/fs;
G0 = 0.5*(s-0.2)/(s^2+s+1.25);
kappa0 = 0.5; alpha_val = 500;
F = kappa0 * 2*alpha_val^2*(s^2+s+1.25)/((s+alpha_val)^2*(s+0.2));
Delta_m = -0.001 * s;  % 论文乘性不确定性

%% ====== 1. G₀ 频率响应分析 ======
omega = logspace(-1, 4, 2000);  % 0.1 ~ 10000 rad/s
[mag_G0, ~] = bode(G0, omega);
mag_G0 = squeeze(mag_G0);
[mag_F, ~] = bode(F, omega);
mag_F = squeeze(mag_F);
[mag_G0F, ~] = bode(G0*F, omega);
mag_G0F = squeeze(mag_G0F);

omega_dist = [70, 187];  % 扰动频率
mag_G0_dist = interp1(omega, mag_G0, omega_dist);
mag_F_dist = interp1(omega, mag_F, omega_dist);

fprintf('===== G₀ 频响分析 =====\n');
for i = 1:length(omega_dist)
    fprintf('  ω=%.0f rad/s (%.1f Hz): |G₀|=%.4f (%.1f dB), |F|=%.1f (%.1f dB), |G₀F|=%.3f (%.1f dB)\n', ...
        omega_dist(i), omega_dist(i)/(2*pi), ...
        mag_G0_dist(i), 20*log10(mag_G0_dist(i)), ...
        mag_F_dist(i), 20*log10(mag_F_dist(i)), ...
        mag_G0F(i)*mag_F_dist(i), 20*log10(mag_G0_dist(i)*mag_F_dist(i)));
end
fprintf('  所需控制增益 (|G₀|⁻¹): %.0f ~ %.0f (%.0f ~ %.0f dB)\n', ...
    1/mag_G0_dist(1), 1/mag_G0_dist(2), ...
    20*log10(1/mag_G0_dist(1)), 20*log10(1/mag_G0_dist(2)));

%% ====== 2. 鲁棒稳定性分析 ======
% 论文定理条件: κ₀ θ̄_M ||F̄ G₀ Δₘ||₁ < 1
% F̄ = F/κ₀ (归一化 F)
F_bar = F / kappa0;
T_rob = F_bar * G0 * Delta_m;

% 计算 L1 范数 (脉冲响应的 L1 范数)
% 对于连续系统，近似为 H∞ 范数或通过频率响应估算
[mag_T, ~] = bode(T_rob, omega);
mag_T = squeeze(mag_T);
Hinf_T = max(mag_T);  % H∞ 上界
L1_est = Hinf_T;       % 对 SISO: ||T||₁ ≥ ||T||∞

fprintf('\n===== 鲁棒稳定性分析 =====\n');
fprintf('  ||F̄ G₀ Δₘ||∞ ≈ %.4f\n', Hinf_T);
fprintf('  ||F̄ G₀ Δₘ||₁ 估计 ≥ %.4f\n', L1_est);

% 论文的鲁棒条件: κ₀ * θ_max * ||F̄ G₀ Δₘ||₁ < 1
% θ_max 来自投影界 (默认 1e3)
theta_bound = 1e3;
robust_LHS = kappa0 * theta_bound * L1_est;
fprintf('  鲁棒条件 LHS = κ₀·θ_max·||F̄G₀Δₘ||₁ ≈ %.2f (需 < 1)\n', robust_LHS);
if robust_LHS < 1
    fprintf('  ✓ 满足鲁棒稳定性条件\n');
else
    fprintf('  ❌ 违反鲁棒条件！需要降低 κ₀ 或收紧 θ 界\n');
end

% 检查不同频率下的环路增益
omega_high = logspace(2, 4, 500);
[mag_Delta, ~] = bode(Delta_m, omega_high);
mag_Delta = squeeze(mag_Delta);
[mag_G0F_hf, ~] = bode(G0*F, omega_high);
mag_G0F_hf = squeeze(mag_G0F_hf);

fprintf('\n===== 高频鲁棒性 (ω > 100 rad/s) =====\n');
loop_hf = mag_G0F_hf .* mag_Delta';
[max_loop, idx] = max(loop_hf);
fprintf('  最大环路增益 |G₀F Δₘ| = %.4f at ω = %.0f rad/s\n', max_loop, omega_high(idx));
if max_loop < 1
    fprintf('  ✓ 高频环路满足小增益条件\n');
else
    fprintf('  ❌ 高频环路违反小增益条件，可能不稳定\n');
end

%% ====== 3. 运行短仿真获取控制频谱 ======
fprintf('\n===== 控制频谱分析 =====\n');

Nparam = 20;
lambda_val = 500;
G1_s = tf(lambda_val, [1 lambda_val]);
G1_d = c2d(G1_s, Ts, 'tustin');
[numG1, denG1] = tfdata(G1_d, 'v');

Lambda_order = zeros(Nparam, 1);
for i = 1:Nparam
    Lambda_order(i) = Nparam - i + 1;
end

Fd = c2d(F, Ts, 'tustin');
[numFd, denFd] = tfdata(Fd, 'v');
Gd = c2d(G0, Ts, 'tustin');
[numGd, denGd] = tfdata(Gd, 'v');
G0d = c2d(G0, Ts, 'tustin');
[numG0d, denG0d] = tfdata(G0d, 'v');

% 扰动
Tsim = 5; Nsim = Tsim * fs; t = (0:Nsim-1) * Ts;
omega_d = [70, 187]; A = [0.6, 0.7]; phi0 = [pi/4, pi/2];
d = zeros(1, Nsim);
for k = 1:length(omega_d)
    d = d + A(k)*sin(omega_d(k)*t + phi0(k));
end
d = d + 0.02*randn(1, Nsim);

% 初始化 (与修复版相同)
hist_len = max([length(denFd), length(numFd), length(denGd), length(numGd)]) + Nparam + 5;
G_u_hist = zeros(1, hist_len); G_y_hist = zeros(1, hist_len);
G0_u_hist = zeros(1, hist_len); G0_y_hist = zeros(1, hist_len);
F_ctrl_u_hist = zeros(1, hist_len); F_ctrl_y_hist = zeros(1, hist_len);
max_order = Nparam;
phi_Lambda_u = zeros(Nparam, max_order, hist_len);
phi_Lambda_y = zeros(Nparam, max_order, hist_len);
phi_F_u = zeros(Nparam, hist_len); phi_F_y = zeros(Nparam, hist_len);
phi_G0_u = zeros(Nparam, hist_len); phi_G0_y = zeros(Nparam, hist_len);

y = zeros(1, Nsim); u = zeros(1, Nsim);
theta = zeros(Nparam, 1); P = eye(Nparam) * 500;
gamma0 = 1.0; theta_max_val = 1e3;

for k = hist_len+1:Nsim
    [y_G_k, G_u_hist, G_y_hist] = IIR_filter(numGd, denGd, u(k-1), G_u_hist, G_y_hist);
    y(k) = y_G_k + d(k);
    [G0_u_k, G0_u_hist, G0_y_hist] = IIR_filter(numG0d, denG0d, u(k-1), G0_u_hist, G0_y_hist);
    z_k = y(k) - G0_u_k;

    phi_reg = zeros(Nparam, 1);
    lambda_vals = zeros(Nparam, 1);
    for i = 1:Nparam
        casc_out = z_k;
        for stage = 1:Lambda_order(i)
            uu = reshape(phi_Lambda_u(i,stage,:), 1, []);
            yy = reshape(phi_Lambda_y(i,stage,:), 1, []);
            [casc_out, uu_n, yy_n] = IIR_filter(numG1, denG1, casc_out, uu, yy);
            phi_Lambda_u(i,stage,:) = uu_n;
            phi_Lambda_y(i,stage,:) = yy_n;
        end
        lambda_vals(i) = casc_out;
        [f_out, phi_F_u(i,:), phi_F_y(i,:)] = ...
            IIR_filter(numFd, denFd, casc_out, phi_F_u(i,:), phi_F_y(i,:));
        [phi_reg(i), phi_G0_u(i,:), phi_G0_y(i,:)] = ...
            IIR_filter(numG0d, denG0d, f_out, phi_G0_u(i,:), phi_G0_y(i,:));
    end

    m_s_sq = 1 + gamma0 * (phi_reg.' * phi_reg);
    pred_err = z_k - theta.' * phi_reg;
    epsilon = pred_err / m_s_sq;
    P = P - Ts * (P * (phi_reg * phi_reg.') * P) / m_s_sq;
    theta_new = theta + Ts * P * epsilon * phi_reg;
    if norm(theta_new) > theta_max_val
        theta = theta_new * (theta_max_val / norm(theta_new));
    else
        theta = theta_new;
    end
    K_z_k = theta.' * lambda_vals;
    [u(k), F_ctrl_u_hist, F_ctrl_y_hist] = ...
        IIR_filter(numFd, denFd, -K_z_k, F_ctrl_u_hist, F_ctrl_y_hist);
end

% 稳态频谱分析
final_idx = round(0.5*Nsim):Nsim;
y_ss = y(final_idx); u_ss = u(final_idx); d_ss = d(final_idx);

[Pyy, f_pyy] = pwelch(y_ss, hann(2048), 1024, 2048, fs);
[Puu, ~] = pwelch(u_ss, hann(2048), 1024, 2048, fs);
[Pdd, ~] = pwelch(d_ss, hann(2048), 1024, 2048, fs);

% 计算各频段能量分布
bands = {[0 5], [5 15], [15 30], [30 50], [50 100], [100 500], [500 fs/2]};
fprintf('\n频段能量分布:\n');
fprintf('  %-20s %12s %12s %12s\n', '频段 (Hz)', '扰动 PSD', '输出 PSD', '抑制 dB');
for b = 1:length(bands)
    idx_b = f_pyy >= bands{b}(1) & f_pyy < bands{b}(2);
    if any(idx_b)
        p_d = mean(Pdd(idx_b));
        p_y = mean(Pyy(idx_b));
        att = 10*log10(max(p_y, 1e-20) / max(p_d, 1e-20));
        fprintf('  [%3.0f, %3.0f] Hz     %12.4f %12.4f %12.1f\n', ...
            bands{b}(1), bands{b}(2), p_d, p_y, att);
    end
end

% 控制 RMS 分解
rms_u_total = sqrt(mean(u_ss.^2));
fprintf('\n控制信号总 RMS: %.1f\n', rms_u_total);

%% ====== 4. 绘制综合分析图 ======
out_dir = 'E:/Code/MATLAB_ANC/ANC-Classic-Control-Simulations/demo3_Robust/demos/';

% 图1: G₀, F, G₀F Bode图 + 扰动频率标记
fig1 = figure('Position', [100 100 1000 700]);

subplot(2,1,1);
semilogx(omega/(2*pi), 20*log10(mag_G0), 'b-', 'LineWidth', 1.5); hold on;
semilogx(omega/(2*pi), 20*log10(mag_F), 'g-', 'LineWidth', 1.5);
semilogx(omega/(2*pi), 20*log10(mag_G0F), 'r-', 'LineWidth', 1.5);
xline(omega_dist(1)/(2*pi), 'k--', 'LineWidth', 1);
xline(omega_dist(2)/(2*pi), 'k--', 'LineWidth', 1);
text(omega_dist(1)/(2*pi)*1.1, -20, sprintf('ω₁=%.0f\\n|G₀|=-%.0fdB', omega_dist(1), -20*log10(mag_G0_dist(1))), 'FontSize', 9);
text(omega_dist(2)/(2*pi)*1.1, -20, sprintf('ω₂=%.0f\\n|G₀|=-%.0fdB', omega_dist(2), -20*log10(mag_G0_dist(2))), 'FontSize', 9);
grid on; legend('G₀(s)', 'F(s)', 'G₀F(s)', 'Location', 'best');
title('传递函数幅频特性'); xlabel('频率 (Hz)'); ylabel('幅值 (dB)');
ylim([-80, 40]);

subplot(2,1,2);
semilogx(omega/(2*pi), 20*log10(mag_G0), 'b-', 'LineWidth', 1.5); hold on;
semilogx(omega/(2*pi), 20*log10(1./mag_G0), 'r--', 'LineWidth', 1.5);
xline(omega_dist(1)/(2*pi), 'k--'); xline(omega_dist(2)/(2*pi), 'k--');
grid on; legend('|G₀|', '|G₀|^{-1} (所需控制增益)', 'Location', 'best');
title('G₀ 与控制需求'); xlabel('频率 (Hz)'); ylabel('幅值 (dB)');
ylim([-80, 80]);
saveas(fig1, [out_dir 'robust_G0_analysis.png']);

% 图2: 鲁棒性分析 + 高频环路增益
fig2 = figure('Position', [100 100 1000 500]);
subplot(2,1,1);
semilogx(omega_high/(2*pi), 20*log10(loop_hf), 'b-', 'LineWidth', 1.5); hold on;
yline(0, 'r--', '0 dB (稳定边界)');
grid on;
title(sprintf('高频环路增益 |G₀F Δₘ| (max=%.4f)', max_loop));
xlabel('频率 (Hz)'); ylabel('幅值 (dB)');

subplot(2,1,2);
semilogx(omega_high/(2*pi), mag_Delta, 'r-', 'LineWidth', 1.5); hold on;
semilogx(omega_high/(2*pi), mag_G0F_hf, 'b-', 'LineWidth', 1.5);
grid on; legend('|Δₘ| = 0.001ω', '|G₀F|', 'Location', 'best');
title('未建模动态 Δₘ 与标称环路 |G₀F|'); xlabel('频率 (Hz)');
saveas(fig2, [out_dir 'robust_loop_gain.png']);

% 图3: 控制信号频谱
fig3 = figure('Position', [100 100 1200 500]);
subplot(1,2,1);
semilogy(f_pyy(f_pyy<=200), Pdd(f_pyy<=200), 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5); hold on;
semilogy(f_pyy(f_pyy<=200), Pyy(f_pyy<=200), 'b-', 'LineWidth', 1.5);
xline(omega_dist(1)/(2*pi), 'r--'); xline(omega_dist(2)/(2*pi), 'r--');
grid on; legend('扰动 PSD', '输出 PSD');
title('输出功率谱密度 (0-200 Hz)'); xlabel('频率 (Hz)');

subplot(1,2,2);
semilogy(f_pyy(f_pyy<=fs/2), Puu(f_pyy<=fs/2), 'g-', 'LineWidth', 1);
hold on;
% 在扰动频率处标记
for i = 1:length(omega_dist)
    f_hz = omega_dist(i)/(2*pi);
    xline(f_hz, 'r--', 'LineWidth', 1);
    text(f_hz*1.1, max(Puu)/10, sprintf('%.1f Hz', f_hz), 'FontSize', 9);
end
grid on; title('控制信号功率谱密度'); xlabel('频率 (Hz)');
saveas(fig3, [out_dir 'robust_spectrum.png']);

% 图4: κ₀ 对鲁棒性的影响
kappa_vec = logspace(-2, 1, 50);
robust_margin = zeros(size(kappa_vec));
for i = 1:length(kappa_vec)
    F_test = kappa_vec(i) * F_bar;
    T_test = F_bar * G0 * Delta_m;  % F_bar 不含 κ₀
    [mag_T, ~] = bode(T_test, omega);
    robust_margin(i) = kappa_vec(i) * theta_bound * max(squeeze(mag_T));
end

fig4 = figure('Position', [100 100 600 400]);
semilogx(kappa_vec, robust_margin, 'b-', 'LineWidth', 2); hold on;
yline(1, 'r--', '稳定边界');
xline(kappa0, 'g--', sprintf('κ₀=%.1f', kappa0));
grid on;
title('κ₀ 对鲁棒裕度的影响'); xlabel('κ₀'); ylabel('κ₀·θ̄ₘ·||F̄G₀Δₘ||₁');

% 找出安全范围
safe_kappa = kappa_vec(robust_margin < 1);
if ~isempty(safe_kappa)
    xline(max(safe_kappa), 'g-', sprintf('κ₀,max=%.2f', max(safe_kappa)));
    fprintf('\n===== κ₀ 鲁棒范围 =====\n');
    fprintf('  当前 κ₀ = %.2f, 鲁棒裕度 = %.4f\n', kappa0, interp1(kappa_vec, robust_margin, kappa0));
    fprintf('  最大安全 κ₀ ≈ %.4f\n', max(safe_kappa));
end
saveas(fig4, [out_dir 'robust_kappa_sweep.png']);

fprintf('\n分析完成。图片保存到 %s\n', out_dir);
