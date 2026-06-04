%% Jafari 2017 JVC 自适应 CC 控制器 — 修正版
% 修复: (1)频率rad/s (2)N=20 (3)控制律u=-F[θᵀΛ[z]] (4)级联Λ
% 保留原版全部分析工具 + 新增PSD/G₀F分析
clear; close all; clc;
addpath('../../functions');

%% ====== 论文参数 (2017 JVC Section 6) ======
fs = 10000;  Ts = 1/fs;  Tsim = 10;  Nsim = Tsim*fs;  t = (0:Nsim-1)*Ts;
s = tf('s');
G0 = 0.5*(s-0.2)/(s^2+s+1.25);
G = G0;

% 扰动: 论文 ω=[70,187] rad/s
omega = [70, 187];  A = [0.6 0.7];  phi0 = [pi/4 pi/2];
d = zeros(1,Nsim);
for k=1:length(omega), d = d + A(k)*sin(omega(k)*t + phi0(k)); end
rng(42); d = d + 0.02*randn(1,Nsim);

% F(s) 滤波器
kappa0 = 0.5;  alpha_val = 500;
F = kappa0 * 2*alpha_val^2*(s^2+s+1.25)/((s+alpha_val)^2*(s+0.2));

%% ====== Λ(s) 基函数 — 级联一阶节 ======
Nparam = 20;  lambda_val = 500;
G1_d = c2d(tf(lambda_val,[1 lambda_val]), Ts, 'tustin');
[numG1, denG1] = tfdata(G1_d, 'v');

% 离散化 F, G
Fd = c2d(F, Ts, 'tustin');
[numF, denF] = tfdata(Fd, 'v');
Gd = c2d(G, Ts, 'tustin');
[numG, denG] = tfdata(Gd, 'v');
G0d = c2d(G0, Ts, 'tustin');
[numG0, denG0] = tfdata(G0d, 'v');

%% ====== 自适应参数 ======
P = eye(Nparam) * 500;
gamma0 = 1.0;  theta_max_val = 5;
theta = zeros(Nparam, 1);

%% ====== 滤波器状态 ======
nF = max(length(numF), length(denF))-1;
nG0 = max(length(numG0), length(denG0))-1;
Lambda_zf = zeros(Nparam, 1);
phi_F_zf = zeros(Nparam, nF);
phi_G0_zf = zeros(Nparam, nG0);
G_zf = zeros(nG0, 1);  G0_zf = zeros(nG0, 1);
F_ctrl_zf = zeros(nF, 1);

%% ====== 信号存储 ======
y = zeros(1,Nsim);  u = zeros(1,Nsim);
plant_out = zeros(1,Nsim);  z_hist = zeros(1,Nsim);
theta_hist = zeros(Nparam,Nsim);  P_tr_hist = zeros(1,Nsim);
phi_hist = zeros(Nparam,Nsim);  % 原版: 回归向量历史
pred_err_hist = zeros(1,Nsim);  epsilon_hist = zeros(1,Nsim);
m_sq_hist = zeros(1,Nsim);  theta_norm_hist = zeros(1,Nsim);
Kz_hist = zeros(1,Nsim);  phi_norm_hist = zeros(1,Nsim);

fprintf('===== 2017 JVC CC 修正版 =====\n');
fprintf('N=%d, λ=%.0f, κ₀=%.1f, α=%.0f, θ_max=%.0f\n', Nparam, lambda_val, kappa0, alpha_val, theta_max_val);
fprintf('扰动: ω=[%.0f, %.0f] rad/s (%.1f, %.1f Hz)\n', omega(1), omega(2), omega(1)/(2*pi), omega(2)/(2*pi));

%% ====== 仿真循环 ======
lambda_vals = zeros(Nparam, 1);
phi_reg = zeros(Nparam, 1);

for k = 2:Nsim
    % 1. y(k)=G[u(k-1)]+d(k)
    [y_G_k, G_zf] = filter(numG, denG, u(k-1), G_zf);
    plant_out(k) = y_G_k;  y(k) = y_G_k + d(k);

    % 2. z(k)=y(k)-G₀[u(k-1)]
    [G0_u_k, G0_zf] = filter(numG0, denG0, u(k-1), G0_zf);
    z_k = y(k) - G0_u_k;  z_hist(k) = z_k;

    % 3. φ=G₀FΛ[z], λ=Λ[z] (共享级联链)
    y_stage = z_k;
    for stage = 1:Nparam
        [y_stage, Lambda_zf(stage)] = filter(numG1, denG1, y_stage, Lambda_zf(stage));
        lambda_vals(Nparam - stage + 1) = y_stage;
    end
    for i = 1:Nparam
        [f_out, phi_F_zf(i,:)] = filter(numF, denF, lambda_vals(i), phi_F_zf(i,:)');
        [phi_reg(i), phi_G0_zf(i,:)] = filter(numG0, denG0, f_out, phi_G0_zf(i,:)');
    end
    phi_hist(:,k) = phi_reg;

    % 4. 自适应律 (Euler离散)
    m_s_sq = 1 + gamma0*(phi_reg.'*phi_reg);
    pred_err = z_k - theta.'*phi_reg;
    epsilon = pred_err / m_s_sq;
    P_phi = P * phi_reg;
    P = P - Ts*(P_phi*P_phi.')/m_s_sq;
    theta_new = theta + Ts*P*epsilon*phi_reg;
    th_norm = norm(theta_new);
    theta = theta_new * min(1, theta_max_val/th_norm);
    theta_hist(:,k) = theta;  P_tr_hist(k) = trace(P);

    % 5. u=-F[θᵀΛ[z]]
    Kz = theta.'*lambda_vals;
    [u(k), F_ctrl_zf] = filter(numF, denF, -Kz, F_ctrl_zf);

    % 监控
    pred_err_hist(k) = pred_err;  epsilon_hist(k) = epsilon;
    m_sq_hist(k) = m_s_sq;  theta_norm_hist(k) = th_norm;
    Kz_hist(k) = Kz;  phi_norm_hist(k) = norm(phi_reg);

    if mod(k, round(Nsim/8))==0
        fprintf('  t=%.1fs: ||θ||=%.3f, tr(P)=%.2e, z=%.4f\n', t(k), th_norm, trace(P), z_k);
    end
end

%% ========================================================================
%                         原版分析工具 (保留)
%% ========================================================================

% ---- 图1: Λ_i(s) 与 K(s) Bode对比 (原版) ----
theta_final = theta;
f_plot = logspace(0, log10(fs/2), 800);
w_plot = 2*pi*f_plot;

% 构造 Λ_i(s) 连续传递函数用于Bode图
Lambda_tf = cell(Nparam,1);
for i = 1:Nparam
    order = Nparam - i + 1;  den_poly = 1;
    for j = 1:order, den_poly = conv(den_poly, [1 lambda_val]); end
    Lambda_tf{i} = tf(lambda_val^order, den_poly);
end
K_tf = 0;
for i = 1:Nparam, K_tf = K_tf + theta_final(i)*Lambda_tf{i}; end

magL = zeros(Nparam, length(w_plot));  phL = zeros(Nparam, length(w_plot));
for i = 1:Nparam
    [m, p] = bode(Lambda_tf{i}, w_plot);
    magL(i,:) = squeeze(m);  phL(i,:) = squeeze(p);
end
[magK, phaseK] = bode(K_tf, w_plot);
magK = squeeze(magK);  phaseK = squeeze(phaseK);

figure('Name','Λ Bode分析','Position',[50 50 900 600]);
subplot(2,1,1);
semilogx(f_plot, 20*log10(magK), 'k-', 'LineWidth',2.2); hold on;
for i = 1:5:Nparam
    semilogx(f_plot, 20*log10(magL(i,:)), '--', 'LineWidth',1);
end
hold off; grid on;
legend(['K(s)'], 'Location','best');
xlabel('频率 (Hz)'); ylabel('幅值 (dB)');
title('Λ_i(s) 与 K(s) 幅频响应对比');

subplot(2,1,2);
semilogx(f_plot, phaseK, 'k-', 'LineWidth',2.2); hold on;
for i = 1:5:Nparam
    semilogx(f_plot, phL(i,:), '--', 'LineWidth',1);
end
hold off; grid on; xlabel('频率 (Hz)'); ylabel('相位 (deg)');
title('Λ_i(s) 与 K(s) 相频响应对比');

% ---- 图2: 关键时域对比 (原版3×2) ----
figure('Name','时域对比','Position',[100 100 1000 650]);
subplot(3,2,1);
plot(t, d, 'Color',[0.7 0.7 0.7]); hold on;
plot(t, plant_out, 'r--', 'LineWidth',1.2);
plot(t, y, 'b-', 'LineWidth',1);
legend('扰动 d(t)', 'G输出 y_G(t)', '总输出 y(t)');  title('输出信号对比'); grid on;

subplot(3,2,2);
plot(t, u, 'g-', 'LineWidth',1); title('控制信号 u(t)'); grid on;

subplot(3,2,3);
plot(t, z_hist, 'k-', 'LineWidth',1); title('辅佐信号 z(t)'); grid on;

subplot(3,2,4);
plot(t, theta_hist(1:5:end,:)'); title('参数演化 θ(t) (每5个)'); grid on;

subplot(3,2,5);
semilogy(t, max(P_tr_hist,1e-10), 'LineWidth',1); title('P矩阵迹 trace(P)'); grid on;

subplot(3,2,6);
plot(t, phi_norm_hist, 'LineWidth',1); title('回归向量范数 ||φ(t)||'); grid on;

% ---- 图3: 误差与收敛 (原版2×2) ----
figure('Name','误差与收敛','Position',[100 100 900 600]);
subplot(2,2,1);
plot(t, abs(pred_err_hist), 'r-', 'LineWidth',1); hold on;
plot(t, abs(epsilon_hist), 'b-', 'LineWidth',1);
legend('|预测误差|', '|归一化ε|'); title('误差对比'); grid on;

subplot(2,2,2);
plot(t, m_sq_hist, 'LineWidth',1); title('归一化因子 m_s^2'); grid on;

subplot(2,2,3);
plot(t, theta_norm_hist, 'LineWidth',1.2); title('参数范数 ||θ(t)||'); grid on;

subplot(2,2,4);
yyaxis left; plot(t, Kz_hist, 'b-', 'LineWidth',1); ylabel('K(s,θ)z');
yyaxis right; plot(t, phi_norm_hist, 'r-', 'LineWidth',1); ylabel('||φ||');
title('K(s,θ)z 与 φ范数'); grid on;

%% ========================================================================
%                         新增分析工具
%% ========================================================================

final_idx = round(0.5*Nsim):Nsim;
rms_output = sqrt(mean(y(final_idx).^2));
rms_dist = sqrt(mean(d(final_idx).^2));
rms_yG = sqrt(mean(plant_out(final_idx).^2));
supp_db = 20*log10(rms_output/rms_dist);
rms_u_val = sqrt(mean(u(final_idx).^2));

% ---- 图4: 性能指标与频谱 (原版2×2 + 新增PSD) ----
figure('Name','性能与频谱','Position',[100 100 1100 700]);
subplot(2,3,1);
bar([rms_dist, rms_yG, rms_output]);
set(gca, 'XTickLabel', {'扰动RMS','G输出RMS','总输出RMS'});
title(sprintf('抑制: %.1f dB', supp_db)); grid on;

subplot(2,3,2);
vars = [mean(abs(epsilon_hist(final_idx))), mean(abs(pred_err_hist(final_idx))), ...
        mean(m_sq_hist(final_idx)), theta_norm_hist(end)];
bar(vars);
set(gca, 'XTickLabel', {'|ε|均值','|预误|均值','m_s²均值','||θ||终值'});
title('自适应关键指标'); grid on;

subplot(2,3,3);
cw = round(0.3*Nsim):Nsim;
plot(t(cw), theta_norm_hist(cw), 'b-', 'LineWidth',2); hold on;
plot(t(cw), P_tr_hist(cw)/1000, 'r-', 'LineWidth',2);
title('收敛分析（后半程）'); legend('||θ||','tr(P)/1000'); grid on;

subplot(2,3,4);
Nfft = 2048;  fs_idx = final_idx(end-Nfft+1:end);
[Pyy, f_p] = pwelch(y(fs_idx).*hann(Nfft)', hann(512), 256, 512, fs);
[Pdd, ~]  = pwelch(d(fs_idx).*hann(Nfft)', hann(512), 256, 512, fs);
semilogy(f_p(f_p<=100), Pdd(f_p<=100), 'Color',[0.7 0.7 0.7]); hold on;
semilogy(f_p(f_p<=100), Pyy(f_p<=100), 'b-', 'LineWidth',1.2);
xline(omega(1)/(2*pi),'r--'); xline(omega(2)/(2*pi),'r--');
title('稳态PSD频谱'); xlabel('Hz'); grid on; xlim([0 50]);
legend('扰动','输出');

subplot(2,3,5);
[Puu, ~] = pwelch(u(fs_idx).*hann(Nfft)', hann(512), 256, 512, fs);
semilogy(f_p(f_p<=500), Puu(f_p<=500), 'g-', 'LineWidth',1);
xline(omega(1)/(2*pi),'r--'); xline(omega(2)/(2*pi),'r--');
title('控制信号PSD'); xlabel('Hz'); grid on;

subplot(2,3,6);
w_plot2 = logspace(0, log10(fs/2*pi), 500);
[m_g0f, ~] = bode(G0*F, w_plot2); m_g0f = squeeze(m_g0f);
semilogx(w_plot2/(2*pi), 20*log10(m_g0f), 'b-', 'LineWidth',1.5); hold on;
xline(omega(1)/(2*pi),'r--'); xline(omega(2)/(2*pi),'r--');
title('G₀(s)F(s) 频响'); xlabel('Hz'); ylabel('dB'); grid on;

%% ====== 详细性能输出 (原版) ======
fprintf('\n=== 自适应控制性能评估 ===\n');
fprintf('扰动RMS: %.4f\n', rms_dist);
fprintf('G输出RMS: %.4f\n', rms_yG);
fprintf('总输出RMS: %.4f\n', rms_output);
fprintf('抑制效果: %.1f dB\n', supp_db);
fprintf('控制RMS: %.1f\n', rms_u_val);
fprintf('最终参数范数: %.4f\n', norm(theta));
fprintf('最终P矩阵迹: %.4e\n', trace(P));
fprintf('稳态预测误差RMS: %.6f\n', sqrt(mean(pred_err_hist(final_idx).^2)));
fprintf('稳态ε RMS: %.6f\n', sqrt(mean(epsilon_hist(final_idx).^2)));
fprintf('稳态φ范数均值: %.4f\n', mean(phi_norm_hist(final_idx)));
fprintf('调试完成\n');
