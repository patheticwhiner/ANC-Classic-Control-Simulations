%% Jafari 2017 JVC 离散域 CD 自适应控制器 — 统一对比脚本
% 按顺序执行 2 种基函数解法，最后并列对比。
%
% Section A: α(z) 时移基 —— φ_i = G0·F·z(k-i), i=1..N, 离散 RLS
% Section B: Λ(z) 级联基 —— Λ(z) 一阶节级联, Euler RLS + 投影
% Section C: 并列对比       —— 性能表 + 统一绘图 + K(z) 频响

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
run(fullfile(projectRoot, 'project_init.m'));
addpath(fileparts(scriptDir));

%% ====== 配置区 ======
cfg.Nparam     = 20;
cfg.lambda_val = 500;
cfg.kappa0     = 0.5;
cfg.alpha_val  = 500;
cfg.fs         = 10000;
cfg.Tsim       = 10;
cfg.gamma0     = 1.0;
cfg.P0         = 500;
cfg.thetaMax   = 5;
cfg.rng_seed   = 42;

%% ====== Section 0: 共享题设 ======
modelFile = fullfile(projectRoot, 'dataset', 'syn_JVC2017_3rd.mat');
load(modelFile, 'model');
G0 = model.G0;
fprintf('===== 共享题设: G0=0.5(s-0.2)/(s^2+s+1.25) [离散域 CD] =====\n');
s = tf('s');
Ts = 1/cfg.fs;
F = cfg.kappa0 * 2*cfg.alpha_val^2*(s^2+s+1.25)/((s+cfg.alpha_val)^2*(s+0.2));

omega = [70, 187];   % rad/s
A_m   = [0.6, 0.7];
phi0  = [pi/4, pi/2];

fprintf('扰动: ω=[%.0f, %.0f] rad/s\n', omega(1), omega(2));
fprintf('F(s) = %.1f·2·%d²(s²+s+1.25)/((s+%d)²(s+0.2))\n', ...
    cfg.kappa0, cfg.alpha_val, cfg.alpha_val);

%% ====== 共享离散化 (Tustin) ======
Fd_tustin = c2d(F, Ts, 'tustin');
[numF_t, denF_t] = tfdata(Fd_tustin, 'v');
G0d_tustin = c2d(G0, Ts, 'tustin');
[numG0_t, denG0_t] = tfdata(G0d_tustin, 'v');

% G1(z) — 一阶 Laguerre 节 (Section B 用)
G1_d = c2d(tf(cfg.lambda_val, [1 cfg.lambda_val]), Ts, 'tustin');
[numG1, denG1] = tfdata(G1_d, 'v');

nF_t = max(length(numF_t), length(denF_t)) - 1;
nG0_t = max(length(numG0_t), length(denG0_t)) - 1;

% 生成扰动信号 (共用)
Nsim = cfg.Tsim * cfg.fs;
t = (0:Nsim-1) * Ts;
d = zeros(1, Nsim);
for k = 1:length(omega)
    d = d + A_m(k) * sin(omega(k)*t + phi0(k));
end
rng(cfg.rng_seed);
d = d + 0.02*randn(1, Nsim);

N = cfg.Nparam;


%% ========================================================================
%%  Section A: α(z) 时移基 — φ_i = G0·F·z(k-i)
%% ========================================================================
fprintf('\n========== Section A: α(z) 时移基 ==========\n');
fprintf('N=%d, 基函数: α(z)=[z^{-1}, z^{-2}, ..., z^{-%d}]^T, 离散 RLS\n', N, N);

tic_A = tic;

% ZOH 离散化 (匹配原版 α 版本行为)
Fd_zoh = c2d(F, Ts, 'zoh');
[numFz, denFz] = tfdata(Fd_zoh, 'v');
G0d_zoh = c2d(G0, Ts, 'zoh');
[numG0z, denG0z] = tfdata(G0d_zoh, 'v');
Gd_zoh = c2d(G0, Ts, 'zoh');  % 控制路径
[numGz, denGz] = tfdata(Gd_zoh, 'v');

hist_len = N + 1;

% IIR_filter 缓冲区
F_u_hist = zeros(1, hist_len);  F_y_hist = zeros(1, hist_len);
G_u_hist = zeros(1, hist_len);  G_y_hist = zeros(1, hist_len);
G0_u_hist = zeros(1, hist_len); G0_y_hist = zeros(1, hist_len);
phi_G0_u_hist = zeros(N, hist_len); phi_G0_y_hist = zeros(N, hist_len);
phi_F_u_hist = zeros(N, hist_len);  phi_F_y_hist = zeros(N, hist_len);

% 信号存储
y_A = zeros(1, Nsim);  u_A = zeros(1, Nsim);
plant_out_A = zeros(1, Nsim);  z_hist_A = zeros(1, Nsim);
theta_hist_A = zeros(N, Nsim);  P_tr_hist_A = zeros(1, Nsim);
phi_hist_A = zeros(N, Nsim);
pred_err_hist_A = zeros(1, Nsim);  epsilon_hist_A = zeros(1, Nsim);
m_sq_hist_A = zeros(1, Nsim);  theta_norm_hist_A = zeros(1, Nsim);
Kz_hist_A = zeros(1, Nsim);  phi_norm_hist_A = zeros(1, Nsim);

theta_A = zeros(N, 1);
P_A = eye(N) * cfg.P0;
phi_reg_A = zeros(N, 1);

for k = hist_len+1:Nsim
    % 系统输出
    [y_G_k, G_u_hist, G_y_hist] = IIR_filter(numGz, denGz, u_A(k-1), G_u_hist, G_y_hist);
    plant_out_A(k) = y_G_k;
    y_A(k) = y_G_k + d(k);

    % z(k)
    [G0_u_k, G0_u_hist, G0_y_hist] = IIR_filter(numG0z, denG0z, u_A(k-1), G0_u_hist, G0_y_hist);
    z_k = y_A(k) - G0_u_k;
    z_hist_A(k) = z_k;

    % φ_i = G0·F·z(k-i), i=1..N
    for i = 1:N
        z_delay = 0;
        if k - i >= 1
            z_delay = z_hist_A(k - i);
        end
        [f_out, phi_F_u_hist(i,:), phi_F_y_hist(i,:)] = ...
            IIR_filter(numFz, denFz, z_delay, phi_F_u_hist(i,:), phi_F_y_hist(i,:));
        [phi_reg_A(i), phi_G0_u_hist(i,:), phi_G0_y_hist(i,:)] = ...
            IIR_filter(numG0z, denG0z, f_out, phi_G0_u_hist(i,:), phi_G0_y_hist(i,:));
    end
    phi_hist_A(:, k) = phi_reg_A;

    % 离散 RLS (无 Ts 因子)
    m_s_sq = 1 + cfg.gamma0*(phi_reg_A.'*phi_reg_A);
    pred_err = z_k - theta_A.'*phi_reg_A;
    epsilon = pred_err / m_s_sq;
    denom_rls = m_s_sq + (phi_reg_A.' * P_A * phi_reg_A);
    if denom_rls <= 0, denom_rls = 1e-12; end
    K_gain = (P_A * phi_reg_A) / m_s_sq;
    theta_new = theta_A + K_gain * epsilon;
    P_A = P_A - (P_A * (phi_reg_A * phi_reg_A.') * P_A) / denom_rls;

    % 投影
    theta_norm_new = norm(theta_new);
    if theta_norm_new > cfg.thetaMax
        theta_A = (theta_new / theta_norm_new) * cfg.thetaMax;
    else
        theta_A = theta_new;
    end

    theta_hist_A(:, k) = theta_A;
    P_tr_hist_A(k) = trace(P_A);

    % 控制律 u = -F[θᵀα[z]] = -F[Σ θ_i·z(k-i)]
    Kz_k = theta_A.' * phi_reg_A;
    [u_A(k), F_u_hist, F_y_hist] = IIR_filter(numFz, denFz, -Kz_k, F_u_hist, F_y_hist);

    pred_err_hist_A(k) = pred_err;  epsilon_hist_A(k) = epsilon;
    m_sq_hist_A(k) = m_s_sq;  theta_norm_hist_A(k) = norm(theta_A);
    Kz_hist_A(k) = Kz_k;  phi_norm_hist_A(k) = norm(phi_reg_A);

    if any(abs(theta_A) > 1000) || trace(P_A) > 1e6
        fprintf('  [α(z)] 发散! k=%d\n', k); break;
    end
end
time_A = toc(tic_A);

% --- A 结果封装 ---
final_idx = round(0.5*Nsim):Nsim;
% K(z) FIR 系数
b_K_A = [0; flip(theta_A(:))];
result_A = struct(...
    'name',       '\alpha(z) 时移基', ...
    't',          t, ...
    'y',          y_A, ...
    'u',          u_A, ...
    'd',          d, ...
    'plant_out',  plant_out_A, ...
    'z_hist',     z_hist_A, ...
    'theta_hist', theta_hist_A, ...
    'P_tr_hist',  P_tr_hist_A, ...
    'phi_norm_hist', phi_norm_hist_A, ...
    'pred_err_hist', pred_err_hist_A, ...
    'epsilon_hist',  epsilon_hist_A, ...
    'm_sq_hist',     m_sq_hist_A, ...
    'theta_norm_hist', theta_norm_hist_A, ...
    'Kz_hist',    Kz_hist_A, ...
    'supp_db',    20*log10(rms(y_A(final_idx))/rms(d(final_idx))), ...
    'rms_dist',   rms(d(final_idx)), ...
    'rms_y',      rms(y_A(final_idx)), ...
    'rms_u',      rms(u_A(final_idx)), ...
    'theta_final', theta_A, ...
    'P_final',     P_A, ...
    'sim_time',    time_A, ...
    'final_idx',   final_idx, ...
    'Lambda_tf',   {{}}, ...
    'K_tf',        [], ...
    'b_K',         b_K_A, ...
    'Hk',          [], ...
    'fk',          []);

fprintf('[α(z)] 抑制: %.1f dB, ||θ||=%.3f, 耗时: %.1f s\n', ...
    result_A.supp_db, norm(theta_A), time_A);


%% ========================================================================
%%  Section B: Λ(z) 级联基 — 级联一阶节 + Euler RLS
%% ========================================================================
fprintf('\n========== Section B: Λ(z) 级联基 ==========\n');
fprintf('N=%d, 基函数: Λ_i(z)=[λ/(s+λ)]^{N-i+1} → Tustin 离散, Euler RLS\n', N);

tic_B = tic;

% 滤波器状态
Lambda_zf_B = zeros(N, 1);
phi_F_zf_B  = zeros(N, nF_t);
phi_G0_zf_B = zeros(N, nG0_t);
G_zf_B = zeros(nG0_t, 1);  G0_zf_B = zeros(nG0_t, 1);
F_ctrl_zf_B = zeros(nF_t, 1);

% 信号存储
y_B = zeros(1, Nsim);  u_B = zeros(1, Nsim);
plant_out_B = zeros(1, Nsim);  z_hist_B = zeros(1, Nsim);
theta_hist_B = zeros(N, Nsim);  P_tr_hist_B = zeros(1, Nsim);
phi_hist_B = zeros(N, Nsim);
pred_err_hist_B = zeros(1, Nsim);  epsilon_hist_B = zeros(1, Nsim);
m_sq_hist_B = zeros(1, Nsim);  theta_norm_hist_B = zeros(1, Nsim);
Kz_hist_B = zeros(1, Nsim);  phi_norm_hist_B = zeros(1, Nsim);

theta_B = zeros(N, 1);
P_B = eye(N) * cfg.P0;
lambda_vals_B = zeros(N, 1);
phi_reg_B = zeros(N, 1);

for k = 2:Nsim
    [y_G_k, G_zf_B] = filter(numG0_t, denG0_t, u_B(k-1), G_zf_B);
    plant_out_B(k) = y_G_k;  y_B(k) = y_G_k + d(k);

    [G0_u_k, G0_zf_B] = filter(numG0_t, denG0_t, u_B(k-1), G0_zf_B);
    z_k = y_B(k) - G0_u_k;  z_hist_B(k) = z_k;

    % 共享级联链
    y_stage = z_k;
    for stage = 1:N
        [y_stage, Lambda_zf_B(stage)] = filter(numG1, denG1, y_stage, Lambda_zf_B(stage));
        lambda_vals_B(N - stage + 1) = y_stage;
    end

    for i = 1:N
        [f_out, phi_F_zf_B(i,:)] = filter(numF_t, denF_t, lambda_vals_B(i), phi_F_zf_B(i,:)');
        [phi_reg_B(i), phi_G0_zf_B(i,:)] = filter(numG0_t, denG0_t, f_out, phi_G0_zf_B(i,:)');
    end
    phi_hist_B(:, k) = phi_reg_B;

    % Euler RLS + 投影
    m_s_sq = 1 + cfg.gamma0*(phi_reg_B.'*phi_reg_B);
    pred_err = z_k - theta_B.'*phi_reg_B;
    epsilon = pred_err / m_s_sq;
    P_phi = P_B * phi_reg_B;
    P_B = P_B - Ts*(P_phi*P_phi.')/m_s_sq;
    theta_new = theta_B + Ts*P_B*epsilon*phi_reg_B;
    th_norm = norm(theta_new);
    theta_B = theta_new * min(1, cfg.thetaMax/th_norm);
    theta_hist_B(:, k) = theta_B;  P_tr_hist_B(k) = trace(P_B);

    Kz = theta_B.' * lambda_vals_B;
    [u_B(k), F_ctrl_zf_B] = filter(numF_t, denF_t, -Kz, F_ctrl_zf_B);

    pred_err_hist_B(k) = pred_err;  epsilon_hist_B(k) = epsilon;
    m_sq_hist_B(k) = m_s_sq;  theta_norm_hist_B(k) = th_norm;
    Kz_hist_B(k) = Kz;  phi_norm_hist_B(k) = norm(phi_reg_B);
end
time_B = toc(tic_B);

% K(z) freqz 计算
nfft_bode = 1024;
w_freqz = linspace(0, pi, nfft_bode);
Hk_B = zeros(1, nfft_bode);
for i = 1:N
    order = N - i + 1;
    [H_lam, ~] = freqz(numG1, denG1, w_freqz);
    Hk_B = Hk_B + theta_B(i) * H_lam.^order;
end
fk_B = w_freqz * cfg.fs / (2*pi);

result_B = struct(...
    'name',       '\Lambda(z) 级联基', ...
    't',          t, ...
    'y',          y_B, ...
    'u',          u_B, ...
    'd',          d, ...
    'plant_out',  plant_out_B, ...
    'z_hist',     z_hist_B, ...
    'theta_hist', theta_hist_B, ...
    'P_tr_hist',  P_tr_hist_B, ...
    'phi_norm_hist', phi_norm_hist_B, ...
    'pred_err_hist', pred_err_hist_B, ...
    'epsilon_hist',  epsilon_hist_B, ...
    'm_sq_hist',     m_sq_hist_B, ...
    'theta_norm_hist', theta_norm_hist_B, ...
    'Kz_hist',    Kz_hist_B, ...
    'supp_db',    20*log10(rms(y_B(final_idx))/rms(d(final_idx))), ...
    'rms_dist',   rms(d(final_idx)), ...
    'rms_y',      rms(y_B(final_idx)), ...
    'rms_u',      rms(u_B(final_idx)), ...
    'theta_final', theta_B, ...
    'P_final',     P_B, ...
    'sim_time',    time_B, ...
    'final_idx',   final_idx, ...
    'Lambda_tf',   {{}}, ...
    'K_tf',        [], ...
    'b_K',         [], ...
    'Hk',          Hk_B, ...
    'fk',          fk_B);

fprintf('[Λ(z)] 抑制: %.1f dB, ||θ||=%.3f, 耗时: %.1f s\n', ...
    result_B.supp_db, norm(theta_B), time_B);


%% ========================================================================
%%  Section C: 并列对比
%% ========================================================================
fprintf('\n========== Section C: 并列对比 ==========\n');

results = [result_A, result_B];

% 性能对比表
plotComparisonTable(results, ...
    {{'抑制(dB)', 'supp_db', '%.1f'}, ...
     {'||θ||',   'theta_final', '%.3f'}, ...
     {'输出RMS', 'rms_y', '%.4f'}, ...
     {'控制RMS', 'rms_u', '%.2f'}, ...
     {'耗时(s)', 'sim_time', '%.1f'}}, ...
    'CD 系列 - 两种基函数性能对比');

% 统一绘图
cfg_plot = struct('fs', cfg.fs, 'omega', omega, ...
    'lambda_val', cfg.lambda_val, 'Nparam', N);
plotAdaptiveResults(results, cfg_plot);

% ---- K(z) 频响对比 (CD 专用) ----
figure('Name','K(z) 频率响应对比','Position',[50 50 1000 600]);
colors = lines(2);

% α(z) 基: FIR freqz
subplot(2,1,1);
[Hk_A, fk] = freqz(result_A.b_K, 1, 1024, cfg.fs);
semilogx(fk, 20*log10(abs(Hk_A)), 'Color', colors(1,:), 'LineWidth', 2); hold on;
semilogx(result_B.fk, 20*log10(abs(result_B.Hk)), 'Color', colors(2,:), 'LineWidth', 2);
hold off; grid on;
legend({results(1).name, results(2).name}, 'Location','best');
xlabel('频率 (Hz)'); ylabel('幅值 (dB)');
title('K(z) 幅频响应对比');

subplot(2,1,2);
semilogx(fk, angle(Hk_A)*180/pi, 'Color', colors(1,:), 'LineWidth', 2); hold on;
semilogx(result_B.fk, angle(result_B.Hk)*180/pi, 'Color', colors(2,:), 'LineWidth', 2);
hold off; grid on;
legend({results(1).name, results(2).name}, 'Location','best');
xlabel('频率 (Hz)'); ylabel('相位 (deg)');
title('K(z) 相频响应对比');

fprintf('\n===== CD 系列统一对比完成 =====\n');
fprintf('总结: α(z) 用延迟线参数化, Λ(z) 用 Laguerre 级联基, 数学不等价\n');
