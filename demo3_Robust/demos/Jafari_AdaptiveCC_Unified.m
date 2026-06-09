%% Jafari 2017 JVC 连续时间 CC 自适应控制器 — 统一对比脚本
% 按顺序执行 3 种解法，最后并列对比。
%
% Section A: 原版 (Basic)    — N=6, T=2s, Hz频率, 单高阶Λ, IIR_filter
% Section B: 修正版 (Fixed)  — N=20, T=10s, rad/s, 独立级联Λ, filter()
% Section C: 极速版 (VFast)  — N=20, T=10s, rad/s, 共享级联链Λ, filter()
% Section D: 并列对比         — 性能表 + 统一绘图

clear; close all; clc;
addpath('..');

%% ====== 配置区 ======
% 共享参数（Section B/C 使用；Section A 有独立覆盖）
cfg.Nparam     = 20;       % Λ 基函数个数
cfg.lambda_val = 500;      % Laguerre 基极点 (rad/s)
cfg.kappa0     = 0.5;      % F(s) 增益
cfg.alpha_val  = 500;      % F(s) 截止频率 (rad/s)
cfg.fs         = 10000;    % 采样频率 (Hz)
cfg.Tsim       = 10;       % 仿真时长 (s)
cfg.gamma0     = 1.0;      % RLS 归一化增益
cfg.P0         = 500;      % 协方差矩阵初值缩放
cfg.thetaMax   = 5;        % 参数投影界
cfg.rng_seed   = 42;       % 随机数种子

%% ====== Section 0: 共享题设 ======
fprintf('===== 共享题设: G0=0.5(s-0.2)/(s²+s+1.25) =====\n');
s = tf('s');
G0 = 0.5*(s-0.2)/(s^2+s+1.25);

% F(s) 设计: RHP 零点镜像反射 + 低通
F = cfg.kappa0 * 2*cfg.alpha_val^2*(s^2+s+1.25)/((s+cfg.alpha_val)^2*(s+0.2));

% 扰动定义（统一使用 rad/s）
omega = [70, 187];   % rad/s
A_m   = [0.6, 0.7];
phi0  = [pi/4, pi/2];

fprintf('扰动: ω=[%.0f, %.0f] rad/s (%.1f, %.1f Hz)\n', ...
    omega(1), omega(2), omega(1)/(2*pi), omega(2)/(2*pi));
fprintf('F(s) = %.1f·2·%d²(s²+s+1.25)/((s+%d)²(s+0.2))\n', ...
    cfg.kappa0, cfg.alpha_val, cfg.alpha_val);

%% ====== 预计算离散化滤波器系数 (Section B/C 共用) ======
Ts = 1/cfg.fs;

% G1(z) — 一阶 Laguerre 节
G1_d = c2d(tf(cfg.lambda_val, [1 cfg.lambda_val]), Ts, 'tustin');
[numG1, denG1] = tfdata(G1_d, 'v');

% F_d, G0_d
Fd = c2d(F, Ts, 'tustin');
[numF, denF] = tfdata(Fd, 'v');
G0d = c2d(G0, Ts, 'tustin');
[numG0, denG0] = tfdata(G0d, 'v');

nF = max(length(numF), length(denF)) - 1;
nG0 = max(length(numG0), length(denG0)) - 1;


%% ========================================================================
%%  Section A: 原版 (Basic) — 含已知缺陷, 作为"修复前"基线
%% ========================================================================
fprintf('\n========== Section A: 原版 (Basic) ==========\n');
fprintf('N=6, T=2s, 频率作为 Hz, 单高阶Λ, IIR_filter\n');

% --- A.1 覆盖参数 ---
N_A = 6;
Tsim_A = 2;
Nsim_A = Tsim_A * cfg.fs;
t_A = (0:Nsim_A-1) * Ts;

% 扰动: Hz 频率 (已知缺陷 #1)
f_Hz = [70, 187];
d_A = zeros(1, Nsim_A);
for k = 1:length(A_m)
    d_A = d_A + A_m(k) * sin(2*pi*f_Hz(k)*t_A + phi0(k));
end
rng(cfg.rng_seed);
d_A = d_A + 0.02*randn(1, Nsim_A);

% Lambda: 单高阶 TF (已知缺陷 #4)
Lambda_num_coeffs_A = cell(N_A, 1);
Lambda_den_coeffs_A = cell(N_A, 1);
for i = 1:N_A
    order = i;
    den_poly = 1;
    for j = 1:order
        den_poly = conv(den_poly, [1 cfg.lambda_val]);
    end
    Lambda_s = tf(cfg.lambda_val^order, den_poly);
    Lambda_d = c2d(Lambda_s, Ts, 'tustin');
    [n_tmp, d_tmp] = tfdata(Lambda_d, 'v');
    Lambda_num_coeffs_A{i} = n_tmp;
    Lambda_den_coeffs_A{i} = d_tmp;
end

% 离散化 (IIR_filter 用)
Fd_A = c2d(F, Ts, 'tustin');
[numFd_A, denFd_A] = tfdata(Fd_A, 'v');
Gd_A = c2d(G0, Ts, 'tustin');
[numGd_A, denGd_A] = tfdata(Gd_A, 'v');
[numG0d_A, denG0d_A] = tfdata(G0d, 'v');

% --- A.2 仿真 ---
tic_A = tic;

P_A = eye(N_A) * cfg.P0;
theta_A = zeros(N_A, 1);

% 信号存储
y_A = zeros(1, Nsim_A);  u_A = zeros(1, Nsim_A);
plant_out_A = zeros(1, Nsim_A);  z_hist_A = zeros(1, Nsim_A);
theta_hist_A = zeros(N_A, Nsim_A);  P_tr_hist_A = zeros(1, Nsim_A);
phi_hist_A = zeros(N_A, Nsim_A);
pred_err_hist_A = zeros(1, Nsim_A);  epsilon_hist_A = zeros(1, Nsim_A);
m_sq_hist_A = zeros(1, Nsim_A);  theta_norm_hist_A = zeros(1, Nsim_A);
Kz_hist_A = zeros(1, Nsim_A);  phi_norm_hist_A = zeros(1, Nsim_A);

% IIR_filter 历史缓冲区
hist_len_A = N_A + 1;
F_u_hist_A = zeros(1, hist_len_A);  F_y_hist_A = zeros(1, hist_len_A);
G_u_hist_A = zeros(1, hist_len_A);  G_y_hist_A = zeros(1, hist_len_A);
G0_u_hist_A = zeros(1, hist_len_A); G0_y_hist_A = zeros(1, hist_len_A);
phi_G0_u_hist_A = zeros(1, hist_len_A); phi_G0_y_hist_A = zeros(1, hist_len_A);
phi_F_u_hist_A = zeros(1, hist_len_A); phi_F_y_hist_A = zeros(1, hist_len_A);
Lambda_u_hist_A = zeros(N_A, hist_len_A);
Lambda_y_hist_A = zeros(N_A, hist_len_A);

phi_reg_A = zeros(N_A, 1);

for k = hist_len_A+1:Nsim_A
    % 系统输出
    [y_G_k, G_u_hist_A, G_y_hist_A] = IIR_filter(numGd_A, denGd_A, u_A(k-1), G_u_hist_A, G_y_hist_A);
    plant_out_A(k) = y_G_k;
    y_A(k) = y_G_k + d_A(k);

    % z(k)
    [G0_u_k, G0_u_hist_A, G0_y_hist_A] = IIR_filter(numG0d_A, denG0d_A, u_A(k-1), G0_u_hist_A, G0_y_hist_A);
    z_k = y_A(k) - G0_u_k;
    z_hist_A(k) = z_k;

    % phi = G0·F·Λ[z] (控制律混淆 — 缺陷 #3 保留)
    for i = 1:N_A
        [lam_out, Lambda_u_hist_A(i,:), Lambda_y_hist_A(i,:)] = ...
            IIR_filter(Lambda_num_coeffs_A{i}, Lambda_den_coeffs_A{i}, z_k, ...
                      Lambda_u_hist_A(i,:), Lambda_y_hist_A(i,:));
        [f_out, phi_F_u_hist_A, phi_F_y_hist_A] = IIR_filter(numFd_A, denFd_A, lam_out, phi_F_u_hist_A, phi_F_y_hist_A);
        [phi_reg_A(i), phi_G0_u_hist_A, phi_G0_y_hist_A] = IIR_filter(numG0d_A, denG0d_A, f_out, phi_G0_u_hist_A, phi_G0_y_hist_A);
    end
    phi_hist_A(:, k) = phi_reg_A;

    % RLS (Euler 离散)
    m_s_sq = 1 + cfg.gamma0*(phi_reg_A.'*phi_reg_A);
    pred_err = z_k - theta_A.'*phi_reg_A;
    epsilon = pred_err / m_s_sq;
    theta_A = theta_A + Ts * P_A * epsilon * phi_reg_A;
    P_A = P_A - Ts * P_A * (phi_reg_A * phi_reg_A.') * P_A / m_s_sq;

    % 控制律 u = -F[θᵀΛ[z]]
    Kz_k = theta_A.' * phi_reg_A;
    [u_A(k), F_u_hist_A, F_y_hist_A] = IIR_filter(numFd_A, denFd_A, -Kz_k, F_u_hist_A, F_y_hist_A);

    theta_hist_A(:, k) = theta_A;  P_tr_hist_A(k) = trace(P_A);
    pred_err_hist_A(k) = pred_err;  epsilon_hist_A(k) = epsilon;
    m_sq_hist_A(k) = m_s_sq;  theta_norm_hist_A(k) = norm(theta_A);
    Kz_hist_A(k) = Kz_k;  phi_norm_hist_A(k) = norm(phi_reg_A);

    if any(abs(theta_A) > 100) || trace(P_A) > 1e6
        fprintf('  [Basic] 发散! k=%d\n', k); break;
    end
end
time_A = toc(tic_A);

% --- A.3 结果封装 ---
final_idx_A = round(0.5*Nsim_A):Nsim_A;
result_A = struct(...
    'name',       'Basic', ...
    't',          t_A, ...
    'y',          y_A, ...
    'u',          u_A, ...
    'd',          d_A, ...
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
    'supp_db',    20*log10(rms(y_A(final_idx_A))/rms(d_A(final_idx_A))), ...
    'rms_dist',   rms(d_A(final_idx_A)), ...
    'rms_y',      rms(y_A(final_idx_A)), ...
    'rms_u',      rms(u_A(final_idx_A)), ...
    'theta_final', theta_A, ...
    'P_final',     P_A, ...
    'sim_time',    time_A, ...
    'final_idx',   final_idx_A, ...
    'Lambda_tf',   {{}}, ...    % Basic 不参与 Bode 对比
    'K_tf',        []);

fprintf('[Basic] 抑制: %.1f dB, ||θ||=%.3f, 耗时: %.1f s\n', ...
    result_A.supp_db, norm(theta_A), time_A);


%% ========================================================================
%%  Section B: 修正版 (Fixed) — N=20, rad/s, 独立级联Λ
%% ========================================================================
fprintf('\n========== Section B: 修正版 (Fixed) ==========\n');
fprintf('N=%d, T=%ds, rad/s, 独立级联Λ, filter()\n', cfg.Nparam, cfg.Tsim);

% --- B.1 参数 ---
N_B = cfg.Nparam;
Tsim_B = cfg.Tsim;
Nsim_B = Tsim_B * cfg.fs;
t_B = (0:Nsim_B-1) * Ts;

% 扰动 (rad/s)
d_B = zeros(1, Nsim_B);
for k = 1:length(omega)
    d_B = d_B + A_m(k) * sin(omega(k)*t_B + phi0(k));
end
rng(cfg.rng_seed);
d_B = d_B + 0.02*randn(1, Nsim_B);

% --- B.2 仿真 ---
tic_B = tic;

P_B = eye(N_B) * cfg.P0;
theta_B = zeros(N_B, 1);

% 滤波器状态
Lambda_zf_B = zeros(N_B, 1);
phi_F_zf_B  = zeros(N_B, nF);
phi_G0_zf_B = zeros(N_B, nG0);
G_zf_B = zeros(nG0, 1);  G0_zf_B = zeros(nG0, 1);
F_ctrl_zf_B = zeros(nF, 1);

% 信号存储
y_B = zeros(1, Nsim_B);  u_B = zeros(1, Nsim_B);
plant_out_B = zeros(1, Nsim_B);  z_hist_B = zeros(1, Nsim_B);
theta_hist_B = zeros(N_B, Nsim_B);  P_tr_hist_B = zeros(1, Nsim_B);
phi_hist_B = zeros(N_B, Nsim_B);
pred_err_hist_B = zeros(1, Nsim_B);  epsilon_hist_B = zeros(1, Nsim_B);
m_sq_hist_B = zeros(1, Nsim_B);  theta_norm_hist_B = zeros(1, Nsim_B);
Kz_hist_B = zeros(1, Nsim_B);  phi_norm_hist_B = zeros(1, Nsim_B);

lambda_vals_B = zeros(N_B, 1);
phi_reg_B = zeros(N_B, 1);

for k = 2:Nsim_B
    % 系统输出
    [y_G_k, G_zf_B] = filter(numG0, denG0, u_B(k-1), G_zf_B);
    plant_out_B(k) = y_G_k;  y_B(k) = y_G_k + d_B(k);

    % z(k)
    [G0_u_k, G0_zf_B] = filter(numG0, denG0, u_B(k-1), G0_zf_B);
    z_k = y_B(k) - G0_u_k;  z_hist_B(k) = z_k;

    % 独立级联: 每个 Λ_i 独立从 z 开始滤波
    for i = 1:N_B
        order = N_B - i + 1;
        y_tmp = z_k;
        for s = 1:order
            [y_tmp, Lambda_zf_B(i)] = filter(numG1, denG1, y_tmp, Lambda_zf_B(i));
        end
        lambda_vals_B(i) = y_tmp;
    end

    % φ = G0·F·Λ[z]
    for i = 1:N_B
        [f_out, phi_F_zf_B(i,:)] = filter(numF, denF, lambda_vals_B(i), phi_F_zf_B(i,:)');
        [phi_reg_B(i), phi_G0_zf_B(i,:)] = filter(numG0, denG0, f_out, phi_G0_zf_B(i,:)');
    end
    phi_hist_B(:, k) = phi_reg_B;

    % RLS (Euler + 投影)
    m_s_sq = 1 + cfg.gamma0*(phi_reg_B.'*phi_reg_B);
    pred_err = z_k - theta_B.'*phi_reg_B;
    epsilon = pred_err / m_s_sq;
    P_phi = P_B * phi_reg_B;
    P_B = P_B - Ts*(P_phi*P_phi.')/m_s_sq;
    theta_new = theta_B + Ts*P_B*epsilon*phi_reg_B;
    th_norm = norm(theta_new);
    theta_B = theta_new * min(1, cfg.thetaMax/th_norm);
    theta_hist_B(:, k) = theta_B;  P_tr_hist_B(k) = trace(P_B);

    % 控制律 u = -F[θᵀΛ[z]]
    Kz = theta_B.' * lambda_vals_B;
    [u_B(k), F_ctrl_zf_B] = filter(numF, denF, -Kz, F_ctrl_zf_B);

    pred_err_hist_B(k) = pred_err;  epsilon_hist_B(k) = epsilon;
    m_sq_hist_B(k) = m_s_sq;  theta_norm_hist_B(k) = th_norm;
    Kz_hist_B(k) = Kz;  phi_norm_hist_B(k) = norm(phi_reg_B);
end
time_B = toc(tic_B);

% --- B.3 构建 K_tf 用于 Bode 图 ---
Lambda_tf_B = cell(N_B, 1);
for i = 1:N_B
    order = N_B - i + 1;
    den_poly = 1;
    for j = 1:order, den_poly = conv(den_poly, [1 cfg.lambda_val]); end
    Lambda_tf_B{i} = tf(cfg.lambda_val^order, den_poly);
end
K_tf_B = 0;
for i = 1:N_B, K_tf_B = K_tf_B + theta_B(i)*Lambda_tf_B{i}; end

% --- B.4 结果封装 ---
final_idx_B = round(0.5*Nsim_B):Nsim_B;
result_B = struct(...
    'name',       'Fixed', ...
    't',          t_B, ...
    'y',          y_B, ...
    'u',          u_B, ...
    'd',          d_B, ...
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
    'supp_db',    20*log10(rms(y_B(final_idx_B))/rms(d_B(final_idx_B))), ...
    'rms_dist',   rms(d_B(final_idx_B)), ...
    'rms_y',      rms(y_B(final_idx_B)), ...
    'rms_u',      rms(u_B(final_idx_B)), ...
    'theta_final', theta_B, ...
    'P_final',     P_B, ...
    'sim_time',    time_B, ...
    'final_idx',   final_idx_B, ...
    'Lambda_tf',   {Lambda_tf_B}, ...
    'K_tf',        K_tf_B);

fprintf('[Fixed] 抑制: %.1f dB, ||θ||=%.3f, 耗时: %.1f s\n', ...
    result_B.supp_db, norm(theta_B), time_B);


%% ========================================================================
%%  Section C: 极速版 (VFast) — N=20, rad/s, 共享级联链Λ
%% ========================================================================
fprintf('\n========== Section C: 极速版 (VFast) ==========\n');
fprintf('N=%d, T=%ds, rad/s, 共享级联链Λ (%.0f× 理论加速)\n', ...
    cfg.Nparam, cfg.Tsim, 250/60);

% --- C.1 参数 (与 Fixed 完全相同) ---
N_C = cfg.Nparam;
Tsim_C = cfg.Tsim;
Nsim_C = Nsim_B;
t_C = t_B;
d_C = d_B;

% --- C.2 仿真 ---
tic_C = tic;

P_C = eye(N_C) * cfg.P0;
theta_C = zeros(N_C, 1);

Lambda_zf_C = zeros(N_C, 1);
phi_F_zf_C  = zeros(N_C, nF);
phi_G0_zf_C = zeros(N_C, nG0);
G_zf_C = zeros(nG0, 1);  G0_zf_C = zeros(nG0, 1);
F_ctrl_zf_C = zeros(nF, 1);

y_C = zeros(1, Nsim_C);  u_C = zeros(1, Nsim_C);
plant_out_C = zeros(1, Nsim_C);  z_hist_C = zeros(1, Nsim_C);
theta_hist_C = zeros(N_C, Nsim_C);  P_tr_hist_C = zeros(1, Nsim_C);
phi_hist_C = zeros(N_C, Nsim_C);
pred_err_hist_C = zeros(1, Nsim_C);  epsilon_hist_C = zeros(1, Nsim_C);
m_sq_hist_C = zeros(1, Nsim_C);  theta_norm_hist_C = zeros(1, Nsim_C);
Kz_hist_C = zeros(1, Nsim_C);  phi_norm_hist_C = zeros(1, Nsim_C);

lambda_vals_C = zeros(N_C, 1);
phi_reg_C = zeros(N_C, 1);

for k = 2:Nsim_C
    [y_G_k, G_zf_C] = filter(numG0, denG0, u_C(k-1), G_zf_C);
    plant_out_C(k) = y_G_k;  y_C(k) = y_G_k + d_C(k);

    [G0_u_k, G0_zf_C] = filter(numG0, denG0, u_C(k-1), G0_zf_C);
    z_k = y_C(k) - G0_u_k;  z_hist_C(k) = z_k;

    % ★ 共享级联链: 一遍扫描算出全部 Λ_i[z]
    y_stage = z_k;
    for stage = 1:N_C
        [y_stage, Lambda_zf_C(stage)] = filter(numG1, denG1, y_stage, Lambda_zf_C(stage));
        lambda_vals_C(N_C - stage + 1) = y_stage;
    end

    % φ = G0·F·Λ[z]
    for i = 1:N_C
        [f_out, phi_F_zf_C(i,:)] = filter(numF, denF, lambda_vals_C(i), phi_F_zf_C(i,:)');
        [phi_reg_C(i), phi_G0_zf_C(i,:)] = filter(numG0, denG0, f_out, phi_G0_zf_C(i,:)');
    end
    phi_hist_C(:, k) = phi_reg_C;

    % RLS (与 Fixed 相同)
    m_s_sq = 1 + cfg.gamma0*(phi_reg_C.'*phi_reg_C);
    pred_err = z_k - theta_C.'*phi_reg_C;
    epsilon = pred_err / m_s_sq;
    P_phi = P_C * phi_reg_C;
    P_C = P_C - Ts*(P_phi*P_phi.')/m_s_sq;
    theta_new = theta_C + Ts*P_C*epsilon*phi_reg_C;
    th_norm = norm(theta_new);
    theta_C = theta_new * min(1, cfg.thetaMax/th_norm);
    theta_hist_C(:, k) = theta_C;  P_tr_hist_C(k) = trace(P_C);

    Kz = theta_C.' * lambda_vals_C;
    [u_C(k), F_ctrl_zf_C] = filter(numF, denF, -Kz, F_ctrl_zf_C);

    pred_err_hist_C(k) = pred_err;  epsilon_hist_C(k) = epsilon;
    m_sq_hist_C(k) = m_s_sq;  theta_norm_hist_C(k) = th_norm;
    Kz_hist_C(k) = Kz;  phi_norm_hist_C(k) = norm(phi_reg_C);
end
time_C = toc(tic_C);

% --- C.3 K_tf ---
K_tf_C = 0;
for i = 1:N_C, K_tf_C = K_tf_C + theta_C(i)*Lambda_tf_B{i}; end

% --- C.4 结果封装 ---
final_idx_C = round(0.5*Nsim_C):Nsim_C;
result_C = struct(...
    'name',       'VFast', ...
    't',          t_C, ...
    'y',          y_C, ...
    'u',          u_C, ...
    'd',          d_C, ...
    'plant_out',  plant_out_C, ...
    'z_hist',     z_hist_C, ...
    'theta_hist', theta_hist_C, ...
    'P_tr_hist',  P_tr_hist_C, ...
    'phi_norm_hist', phi_norm_hist_C, ...
    'pred_err_hist', pred_err_hist_C, ...
    'epsilon_hist',  epsilon_hist_C, ...
    'm_sq_hist',     m_sq_hist_C, ...
    'theta_norm_hist', theta_norm_hist_C, ...
    'Kz_hist',    Kz_hist_C, ...
    'supp_db',    20*log10(rms(y_C(final_idx_C))/rms(d_C(final_idx_C))), ...
    'rms_dist',   rms(d_C(final_idx_C)), ...
    'rms_y',      rms(y_C(final_idx_C)), ...
    'rms_u',      rms(u_C(final_idx_C)), ...
    'theta_final', theta_C, ...
    'P_final',     P_C, ...
    'sim_time',    time_C, ...
    'final_idx',   final_idx_C, ...
    'Lambda_tf',   {Lambda_tf_B}, ...
    'K_tf',        K_tf_C);

fprintf('[VFast] 抑制: %.1f dB, ||θ||=%.3f, 耗时: %.1f s (vs Fixed %.1f×)\n', ...
    result_C.supp_db, norm(theta_C), time_C, time_B/time_C);


%% ========================================================================
%%  Section D: 并列对比
%% ========================================================================
fprintf('\n========== Section D: 并列对比 ==========\n');

results = [result_A, result_B, result_C];

% 性能对比表
plotComparisonTable(results, ...
    {{'抑制(dB)', 'supp_db', '%.1f'}, ...
     {'||θ||',   'theta_final', '%.3f'}, ...
     {'输出RMS', 'rms_y', '%.4f'}, ...
     {'控制RMS', 'rms_u', '%.2f'}, ...
     {'耗时(s)', 'sim_time', '%.1f'}}, ...
    'CC 系列 - 三种解法性能对比');

% 统一绘图 (时域 + 误差 + 性能 + Bode)
cfg_plot = struct('fs', cfg.fs, 'omega', omega, ...
    'lambda_val', cfg.lambda_val, 'Nparam', cfg.Nparam);
plotAdaptiveResults(results, cfg_plot);

fprintf('\n===== CC 系列统一对比完成 =====\n');
fprintf('总结: Basic 含已知缺陷(Hz/N=6/控制律混淆), Fixed 修正, VFast %s\n', ...
    '优化级联链');
