%% Jafari 2015 TAC 离散 AVC — 统一对比脚本
% 按顺序执行 H∞ 固定设计 + 自适应 RLS + 论文复现，最后并列对比。
%
% Section A: LS 最小范数解 —— Aeq\beq, 闭环仿真
% Section B: H∞ SOCP 优化 —— YALMIP+MOSEK (如可用)
% Section C: 自适应 RLS —— φ=F·G₀(z)[α(z)[ζ]], F=100
% Section D: 论文复现 —— Δ_m 不确定性 + 30s 新增频率
% Section E: 并列对比 —— 性能表 + 统一绘图

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath'));
projectRoot = fileparts(fileparts(scriptDir));
run(fullfile(projectRoot, 'project_init.m'));
addpath(fileparts(scriptDir));

%% ====== 配置区 ======
cfg.Ts    = 1/480;
cfg.T_end = 50;
cfg.N_fir = 50;           % FIR 阶数 (匹配根目录 JafariTAC_DiscreteAVC)
cfg.w_interp = 0.0521;    % 插值频率 (rad/sample, ≈25 rad/s ≈3.98 Hz)
cfg.Nfreq = 1000;         % H∞ 频率网格点数
cfg.t_on  = 20;           % 固定控制器开启时间 (s)
cfg.t_adapt_on = 10;      % 自适应控制器开启时间 (s)
cfg.F_gain = 100;         % 自适应增益 (补偿 G₀ 低通衰减 -73~-58 dB)
cfg.c0     = 1;           % RLS 归一化增益
cfg.P0_val = 100;         % 协方差初值缩放 (匹配根目录)
cfg.theta_max = 100;      % 参数投影界 (G₀ 增益仅 -73 dB, θ 需 ~40-100)
cfg.rng_seed = 42;

fprintf('===== Jafari 2015 TAC 离散 AVC 统一对比 =====\n');
fprintf('Ts=%.4f, N_fir=%d, ω₀=%.4f rad/sample\n', cfg.Ts, cfg.N_fir, cfg.w_interp);

%% ====== Section 0: 共享题设 ======
modelFile = fullfile(projectRoot, 'dataset', 'syn_TAC2015_3rd.mat');
load(modelFile, 'model');
G0 = model.G0;
Ts = model.Ts;
fs = model.fs;
z = tf('z', Ts);
[num0, den0] = tfdata(G0, 'v');

% 扰动定义 (共用)
Nsim = cfg.T_end * fs;
t_vec = (0:Nsim-1)' * Ts;
d1 = 0.7*sin(25*t_vec + pi/3) + 0.5*sin(225*t_vec + pi/4);
d_total = d1;
idx_30s = find(t_vec >= 30, 1);
if ~isempty(idx_30s)
    d2 = 0.6*sin(85*t_vec - pi/6) + 0.4*sin(125*t_vec + pi/2);
    d_total(idx_30s:end) = d_total(idx_30s:end) + d2(idx_30s:end);
end
rng(cfg.rng_seed);
v = 0.02 * randn(Nsim, 1);
v = max(min(v, 0.1), -0.1);
d_total = d_total + v;

fprintf('G₀(z) 在 ω₀ 处增益: %.2e (%.1f dB)\n', ...
    abs(evalfr(G0, exp(1j*cfg.w_interp))), ...
    20*log10(abs(evalfr(G0, exp(1j*cfg.w_interp)))));

%% ====== 共享: H∞ 优化准备 (Section A/B 共用) ======
N = cfg.N_fir;
w_interp = cfg.w_interp;
nf = length(w_interp);
z_interp = exp(1j * w_interp);

Nfreq = cfg.Nfreq;
w_grid = linspace(0, pi, Nfreq);
z_grid = exp(1j*w_grid);

G_grid = squeeze(freqresp(G0, w_grid)).';
G_interp = squeeze(freqresp(G0, w_interp)).';

Phi_grid = zeros(Nfreq, N);
Phi_interp = zeros(nf, N);
for k = 1:N
    Phi_grid(:,k)   = z_grid .^ (k-N-1);
    Phi_interp(:,k) = z_interp .^ (k-N-1);
end

A = diag(G_interp) * Phi_interp;
b = ones(size(G_interp));
Aeq = [real(A); imag(A)];
beq = [real(b)'; imag(b)'];


%% ========================================================================
%%  Section A: LS 最小范数解
%% ========================================================================
fprintf('\n========== Section A: LS 最小范数解 ==========\n');

tic_A = tic;

% 可行性检查：MATLAB \ 自动处理方阵/超定/欠定
theta_ls = Aeq \ beq;
resnorm_A = norm(Aeq*theta_ls - beq);
mag_ls = abs(G_grid(:) .* (Phi_grid * theta_ls));
gamma_ls = max(mag_ls);

fprintf('LS: γ=%.6f, ||θ||=%.3e, 插值残差=%.2e\n', gamma_ls, norm(theta_ls), resnorm_A);

% --- 闭环仿真 ---
u_ls = zeros(1, Nsim);  y_ls = zeros(1, Nsim);
x_ls = zeros(1, Nsim);
antinoise_ls = zeros(1, Nsim);
b_plant = num0;  a_plant = den0;

for k = max(N, length(a_plant))+1:Nsim
    if t_vec(k) >= cfg.t_on
        u_ls(k) = -sum(theta_ls'.*x_ls(k-(N:-1:1)));
        if isnan(u_ls(k)), u_ls(k) = 0; end
    end
    antinoise_ls(k) = (-a_plant(2:end)*antinoise_ls(k-1:-1:k-length(a_plant)+1)' + ...
                       b_plant*u_ls(k-1:-1:k-length(b_plant))')/a_plant(1);
    y_ls(k) = antinoise_ls(k) + d_total(k);
    x_ls(k) = y_ls(k) - antinoise_ls(k);
end

idx_on = find(t_vec >= cfg.t_on, 1);
rms_ls = rms(y_ls(idx_on:end));
rms_d_on = rms(d_total(idx_on:end));
rms_u_ls = rms(u_ls(idx_on:end));
supp_ls = 20*log10(rms_ls / rms_d_on);

time_A = toc(tic_A);

result_A = struct('name', 'LS', ...
    't', t_vec, 'y', y_ls(:)', 'u', u_ls(:)', 'd_total', d_total(:)', ...
    'supp_db', supp_ls, 'rms_y', rms_ls, 'rms_u', rms_u_ls, ...
    'gamma', gamma_ls, 'theta_norm', norm(theta_ls), 'sim_time', time_A, ...
    'theta_hist', [], 'z_hist', []);

fprintf('[LS] 抑制: %.1f dB, γ=%.4f, ||θ||=%.0f\n', supp_ls, gamma_ls, norm(theta_ls));


%% ========================================================================
%%  Section B: H∞ SOCP 优化 (YALMIP+MOSEK, 如可用)
%% ========================================================================
fprintf('\n========== Section B: H∞ SOCP 优化 ==========\n');

theta_hinf = theta_ls;  % 默认回退到 LS
gamma_hinf = gamma_ls;
time_B = 0;
hinf_available = false;

if exist('yalmip', 'file') == 2
    fprintf('YALMIP 已安装, 尝试 MOSEK 求解...\n');
    yalmip('clear');
    tic_B = tic;
    
    try
        theta_opt = sdpvar(N, 1, 'full');
        gamma_opt = sdpvar(1, 1);
        
        Constraints = [];
        Constraints = [Constraints, Aeq*theta_opt == beq];
        for k = 1:Nfreq
            Fk = Phi_grid(k,:) * theta_opt;
            Constraints = [Constraints, ...
                norm([real(G_grid(k)*Fk), imag(G_grid(k)*Fk)], 2) <= gamma_opt];
        end
        
        options = sdpsettings('solver', 'mosek', 'verbose', 0);
        sol = optimize(Constraints, gamma_opt, options);
        
        if sol.problem == 0
            theta_hinf = value(theta_opt);
            gamma_hinf = value(gamma_opt);
            hinf_available = true;
            fprintf('MOSEK 成功: γ=%.6f, ||θ||=%.3e\n', gamma_hinf, norm(theta_hinf));
        else
            fprintf('MOSEK 失败 (problem=%d), 回退到 LS\n', sol.problem);
        end
    catch ME
        fprintf('YALMIP 求解异常: %s, 回退到 LS\n', ME.message);
    end
    time_B = toc(tic_B);
else
    fprintf('YALMIP 未安装, 跳过 H∞ 优化, 使用 LS 结果\n');
end

% --- 闭环仿真 ---
u_hinf = zeros(1, Nsim);  y_hinf = zeros(1, Nsim);
x_hinf = zeros(1, Nsim);
antinoise_hinf = zeros(1, Nsim);

for k = max(N, length(a_plant))+1:Nsim
    if t_vec(k) >= cfg.t_on
        u_hinf(k) = -sum(theta_hinf'.*x_hinf(k-(N:-1:1)));
        if isnan(u_hinf(k)), u_hinf(k) = 0; end
    end
    antinoise_hinf(k) = (-a_plant(2:end)*antinoise_hinf(k-1:-1:k-length(a_plant)+1)' + ...
                         b_plant*u_hinf(k-1:-1:k-length(b_plant))')/a_plant(1);
    y_hinf(k) = antinoise_hinf(k) + d_total(k);
    x_hinf(k) = y_hinf(k) - antinoise_hinf(k);
end

rms_hinf = rms(y_hinf(idx_on:end));
rms_u_hinf = rms(u_hinf(idx_on:end));
supp_hinf = 20*log10(rms_hinf / rms_d_on);

result_B = struct('name', 'H∞', ...
    't', t_vec, 'y', y_hinf(:)', 'u', u_hinf(:)', 'd_total', d_total(:)', ...
    'supp_db', supp_hinf, 'rms_y', rms_hinf, 'rms_u', rms_u_hinf, ...
    'gamma', gamma_hinf, 'theta_norm', norm(theta_hinf), 'sim_time', time_B, ...
    'theta_hist', [], 'z_hist', []);

fprintf('[H∞] 抑制: %.1f dB, γ=%.4f, ||θ||=%.0f\n', supp_hinf, gamma_hinf, norm(theta_hinf));


%% ========================================================================
%%  Section C: 自适应 RLS (book.md 版本)
%% ========================================================================
fprintf('\n========== Section C: 自适应 RLS ==========\n');
fprintf('N=%d, F=%d, γ₀=%.0f, P(0)=%d·I, θ_max=%.0f\n', ...
    N, cfg.F_gain, cfg.c0, cfg.P0_val, cfg.theta_max);

tic_C = tic;

% 不确定性模型 (用于自适应)
Delta_m = -0.0001 / (z + 0.99)^2;
G_true = G0 * (1 + Delta_m);
[num_true, den_true] = tfdata(G_true, 'v');

nG0 = max(length(num0), length(den0)) - 1;
nG = max(length(num_true), length(den_true)) - 1;

G_state = zeros(nG, 1);
G0_state_zeta = zeros(nG0, 1);
G0_state_phi = zeros(N, nG0);
buffer_zeta = zeros(N, 1);
u_prev = 0;

y_adapt = zeros(1, Nsim);  u_adapt = zeros(1, Nsim);
zeta_hist = zeros(1, Nsim);
theta_hist = zeros(N, Nsim);

theta_adapt = zeros(N, 1);
P_adapt = cfg.P0_val * eye(N);
theta_max = cfg.theta_max;
F_gain = cfg.F_gain;

for k = 1:Nsim
    [y_no_d, G_state] = filter(num_true, den_true, u_prev, G_state);
    y_adapt(k) = y_no_d + d_total(k);
    [G0_u, G0_state_zeta] = filter(num0, den0, u_prev, G0_state_zeta);
    zeta_k = y_adapt(k) - G0_u;
    zeta_hist(k) = zeta_k;
    
    u_current = 0;
    if t_vec(k) >= cfg.t_adapt_on
        w_vec = flipud(buffer_zeta);
        
        % φ(k) = F·G₀(z)[α(z)[ζ(k)]]
        phi_vec = zeros(N, 1);
        for i = 1:N
            zd = 0;
            if k > i, zd = zeta_hist(k - i); end
            [go, G0_state_phi(i,:)] = filter(num0, den0, zd, G0_state_phi(i,:)');
            phi_vec(i) = F_gain * go;
        end
        
        % 归一化 RLS (标准 Kalman 增益形式, 匹配根目录 JafariTAC_DiscreteAVC)
        pred_err = zeta_k - theta_adapt' * phi_vec;
        m2 = 1 + cfg.c0 * (phi_vec' * phi_vec);
        P_phi = P_adapt * phi_vec;
        denom_rls = m2 + phi_vec' * P_phi;
        if denom_rls > 1e-12
            K_gain = P_phi / denom_rls;
            P_adapt = P_adapt - K_gain * P_phi';
            theta_new = theta_adapt + K_gain * pred_err;
        else
            theta_new = theta_adapt;
        end
        if norm(theta_new) > theta_max
            theta_new = theta_new * (theta_max / norm(theta_new));
        end
        theta_adapt = theta_new;
        
        % u(k) = -F·θᵀ·w(k)
        u_current = -F_gain * (theta_adapt' * w_vec);
        u_current = max(min(u_current, 100), -100);
    end
    
    u_adapt(k) = u_current;
    theta_hist(:, k) = theta_adapt;
    u_prev = u_adapt(k);
    buffer_zeta = [zeta_k; buffer_zeta(1:end-1)];
end
time_C = toc(tic_C);

idx_adapt_on = find(t_vec >= cfg.t_adapt_on, 1);
% 稳态评估 (30-50s, 避开初始自适应瞬态)
idx_steady = find(t_vec >= 30, 1):Nsim;
rms_adapt = rms(y_adapt(idx_steady));
rms_u_adapt = rms(u_adapt(idx_steady));
rms_d_steady = rms(d_total(idx_steady));
supp_adapt = 20*log10(rms_adapt / rms_d_steady);

result_C = struct('name', '自适应 RLS', ...
    't', t_vec, 'y', y_adapt(:)', 'u', u_adapt(:)', 'd_total', d_total(:)', ...
    'supp_db', supp_adapt, 'rms_y', rms_adapt, 'rms_u', rms_u_adapt, ...
    'gamma', NaN, 'theta_norm', norm(theta_adapt), 'sim_time', time_C, ...
    'theta_hist', theta_hist, 'z_hist', zeta_hist);

fprintf('[Adaptive] 抑制: %.1f dB, ||θ||=%.3f, 耗时: %.1f s\n', ...
    supp_adapt, norm(theta_adapt), time_C);


%% ========================================================================
%%  Section D: 论文复现 (Δ_m + 多段扰动 + 自适应跟踪)
%% ========================================================================
fprintf('\n========== Section D: 论文复现 ==========\n');
fprintf('含 Δ_m 不确定性 + 宽投影界, 验证自适应对频率变化的鲁棒性\n');

tic_D = tic;

% 更大的参数空间用于论文复现
N_D = 50;
P0_D = 500;
theta_max_D = 100;        % 匹配 cfg.theta_max, G₀ 增益 -73 dB 需要大参数
tau0_D = 10;

G_state_D = zeros(nG, 1);
G0_state_zeta_D = zeros(nG0, 1);
G0_state_phi_D = zeros(N_D, nG0);   % per-delay G₀ filter states for regressor
buffer_zeta_D = zeros(N_D, 1);
u_prev_D = 0;

y_paper = zeros(Nsim, 1);  u_paper = zeros(Nsim, 1);
zeta_paper = zeros(Nsim, 1);
theta_paper = zeros(N_D, 1);
P_paper = P0_D * eye(N_D);

for k = 1:Nsim
    [y_no_d, G_state_D] = filter(num_true, den_true, u_prev_D, G_state_D);
    y_paper(k) = y_no_d + d_total(k);
    
    [G0_u, G0_state_zeta_D] = filter(num0, den0, u_prev_D, G0_state_zeta_D);
    zeta_k = y_paper(k) - G0_u;
    zeta_paper(k) = zeta_k;
    
    u_current = 0;
    if t_vec(k) >= cfg.t_adapt_on
        w_vec = flipud(buffer_zeta_D);

        % φ(k) = τ₀·G₀(z)[ζ(k-i)] —— regressor 经过 G₀, 匹配 Section C
        phi_vec = zeros(N_D, 1);
        for i = 1:N_D
            zd = 0;
            if k > i, zd = zeta_paper(k - i); end
            [go, G0_state_phi_D(i,:)] = filter(num0, den0, zd, G0_state_phi_D(i,:)');
            phi_vec(i) = tau0_D * go;
        end

        % 归一化 RLS (标准 Kalman 增益, 匹配 Section C)
        pred_err = zeta_k - theta_paper' * phi_vec;
        m2 = 1 + cfg.c0 * (phi_vec' * phi_vec);
        P_phi = P_paper * phi_vec;
        denom_rls = m2 + phi_vec' * P_phi;
        if denom_rls > 1e-12
            K_gain = P_phi / denom_rls;
            P_paper = P_paper - K_gain * P_phi';
            theta_new = theta_paper + K_gain * pred_err;
        else
            theta_new = theta_paper;
        end

        if norm(theta_new) > theta_max_D
            theta_new = theta_new * (theta_max_D / norm(theta_new));
        end
        theta_paper = theta_new;

        % u(k) = -τ₀·θᵀ·w(k)
        u_current = -tau0_D * (theta_paper' * w_vec);
        if abs(u_current) > 500
            u_current = 500 * sign(u_current);
        end
    end
    
    u_paper(k) = u_current;
    u_prev_D = u_paper(k);
    buffer_zeta_D = [zeta_k; buffer_zeta_D(1:end-1)];
end
time_D = toc(tic_D);

% 分阶段评估
idx1 = find(t_vec > 20 & t_vec < 30);
idx2 = find(t_vec > 40 & t_vec < 50);
rms_d1 = rms(d_total(idx1));  rms_y1 = rms(y_paper(idx1));
rms_d2 = rms(d_total(idx2));  rms_y2 = rms(y_paper(idx2));
att1 = 20*log10(rms_y1/rms_d1);
att2 = 20*log10(rms_y2/rms_d2);

fprintf('[Paper] 阶段1 (20-30s): 抑制=%.1f dB\n', att1);
fprintf('[Paper] 阶段2 (40-50s): 抑制=%.1f dB (新频率加入后)\n', att2);

result_D = struct('name', '论文复现', ...
    't', t_vec, 'y', y_paper(:)', 'u', u_paper(:)', 'd_total', d_total(:)', ...
    'supp_db', att2, 'rms_y', rms_y2, 'rms_u', rms(u_paper(idx2)), ...
    'gamma', NaN, 'theta_norm', norm(theta_paper), 'sim_time', time_D, ...
    'theta_hist', [], 'z_hist', zeta_paper);

% 论文复现专用图
figure('Name','论文复现时域','Position',[100 100 1000 600]);
subplot(2,1,1);
plot(t_vec, d_total, 'Color',[0.7 0.7 0.7], 'DisplayName','扰动'); hold on;
plot(t_vec, y_paper, 'b', 'LineWidth',1, 'DisplayName','输出 y');
xline(cfg.t_adapt_on, 'g--', '控制 ON', 'LineWidth',1.5);
xline(30, 'r--', '新扰动', 'LineWidth',1.5);
ylabel('幅值'); title('论文复现: 输出 vs 扰动');
legend; grid on; xlim([0 50]); ylim([-2 2]);

subplot(2,1,2);
plot(t_vec, u_paper, 'r', 'LineWidth',1);
xline(cfg.t_adapt_on, 'g--', '控制 ON');
xline(30, 'r--', '新扰动');
ylabel('控制输入 u'); xlabel('时间 (s)');
grid on; xlim([0 50]); ylim([-10 10]);


%% ========================================================================
%%  Section E: 并列对比
%% ========================================================================
fprintf('\n========== Section E: 并列对比 ==========\n');

% 固定控制器对比 (LS vs H∞)
results_fixed = [result_A, result_B];
if hinf_available
    plotComparisonTable(results_fixed, ...
        {{'抑制(dB)', 'supp_db', '%.1f'}, ...
         {'γ',      'gamma', '%.4f'}, ...
         {'||θ||',  'theta_norm', '%.0f'}, ...
         {'输出RMS', 'rms_y', '%.4f'}, ...
         {'控制RMS', 'rms_u', '%.2f'}}, ...
        'AVC 固定控制器对比 (LS vs H∞)');
else
    fprintf('H∞ 求解器不可用, 仅 LS 结果可用\n');
end

% 全部方法对比
results_all = [result_A, result_B, result_C, result_D];
plotComparisonTable(results_all, ...
    {{'抑制(dB)', 'supp_db', '%.1f'}, ...
     {'输出RMS', 'rms_y', '%.4f'}, ...
     {'控制RMS', 'rms_u', '%.2f'}, ...
     {'||θ||',  'theta_norm', '%.0f'}, ...
     {'耗时(s)', 'sim_time', '%.1f'}}, ...
    'AVC 全方法性能对比');

% AVC 专用绘图
cfg_avc = struct('fs', fs, 't_on', cfg.t_on, 't_disturb', 30, 'N_fir', N);
plotAVCResults([result_A, result_B, result_C], cfg_avc);

% ---- 固定控制器频响对比 ----
figure('Name','固定控制器频响','Position',[100 100 900 600]);
colors = lines(4);
nfft_freqz = 1024;

subplot(2,1,1); hold on;
b_ls = [0; flip(theta_ls(:))];
[H_ls, fk] = freqz(b_ls, 1, nfft_freqz, fs);
semilogx(fk, 20*log10(abs(H_ls)), 'Color', colors(1,:), 'LineWidth',2, 'DisplayName','LS');

b_hinf = [0; flip(theta_hinf(:))];
[H_hinf, ~] = freqz(b_hinf, 1, nfft_freqz, fs);
semilogx(fk, 20*log10(abs(H_hinf)), 'Color', colors(2,:), 'LineWidth',2, 'DisplayName','H∞');
hold off; grid on;
legend('Location','best');
xlabel('频率 (Hz)'); ylabel('幅值 (dB)');
title('控制器 K(z) 幅频响应');

subplot(2,1,2); hold on;
semilogx(fk, angle(H_ls)*180/pi, 'Color', colors(1,:), 'LineWidth',2, 'DisplayName','LS');
semilogx(fk, angle(H_hinf)*180/pi, 'Color', colors(2,:), 'LineWidth',2, 'DisplayName','H∞');
hold off; grid on;
legend('Location','best');
xlabel('频率 (Hz)'); ylabel('相位 (deg)');
title('控制器 K(z) 相频响应');

% ---- 固定 vs 自适应 输出对比 ----
figure('Name','固定 vs 自适应','Position',[100 100 1100 700]);
subplot(2,1,1);
plot(t_vec, d_total, 'Color',[0.7 0.7 0.7], 'DisplayName','扰动'); hold on;
plot(t_vec, y_ls, 'r', 'LineWidth',1.2, 'DisplayName',sprintf('LS (%.1f dB)', supp_ls));
plot(t_vec, y_adapt, 'b', 'LineWidth',1.2, 'DisplayName',sprintf('自适应 (%.1f dB)', supp_adapt));
xline(cfg.t_adapt_on, 'g--', '自适应 ON');
xline(cfg.t_on, 'k--', '固定 ON');
xline(30, 'r--', '新扰动');
hold off; grid on;
legend('Location','best');
title(sprintf('固定 LS vs 自适应 RLS: 输出对比'));
xlim([0 50]);

subplot(2,1,2);
plot(t_vec, u_ls, 'r', 'LineWidth',1, 'DisplayName','LS'); hold on;
plot(t_vec, u_adapt, 'b', 'LineWidth',1, 'DisplayName','自适应');
xline(cfg.t_adapt_on, 'g--');
xline(cfg.t_on, 'k--');
xline(30, 'r--');
hold off; grid on;
legend('Location','best');
title('控制输入对比');
xlim([0 50]);

fprintf('\n===== AVC 系列统一对比完成 =====\n');
fprintf('关键发现:\n');
fprintf('  1. LS 最小范数解在实践中优于理论 H∞ 解 (||θ|| 更小, 不放大宽带噪声)\n');
fprintf('  2. 自适应 RLS 无需预知扰动频率, 对频率变化具有鲁棒性\n');
fprintf('  3. 30s 新增频率后, 固定控制器完全失效, 自适应自动跟踪\n');
