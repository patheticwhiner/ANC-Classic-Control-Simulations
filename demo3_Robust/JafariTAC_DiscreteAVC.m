%[text] # 未知周期干扰的主动悬架鲁棒性与性能分析（Jafari, Ioannou,et al. 2015）
%%
%[text] ## 1 被控系统建模
%[text] 考虑以下三阶模型：$G\_0(z) = \\frac{-0.00146(z - 0.1438)(z - 1)}{(z - 0.7096)(z^2 - 0.04369z + 0.01392)}$
%[text] 离线设计控制器理论上能够完全拒绝干扰的周期性部分：
%[text] -  $d\_s${"editStyle":"visual"} 是幅度为1的单频正弦信号，$v${"editStyle":"visual"} 是均值为零、标准差为0.02且 $|v| \\leq 0.1$ 的高斯噪声。
%[text] - 采样周期 $T\_s = 1/480 \\text{sec}\n$ \
clear; close all; clc;
addpath('..\functions'); %[output:20536631]
% 离散系统G0(z)
modelFile = fullfile('..', 'dataset', 'syn_TAC2015_3rd.mat');
load(modelFile, 'model');
G0 = model.G0;
Ts = model.Ts;
fs = model.fs;
z = tf('z', Ts);

% 计算灵敏度函数和互补灵敏度函数
% 假设单位反馈: S = 1/(1+G0), T = G0/(1+G0)
S = feedback(1, G0);  % 灵敏度函数 S = 1/(1+G0)
T = feedback(G0, 1);  % 互补灵敏度函数 T = G0/(1+G0)

% 计算频率响应
[H_G0, f_bode] = freqz(G0.num{1}, G0.den{1}, 1000, fs);
[H_S, ~] = freqz(S.num{1}, S.den{1}, 1000, fs);
[H_T, ~] = freqz(T.num{1}, T.den{1}, 1000, fs);

% 绘制三个传递函数的幅频响应
figure('Name', '开环与灵敏度函数幅频响应'); hold on; %[output:54280bf2]
semilogx(f_bode, 20*log10(abs(H_G0)), 'b-', 'LineWidth', 1.5, 'DisplayName', 'G_0(z) 开环传函'); %[output:54280bf2]
semilogx(f_bode, 20*log10(abs(H_S)), 'r-', 'LineWidth', 1.5, 'DisplayName', 'S(z) 灵敏度函数'); %[output:54280bf2]
semilogx(f_bode, 20*log10(abs(H_T)), 'g--', 'LineWidth', 1.5, 'DisplayName', 'T(z) 互补灵敏度函数'); %[output:54280bf2]
xlabel('频率 (Hz)'); %[output:54280bf2]
ylabel('幅值 (dB)'); %[output:54280bf2]
title('开环传函G_0(z)与灵敏度函数S(z), T(z)的幅频响应'); %[output:54280bf2]
legend('Location', 'best'); grid on; %[output:54280bf2]
xlim([f_bode(2), f_bode(end)]); % 避免DC点的显示问题 %[output:54280bf2]
%%
%[text] ## 2 求解准备
w_interp = 0.0521; % interpolation frequencies (rad)
N = 15;
nf = length(w_interp);
z_interp = exp(1j * w_interp);

% frequency grid for Hinf constraints
Nfreq = 1000;                     
w_grid = linspace(0, pi, Nfreq);
z_grid = exp(1j*w_grid);
% precompute frequency responses
G_grid = squeeze(freqresp(G0, w_grid * fs)).';        % rad/sample → rad/s
G_interp = squeeze(freqresp(G0, w_interp * fs)).';   % rad/sample → rad/s
% build Phi matrices
Phi_grid = zeros(Nfreq, N);
Phi_interp = zeros(nf, N);
for k = 1:N
    Phi_grid(:,k)   = z_grid .^ (-(k-1));      % 修正：标准 z^{-(k-1)} 基
    Phi_interp(:,k) = z_interp .^ (-(k-1));   % 修正：标准 z^{-(k-1)} 基
end

% build interpolation A * theta = b (complex)
A = diag(G_interp) * Phi_interp; % size N x N (complex)
b = ones(size(G_interp));        % complex ones

% if theta real, convert to real linear eqns
Aeq = [ real(A); imag(A) ];
beq = [ real(b)'; imag(b)' ];

% ---------- 可行性检查（先看是否存在满足等式的 theta） ----------
feasible_flag = true;
% 线性可行性：是否有 theta s.t. Aeq * theta = beq
% 简单尝试解方程并测残差
theta_test = Aeq \ beq;
resnorm = norm(Aeq*theta_test - beq);
fprintf('线性可行性检验: 等式最小残差 = %.3e\n', resnorm); %[output:4bd508b1]
if resnorm > 1e-8
    warning('插值等式可能无法严格满足 (残差 %.3e). 若要强制满足，可考虑增大 N 或放宽约束。', resnorm);
    % 依使用场景，可视为不可行 or 允许松弛
end
%%
%[text] ## 3 求解过程
%%
%[text] ### （1）最小二乘法求解优化问题
% ---------- 最小二乘基准解（仅供对比） ----------
theta_ls = theta_test;  % 使用前面计算的最小二乘解
fprintf('LS解的插值残差: %.3e\n', resnorm); %[output:09f458ac]

% 计算最小二乘解的H∞性能
mag_ls = abs(G_grid(:) .* (Phi_grid * theta_ls));
gamma_ls = max(mag_ls);
fprintf('LS解的H∞范数: %.6f\n', gamma_ls); %[output:1944086b]
fprintf('LS解的参数范数: %.6f\n', norm(theta_ls)); %[output:133cb305]
%%
%[text] ### （2）cvx求解优化问题
% ---------- CVX: 直接最小化 gammma，等式作为硬约束 ----------
if exist('cvx_begin','file') ~= 2 %[output:group:47dbbc57]
    fprintf('CVX 未安装，使用 LS 解。\n'); %[output:035dd64a]
    theta_opt = theta_ls;
    gammma = gamma_ls;
else
    r0 = 1e4;                     % theta norm bound
    cvx_begin quiet
        cvx_precision high
        variable theta(N,1)
        variable gammma nonnegative
        minimize( gammma )
        subject to
            for i = 1:Nfreq
                abs( G_grid(i) * (Phi_grid(i,:) * theta) ) <= gammma;
            end
            Aeq * theta == beq;
            norm(theta,2) <= r0;
    cvx_end
    theta_opt = theta;
end %[output:group:47dbbc57]
%%
%[text] ### （3）最优化问题总结分析
if exist('cvx_begin','file') == 2 %[output:group:8c48a9ac]
    theta_opt = theta;
    compare_CVX_LS(gamma_ls, theta_ls, resnorm, gammma, theta_opt, Aeq, beq);
else
    theta_opt = theta_ls;
    fprintf('跳过 CVX 对比（CVX 未安装）。\n'); %[output:8d534122]
end %[output:group:8c48a9ac]
% dense check validation
Mchk = 3000;
w_chk = linspace(0, pi, Mchk);
G_chk = squeeze(freqresp(G0, w_chk * fs)).';
Phi_chk = zeros(Mchk, N);
for k=1:N, Phi_chk(:,k) = exp(-1j*w_chk*(k-1)); end
mag = abs( G_chk(:) .* (Phi_chk * theta_opt) );
gammma_check = max(mag);
fprintf('CVX gammma = %.6e, dense-grid gammma_check = %.6e, residual ||A theta - b|| = %.3e\n', gammma, gammma_check, norm(A*theta_opt - b)); %[output:01378f4d]
%%
%[text] ## 4 抗干扰仿真分析
% 扰动参数 (Ts, fs 已从模型文件加载)
omega0 = w_interp(:); % rad/sample, 25rad/sec
Nsim = 480*50;   % 50秒
t = (0:Nsim-1)*Ts;
d = sum(sin(omega0*t*fs) + 0.02*randn(length(omega0),Nsim),1);
figure('Name', '扰动信号'); %[output:9ce7464c]
plot(t, d); grid on; %[output:9ce7464c]
xlabel('时间 (s)'); ylabel('幅值'); %[output:9ce7464c]
title('扰动信号 d(t) 时域波形'); %[output:9ce7464c]
xlim([0,5]); %[output:9ce7464c]
% 构造 FIR 系数（z^-1 顺序）
b_K_cvx = [0; flip(theta_opt(:))]; % CVX优化解
b_K_ls = [0; flip(theta_ls(:))];   % 最小二乘解

% 绘制频率响应对比
nfft = 1024;
[Hk_cvx, fk] = freqz(b_K_cvx, 1, nfft, fs);
[Hk_ls, ~] = freqz(b_K_ls, 1, nfft, fs);

figure('Name', '控制器频率响应对比: CVX vs 最小二乘'); %[output:65ac8863]
subplot(2,1,1);  %[output:65ac8863]
semilogx(fk, 20*log10(abs(Hk_cvx)), 'r-', 'LineWidth', 2, 'DisplayName', 'CVX优化解'); %[output:65ac8863]
hold on; %[output:65ac8863]
semilogx(fk, 20*log10(abs(Hk_ls)), 'b--', 'LineWidth', 2, 'DisplayName', '最小二乘解'); %[output:65ac8863]
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on; %[output:65ac8863]
title('控制器幅频响应对比'); %[output:65ac8863]
legend('Location', 'best'); %[output:65ac8863]
subplot(2,1,2);  %[output:65ac8863]
semilogx(fk, angle(Hk_cvx)*180/pi, 'r-', 'LineWidth', 2, 'DisplayName', 'CVX优化解'); %[output:65ac8863]
hold on; %[output:65ac8863]
semilogx(fk, angle(Hk_ls)*180/pi, 'b--', 'LineWidth', 2, 'DisplayName', '最小二乘解'); %[output:65ac8863]
xlabel('Frequency (Hz)'); ylabel('Phase (deg)'); grid on; %[output:65ac8863]
title('控制器相频响应对比'); %[output:65ac8863]
legend('Location', 'best'); %[output:65ac8863]
sgtitle('Jafari TAC 离散 — 固定控制器频率响应: CVX vs 最小二乘'); %[output:65ac8863]
% 构造K(z,theta)
% Kz = tf(theta', [1 zeros(1,N-1)], Ts); % 可选
% 闭环仿真 - CVX解
u_cvx = zeros(1,Nsim);
y_cvx = zeros(1,Nsim);
x_cvx = zeros(1,Nsim);
antinoise_cvx = zeros(1,Nsim);

% 闭环仿真 - LS解
u_ls = zeros(1,Nsim);
y_ls = zeros(1,Nsim);
x_ls = zeros(1,Nsim);
antinoise_ls = zeros(1,Nsim);

% 获取G0(z)分子分母系数（z形式）
b = G0.num{1};
a = G0.den{1};

% 控制器仿真
for k = max(N,length(a))+1:Nsim
    if t(k) >= 20
        u_cvx(k) = -sum(theta_opt'.*x_cvx(k-(N:-1:1)));
        u_ls(k) = -sum(theta_ls'.*x_ls(k-(N:-1:1)));
        if isnan(u_cvx(k)), u_cvx(k) = 0; end
        if isnan(u_ls(k)), u_ls(k) = 0; end
    end
    % 离散系统差分方程递推y(k)
    antinoise_cvx(k) = (-a(2:end)*antinoise_cvx(k-1:-1:k-length(a)+1)' + b*u_cvx(k-1:-1:k-length(b))')/a(1);
    antinoise_ls(k) = (-a(2:end)*antinoise_ls(k-1:-1:k-length(a)+1)' + b*u_ls(k-1:-1:k-length(b))')/a(1);
    y_cvx(k) = antinoise_cvx(k) + d(k);
    y_ls(k) = antinoise_ls(k) + d(k);
    x_cvx(k) = y_cvx(k) - antinoise_cvx(k);
    x_ls(k) = y_ls(k) - antinoise_ls(k);
end

% 绘图对比
figure('Name', '控制输入与输出信号: CVX vs LS'); %[output:1af9a3e1]
subplot(2,1,1); %[output:1af9a3e1]
plot(t, u_cvx, 'r-', 'LineWidth', 1.5, 'DisplayName', 'CVX控制输入'); %[output:1af9a3e1]
hold on; %[output:1af9a3e1]
plot(t, u_ls, 'b--', 'LineWidth', 1.5, 'DisplayName', 'LS控制输入'); %[output:1af9a3e1]
xline(20, 'k:', 'Controller ON', 'HandleVisibility', 'off'); %[output:1af9a3e1]
xlabel('时间 (s)'); ylabel('控制输入 u'); %[output:1af9a3e1]
title(sprintf('控制输入 (CVX RMS=%.2f, LS RMS=%.2f)', rms(u_cvx(20*fs:end)), rms(u_ls(20*fs:end)))); %[output:1af9a3e1]
legend show; grid on; %[output:1af9a3e1]
subplot(2,1,2); %[output:1af9a3e1]
plot(t, d, 'k:', 'LineWidth', 1.0, 'DisplayName', '扰动 d'); %[output:1af9a3e1]
hold on; %[output:1af9a3e1]
plot(t, y_cvx, 'r-', 'LineWidth', 1.5, 'DisplayName', ['CVX (γ=' num2str(gammma, '%.3f') ')']); %[output:1af9a3e1]
plot(t, y_ls, 'b--', 'LineWidth', 1.5, 'DisplayName', ['LS (γ=' num2str(gamma_ls(end), '%.3f') ')']); %[output:1af9a3e1]
xline(20, 'k--', 'Controller ON', 'HandleVisibility', 'off'); %[output:1af9a3e1]
xlabel('时间 (s)'); ylabel('输出 y'); %[output:1af9a3e1]
title('闭环输出性能对比: CVX优化 vs 最小二乘'); %[output:1af9a3e1]
legend show; grid on; ylim([-2,2]); %[output:1af9a3e1]
sgtitle('Jafari TAC 离散 — 固定控制器闭环仿真'); %[output:1af9a3e1]
%%
%[text] ## 5 自适应仿真分析
%[text] 起初，我期望使用2016年专著中书写的成熟的控制方法来进行设计，但是在该书中没有提供完整的控制器设计策略。应对论文中的控制器设计方法，最好还是先从论文复现出发。
%[text] 在TAC2015论文中还没有提出引入一个固定控制器的设计方法，而是直接引入一个参数$\\tau\_{0\\;}${"editStyle":"visual"}作为增益。
rng(42);  % 固定随机种子
d = 0.7*sin(25*t + pi/3) + 0.5*sin(225*t + pi/4);
idx_add = t >= 30;  % step activation
d(idx_add) = d(idx_add) + ...
             0.6*sin(85*t(idx_add) - pi/6) + ...
             0.4*sin(125*t(idx_add) + pi/2);
v = 0.02*randn(size(t));         % std = 0.02
v = max(min(v,0.1),-0.1);        % |v| ≤ 0.1
d = d + v;
N = 50;
% ---------- 1. F(z) 滤波器 ----------
% 注: G0(z) 在扰动频段增益仅 -73~-58dB, 纯延迟基自适应控制器在此条件下
% 无法获得有效抑制 (详见 DebugReport_JafariJVC_Discrete.md §5.2).
% 此处保留 F=100 静态增益作教学对照, 正确方案 (Lambda基+F(s)频谱平坦化)
% 参考 demo3_Robust\demos\Jafari_AdaptiveCC_Fixed.m
F_gain = 100;
b_F = F_gain;
a_F = 1;

% ---------- 2. 初始化自适应仿真变量 ----------
% 关键信号
kz_adapt = zeros(1,Nsim);
antinoise_term = zeros(1,Nsim);
u_adapt = zeros(1,Nsim);
y_adapt = zeros(1,Nsim);
z_adapt = zeros(1,Nsim); % 观测信号 z = y - G0*u
ep_adapt = zeros(1, Nsim);
phi_hist = zeros(N,1);    % 回归向量

% 自适应参数
theta_adapt = zeros(N,1); % 参数初值为0
P = eye(N) * 100;         % 协方差矩阵 P(0)=100I (适度增大以适应参数尺度)
gamma0 = 1;               % 归一化增益
theta_max = 100;           % 投影界 (F_gain=100 补偿 G0 低增益，θ 约需 12–50)

% 历史缓冲区 (用于IIR滤波器)
hist_len = max([length(a), length(b)]) + N + 1;
phi_G0_u_hist = zeros(1, length(b));
phi_G0_y_hist = zeros(1, length(a));
phi_N = zeros(1, Nsim);
F_u_hist = zeros(1, length(b_F));
F_y_hist = zeros(1, length(a_F));

% 监控变量
theta_history = zeros(N, Nsim);
rms_history = zeros(1, Nsim);

% ---------- 3. 自适应控制仿真循环 ----------
% phi_adpat   N×1
% theta_adpat N×1
for k = hist_len:Nsim %[output:group:7a971160]
    % --- 步骤1: 计算系统输出 y(k) ---
    % y(k) = G(u(k-1)) + d(k)，这里假设 G=G0
    % 使用差分方程直接计算
    antinoise_term(k) = (-a(2:end)*antinoise_term(k-1:-1:k-length(a)+1)' + b*u_adapt(k-1:-1:k-length(b))')/a(1);
    y_adapt(k) = antinoise_term(k) + d(k);
    
    if t(k) >= 20  % 20秒后开启自适应控制
        % --- 步骤2: 计算观测信号 z(k) ---
        % z(k) = y(k) - G0(u(k-1))
        % G0(u(k-1))与上面计算的antinoise_term是同一个值
        z_adapt(k) = y_adapt(k) - antinoise_term(k);

        % --- 步骤3: 回归量 phi = F * G0[z_delays] (F=100 静态增益)
        [phi_N(k), phi_G0_u_hist, phi_G0_y_hist] = IIR_filter(b, a, z_adapt(k), phi_G0_u_hist, phi_G0_y_hist); %[output:71f946c6]
        phi_hist = F_gain * phi_N(k-(N:-1:1)+1)';

        % --- 步骤4: 更新自适应参数 theta(k) ---
        % 归一化RLS + 参数投影
        m_s_squared = 1 + gamma0 * (phi_hist' * phi_hist);
        epsilon = (z_adapt(k) - theta_adapt' * phi_hist) / m_s_squared;
        ep_adapt(k) = epsilon;
        P_phi = P * phi_hist;
        denom = m_s_squared + phi_hist' * P_phi;
        K_rls = P_phi / denom; % Kalman增益
        P = P - (K_rls * P_phi');
        theta_new = theta_adapt + K_rls * epsilon;  % 修复: 用Kalman增益
        th_n = norm(theta_new);
        if th_n > theta_max, theta_adapt = theta_new * theta_max/th_n;
        else, theta_adapt = theta_new; end

        % --- 步骤5: 控制 u = -F * theta' * z_delays (F=100 静态增益)
        z_delays = z_adapt(k-(N:-1:1)+1)';
        K_z_k = theta_adapt' * z_delays;
        kz_adapt(k) = K_z_k;
        u_adapt(k) = -F_gain * K_z_k;
    else
        % 控制器关闭时，系统自由响应
        z_adapt(k) = y_adapt(k); % z=y
        u_adapt(k) = 0;
    end
    
    % 记录历史数据
    theta_history(:,k) = theta_adapt;
    if k > fs
        window_start = max(1, k - fs + 1);
        rms_history(k) = rms(y_adapt(window_start:k));
    end
end %[output:group:7a971160]

% ---------- 4. 性能分析与可视化 ----------
% 自适应控制性能对比
figure('Name', '自适应(Jafari结构) vs 固定MOSEK控制器');
subplot(3,1,1);
hold on;
plot(t, y_adapt, 'b', 'LineWidth', 1.5, 'DisplayName', '自适应 (Jafari)');
plot(t, d, 'k:', 'LineWidth', 1.0, 'DisplayName', '扰动信号');
xline(20, 'g--', 'Controller ON');
xlabel('时间 (s)'); ylabel('输出 y');
title('输出性能对比');
legend show; grid on; ylim([-2,2]);

subplot(3,1,2);
hold on;
plot(t, u_adapt, 'b--', 'LineWidth', 1.5, 'DisplayName', '自适应 (Jafari)');
xline(20, 'g--', 'Controller ON');
xlabel('时间 (s)'); ylabel('控制输入 u');
title('控制输入对比');
legend show; grid on;

subplot(3,1,3);
plot(t, rms_history, 'b-', 'LineWidth', 1.5);
xline(20, 'g--', 'Controller ON');
xlabel('Time (s)'); ylabel('Rolling RMS');
title('自适应控制 RMS 收敛轨迹');
grid on;
sgtitle('Jafari TAC 离散 — 自适应 RLS 控制器时域仿真');

% 参数演化对比
figure('Name', '控制器参数演化对比: MOSEK vs 自适应(Jafari)');
param_indices = [1, 5, 10, 15];
colors = {'r-', 'g-', 'b-', 'm-'};
for i = 1:length(param_indices)
    subplot(2,2,i);
    idx = param_indices(i);
    plot(t, repmat(theta_opt(idx), 1, Nsim), 'k--', 'LineWidth', 2, 'DisplayName', sprintf('MOSEK θ_%d', idx));
    hold on;
    plot(t, theta_history(idx,:), colors{i}, 'LineWidth', 1.5, 'DisplayName', sprintf('自适应 θ_%d', idx));
    xline(20, 'g--', 'Controller ON', 'Alpha', 0.7);
    xlabel('Time (s)'); ylabel(sprintf('θ_%d', idx));
    title(sprintf('参数 θ_%d 演化', idx));
    legend('Location', 'best'); grid on;
end

% 所有参数的范数演化
figure('Name', '控制器参数范数演化 (Jafari)');
theta_norm_history = sqrt(sum(theta_history.^2, 1));
plot(t, theta_norm_history, 'b-', 'LineWidth', 1.5, 'DisplayName', '自适应参数范数');
hold on;
plot(t, repmat(norm(theta_opt), 1, Nsim), 'r--', 'LineWidth', 2, ...
     'DisplayName', sprintf('MOSEK参数范数 = %.4f', norm(theta_opt)));
xline(20, 'g--', 'Controller ON');
xlabel('Time (s)'); 
ylabel('||θ||_2');
title(sprintf('控制器参数范数演化 (终值=%.2f)', norm(theta_adapt)));
legend('Location', 'best'); 
grid on;
figure('Name', '输出 RMS 收敛');
plot(t, rms_history, 'LineWidth',1.5); grid on;
xlabel('时间 (s)'); ylabel('RMS(y)');
title(sprintf('输出滚动 RMS (终值=%.4f)', rms_history(end)));

% ---- 性能统计 ----
rms_d_s5 = rms(d(20*fs:end));
rms_adapt_s5 = rms(y_adapt(20*fs:end));
fprintf('自适应: 抑制=%.1f dB, ||θ||=%.1f\n', 20*log10(rms_d_s5/rms_adapt_s5), norm(theta_adapt));
%%
%[text] ## 6 结构鲁棒性测试：从回归量中移除 G₀(z)
%[text] **对照实验**：使用 φ = τ₀·w(k) 作为回归量（不含 G₀(z) 滤波），
%[text] 验证 G₀ 滤波的不可或缺性。**预期结果：抑制 ≈ 0 dB**——
%[text] 证明移除 G₀ 后 RLS 无法学习逆对象动态，自适应算法失效。
%% 1. 系统定义 (含不确定性)
clear G0 z S T
Ts = 1/480;  z = tf('z', Ts);  fs = 1/Ts;
G0_tf = (-0.00146*(z-0.1438)*(z+1)) / ((z - 0.7096)*(z^2 - 0.04369*z + 0.01392));
[num0, den0] = tfdata(G0_tf, 'v');
Delta_m = -0.0001 / (z + 0.99)^2;
G_true_tf = G0_tf * (1 + Delta_m);
[num_true, den_true] = tfdata(G_true_tf, 'v');

%% 2. 仿真参数 (专著 Table 4.1)
T_end = 50;  Nsim6 = T_end*fs;  t6 = (0:Nsim6-1)'*Ts;
% 多频扰动 + 30s阶跃
d1 = 0.7*sin(25*t6 + pi/3) + 0.5*sin(225*t6 + pi/4);
d2 = 0.6*sin(85*t6 - pi/6) + 0.4*sin(125*t6 + pi/2);
idx_30s6 = find(t6 >= 30, 1);
d_total6 = d1;  d_total6(idx_30s6:end) = d_total6(idx_30s6:end) + d2(idx_30s6:end);
rng(42);  v6 = 0.02*randn(Nsim6, 1);  v6 = max(min(v6, 0.1), -0.1);
d_total6 = d_total6 + v6;

N6 = 50;  tau0 = 10;  c0 = 1;
P0_val = 500;  P6 = P0_val * eye(N6);
theta_max = 10;  theta6 = zeros(N6, 1);

%% 3. 初始化
y6 = zeros(Nsim6,1);  u6 = zeros(Nsim6,1);
zeta_hist6 = zeros(Nsim6,1);
nG0 = max(length(num0), length(den0))-1;
nG = max(length(num_true), length(den_true))-1;
G_state6 = zeros(nG,1);  G0_state_zeta6 = zeros(nG0,1);
buffer_zeta6 = zeros(N6,1);  u_prev6 = 0;

% 结构测试: phi = tau0*w (无 G0 滤波), 预期抑制≈0 dB

%% 4. 仿真循环
for k = 1:Nsim6
    [y_no_d6, G_state6] = filter(num_true, den_true, u_prev6, G_state6);
    y6(k) = y_no_d6 + d_total6(k);
    [G0_u6, G0_state_zeta6] = filter(num0, den0, u_prev6, G0_state_zeta6);
    zeta_k6 = y6(k) - G0_u6;  zeta_hist6(k) = zeta_k6;

    u_current6 = 0;
    if t6(k) >= 10
        w_vec6 = flipud(buffer_zeta6);
        phi_vec6 = tau0 * w_vec6;  % phi = tau0*w (不含 G0 滤波)

        pred_err6 = zeta_k6 - theta6' * phi_vec6;
        m2_6 = 1 + c0 * (phi_vec6' * phi_vec6);
        epsilon6 = pred_err6 / m2_6;
        P_phi6 = P6 * phi_vec6;
        denom6 = m2_6 + phi_vec6' * P_phi6;
        if denom6 > 1e-12
            Kgain6 = P_phi6 / denom6;  % Kalman 增益
            P6 = P6 - (Kgain6 * P_phi6');
            theta_new6 = theta6 + Kgain6 * pred_err6;  % 使用 Kalman 增益 × 预测误差
        else
            theta_new6 = theta6;
        end
        if norm(theta_new6) > theta_max
            theta_new6 = theta_new6 * (theta_max / norm(theta_new6));
        end
        theta6 = theta_new6;

        u_current6 = -tau0 * (theta6' * w_vec6);
        u_current6 = max(min(u_current6, 100), -100);
    end

    u6(k) = u_current6;  u_prev6 = u6(k);
    buffer_zeta6 = [zeta_k6; buffer_zeta6(1:end-1)];
end

%% 5. 性能报告
idx1_6 = (t6>20 & t6<30);  idx2_6 = (t6>40 & t6<50);
rms_d1_6 = rms(d_total6(idx1_6));  rms_y1_6 = rms(y6(idx1_6));
supp1_6 = 20*log10(rms_d1_6/rms_y1_6);  % 正值=降噪
rms_d2_6 = rms(d_total6(idx2_6));  rms_y2_6 = rms(y6(idx2_6));
supp2_6 = 20*log10(rms_d2_6/rms_y2_6);

fprintf('结构测试: phi=tau0*w, S1=%.1f dB, S2=%.1f dB (预期≈0 dB)\n', supp1_6, supp2_6);

% ---- 图: 结构鲁棒性测试 ----
figure('Name','结构鲁棒性测试');
subplot(3,1,1);
plot(t6, d_total6, 'Color',[0.7 0.7 0.7]); hold on;
plot(t6, y6, 'b-', 'LineWidth',1);
xline(10,'r--','控制 ON'); xline(30,'r--','新扰动加入');
ylabel('幅值'); title(sprintf('输出信号 (S1: %.1f dB, S2: %.1f dB)', supp1_6, supp2_6)); grid on;
legend('d(t)', 'y(t)');

subplot(3,1,2);
plot(t6, u6, 'g-', 'LineWidth',1);
xline(10,'r--'); xline(30,'r--');
ylabel('控制 u'); title(sprintf('控制信号 (RMS=%.2f)', rms(u6))); grid on;

subplot(3,1,3);
plot(t6, zeta_hist6, 'k-', 'LineWidth',1);
xline(10,'r--'); xline(30,'r--');
xlabel('时间 (s)'); ylabel('ζ(t)'); title('观测信号 ζ(t) = y - G_0[u]'); grid on;
sgtitle('Jafari TAC 离散 — 结构鲁棒性测试: 移除 G_0(z) 滤波后自适应失效');
%%
%[text] ## 附录
function compare_CVX_LS(gamma_ls, theta_ls, resnorm, gammma, theta_opt, Aeq, beq)
    fprintf('CVX vs LS 对比:\n');
    fprintf('方法          | H∞范数    | 参数范数   | 插值残差\n');
    fprintf('-------------|----------|----------|----------\n');
    fprintf('最小二乘(LS) | %.6f | %.6f | %.3e\n', gamma_ls, norm(theta_ls), resnorm);
    fprintf('CVX优化     | %.6f | %.6f | %.3e\n', gammma, norm(theta_opt), norm(Aeq*theta_opt - beq));

    % 性能改善分析
    improvement_ratio = gamma_ls / gammma;
    fprintf('\nCVX相对LS的H∞性能改善: %.2f倍\n', improvement_ratio);
    if improvement_ratio > 1.2
        fprintf('✓ CVX优化显著改善了H∞性能\n');
    elseif improvement_ratio > 1.05
        fprintf('△ CVX优化略有改善H∞性能\n'); 
    else
        fprintf('○ CVX与LS性能接近，可能H∞约束不起作用\n');
    end
end

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":40}
%---
%[output:20536631]
%   data: {"dataType":"warning","outputData":{"text":"Warning: Name is nonexistent or not a directory: \/home\/dcol\/Projects\/MATLAB\/ANC-Classic-Control-Simulations\/..\\functions"}}
%---
%[output:54280bf2]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAjAAAAFRCAYAAABqsZcNAAAAAXNSR0IArs4c6QAAIABJREFUeF7snQm8ncP9\/79XBJFNlNqCxBZapH6KWEKCqlIqLWopEiqotWiUWmIPUaq1RiWhVUuL2IoWiVL7EktJrKnYSiQVQRKS+\/+\/J+b0uSfnnjvnnDmTc879jFdekXO+z8zzvGfOzOf5zndmmpqbm5tNSQREQAREQAREQATqiECTBEwd1ZZuVQREQAREQAREwBGQgFFDqAkC77zzjr322mu2+eab2xJLLOHu6YknnrA5c+bYZpttlvusJm62Sjcxbdo0e++996xXr17WtWvXKpXSMtu\/\/OUv1rNnT+vXr99C5Y0bN87mz59vP\/zhD4Puhbr65z\/\/6Wy33XbboGvKMfr8889dW7n\/\/vvtRz\/6kXXo0MFWXHFFW2yxxezBBx+0ZZZZxvr27VtO1g19zfTp0+1vf\/ub\/d\/\/\/Z+ts846uWe98cYbXZvjd9ZauuGGG9xXe+21V1mMuH7mzJk2dOjQsq7Pv+j11183\/nzrW9+yr3\/960F53nfffa59\/vKXv7Qll1wy6BoZ1TYBCZjarp9kdzd8+HCbNWtWSeWdcMIJ9ve\/\/92eeeaZFtctt9xyxmBcLA0cONB23nnnnMmvfvUru\/766+3SSy+1733ve+7z73\/\/+66Tev75561jx44Fs7v99tvdPRRLiIGTTjrJunTpUtLzpTYePXq0nXXWWXbxxRfbLrvsElT8lClT7Le\/\/W2QLcJwxIgROdtPPvnEDWbf+MY37Lbbblsoj6222sp99vDDDwflP2PGDNtiiy2c+GKgaK3OyKycemO2GzFF+xowYICdcsopdvXVVxttd9lll7WbbrrJCZdvfvOb9uc\/\/znontuT0QMPPGA\/\/elPbb\/99rPTTz\/dPfqnn35qG2+8sROxDPCtpR\/84AfuK9rJyy+\/7LgXS\/QNX\/va12zSpEm20kor2SGHHGLvvvturi0NGzbMieNCiXrt3r170fzPOecc+\/3vf28nn3yyHXjggUHV+Otf\/9r1LyeeeKIdfPDBQdfIqLYJSMDUdv0kuzs6MQagUhIdIoMnb+rZtPrqq9u\/\/\/3volnRodHJkT788ENjsGTg69+\/v3ujJt1yyy3WrVs323777XN58ZZ9\/vnn5\/593nnn2ZVXXlm0rKWXXtpGjRrlBtdaTj\/\/+c\/dAIG469GjR6u3usMOOxh\/SM8++6zzQoQkOLz44os50zvuuMOOPvpoo76o\/2zaZJNNXEdPOVdcccVC2T\/yyCM2duzYhT5\/8skn7eOPP3Zv8\/lepG222cb23Xdfd0259eYHoX322cf+9Kc\/2Xe+8x0nYBEx3DOi+Cc\/+YmdccYZIUjalc0ll1xiF154oY0ZM8aoCxKiBa\/I3nvvbWeffXarPLICBi\/X4YcfXpTdX\/\/6VydgNthgAzvqqKOM9oLw8W312muvbfV6xO8999xjTz31VM4GjxFt1acdd9zRXn31VXvsscds+eWXXyivU0891T777LMWn9Mu8dohdhHA+emCCy5oV+2hER5WAqYRajHCMyAifBo0aJDNnj3b7r777txnl19+uRuwLrvsstxgR0cwd+5c+\/LLL22nnXZyHg7egvk315N4a6ZzZGrkzjvvzA1qnTt3duKERMdBvnR0bXkTEDd0XMUSzzJy5EhjemTllVe2q666ytZbb72FLkGw\/e53v3Od2vvvv+\/ujcGcARAGIQmv1Xbbbee8Rbw5FksMqnTsDBqFvEHk8+abb7ZZLJyOOeaYnB2dOF6qYum6666zjz76qIWA+dnPfuYGikIJT8Zzzz1n6667rvPQ+LThhhva\/vvv79jyFl1KamuQbKveGDiZ7uCem5qaXNvy6aCDDrJ58+a5Nsq9r7nmmrnvEK54bsj\/N7\/5jY0fP955CBF0CDd45k85IXgZ7PE+4fFpLZHPlltu6VhwDyFp8uTJOS9ja\/b8JvzUXaX3AjfEKlOyb7\/9tmureOP23HNPJ5jxfOK1YhrJJ34Hxx9\/vCF68PLxGyHRRvGoIFCzAiP7HPQFCA68O1kBk7WnbcGf3wPChJcS\/5vjN37ssce2eDGCMWVyL9Qz983vFQHrE\/0RnlYS7bRUj\/Ibb7wRUn2yqSECEjA1VBm1cit4APjx8zbtk3\/z5a23ULwE13Tq1Ml1SNlER0XHRMeHkMhPdKh0eCTe0vygdOaZZ7pB8o9\/\/KPrBLPJv9lzHa5v3MmIITplyuMtn86Tt3QGW7w2f\/jDH5xwYOqCAYlyd999d\/vggw9c1ng8iOHwb2242n1nWKxe8CTgocAb1bt376JViDghNgT3dX7eTOfQ6dJJI6qKJTp+BKBP3p3eVvvJemBeeOEF23XXXR2\/H\/\/4x+7ZEZ78Ia4AVoWmlZjeww3vBQwDDxwRuIhGBAaeEBLxEk8\/\/fRCgrPcervmmmtc2fmJeAYGNZ6POIv8hHBCPHLvxM4stdRSttpqq7kpDdo5YpJ24QdwRAlv6AyOF110UVtYDc8ZvGgDhbwB+RnQDvBA+sQ9kRANftoNQYTnMca9IICyU4e+XAQK3kvaXn7aaKON7Oabb3YiJ1+oIE74vbfm+fRToPkCJjuFhHeV3\/u\/\/vUvQ9DxrPyGEctMCSFgEF28rFAW7ZQ22ppo4v55WfHTnfyWqAvaqU+0bTyWfOdjYPg37aZPnz5FvZ5tNgIZLBICEjCLBHttFxpTwBxxxBFO1DA1RKedncbBe0MHifeAwZ9BAzcv6R\/\/+IfrvHFd+yklT817AdZff30XsElae+21nZDA\/Yz4oAOkbATMo48+6gYsBo7HH3\/cdWy4sn1nxpv2Gmus4QZvOnvefhFSvKln30rza43BkqkvyqazD0m8VdMpc5\/eA8V1lMUbPO58ggxjJ8QE00d+Cglxh+emUGIw4X4QBcSS8HyIEsQVIpTO3wsY\/xbMIMRgtPXWW+c8Frz5U4d+2mC33XZzvCqpN+6JqY7s9BVeImIhDj30UCdIqEOSZ0nMB4INBgxyeBa5bwZEPkP0ZgUrogURiWCi3baVHnroITvggAPsyCOPdGKm1IS3iPbGfTGQZlOMe\/ECBkHKywCDPL8Ln\/I9eogIPFOIfn6fiy++eE7kILIQXF7A+N8TefkYGwTMt7\/9bedVoUzqB+8XwddM\/ZAmTpxoeDwIvvYJLy0vOsS2tCZgiIlDFOcnrqM+swKGuiauDKGEyOa++ENA8aabbuqy8F5IXpRqfYq51HbVHuwlYNpDLZf4jMyP44kI8cAwSDFtxIDH2yMDFFMpdEz58RkE5jGNtMoqq7g7Ou2001wnSULA0MExldNW8l4ABkI6f4QPAw6dIwKBgYsBjcHu3HPPtSFDhjhPy1133eUEDG+CeIXomHlzRtxkE4M298i0QtbTkX9ffhBHcFBmsXgUBA5vtd5Tgdcku6LDT6Ph\/UBMFUu8AfvVHIgyBFhbCbECXy9g8JQQJ\/DFF1+4mABECZ04A0E2wY\/Ons\/p7HlrJ2WnkPBgMHiQF2+zeLxIlEcd+OkyvE48c7n1hkhCvPA3K0+oU4QHz5L1IhC7gzeOQFXvEULAIIapT6YpaePcJ\/wYFBG1fhqJNkx+vO3TphEniJT8hLhB5CDEGfDJuzVRWKx+igmYGPdSTMDgjeR3yG\/K1xv3g6cQsYqAoW14zxa\/Fe8VwQNDzJH3ZvDbRSwgEhBiXqwUenb\/UsLUJh4XpvaYFsW725aAwdOX9fRSl1xfSMBQr0yd0ZfxooRw80KTtkkQO+0WUYQHWam+CEjA1Fd9VeVucdVmpwv8W4xfhUKhdPIE5jLAMuCR6MCJRSCAMpt8jAdiArFA8gMOb3cM5gwMvK0zUODSZ5DDrczKBAZHBk5WJhVawsu12PuBkHvnjY1Ojc4QT8yqq67qpol4BlZM8GbsBQx2DG68GeZPeZUCmI6QPIkhwDvBG2U2EJJBEw8EgymBptwTHTwxNkyjEF\/gE9M4WcFY7D6y11YSxIvo4w2fAYc4BAIz\/RQNAZhMASD0eDPlWbMrxBAGCAXEE54VL8yw582XdOutt9p\/\/vOfhaaQyq03gsUp7xe\/+IVbDszAhNhkAGU6DM6ICAYlxDLtiTrgMwYrppL8FARThizZR\/wyeHuhSjtn6oiBlOlSEoy88IMJ5ZIQYwhRkq+\/e++913msSkmtCZhY9+IFDJ416hVPBx5O7pkXFX7\/3guDGON3gTChfbQmYOCN4C2UmK5D8JIXwoCYGf7Nb5nfNOKF3ze\/cwKKERXcY6iAQbTi3UHov\/XWW+4FhRidrIAh5o2XIp6R3zr2PBMilT+0G\/\/b8S9EpdSZbGuDgARMbdTDIr0LP5CVehN0DnREpGwMDG58OveXXnrJeWJ4M6OD5I0Z961319OB8YZEHAmdjw8UZCoFj07+ahw6RAZIPidlB0ICcplCoFNiKsgn7\/nwYgMPDB0vz8zA7DthP1hkGTCgH3fcccYeNXSufq6eZ2bQZEBAtGTd6P56hMrgwYPdvfhBnu8YPAhIpHP1z8vnCLnsfReqC54tfxmsv5YBifvhvhi48ZRkpwYYgHnTxAtEwnPA2y\/TAXTsCDGfeCtnACAQm\/tlkMdTxjMxGJFYRYJw9aLGC5j8+y4UdF1uvVH\/DI4IJO4dkcJyc0QvdYCwe+WVVxxXBAlxE4gon2gjxHIR75IN8KQNMiXFM8OBtlIo4BgRhEDGG8MgzyDohQ9iDlGHJ4DBvZTUmoCJdS+txcDQjpne22OPPdw0GoINBrQRplwY9AsJGLyufjVZoedESPLCwG+aqSReepgypb+gPyBgn8\/YfgEvDm2a9njYYYcFeWCwp93RpvkNE5ScFTD8jvBSMtVMWfze8NJQx4htVkMx3UmbRZyWU2el1K9sq0dAAqZ6bOsyZ1zFdC54LrJLHUODeHmbQ4AwoDINw3QSwoLBlBiY7373u84z4r0WQMLVzOBKJ8RATqfHvwsl3uL8ckc\/EDKQsiSbwZgyERy8RXMPeGoYWOio8j0wWRExderU3CoSBm3etL2A4f4ZGJm2YdqFKSMGMjjBiw4xuzEW0xJ0lDyLz8M\/CwMnHSsxGAwapSRiJMgPccJ0jk+4+hmMvFDydZVdyeI9BDDmuRmk8T5lV\/Iw8CAS8HAhFFdYYQUn9PicZ\/WbmVGuF70MdD4OqbVnYVoJr5NfzVNuvSGE\/XLe\/LL80mk\/6PM97SIb6+GvoW0xeBIbRTskPoPnQ\/gxkPLGnl3m76\/zU554MfAG+alQvveB1LR3PAKlpNYETKx7aW0KiTZDu0C88xskvoW2hLcE8U1bLyRgsmK30HP6YHFfF7Bh+pCy8Mjw20JQIABps3j\/eNEgVilkCqktAUMcHSLMe9EQOPQN\/N6oW9oxvyUCpdmnBs9nNh6tlLqT7aIlIAGzaPnXXOm8XdEB+w7M32CogCFojvgT5puJjWFwoOPy3gDeXhkIePPyrnb2GuHtlsGVv\/P3b6Ajp6OjU2LA8Tv1MhB6WwZdBnYGZx8Dw\/cMYLyp8+bHdAJv6kzr+E306GTzl1hzPdNgXnwgOAh6RbgwsDNg4\/rGFY1I4W3Oz5\/TeSKyEDEMCrwpMs2WL2CYwqAzx9OEN6GthNeJN0s6YOrCLzllUPdeCDwkCCm+o7Omg0Z8+aBNRFU2XgY7pvioB56BemEqzg\/8PIOPifGDnb9PBje8Vyy\/RfwVWv3jbakvBgy\/SqvcesNTwv3CgrgcvEJ4ZHiz9p423vqpL0RK1jNGG8LzQp1kNzFjOghvIfFQiFTqlRVw+avQvMeJdoYgRvBkkxcwrYmmYvXbmoCJdS9tCRjuzW8kSRtnxZL\/vWYFDL8f\/4LANfzO8cTw+\/H9BfFCtGv+8BJBHBTTaryU4Dmj\/uhf8O7RPzAdSHm8bFA3WQGD9wyvGr8RYoyoW8QWdYPXBO8ovy\/6l6wHhhckpqcRu7y4UJ+wpL3Spv1LFVPN+b+Jtn6H+r62CEjA1FZ9LPK78YN3\/nLpYgKGTsov+WSQZSBk4KOz4+0mK2B4QDqirMeCa\/EE+CkVRAydDLE2xDMw4DB1QYeWXabKQMi\/+UPHySoIYi4YyOm06FzZM4Tdbf29+EGNQYvBi8GDGBBic7h3BlquReQgYBAArNygo8XFjXcJ0YDw8VNIiD7iWxg0EUyUgbiBYX5gIN4hYodwcXNviJe29o+BmV8OzmDNm6vfeA7RhueAaTimAvyKGK7B+8TzML2BxwFRll31wcoixBVvwwgAeLPdPPcHOwYTPHF4jbxA8A0UfuTJoMYUDlMRvOVyHwgjpmt4s0Y8MThm90ippN4YuAjCRUggFBkEuUcECol4Jx+H5IOr+dzH7CBm\/VQJ7YyAUAZ4vwTXew14Bto8ifrF40YbzcYBZX+stA8Y+e\/xMpA\/9V9sR2LyaGsKqdJ7CREwTO8yuPNCQF3xG+Y3lS9gEKO0D36rtHFEKd4onoF\/U9+0F9o2zMgDL51\/SaBNEr+EB4+6Q7TQLhHc\/NazAia7QSbTQFkBQ735hKDJChgvSr0Xjfz5f9om+SC2aS94YrK\/pUXe+eoGSiYgAVMyssa9gCkFlh7TKfGjz3oOWhMwDCisGvCeEAYTBi4\/wBYSMPkE\/Vu+FzAE5jElwJsUiTgDOrt8T0k2loLBhUBaNmtjBQSucO\/ip0NlEPdTSAgehAReFbwHJKZ0ECB+Az46Wt7Q6JTzBQwiidVUfrMtPx2WHTxZfp1dwYQIIzbCC4zsVFhoi8LFTiwBb7d+p16mfKg33iSz3Hj7ZABgGggBh8jwu\/dSnp8qbK1sBh42DWPwJsGU\/BCesIYLos1PJ2BLnfGGTRtA2OCN4p4ZDP0KF\/KqpN64HoGKd8kHiPvBlqW5iGWCwrlH6hOvCx4mBj+YM91HQhDTLvzqJS94eQ5W4GSDeP2maOSb3SCP+vXTan5ZPvVLu\/Mrlxios0dmFOLdmoCJdS9ewPBceA8RKwiNfK8a9+a9PlyDQEfU0A\/4XXAJmPcruuBBe6CefWA9PBE5xKfQNomBQswxZem9r7ChDnkhYcoZO8Q0QprfOB627DJqzwwx5T0wPgbGf5ddRu1\/h5499cw9+v6A3z4saMP0K0r1S0ACpn7rLuqd07ngYUCIZDcj84W0JmCIs2Bg4y0LD0A2LoBryxEwXEdHiIsYIUXMAS78QgKGDtKvluKNN3\/PGH\/\/BMDiOclOK1AGnR1vi4glOuq11lrLub7xQPjzWBhwGYiJs2GgYhDgOr+UmIBFOmjvHShUMf7MFh\/smg3sDalIvFyIQuJIcLdnE8KFjph7YqBmwMRLwXQeUxsM5nTYDABMw7W22ZqfBsGjQ+eP6563bcQdIgnOfE4wMx4YAl2pI1anMf2CTXa5MXEFTEVQJiLCL6eGZzn1hoeIuBWmAniLR8Qx+NFmqQMEJ5zwUjF1iXeKe2ew497xCsEfoUx9M9ByXwiv7Go3nhkPHF4BbFpb1u7jmPzmf4hFL8K9gCGItNhyYuqx2DLqGPdSLIgXwcf9M4VGwtOHuMlfTZU9SsALGAJo8VzyAoPI598Iezj7VVsIFQQmwiMrZhFGlIWYoL3yG2SlIm0WsUz90qZ4ocoXMHhYuOesFzcrYJjO4jeCB5V4tvzEvVC3+Ts204by+6+Q36ZsFh0BCZhFx75mSkaAMBDxls5AW2hL9NYEDG83vIlmvTXZBwsRMAyuXJ9dlePz8IMqgyFTRH4JN99nYylCYRZaMdTWtXhP6Ez5mxgUAhKZx+fN0nfieGFCkl+dxX0U23eCAZi3UNgSFIm3gY6bgZEBgcRgzLSF3wIdoYcI5d58fXANb9rZgGEEBSLVD1rEkeC5YrBB1HkPFXEmDE4kxA8Bj8S88AfBhJjLbuCHUGIagMEEb8iECRNarKzyUz0MdvlxTm2xgxdTAIgrEt4Bpopof3Di+fACEE\/h44N8m8XbghjzsVNtleXPDMpueFbsGvZ+wZvnN0\/0ttwj95x\/zlRb5We\/j3EvXsAwXcJvxifqHyHpB3oGb4QLvzXaTdZrVkjAtPYctFeEJsIQsUM5iAPaJsKQ+CleTHjhYHqRqSZEE99TJtOq2fvMFzDZjeyoI+6V\/Gi3xEXRzvwRGL69hDAvFA8Xcp1sFh0BCZhFx75mSmYQoBNjkGrN3U0niBeEN3C\/FDfkAXgT402Yzow\/hRIDLp1coaMGsGcgYvDBy1FuYtBHIDGlkBVBIfnxlsjgxLJwvDN+GTXXMkDBhg6bN8hiCQ+Q35HYbwhXzN7vYowN4gDuBD97DwpeJWKWyJNOO7upWH6+DPAM4nT2rGDK7vEDd+qVN1oYkz9B2Nyj92jRRhiEiDdA5PDMxBEwlcN0mV9+65dZUz6DFW\/CtAGChBF++WcOtcU\/W2+ILIQiIpAYIp8Qe4hLbPPbJkKN+2N6MTSRH6KaQZvpt7YSop+pENqXX2mF0KW9cL9ZIdBWXvnfx7gXxDf1ibAsdEgoIsAfbYFXienW7CaL3BPtn4TnEO8J4hGb7EGrfE+7QFAjJJmWQ6jSjrOeDfLCc0Z8WVaoIEz8LsuFVgUhEsk7u1mgP4wSbylbHfAb4GWM+kPI4pEMTbDxoj70GtktWgISMIuWv0qvcwIM7HSWuOLpOIsl9tXgLRghUegwx0LXEsdBZ0\/HmnWZ1zm2mr99pqPwfiBM80\/Vzt48sRVsiIe4y65uQjQwNdTalGYpACq9l5CymJLBA4sXqzVvakg+tWDj459CPW61cM+6h\/IISMCUx01XiYAIiIAIiIAILEICEjCLEL6KFgEREAEREAERKI9AwwoY4haYh2bOHzcuS\/TYb8PPAeOWJeKdKQDmcQkWDXXrl4daV4mACIiACIiACMQi0LAChhUSLJEkwI4YAgLHiI5nWSOrbtgIieWWLH0kWJHVHuzaqCQCIiACIiACIlD7BBpSwBB5zv4TeGD8zqNEq7Nkj+WY7EfA8leW7ZHYR4MIfU6q9Xt\/1H7V6Q5FQAREQAREoP0SaEgBU6g6WRrINuksQ2WXTZZHZpcAsm8A+w9UsmdD+21GenIREAEREAERSEugXQgY4l3Yw4QlrOxHwF4UbAKGF8YnPvM7duZXATuO4rnxiR0k2VBMKS4BtrfXqbBxmRbKTZyrz9iXINZpWItzcc7sHl5oD6A0tVO9UhpawDCVRAwMm2hlz8ZhUys2I2P\/Dp\/YvpoNnQpt0sZW4n630+pVhXJmgyx\/YrFoVI+AOFePbX7OYp2GtTin4VxrpTSsgGGDKXZ4JQaGgN3swXrEwLB7KbuYkoiB4QAxtj8vtGmVBEyaZqtOSJzTEEhXitp0GtbinIZzrZXSkAKGnRjZ5prpiOHDh+d2lmRLbwQN52OwJTYnr7L1Oduhs\/skU0iFkgRMmmarTkic0xBIV4radBrW4pyGc62V0pAChr1f8s\/yADzHv3M+Ccmf7cPcKefCcHZHa\/EXEjBpmq06IXFOQyBdKWrTaViLcxrOtVZKQwqY2JAlYGITLZyfOiFxTkMgXSlq02lYi3MazrVWigRMQI1IwARAimCiTigCxIAsxDkAUiQTsY4Eso1sxDkN51orRQImoEYkYAIgRTBRJxQBYkAW4hwAKZKJWEcCKQGTBmSdlSIBE1BhEjABkCKYqLOPADEgC3EOgBTJRKwjgZSASQOyzkqRgAmoMAmYAEgRTNTZR4AYkIU4B0CKZCLWkUBKwKQBWWelSMAEVJgETACkCCbq7CNADMhCnAMgRTIR60ggJWDSgKyzUiRgAipMAiYAUgQTdfYRIAZkIc4BkCKZiHUkkBIwaUDWWSkSMAEVJgETACmCiTr7CBADshDnAEiRTMQ6EkgJmDQg66wUCZiACpOACYAUwUSdfQSIAVmIcwCkSCZiHQmkBEwakHVWigRMQIVJwARAimCizj4CxIAsxDkAUiQTsY4EUgImDcg6K0UCJqDCJGACIEUwUWcfAWJAFuIcACmSiVhHAikBkwZknZUiARNQYRIwAZAimKizjwAxIAtxDoAUyUSsI4GUgEkDss5KkYAJqDAJmABIEUzU2UeAGJCFOAdAimQi1pFASsCkAVlnpUjABFSYBEwApAgm6uwjQAzIQpwDIEUyEetIICVg0oCss1IkYAIqTAImAFIEE3X2ESAGZCHOAZAimYh1JJASMGlA1lkpEjABFSYBEwApgok6+wgQA7IQ5wBIkUzEOhJICZg0IOusFAmYgAqTgAmAFMFEnX0EiAFZiHMApEgmYh0JpARMGpB1VooETECFScAEQIpgos4+AsSALMQ5AFIkE7GOBFICJg3IOitFAiagwiRgAiBFMFFnHwFiQBbiHAApkolYRwIpAZMGZJ2VIgETUGESMAGQIpios48AMSALcQ6AFMlErCOBlIBJA7LOSpGACagwCZgASBFM1NlHgBiQhTgHQIpkItaRQErApAFZZ6VIwARUmARMAKQIJursI0AMyEKcAyBFMhHrSCAlYNKArLNSJGACKkwCJgBSBBN19hEgBmQhzgGQIpmIdSSQEjBpQNZZKRIwARUmARMAKYKJOvsIEAOyEOcASJFMxDoSSAmYNCDrrBQJmIAKk4AJgBTBRJ19BIgBWYhzAKRIJmIdCaQETBqQdVaKBExAhUnABECKYKLOPgLEgCzEOQBSJBOxjgRSAiYNyDorRQImoMIkYAIgRTBRZx8BYkAW4hwAKZKJWEcCKQGTBmSdlSIBE1BhEjABkCKYqLOPADEgC3EOgBTJRKwjgZSASQOyzkqRgAmoMAmYAEgRTNTZR4AYkIU4B0CKZCLWkUBKwKQBWWelSMAEVJgETACkCCbq7CNADMhCnAMgRTIR60ggJWDSgKyzUiRgAipMAiYAUgQTdfYRIAZkIc4BkCKZiHUkkBIwaUDWWSkSMAEVJgETACmCiTr7CBADshDnAEiRTMQ6EkgJmDQg66wUCZiACpOACYAUwUSdfQSIAVmIcwCkSCZiHQmkBEwakHVWigRMQIVJwARAimCizj4CxIAsxDkAUiQTsY4EUgLZMKioAAAgAElEQVQmDcg6K0UCJqDCJGACIEUwUWcfAWJAFuIcACmSiVhHAikBkwZknZUiARNQYRIwAZAimKizjwAxIAtxDoAUyUSsI4GUgEkDss5KkYAJqDAJmABIEUzU2UeAGJCFOAdAimQi1pFASsCkAVlnpbRbAXPFFVfY6NGjbdasWbb99tvbOeecY126dClYfRIwaVq1OntxTkMgXSlq02lYi3MazrVWSrsUMOPGjbORI0fa2LFjbYUVVrDjjz\/eOnfubBdddJEEzCJsoeqE0sAX5zScKUWs07AW5zSca62Udilg9ttvP9tiiy3ssMMOc\/Xx1ltv2XbbbWdPPfWUde\/efaE62nq11ewf11674PNevRb8PWVKy\/\/33\/F5LdsVuvfsMxW79yrbTe3QwVZdddWF2Va53FbrtEHLrTnO9fLbCf1tZ35jU6dOtVXnzSvcV9TKb9H3eLT30GcsZFeoHhPZvdncbL2bmv7X9yYq1\/Ud2bJqvdxaUyAV3k+7FDD9+vWzs846y00d+bTuuuvaddddZxtvvPHCSP0Po0LYulwEREAEREAEFgmB5uZFUmw1C22XAqZv3752+eWXOy+MT3x26aWX2lZbbbUQ75u7drUf7b77wvWw5JJmc+Ys+ByR01oDydottpjZ\/PmF67RTJ7PPP1\/wXYcOZry5FUqdO5t9+umCbxZf3OzLLwvbEdMza9aC7zp2NPvii8L32rWr2SefLPhuiSXM5s4tbNetm9nMmeXb+WfyrPLym\/XRR9aFe2nDLnd\/pdh97WtmH364oI5ae17qBrtp0xbYtcavmF22PnhO8vvoowX5Fas37KZPX2C39NJmn31WuB0su+wCO1K2veS3q2WWMfvvfxfY5bXTWTNnLuBMwuP48cel1Wm2LeW3vCyzYm0zy6JYWy\/2jNmyl1rKbPbswu02a5dt34V\/NQv\/XorZFbt3M5v1yScLWBf73WfzL9aPlGNX7N79d95rEGJbAYvcpVVgMfuzz2wp2opS6wTGj284Ou1SwGy22WZ27rnn2rbbbpur0D59+tgNN9xgG2200UKVrCDeNO1e89jinIZAulLUptOwFuc0nGutlHYpYIiB6d+\/vw0dOtTVBzEwAwcOtIkTJ1pX\/2aaqSkJmDTNVp2QOKchkK4Utek0rMU5DedaK6VdCphbbrnFLrzwQhszZoyttNJKNmzYMGtqanJTSIWSBEyaZqtOSJzTEEhXitp0GtbinIZzrZXSLgUMlTBq1Ci79tprbebMmS7uZcSIEdaNmIwCSQImTbNVJyTOaQikK0VtOg1rcU7DudZKabcCppSKkIAphVb5tuqEymdXypXiXAqtymzFujJ+oVeLcyipxrKTgAmoTwmYAEgRTNQJRYAYkIU4B0CKZCLWkUC2kY04p+Fca6VIwATUiARMAKQIJuqEIkAMyEKcAyBFMhHrSCAlYNKArLNSJGACKkwCJgBSBBN19hEgBmQhzgGQIpmIdSSQEjBpQNZZKRIwARUmARMAKYKJOvsIEAOyEOcASJFMxDoSSAmYNCDrrBQJmIAKk4AJgBTBRJ19BIgBWYhzAKRIJmIdCaQETBqQdVaKBExAhUnABECKYKLOPgLEgCzEOQBSJBOxjgRSAiYNyDorRQImoMIkYAIgRTBRZx8BYkAW4hwAKZKJWEcCKQGTBmSdlSIBE1BhEjABkCKYqLOPADEgC3EOgBTJRKwjgZSASQOyzkqRgAmoMAmYAEgRTNTZR4AYkIU4B0CKZCLWkUBKwKQBWWelSMAEVJgETACkCCbq7CNADMhCnAMgRTIR60ggJWDSgKyzUiRgAipMAiYAUgQTdfYRIAZkIc4BkCKZiHUkkBIwaUDWWSkSMAEVJgETACmCiTr7CBADshDnAEiRTMQ6EkgJmDQg66wUCZiACpOACYAUwUSdfQSIAVmIcwCkSCZiHQmkBEwakHVWigRMQIVJwARAimCizj4CxIAsxDkAUiQTsY4EUgImDcg6K0UCJqDCJGACIEUwUWcfAWJAFuIcACmSiVhHAikBkwZknZUiARNQYRIwAZAimKizjwAxIAtxDoAUyUSsI4GUgEkDss5KkYAJqDAJmABIEUzU2UeAGJCFOAdAimQi1pFASsCkAVlnpUjABFSYBEwApAgm6uwjQAzIQpwDIEUyEetIICVg0oCss1IkYAIqTAImAFIEE3X2ESAGZCHOAZAimYh1JJASMGlA1lkpEjABFSYBEwApgok6+wgQA7IQ5wBIkUzEOhJICZg0IOusFAmYgAqTgAmAFMFEnX0EiAFZiHMApEgmYh0JpARMGpB1VooETECFScAEQIpgos4+AsSALMQ5AFIkE7GOBFICJg3IOitFAiagwiRgAiBFMFFnHwFiQBbiHAApkolYRwIpAZMGZJ2VIgETUGESMAGQIpios48AMSALcQ6AFMlErCOBlIBJA7LOSpGACagwCZgASBFM1NlHgBiQhTgHQIpkItaRQErApAFZZ6VIwARUmARMAKQIJursI0AMyEKcAyBFMhHrSCAlYNKArLNSJGACKkwCJgBSBBN19hEgBmQhzmZ77723Pf744wG0ZCICcQlsttlmdv3118fNtJ3mJgETUPESMAGQIphoYI0AMSALcTbTbzqgocikKgTU9uJhlYAJYKkGFwApgokG1ggQA7IQZwmYgGZSEyZnn322bbLJJrbDDju4+9ljjz3s4osvtpVXXrmk+5s7d64tscQSJV1TLWONJ\/HISsAEsFSDC4AUwUQDawSIAVmIswRMQDMpavLnP\/\/ZCQnSu+++20JQTJs2zbp16+YEw4EHHmgHHHCAEyAffvih+3yfffaxxx57zF577TV3Pfbdu3e3jh072lprrWXXXnut+3zGjBm222672SmnnGIXXHCB++z111+3NddcM3dvxx57bE7cPPTQQ4bgyU9vvfWWDRo0yM4880xbbLHFKn30iq\/XeFIxwlwGEjABLNXgAiBFMNHAGgFiQBbiLAET0EyCTTbeeGN7+umnc\/YIluOPP9422GCDFnn86le\/sl122cX69euX+xyRsv3229s\/\/vEP69y5cwv73\/zmN9ahQwc78sgj3ef33Xef+zNixIige5s1a5YTPiussIIddthhQdekMNJ4Eo+yBEwASzW4AEgRTDSwRoAYkIU4N46AefLJJ+3SSy+15557zubPn289e\/a0Qw45xHbdddeiLeFvf\/ubvfzyy3b00Ue3aodwwMPRu3fvonlVImAuueQSe+mll+yyyy5rUQZeHbw2PNs222zjvkPQLLfccvaTn\/xkofsZN26c8\/hsueWWzptD2nHHHe3EE09017\/yyis2adIkYypp9913D\/iVVM9E40k8thIwASzV4AIgRTDRwBoBYkAW4twYAoYpE0TGGWecYd\/5zndsySWXtCeeeMJ9xlSKjxvJbxIzZ860vfbay2699VZ3TWuJdnLCCSfYTTfdVBUB88UXX9hWW21l3\/jGN2ynnXZy8S2kL7\/80vbff3\/717\/+ZRdddJGb+iG9\/fbbttJKKzmvjE9LL7203XXXXXb55Zfnnn\/nnXd2nh7yfvjhh933CLa+ffvanDlzFrk3RuNJQCcVaCIBEwBKDS4AUgQTDawRIAZkIc6NIWAY9IknyfdIjB8\/3ubNm+emZgolPBl8f9xxx7mpn+HDh+fMpk6dan369LEbb7zRfUbeBx98cM4L4g333HNPF9NC+ve\/\/22rr756Lo\/333\/fevTo4QQF0zc33HCDDR482J555hkXA8MUEOJi9OjR9rvf\/c4GDhzovCMnnXSS+\/y9995zgoR7RIThQZk+fboTOPfff3\/BZ\/IChuXx3PPYsWNtu+22s9NPP90tWb7qqqvsuuuuk4AJ6B\/qyUQCJqC2JGACIEUw0cAaAWJAFuK8sICZMiUA3CI06dWrZeHEd2y44YaGF2aVVVZp9c7++te\/2h\/+8Afr0qWLnX\/++U5YfO9737PTTjutRSwKGTz77LPOe4O9nza64oorbMqUKa3GneAVIb7kjjvuyN1DSAwMU14II0TU888\/bz\/96U9dsO\/NN99syy+\/vMuLmBkvYBAf5513XotgYTwxf\/zjH90zeQFD0PAHH3zgPFB4XciLaSSeGcEkD8wibMRVKLphBcw777zj1DcuVRo6Kp5odho7iR8mDZqOgDeVc845x\/3ICyUJmCq0vAJZamAV5zQEFhYwTU2pSi6vnObmtgUM00LEtRALg0fkzjvvtB\/\/+Mf2pz\/9yR599FF74YUXnNhYd911XTAs8TI+sbpn3333dTEnxLT4hAAaM2aMseqoULrnnnvs7rvvzq1IwqYtAYOHB2HC9BBeFoQGggORwsokxFW+gMkv+\/DDD7f11lvPjjjiCPeVFzAIllGjRrm+nvidv\/zlL04oMRVGHIwETHntr1avalgB86Mf\/cjWXnttO\/nkk13gFo2Z+VKCxQj4GjlypHMz4uIkYp4IeH5QEjCLrqlKwKRhL84LC5g24lTTVEyRUt58c+EviXthagbhkU0PPvhgThAQH8M0DS9q9IEM7giYCRMm2IorruguY8qG6Rk8HngqsglxcfXVV+emlPLvAgGx7bbb2g9\/+MPcV8UEDKIDAcM+LniOuFcvYJqbm2327NnWqVOnggLmkUcecX06\/TcvpjxL01fKEwHz6aefuiXdCKCuXbva1ltv7WJgmNL6xS9+4f5NfM2iXpGkF+J4P6eGFDCo7KOOOsp5YPyPlDcOfqBsH77ffvvZFltskWvI7BPAfOlTTz2Vi2DPIlaDi9fgiuWkgVWc0xBojBgYBn+mfLzwIOYETwrxHkz7XHjhhc4zwdQRXhkCY5lyIXaGQZ4t7T\/++GPnpSF2BOGRnxAv5ImHOj8RP\/Ozn\/3MxaVkvdf5AgYvC\/dE2auuuqq7Z\/pfUlbA5OefnULiOzw0v\/3tb11ZiBhEik\/ZKST6+HPPPdd5kvDE+KQppFS\/rnTlNKSAKYSPeVACwa688ko393vWWWe1CHLjrYQfSNZ96vORgEnTICVgxDkNgcYQMLAibgWvMgGyeFmYFkKgEFPCixwrePI9MCxd5jviTwjoZfBHWCA0EDpeWPA3omfo0KEuuDabJk+ebEOGDHHiiVU\/2VTIA4OAwbuTnbYqVcBgjxfmmGOOsV\/\/+tfWv3\/\/ggKmtTYkAZPq15WunHYhYIh3YXdH5nFxW7Kcjh+tfwsAN58x\/5v\/Q+U7BEw28aMu9LaSrtoasyQCAvM7uMZ80kX7VOJszuP6xhtvLNqKSFB6NgbGB+myjJopdmJkyl1GjUeGfpG4m\/xEzAn7y6y\/\/vptPiFByHhwsiuh\/EWUMWDAgBb9NN+x9w33jofdJ1Y6IZKKJZ6bIN9amEJqbTVVm8DKNFhmmWVy8Z9lZlGTlzWEgOGNgzcCEudm0JhJvGUQA8N+AjRuH62P6xQXI3O3PjEvy3UbbbTRQhUlD0yatisPjDinIdA4Hpi2eBGES79GrAjeFr+IgcBb+ky8Ga0l4lt+\/vOft9i6v63yavl7YmHwMJV6jlLsZ9J4Eo9oQwgYBj4CdUkE6uIOZW4XFycxMATsZrepJgYG9yOuURIxMOxFMHHixBbzqh6zGly8BlcsJwkYcU5DoP0ImFQ8VU44AY0n4azasmwIAZP\/kESz49pk0yRckz5SnYO8EDS33HKLC3BjeSA7Ow4bNszZMIVUKKnBtdWM4nwvAROHY1u5iLMETFttRN9Xj4DGk3hsG1LAsMSu0NzsUkst5c7dILEEj7gY5kWJe2F3SASPBEy8xlVqThpYSyVWnr04S8CU13J0VQwCEjAxKC7IoyEFTDw8C3JSg4tNtHB+GljFOQ0B\/aZDOLOjLfvMFEt+OTY2TOOzIol9XbzXm89Z1sy5S22dIo3nnL1trrnmmqK7C7d2P5TPgY4+4WknVAAvPCuu2Epj0KBB7iWWTe0INl4USeNJPOoSMAEs1eACIEUwkYCJADEgC3GWgAloJm2asEke3mu\/GII9XQj8za7aZMEE+7WwQ26+V5zjBNhcL5s4X4mVUVlvOIs0ECAkloCzAimbWP5NHCMLM7Ib6rGJHeKJ1WasukLQcD3BzCzayN+0r80HjmSg8SQSSHlgwkCqwYVxqtRKA2ulBMOuF2cJmLCWYm4HXoRAoYQng31gvIBhgzqWKHO6dDaxaAKB4RdSbL755m6Duccee8ztvcU+NUzlszKKZd8kjkTAA8QRMMWSP6eJjfiythwlw\/EBn3\/+uSubhR0dO3Z0e38havDAsOCD9N3vftettkqVNJ7EIy0PTABLNbgASBFMNLBGgBiQhTjXv4BhLxSmZFgKzeKETTfd1B1auOaaa+ZaAJvYMSXDxp2FEqIB0UE8YLH9YFprUlkPDG0KjwYHTGbTBRdc4FZ4cr\/LLrus23\/H74GSFTAIjl\/+8pfuIEkSh0MikBAhPrEvDEKFfLp37+6EyUEHHeS2yOBImI8++sjFNeJtwcNCuv32291ScUQVQol9pmDFNBfi5tRTT3VHyqRMGk\/i0ZaACWCpBhcAKYKJBtYIEAOyEOf6FjDsb8V+VxyASMzJF1984QZxvCVM4yBoWMjAwMwuvcUS10ydOtWdB5ef2oqB4XuOKcADwzQNAoFDJD\/77DM76aST3I7AiCfiXzgCgN3OOeLl+uuvd0UhYH7\/+9+7HYR5JvLgrCQSYgSBgeDYc8893dQQm4fi4WGqCGHEuUacX4cAIxYH0cIZT9iwnxcs2EKDvW+YvmLl6U9+8hNDeCGUiLn5xz\/+4bxIKZPGk3i0JWACWKrBBUCKYKKBNQLEgCzEub4FzCeffOKOPCEWZJ111nE1ThwIQmGXXXZx3hQGbOJRiFEhKJZTmX3Ca3PIIYe4uBJEAlMv9957b5s7tXIm0s033+y2nSiUiDFh\/y0ECffHcS1Mz6y11lpGm\/vmN7\/ppqP8+UTswss1iJAQDwwCBi\/N4osv7s5m4oBGngGvzL\/\/\/W+XL3t7IWDYLgOxwllObOTHNBYBvuwQzIGQCCPyYRrJT1sF\/HSimGg8iYLRZSIBE8BSDS4AUgQTDawRIAZkIc4FBMyQIQHkFqHJmDEtCsdzwnlAxHYgFvB0IEQ6dOhgM2bMcAc1MgXDIJ1N7H3FtMqf\/vSn3KnPrDTaddddWwTA+muIUWF6htgRpnXwauD9IJEXHhc8MBxTgJhg9REBuJRLHAtTWEx1IT4QWxwN4GNmiHkhP74neDc\/oPf73\/9+C3HhBQxemhdeeMFN\/0ybNs0JprvuusvtsOsFDB4YvDp8T9AunhgS5yEhqphu49Rr7tmLwFS1q\/EkHmkJmACWanABkCKYaGCNADEgC3EuIGCamgLILUKT5uaFCmcKhYGc+BLiSpi6wdtAgCqxLfnBtwgQvB033XSTLbfccrn8EALElPiVPv4LBABTVJwijUcDrwdth+kqliwjavDsfP3rX3eX8D3xJU899ZTzlBATg0ghsWSZKSI8PV5UcTYd+TGF9d\/\/\/neh52OKB++KPywSAUM8ywMPPOAEE8KLZ+HgSpZKc1+cj+SnkIiPgQ0iC7EyadIkF8DL\/\/vVSniAUieNJ\/GIS8AEsFSDC4AUwUQDawSIAVmIcwEBkziQM6CaWppk9mNBmBA\/kh+3grcDLwZiAzHCQO8T1\/AZAqd3794t8iYmhKXO+StxiKNhmsdPNW299dZOrCAWmI4ZP358i7OUKBuhQnwJgbRM6XAWEwlvEQIie+AiMTMIJw7HxaP0\/PPPu7geEse6UC4CyZ\/fhIDB04Rw+9nPfmYvvvii2z2dGBhOqaYs4nK8gGFaimsI2EVI8Sx+LxpiZhBhxNOkThpP4hGXgAlgqQYXACmCiQbWCBADshDn+o6Bof7wSiA8WNVDHAeeBgZuPAp4V7bccksnAhAaDPScC8dAX+iwWkQHsTNMIxVKxMkgJhAS2BJbcvHFF7upq2wi7gavDX\/Y9RxP0A9+8AO3UgiRQnAx1\/l7QITwnc\/n8ssvd8HHBCgT38PUT3Z6JzuFxFQQNtxLdu+X7BQS00c889NPP+08PF4MwYX4GTw8iCgfOBzw04liovEkCkaXySIRMMxZMleaPWAx3iPFz0kNLj7TQjlqYBXnNATqW8DACI8DAzxxLqxCIjCV6RWmdUgM7IceeqgL4kVYEI+CsMF7QfrWt77lxARig91q2aU2\/ygVhAurehjkCcTFi8JeKpRJEHDfvn2d6EGQsFGdTwQUP\/roozZu3DgXPMt0Fp4PBAxeGgJr6VNZXk3symuvvebEFh4YpngQXeRHrAzihvyZ9skKmNbaSb6AYfoLVgTtsgqKcvG6IGwQNdwPm9ttsMEGqZqednaPSDqJgCEyHRch87Tsvpg9OZp9A3iLQEV7V2HE54uSlQRMFIxtZiIB0yaiKAbiXP8Cpq2GgCeDJcoM1MUSU0rsp5K\/sohAYLwyiAfiULbYYosW2dCHE8Py97\/\/3YkMRIpPLKH+z3\/+4\/Kk399hhx1ye9H4\/V1Y+UOA8I477ugCbNmnxQciI2D++c9\/uqkfxguEDKuHCLzluYolymUKiz778MMPd8+w++67u+BmNr3D84IQW2211Vw2TMURSJzdb6YttpV+r\/GkUoL\/u76qAgZ1j7olWp03BKLk+ftrX\/uazZ8\/3631R83TKHE38tbAPgGo7VpKanBpakMDqzinIdD4AgaOBLQSF5MvPjxjNrJjqfWVV15ZsM9liiV7plGquqmknFdffdVWWmkl69KlSyXZVPVajSfx8FZVwLCRES5MosVZ4lYsoY7ZRRF33z333BPvCSPkpAYXAWJAFhIwAZAimIhz+xAwEZqKsqgCAY0n8aBWVcAQBb7KKquUdLflXFNSAWUYq8GVAa2MSzSwlgGtjEvEWQKmjGajSyIR0HgSCeSiCuKNd\/tpclKDS8NZA6s4pyEgAVMNzsQ6snSamBZ\/thJ7v7A8+7bbbluoSPakIYamtUS4AYHGPrFKiWXX++yzT4tLWOFEYG5bL8us0GJFkl+qXQoDgp8JQuYPiYUorIL6v\/\/7PxfkzL43xOyw2ouVWoROtHa+lMaTUsgXt62qB4ai3333XRfhThAVOyASCU\/lEolOHMxuu+3mzqIgcKtWkxpcmpqRgBHnNAQkYNriTPAsy7JJ06dPdzEyhWITOR4gf28XDlvcd9993bXFBEyxe+C4BGIis8Jn0KBB7tTqbHwLYokVUexJk59Y+s3qJp9YrfXee+\/lAnj95wQys+qJvDl3KT8Rq8kqJpZuZxeasCqKGCJCH1hRy\/2y7JsFKYib1pLGk7ZaX\/j3VRUwxLOwzh\/1SgUjYvg3qprKnz17tjvYi0ZAlHmtJjW4NDUjASPOaQhIwJTCGY8Hq3cK7SHj86FPZ5ddXkrxUiB2OK6AkIDWPDCIFFYsFUpebHgBQ3AuRxgQcJxN7H\/DWNKrV6\/cx2xux2Z9LB9nd2D6b5Zrs3Tce2kQXXhL8o9ayOZNEDMrmljezcneWQGH1wXRhA3PSHm8hLMrMKufEDw+kUd2mbnGk1JaX3HbqgoYFDMVh6uQzZZo5JxKypI7v+sjS9popJwKWqtJDS5NzUjAiHMaAvUtYJiiYFdZEsuZGYT91Aaf0+cyncHeL8VeDOmH2Qclf2fe\/DrgeID8fWLY64UX0kIpO1WDB4YXVJZMl7JNRr4HBm8KYoVdhn3C0\/HMM8+4sYXdgfGEcF8IE1JWwLBrMUu32SSPxNYd7Nzrp3lYBcu+Neyd4\/fDYWYAbxJlc5AliTJY3UVinxx2HPbPiHBh6ToLUpimgi3TWwinbNJ4Eu9XXlUBw8ZBHGHOaiQS2zyz1p9GwOokEttd02hoQLWa1ODS1IwEjDinIVDfAibLiP6UDe342yemQhi8OZ26tTgMbPm9sf8J5wm1lohrQcBkPQrYEhqAQMpPHMrIiytlf\/rpp84jw268iJcrrrjC7ffiU7EYGDwbeDS8BwaPCffJ\/TKFxLQOggYRgveF704++WS3sR736wUMYxBHGCBAOAaBIxNITJGxBxkb87HZHWcrEWPDIZiczYQgYlM8XrQpj+092IcG7w6ne\/NsfMfqWnb8ZWsQhBRCj2flfhBYY8eOXWiDQI0n8X7lVRUwVBR7wLD\/C4nAJxoNmxt5JcwmQihmVGytJjW4NDUjASPOaQg0toBB0OCJwHNA0CqDsU8Mzn369HF9MIm4DTw6vo\/O549IYNoGz0JWIHGUgfcC+c8RHcTNIFrwcPA9XggGev72m8ix82+hhDeeE6i5v0KJ74hXIS9eirk3BBHxJggFRBsvx16k4B1h2gkPU1seGC9g2PiPDfiIzSRPjmRgaoypKDwu7CiMgME7g6cLO44oIH\/EFGIIOzb5Q0z5s6Cyz6PxJN6vXAImgKUaXACkCCYSMBEgBmQhzgsLmCYrfhp1sy18GrRH3da12BW7nu\/byqO16wt5YDgFmsMKs14ZymCQZtqImAw\/bYRXhCkPf8hhfvNBEBD4mt2pFm85q428CPLXcDYTU0UM7EzD3Hzzze5aHwND+cSIMFWz6qqruiMQmHbxuwUjehAHeHsQAUx\/+e8QDHhj1l9\/ffcizPfEVjJVQ754VxA3eH58YpdfAm\/xkpBPvuDiOzwwng0eGOzIm78RbsTMsBEg00Z4nLyAYZNWvFs8M1NTiK7JkycbHigEFUcTEAOKxyY\/aTwJ6KQCTaouYFDBVBiJ+VriYNia2rsk8bzQ0OWBCayxBjbTwJqmcsW5sQUM0yZMZWSnfZhCYWUQA2r2EEY8BHjJ8RgUSngjNt98czct5BMnWzNFQjyjTywjZvBmFQ5xKngt2MCUM5SyQbxMPflNTREcHDR51llnuRVA3\/72t91hlIgUhMy\/\/\/1vO\/XUU10RiBY89RwLgBefIwCYskE8kIffxf3AAw\/M3RNhCkwrIaoKJRaWsCKWGCLvgUHgIf6Y\/iEuE4HF9BVMEVhewHC6NYG8hHBfdFoAACAASURBVEUwvr3\/\/vvuGs6jgj3ByXi+fPhEtnwJmHj9XFUFTLGo9fxHoAHValKDS1MzGljFOQ2Bxp5CYrDlEEYf8IoXZI899nDTHnhnsolgX7a0yPemeBviQHwAq\/8MsYFwYNrEJ0IBEDqXXXaZ88BwhAEeCMRGa6uQEAt4cvCe4AlipRBBswgqhBOiiG03SHiJmMZB0OANId8111zTCQQ8Jow1TOUgfkh4Z5gW40wlvkfgsFcMAospLu6JU6iJvyR5AcM+M0ypYceKWZ6R\/WiYGuNlnFWziCsS4ghPC8KF5de8nHPGEsvOEWOUzbEG+UnjSbxfeVUFTLzbXLQ5qcGl4S8BI85pCDS2gGFQx4tAQCpBuIgNvN4sN85PiBcECuIgNBGrwtQKMSf5ifgbBEwp+8AwtYNQwCND3AkBskwNZQOQ8bgQSEtfjED47LPPXHAtogRvEGfqefHDPTGdQx5MgZEQcYgxPEkcVQMbYmn8WU\/ZGBi8PZyijXgi0NjvO5OdQiJPFqAw1cRUGHE9eF9I7EuDaFpmmWVaeKk8K40noS2tbTsJmLYZ6fjzAEYxTCRgYlBsOw9xbmwBw7QIMRoE8SIoiO0g5oTAXrwPJDwmJJb5MrD7VaFttx5zQgKPQ6FrShEwtEM8GsTEcB2xJtwfU0p4hrg3PsOj5Pdr4XuChRFIiBcED\/kQz8KUFV4ZxBvxKzwvQbasOGJDO56Za1k1xKaqeElYsYTXJCtgWmOQFTDk8+ijjzqvC4HE5M+U0nnnned2C+YZ8OQwjYfnK7vfjARMSCsLs6mqgGHOMDTdf\/\/9oabJ7dTg0iDXwCrOaQg0joApxIspHqY7CCYtdxk1HplCK2gKlcf0k98gjkEcD4z3zvCbZvkx3pVswouDB4MpLTxDfuWQt8GDwkZ47PPClA3eFRKrqlhZhPBgJRSro9hkjudk7ximhphi4hq8TsT8IOaYcmJaiyXXxOogmljFhPeJqSdiWNivplgiHwQSeTPl9fbbbxsxN4wPCBqmlPieWB08O9gjpvDIZPfa0XgS71deVQGDKvaJYCfmNtkTgI19WOP\/wgsvuOAxlKrf2C7eo8XLSQ0uHstiOUnAiHMaAo0tYGBIkCkiwE9rFOLKgE2\/y6CbOjHg++mb1GWXWx5TW+y6S3xMJUnjSSX0Wl5bVQGTLQpFvssuu+TmRv13qF+ixLPL3+I9Xpyc1ODicGwrFwmYtgjF+V6cG1\/AxGkpyqUaBDSexKOaTMAwj8lbQf621Syfxo3I20KtJjW4NDWjgVWc0xCQgEnFWeUsTEDjSbxWkUzAEBnPsrbshkgEZDFfyWFZ8sDEq9R6zUkCJk3NibO5+IjHH388DXCVIgIZAqyAInhZqXICyQQMmx+xrh4PDEewMwdKdDixMewMmX\/gVeWPFi8HKeZ4LIvlpIFVnNMQSFeK2nQa1uKchnOtlZJMwPDgftkcR6NzeirLz3gT8jsz1hocfz8SMGlqRp2QOKchkK4Utek0rMU5DedaK6WqAoYpIn\/Me+iDl3NNaN7l2knAlEuutOvUCZXGq1xrcS6XXOnXiXXpzMq5QpzLoVb\/11RVwLDrIZseDRw4sE1STCmxFwwbA91xxx1t2qc0kIBJQ1udkDinIZCuFLXpNKzFOQ3nWiulqgKGU0s5z4LtrAni3Xrrrd3hXmyxjKeFLaEnTZpkTzzxhNtYiM2KOCyLMypqKUnApKkNdULinIZAulLUptOwFuc0nGutlKoKGB4Wz8r48ePdls+cZ8EhW\/5zNjLiOHO2i+bsiW233bbW+Lj7kYBJUy3qhMQ5DYF0pahNp2Etzmk411opVRcw2QdmJ0MO3eIsig4dOrgtp9kfplOnTlXlwlkUHM3+xz\/+MVcOW0GzUzAHd22\/\/fbuMDN\/aFf+zUjAVLV6cpmrExLnNATSlaI2nYa1OKfhXGulJBUwi+LhORSMI9Px8ngBw3TVyJEjnVdohRVWcCehdu7c2Z08WihJwKSpOXVC4pyGQLpS1KbTsBbnNJxrrZSGFjAc0LX77rsbh0pyaJcXMExXIWgOO+wwVx9vvfWWs3nqqafcSaLywCyaZqpOKA13cU7DmVLEOg1rcU7DudZKaVgBQ+wNQoXzl5i64uRTL2D69evnjjtn6sgnprI4k2njjTeWgFlErVSdUBrw4pyGswSMOKcj0D5LalgBQ4zLs88+a1deeaWLdckKmL59+9rll1\/uvDA+8RlHr2+11VYFBUz2w\/33398dAa8UlwDH0\/fs2TNupsptIQLinK5RiHUa1uJcnDMrf3v06JGmMhKWUjMCBi8JcSjlJJZo+8MgN9lkEzv55JNd3Av7yVBp+QKGsyjOPffcFque+vTpYzfccINttNFG8sCUUwkRrpFnIALEgCzEOQBSJBOxjgSyjWzEOQ3nWiulqgLm1FNPtTPOOCPomddff3178cUXg2zzjWi8c+fOdR8vvfTSdvXVV7vDspZaain32ezZs92+M3z3zDPP2ODBg61\/\/\/42dOhQ9z0xMGy2N3HiROvatasETFm1UPlF6oQqZxiSgziHUIpjI9ZxOLaVizi3Ragxv6+qgFl77bWNc4+OO+44d2hjsfTII4+4QNsYacaMGS7uxacbb7zR7UHDKiOmKG655Ra78MILbcyYMbbSSivZsGHDjD1pmEIqlLQKKUattJ2HOqG2GcWwEOcYFMPyEOswTpVaiXOlBOvz+iQCBu\/GKaec4gh5r8yZZ57pPvN\/syLITwPFRpk\/hUT+o0aNsmuvvdZmzpzp4l7YAbhbt24SMLHhl5CfOqESYFVgKs4VwCvxUrEuEViZ5uJcJrg6vyyZgGE3XtKOO+5o99xzj\/3gBz+w2267Lfd3JVNI1a4DeWCqTXhB\/uqExDkNgXSlqE2nYS3OaTjXWilVETBHHHGEffjhh\/bkk08aq4GIg3n44YclYGqt9mvsftQJpakQcU7DWaJcnNMRaJ8lVUXAPPfcc\/bGG2\/Yscce6w5wfPTRR42DHeWBaZ+NLPSpNbCGkqrMTpwr41fK1WJdCq3ybcW5fHb1fGVVBIwH4oN4iYHRFFI9N5M0965OSJzTEEhXitp0GtbinIZzrZWSTMD4TePuvPNO+\/73v2\/33nuvffe73839fdNNN7kVS7WYFAOTplbUCYlzGgLpSlGbTsNanNNwrrVSkgiY22+\/3T755JOiz85qpEmTJtUaH3c\/EjBpqkWdkDinIZCuFLXpNKzFOQ3nWiulKgLmqquuso8++sg22GAD23nnnYOe+Yc\/\/KHbn6UWkwRMmlpRJyTOaQikK0VtOg1rcU7DudZKqYqAefrpp90+KxMmTLBBgwbZlltuGfTcbDJXaCv\/oIuraCQBU0W4mazVCYlzGgLpSlGbTsNanNNwrrVSqiJg\/EO+8847ds0117g\/X3zxhW244Ya2yiqrtMrg29\/+tg0ZMqTWGGkKKVGNqBNKA1qc03CmFLFOw1qc03CutVKqKmD8w3LWkN8L5phjjnEHLdZTkgcmTW2pExLnNATSlaI2nYa1OKfhXGulJBEw\/qEfeOABW3755V1sTD0lCZg0taVOSJzTEEhXitp0GtbinIZzrZVSNQHDtNEBBxxQ8HnZoXfOnDkLfff1r3\/d1llnnVpjpCmkRDWiTigNaHFOw1lTSOKcjkD7LKlqAoZ4l+eff95NHXGsAGnZZZe1008\/3dgThoBd0sSJE+1b3\/qWsXsvRxAceeSRNVcT8sCkqRINrOKchkC6UtSm07AW5zSca62UqgsYduEdNmyYe+6zzz7bnYmEgHnkkUfcZ5tssok7Myn7Wa1BkoBJUyPqhMQ5DYF0pahNp2Etzmk411opSQSMP0Zgq622koCptRZQQ\/ejTihNZYhzGs6aQhLndATaZ0mLTMDMnDnTEf\/000+tc+fO7m861lpM8sCkqRUNrOKchkC6UtSm07AW5zSca62URSZgNIVUa01h0d+POqE0dSDOaTjLAyPO6Qi0z5KSCJgdd9zR0b3jjjs0hdQ+21nQU2tgDcJUsZE4V4wwOAOxDkZVkaE4V4Svbi+uuoD5y1\/+Yn66qGvXrrbHHnu4gF2f3n\/\/fVtxxRXtP\/\/5jzt24IILLqg5mJpCSlMl6oTEOQ2BdKWoTadhLc5pONdaKVUTMJMnT7Y+ffrYfffdZ\/Pnz2\/x3NOnT7d58+bZ+uuvbyuvvHLuuyWXXNK6detWa4y0D0yiGlEnlAa0OKfhrCkkcU5HoH2WVDUB43EiYoYOHboQ3ebmZhs3bpybUqr1JA9MmhrSwCrOaQikK0VtOg1rcU7DudZKqYqAefDBB+2zzz5zz3r00UfbxRdfXPC52f+FfWB86t69e4vppVqBJQGTpibUCYlzGgLpSlGbTsNanNNwrrVSqiJgRowYYdOmTXPPeuutt7rYlpDE7rwc9lhrSQImTY2oExLnNATSlaI2nYa1OKfhXGulVEXAZB+SKaTf\/e53rT53v379ajLuJXvDEjBpmq06IXFOQyBdKWrTaViLcxrOtVZK1QXMeeed5wJ2W0sHHnigW4VUy0kCJk3tqBMS5zQE0pWiNp2GtTin4VxrpVRdwNTaA5dzPxIw5VAr\/Rp1QqUzK+cKcS6HWnnXiHV53Eq9SpxLJdYY9hIwAfUoARMAKYKJOqEIEAOyEOcASJFMxDoSyDayEec0nGutFAmYgBqRgAmAFMFEnVAEiAFZiHMApEgmtcp6yhSzXr0WPOQUm+L+7mW9jM\/d\/\/vvWrPLXOPz4HrShAktr8\/l14uSprhy8q9x9zBhQOFyuYcBC24se6\/ZPKZOnWrzXu\/f4voWz9FrQuvlfvXsNmXB8xd89l4TcowK3bt7pq+uL8ivwLN75rn8pgxY8Ix59+DsMs9fDj8bMMEG2IL8GylJwATUpgRMAKQIJrXa2Ud4tJrKohY4Zwey1uBkB1M\/cBUbxAoNQH4waV79TZs8\/3X7YtYSNrn5FZs2Z6Yt0eMz47+PmqbZYuO3s253\/9j++98Fd9Ohg9nrr2fubPzAonU4Zew2ZqcPX8jmyy+\/tMW3f9hszJAW37ln6\/WVWuCbCQPMBo5vvYzmpuJtaEovs96tHIZLOW\/2brsNNjWXX\/7YwWZDxpR\/faXPX+z6ARPM2qg\/2k6r\/HiqtvhjU6v8vnr+ZitSv223jpq0kIAJqBYJmABIEUxqYWCN8Bg1lYUXAQyWXgR0mNrBVl111YXehnu+NsA4JP7FF83eesvs88\/N3nvPrLnZrNN5w+3lrz9Y8Nn823P+AJYre1EPICEDeDUH4MFjFxIwC4EscQD3b+Et8mlNwGBURMDk8iomoNoSABO2saYzhru24j0Y2XubMmZgzgPSYd4SNq\/DXPd10\/zFrHmx+WZTVrdO5w93bW6xxcw6dzZbZhmz6dPNmprMmn9ztH3aYaa7Zsl5S9ucDp9Z0\/wOtuS8zrZ8x+427T+fWN9xF9nHH5tNnWq2zjpmbAD\/\/vtmPdb6yGbvN8qaO861qfaO9Zq7ls3rOMfebppqPeetbk3WZO\/Zu7bb\/ZcY60kmTTInZAcMMHvpJbMPPjBbbOjvbc68L+w9e9\/W6LC6rdq8mr0y\/zX7osNs26xpU3tx3su2+DUH2RprLCiXP4svvoAFQnjSdpfamramTW+ebtObpts6to7NtJn2RfOX9rWmZd2\/p0xo3QP00YCbrbt1t1VsFXvZXrZpNs3WsrXsE\/vE3rf33fVNE\/4nsimX3x9\/Vu813z4ZcKftarvWVN8U42YkYAIoSsAEQIpgUu8CpoXH4KvOI+eO\/sqFzVvQq\/aqTZz3gs3o8JG9b+85cu+9bzZn74XfYHMi4M3eOQHSKupqvgEyAGY9Bvk3Ucob8JReBQc5\/wZcaADMvUHnX\/vl4rbY7KVt5U497MPDh9v\/\/Z\/ZEkuYdelittxymZscsGAKIevJ4Vv\/WXYKwHf+Le7jq2fP5pGfX6EpmIcemmr9+6\/qbqTFVE1mKiXCT6fdZ1HvfUe7r8AyAUjABICTgAmAFMGkWp2Qn1Ne4IP4aq59wIL\/7zS\/s\/1t\/v328csr2qxLBxsr\/l97reXDOBHRloAo5oIO8UC419HSXfgt3sQzb+ALiYAxQxaIhrkd7ZO5023FTzay9aZvaYstO8MWX+5jW7tjrwXTGJmUHcinZGIIWtj4uIlemdiGvDiKCE2jbrOoVpuuWyBVunFxrhLYGs9WAiaggiRgAiBVaIKYaH6z2Xr3\/t9cvXeBujflAVPsrvl321tNU+yDpg9avEl7UTJlYMI5+HxPAP8e0jKGAQHQo4fZjBlmG55wt6283BK2\/LIdnIs8\/+3fiayvghvzRYR\/ey84bVAGd3X2ZUAr8xKxLhNciZeJc4nAGsRcAiagIiVgAiAVMMlNf+A+7zXBXpo32V60F+3lDi86a+9Sz7n1y\/BA5IrNi2HIeiDc\/zOH\/5Xo6DF9TWt+dS3bdZMVrekrMeEExZQBLVYAFBIS5ZGoravU2aerD7FOw1qc03CutVIkYAJqRALmf5By0zBmtrKtbOc2j7ApTW8uWEI4cExuGWYLrCFBlFN6Wc\/+r9riRL59tYzTi5ANNjDb4MgJ1rFjy5iFQjEIAdXZ7k3U2adrAmKdhrU4p+Fca6VIwATUSHsTMH7qhv0cHnzwq2j28W0EccLxKw+KFx7+7803N+t61Gjbe8VtF3hd\/PQLcSOZPSfUCQU0xggm4hwBYmAWYh0IqkIzca4QYJ1eLgETUHGNKGC8J+W2eXfYxA7P2ASbUDyIlCDQwV9t5jRlgK33yaa25wbrLaCXESQFV5AEMMZEnVAgqArNxLlCgCVcLtYlwKrAVJwrgFfHlzasgPniiy9s5MiRduutt9rcuXNts802sxEjRtiyyy7rquuKK66w0aNH26xZs2z77be3c845x7qw9rJAqncBgyfl0c0vtL8teccCodJK6tV7wSqYgQPN1l\/f3JLUfG9KNdu6OqFq0v1f3uKchrNEuTinI9A+S2pYAYMgefLJJ23UqFHWqVMnO\/bYY93fF198sY0bN86Jm7Fjx9oKK6xgxx9\/vHXu3NkuuuiihhAwCJbTT1+wpXcuZXaSJHaE\/3Z4d7BtvvLqNbPFtAbWNJ2QOKfhLAEjzukItM+SGlLAsH33RhttZNdff72tjyvBzKZNm2YffvihrbfeerbffvvZFltsYYcddpj77q233rLtttvOnnrqKevevftCLaHWPTB4VU63021C08JbkeNBYUfJA8Ys2Mgr1lLcavxcNLBWg+rCeYpzGs4SMOKcjkD7LKkhBcwrr7xiu+yyi5smwuMyY8YM23bbbXOeln79+tlZZ53lpo58Wnfdde26666zjTfeuKCAyX64\/\/772wEHHLBIW8w9S95n1yx7tT2+1OP\/u4+mZuvZ80vbY49PbNCgWe5z\/l0v6e2337aePXvWy+3W7X2Kc7qqE+s0rMW5OOdlllnGerApVYOlhhQwTB3ts88+Toz88pe\/tObmZjv11FMNkcLUUd++fe3yyy93Xhif+OzSSy+1rbbaqiY9MKwMmtBrrF1j17SIY8GjwimjB9gBNTMVVO5vRJ6BcsmVdp04l8arEmuxroRe+LXiHM6qkSwbQsDstNNONnnyZFcvm2yyiZ144ok2aNAgu\/3223NTSE888YQNHjzY\/vWvfxkemHPPPdd5ZXzq06eP3XDDDW7qKT8tyikk4liuucZs7Nj\/nYiKaNn7jV\/Z0DW2r+kpoVJ\/KOqESiVWnr04l8etnKvEuhxqpV8jzqUza4QrGkLA0HhZaURaeumlDXcZHpWbb745J0ieffZZ22uvvWzSpEnGFFD\/\/v1t6NCh7hpiYAYOHGgTJ060rl27LnIB4\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\/\/\/Kd16NDBBvx\/t8nw4cOtS5curtauuOIKGz16tM2aNcu23357O+ecc3Lf5VerFzA6fbq6DV6dUHX5+tzFOQ1neWDEOR2B9llSwwqY448\/3ubNm+eEydy5c+3ggw+2DTbYwImacePG2ciRI23s2LG2wgorGLadO3e2iy66qGAryBcwTB8xjaQUl4AG1rg8W8tNnNNwloAR53QE2mdJDStgBg0aZKuuuqqdd955Nnv2bDv00ENt0003teOOO872228\/22KLLeywww5ztf7WW2\/ZdtttZ0899ZR17959oZaAgHnggTesd+8FXzU3t8\/GUu2n1sBabcIL8hfnNJzFWpzTEWifJTWsgLn77rudWEG8kDbZZBO79tprbckll7R+\/frZWWed5aaOfFp33XXtuuuus4033riwgHnjAettva2X9bI37c322Vqq\/NQaWKsM+KvsxTkNZwkYcU5HoH2W1NAC5vnnn3deFuJcjj76aFtnnXXs7LPPtr59+9rll1\/uvDA+8dmll15qW221VUEB82XPL23qP6baUo8tZSfcc4IdcMAB7bPFVPGp3377bevZs2cVS1DWEBDndO1ArNOwFufinJdZZhnr0aNHmspIWEpDCJiddtrJ5bZXgQAAFdZJREFUJk+enPO0jBo1yk0XPfTQQ7b88su7z59++mnba6+97KWXXnIi5dxzz7Vtt902h7pPnz52ww032EYbbVRQwOz\/xv7uBGrtwFu91inPQPXYZnMW5zSc5YER53QE2mdJDSFg6JAJ1CUtvfTS9uWXX7qYFgTMKqus4j5\/7rnnbM8997R\/\/etfNmTIEOvfv78NHTrUfUcMzMCBA23ixInWtWvXggJmmze2sbE21sbYGBtsg9tna6nyU2tgrTLgr7IX5zScJWDEOR2B9llSQwiY\/KqbP3++fec73zHiWs4\/\/3ybM2eOHXPMMc4b8+tf\/9puueUWu\/DCC23MmDG20kor2bBhw6ypqclNIRVKBPGu\/sbqxkGO4228O8hRKT4BDazxmRbKUZzTcJaAEed0BNpnSQ0pYLxX5cwzz3Qrizp27Gg77rijEyp+HximmQjqnTlzpptSGjFihHXr1q1VAdP8RrNNsSkugJdAXqX4BDSwxmcqAZOGaWulqE2n4S\/OaTjXWikNK2Bigl5t69VcAK9WIMWkunBe6oSqy9fnLs5pOMsDI87pCLTPkiRgAupdAiYAUgQTDawRIAZkIc4BkCKZiHUkkG1kI85pONdaKRIwATWy8sp723vvXW8Dxoy18YMVwBuArCwTdUJlYSv5InEuGVnZF4h12ehKulCcS8LVMMYSMAFV+fWv\/8I+\/HCkOz5Ap1AHACvTRJ1QmeBKvEycSwRWgblYVwCvhEvFuQRYDWQqARNQmT16XGz\/\/e\/ROoU6gFUlJuqEKqEXfq04h7Oq1FKsKyUYdr04h3FqNCsJmIAalYAJgBTBRJ1QBIgBWYhzAKRIJmIdCWQb2YhzGs61VooETECNdOr0qM2evbmNGWOmEJgAYGWaqBMqE1yJl4lzicAqMBfrCuCVcKk4lwCrgUwlYAIqs1Onx2z27H4u\/oU4GKXqEFAnVB2u+bmKcxrOlCLWaViLcxrOtVaKBExAjUjABECKYKJOKALEgCzEOQBSJBOxjgSyjWzEOQ3nWitFAiagRjp2fNu+\/LKnvfmmWS9twhtArDwTdULlcSv1KnEulVj59mJdPrtSrhTnUmg1jq0ETEBdNlmTs2q25gBrmZRLQJ1QueRKu06cS+NVibVYV0Iv\/FpxDmfVSJYSMAG1iYDRMQIBoCo0USdUIcDAy8U5EFQEM7GOADEgC3EOgNSAJhIwAZUqARMAKYKJOqEIEAOyEOcASJFMxDoSyDayEec0nGutFAmYgBqRgAmAFMFEnVAEiAFZiHMApEgmYh0JpARMGpB1VooETECFIWAG2AAbb+MDrGVSLgF19uWSK+06cS6NVyXWYl0JvfBrxTmcVSNZSsAE1KYETACkCCbqhCJADMhCnAMgRTIR60gg5YFJA7LOSpGACagwBMxgG2xjbEyAtUzKJaDOvlxypV0nzqXxqsRarCuhF36tOIezaiRLCZiA2pSACYAUwUSdUASIAVmIcwCkSCZiHQmkPDBpQNZZKRIwARWGgDnNTrPh\/\/8\/peoRUGdfPbbZnMU5DWdKEes0rMU5DedaK0UCJqBGmprMTjvNbLj0SwCt8k3UCZXPrpQrxbkUWpXZinVl\/EKvFudQUo1lJwETUJ8IGJ1EHQCqQhN1QhUCDLxcnANBRTAT6wgQA7IQ5wBIDWgiARNQqRIwAZAimKgTigAxIAtxDoAUyUSsI4FsIxtxTsO51kqRgAmoEQmYAEgRTNQJRYAYkIU4B0CKZCLWkUBKwKQBWWelSMAEVJgETACkCCbq7CNADMhCnAMgRTIR60ggJWDSgKyzUiRgAioMATN+vNmAAQHGMimbgDr7stGVdKE4l4SrImOxrghf8MXiHIyqoQwlYAKqUwImAFIEE3VCESAGZCHOAZAimYh1JJDywKQBWWelSMAEVBj7wHAOEuchKVWPgDr76rHN5izOaThTilinYS3OaTjXWikSMAE1IgETACmCiTqhCBADshDnAEiRTMQ6Ekh5YNKArLNSJGACKkwCJgBSBBN19hEgBmQhzgGQIpmIdSSQEjBpQNZZKRIwARWGgHnT3rRe1ivAWiblElBnXy650q4T59J4VWIt1pXQC79WnMNZNZKlBExAbUrABECKYKJOKALEgCzEOQBSJBOxjgRSHpg0IOusFAmYgAqTgAmAFMFEnX0EiAFZiHMApEgmYh0JpARMGpB1VooETECFIWCarTnAUiaVEFBnXwm98GvFOZxVpZZiXSnBsOvFOYxTo1lJwATUqARMAKQIJuqEIkAMyEKcAyBFMhHrSCDlgUkDss5KkYAJqLAeF\/ewGUfPCLCUSSUE1NlXQi\/8WnEOZ1WppVhXSjDsenEO49RoVhIwATW6xhpr2BtvvBFgKZNKCKgTqoRe+LXiHM6qUkuxrpRg2PXiHMap0awkYAJqVAImAFIEk+HDhxt\/lKpLQJyryzebu1inYS3OaTjXWikNIWD++Mc\/2t\/\/\/ne75pprcnw\/\/\/xzO\/nkk+2+++6zxRZbzA466CA74ogjct9fccUVNnr0aJs1a5Ztv\/32ds4551iXLl0K1o8ETJpmK87inIZAulLUptOwFuc0nGutlLoWMJ9++qldcsklduWVV1r\/\/v1bCJhf\/vKX9u6779rvfvc7mz59uhMwQ4cOtb322svGjRtnI0eOtLFjx9oKK6xgxx9\/vHXu3NkuuugiCZhF2ELVCaWBL85pOFOKWKdhLc5pONdaKXUtYPbee2\/nNVl11VXt9ddfzwmYL774wjbccEP705\/+ZBtttJFjfuONN9of\/vAHu\/POO22\/\/fazLbbYwg477DD33VtvvWXbbbedPfXUU9a9e\/eF6kg\/jjTNVpzFOQ2BdKWoTadhLc5pONdaKXUtYN577z1baaWV7KqrrrKHH344J2AI6EKQPP\/887lpoaefftr22Wcfmzx5svXr18\/OOussN3Xk07rrrmvXXXedbbzxxgvVEULp8ccfr7W60\/2IgAiIgAiIQJsEjj76aONPo6W6FjC+MvIFzEsvvWTf\/\/73nVemqanJmb388su288472yuvvOJEyuWXX+68MD717dvXLr30Uttqq60arY71PCIgAiIgAiLQcATqQsDcfffdLQJwmQrKio\/WPDAvvPCCi20hPfvssy7+BQ\/MZpttZueee65tu+22uQrt06eP3XDDDbkpp4araT2QCIiACIiACDQQgboQMDNnzjSmi3zq2bNnTpjwWb6AIQZmgw02sD\/\/+c\/ubxIxMMTE3HbbbS4GhqBfgnpJxMAMHDjQJk6caF27dm2g6tWjiIAIiIAIiEBjEqgLAdMW+nwBg312FRICaPDgwXbAAQfY\/vvvb7fccotdeOGFNmbMGBdDM2zYMDfVxBSSkgiIgAiIgAiIQO0TaFgBw\/4uZ5xxhj3wwAOuFpg+Ou6443IxMaNGjbJrr73WEDfEvYwYMcK6detW+zWmOxQBERABERABEbCGEDCqRxEQAREQAREQgfZFQAKmfdW3nlYEREAEREAEGoKABExDVKMeQgREQAREQATaFwEJmPZV34v8aZubm91SdWKPfGLl13PPPef+ySaEJ554ovv38ssvb6eeemqLDQcX+QPUyQ28+OKLtu++++a4ctttnQ9GvBhngk2dOtW+8Y1vuED33r1718kTL5rbLMTZ7zmVvSM2zSTujiTOpdXVE0884ba9ePXVV91O6XvssYfblI2FF2rTpbFsNGsJmEar0Rp\/HjYS3G233eyee+6xxRdf3N0th22uuOKKNm\/ePPvud79r3\/ve99y+P4899pg77uGuu+7SQBpYrwjEm2++2e00PXv2bJs0aVLuymLng7399tv2ne98x50dxhYDV199tdt2YMKECbl6CryFdmFWjDOLA2iz2bPVOnXqZD169DBxLq15cI4dW1ycdNJJ9sMf\/tC94Bx44IF28MEHu1WlatOl8Ww0awmYRqvRGn8eBkUGR\/byYf8dvDF4XPC2PPnkkzZkyBC36WDHjh3dkxx77LH2ta99zX71q1\/V+JPVxu1dfPHFThzuueeedt555+UETFvngyFcOAsse6L7lltuaaeffro8YAWqtjXOmB511FH24YcfOu\/AZ5995kQ5n9Gmxbm03wm7qo8ePdouuOCC3IVnn322ffDBB+6zYmfeiXVprOvRWgKmHmutju+Zk78RKHQ+Sy+9tJ1\/\/vmuM2KDwZtuusmdR3XHHXfknvCyyy5znhjeapXaJvD++++7E9Y5B+zHP\/5xTsC0dT7Yz3\/+cyciedP1iT2TODfsZz\/7WdsFtzOL1jiDgW0ZmII74YQTnJCB6YABA9x0qDhX1lAQ4jvttJM71w6mxc68E+vKWNfD1RIw9VBLDXSP8+fPNzqhJZdc0j3Vxx9\/bN\/+9rfdMQ7Evdx7771u12Sf2Gzwr3\/9q9tVWSmcACyzAqat88EQKeutt54bYH065JBDbJ111nH7JykVJpDPGau5c+dahw4d3B\/SfffdZ4cffrg7j+3QQw8V5zIb05w5c1zsC3t84ZV57bXXip55pzZdJug6ukwCpo4qq1FvlbOpTjvtNBfYe\/311ztvjE8cuokHJju10agcYj5X\/sDqPTCtnQ\/GwIDnJuuBIcYADwxxSErhAibf0rMnGPXMM88U5zIa07vvvutEIFPPeG95AVKbLgNkg10iAdNgFVrrj8OOyLvvvrv7Q0K0EAeDaCFm4Kc\/\/ak988wzubdXYmCWW265FgNrrT9jLdxfvoBp63ww4gW45ve\/\/33u9pkKYTfr7KGntfBstXQP+ZynTZvmvAIIcb+CCw8MXizOWrvkkkvEucQKxHtIbByCOjudqTZdIsgGNJeAacBKreVH+vWvf+2mhIhpYUkkAbzECTCFxPQSK2GY4z7yyCONN1Zc7kwfMb2hFE6g0NRGsfPBpkyZ4rj\/9re\/tW222cadEzZ27FgbP358brovvPT2Y1mIM1N3BJ4T3zVjxgw76KCDHNtjjjnGxLm0tkF8HOxYecTqRZ\/8qi616dJ4Npq1BEyj1WiNPw\/xAbiAx40b51ZobL311m6lCwGkJNzCp5xyigtCZdkpQb+77LJLjT9V7d1eoYG1rfPBHnzwQbffBkt91157bbcnjIRj8botxJkA3+HDh7upT+Jg2LfkF7\/4Rc6rKM7hvxdWe\/EnP7G0mtWMatPhLBvRUgKmEWtVzyQCIiACIiACDU5AAqbBK1iPJwIiIAIiIAKNSEACphFrVc8kAiIgAiIgAg1OQAKmwStYjycCIiACIiACjUhAAqYRa1XPJAIiIAIiIAINTkACpsErWI8nAuUSYBXS+uuvb7vuumu5Wdgbb7zhlhNzJASHdiqJgAiIQCwCEjCxSCofEWgQAldddZX95Cc\/cZvYfetb33JHEpDY2I49e\/ITZ\/zsu+++BZ+eZfNcx07LO++8c4MQ0mOIgAjUAgEJmFqoBd2DCNQQAQ7I+8tf\/uI8J7feeqstscQSdvDBB7vNBUnsb3LWWWfZnXfe2eKu2Sn1oYceKulJskcb5F\/4yiuv2LBhw9w9sCEcZ2ZdccUVtsMOO7Qw5RRzRNbf\/\/53W3PNNYuW\/84777hnIU9\/HldJNyxjERCBmiEgAVMzVaEbEYHaIJAVMFkPjL87NiJsamoqeMgjm7ituOKKuQfhKAg2xeM8JXZW3nTTTd137Ejbq1evVh+4ubnZfvSjHznRxKZl06dPjyJgKJCpMcSLDqmsjfamuxCBcglIwJRLTteJQIMSaE3AnHfeee7Ih86dO+eenEP2Xn31VbfL7NNPP21Dhw51U09MF\/3tb38zDuO86aab7JZbbnEHct5888121113uR1\/77nnntwOzPkoH3jgATv77LPt\/vvvd1\/FFDDcM16cf\/7zn+44CyUREIH6JCABU5\/1prsWgaoQ4GDNjz\/+2HlYOJuKvxEnq622mhv0OQ34e9\/7nr3++uu2ySabOO\/KpEmTctvks7U+IoZDCzmAj8BdDt3jQE5EDF6YESNG2KhRo4oeU7D\/\/vtb3759c16SUAHz5ZdfuvsrlF5++eXctNGgQYOcHfeqJAIiUJ8EJGDqs9501yJQdQIE39544422yiqruLLwwCBgiDPhME5WFuULGOw4YZxYE6aSnn32WXeWDfExq6++uq2wwgpO0PiTmgs9xKeffupOKL\/uuuucSMp6YMiLKaVsogzyJwaGfDljyyeCjjmLiD8nnHBC7vPf\/OY39uijj7rnUxIBEahPAhIw9VlvumsRqBqB0aNHO5Gw3377uQF+9uzZbvURhxAiYBAiTCXhZckKGIQLwbHXX3+99evXz53IzEniHGRIzAv54vFAmHBC8+DBg52HJz\/5oFympDjQMytgij10fhDvnDlzXHDvMsss48rOLuNmGosYGLxHSiIgAvVJQAKmPutNdy0CVSGAd2TLLbd0sSdMsSBgWFaNN4WpJQQMXpR7773XRo4c2ULAcHJ4t27d3Cof4kyYKlp55ZXdVNQ3v\/lNd3Iwq4mOOeYYd8r4N77xjRZeEf9Ad999tx1++OFumsoLHD+FxHXeK+PtJ0+ebCzlzhcwRx99tLHKiZPPua9sYiXVPvvs404979KlS1VYKlMREIHqEpCAqS5f5S4CdUUA7whTK3hXmEL6wx\/+4LwYCBmCcREwLKvGc3HKKacUnELKPvC8efOczX333WdrrLGGsboIUcLfxKt07NhxIT4E\/x5yyCH22muv5WJrQmNg\/DJqprcIICZoeJ111lmojEf+X3tnjKowFETRtAo2QpahTSoLOwMWAUHEBVi5BC2ziKzEDYQUgp1F9mDjIuQMPAkh+SQ\/iAneqX4wL\/\/lvOYyc2dyvdqsm7\/auAd1cNqsCPwgAQmYHzx0vbII1BGgxINXBCEQhqG3WCys\/IPnxXlgHo+HlYcw6boSEs\/j77ZBJ1JZYOBpoYX6dru9u5TaCBjEEgIoSRJvvV5XbulyuXjn89nL87ztlnW\/CIhATwhIwPTkILQNEegbAQbJYX6lBIM5Fk\/KeDy28s\/xeLSuJH7b7\/eWVcFzUg6yLLPZzFqnycCUg2xO2QeDiTcIAitdLZdLW9JUwJDxQfzsdjvvcDhYlscFZbDRaGSXiLH7\/W5+HYUIiMAwCUjADPPctGsR+CgBjLQMnzudTiYGXND+nKapZWaY1cK17\/u1eymXkJpumq6i+Xzeuo2a8hPenKogi0RZjEDkkJ3Br6MQAREYJgEJmGGem3YtAh8hwMcX8Y9giCVL4UowTNhlEB2zXOg0wpzLNd09+GSYE1MV\/xUwlIHiOPayLKvsVOry8pTEoiiyzx64Lqcuz9NaERCB7xCQgPkOd\/1XEeglATIqz+fTSkTT6fS9R1qhJ5OJZWWKGRfMvRhxt9ttrYDhN0SRmyfT9MU3m419SmC1WjVd0ug+JgXzKYHiXJhGC3WTCIhArwhIwPTqOLQZERABR4DJuZSwyPgUZ7h0IcSAPWbQ8Eznh+nyPK0VARH4HoEX9j9mGUE+HRAAAAAASUVORK5CYII=","height":337,"width":560}}
%---
%[output:4bd508b1]
%   data: {"dataType":"text","outputData":{"text":"线性可行性检验: 等式最小残差 = 3.646e-16\n","truncated":false}}
%---
%[output:09f458ac]
%   data: {"dataType":"text","outputData":{"text":"LS解的插值残差: 3.646e-16\n","truncated":false}}
%---
%[output:1944086b]
%   data: {"dataType":"text","outputData":{"text":"LS解的H∞范数: 10.648813\n","truncated":false}}
%---
%[output:133cb305]
%   data: {"dataType":"text","outputData":{"text":"LS解的参数范数: 4077.548000\n","truncated":false}}
%---
%[output:035dd64a]
%   data: {"dataType":"text","outputData":{"text":"CVX 未安装，使用 LS 解。\n","truncated":false}}
%---
%[output:8d534122]
%   data: {"dataType":"text","outputData":{"text":"跳过 CVX 对比（CVX 未安装）。\n","truncated":false}}
%---
%[output:01378f4d]
%   data: {"dataType":"text","outputData":{"text":"CVX gammma = 1.064881e+01, dense-grid gammma_check = 1.064773e+01, residual ||A theta - b|| = 3.646e-16\n","truncated":false}}
%---
%[output:9ce7464c]
%   data: {"dataType":"image","outputData":{"height":337.44619799139173,"width":560}}
%---
%[output:65ac8863]
%   data: {"dataType":"image","outputData":{"height":337.44619799139173,"width":560}}
%---
%[output:1af9a3e1]
%   data: {"dataType":"image","outputData":{"height":337.44619799139173,"width":560}}
%---
%[output:71f946c6]
%   data: {"dataType":"error","outputData":{"errorType":"runtime","text":"Unrecognized function or variable 'IIR_filter'."}}
%---
