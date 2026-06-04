%% Jafari 2017 JVC CC — 极速版
% 核心优化: Λ_{N}到Λ_{1}共享同一级联链,一轮扫描算出全部N个输出
% Λ运算: 210次filter → 20次 (10.5×), 总: 250次→60次 (4.2×)
clear; close all; clc;

%% ====== 论文参数 ======
tic_total = tic;
fs = 10000;  Ts = 1/fs;  Tsim = 10;  Nsim = Tsim*fs;  t = (0:Nsim-1)*Ts;
s = tf('s');
G0 = 0.5*(s-0.2)/(s^2+s+1.25);
kappa0 = 0.5;  alpha_val = 500;
F = kappa0 * 2*alpha_val^2*(s^2+s+1.25)/((s+alpha_val)^2*(s+0.2));

%% 扰动
omega = [70, 187];  A_m = [0.6 0.7];  phi0 = [pi/4 pi/2];
d = zeros(1,Nsim);
for k = 1:length(omega), d = d + A_m(k)*sin(omega(k)*t + phi0(k)); end
rng(42); d = d + 0.02*randn(1,Nsim);

%% ====== 预计算滤波器系数 ======
Nparam = 20;  lambda_val = 500;

% 一阶节 G1(z) = λ/(s+λ) → Tustin
G1_d = c2d(tf(lambda_val,[1 lambda_val]), Ts, 'tustin');
[numG1, denG1] = tfdata(G1_d, 'v');  % 1阶 → 各2系数

% F, G0 滤波器系数
Fd = c2d(F, Ts, 'tustin');
[numF, denF] = tfdata(Fd, 'v');
G0d = c2d(G0, Ts, 'tustin');
[numG0, denG0] = tfdata(G0d, 'v');

%% ====== 滤波器状态初始化 ======
% ★ 关键: 所有Λ共享同一级联链,仅需Nparam个一阶状态!
Lambda_zf = zeros(Nparam, 1);  % Nparam=20个一阶节的终态

% phi链路F/G0状态 (每Λ_i独立)
nF = max(length(numF), length(denF)) - 1;
nG0 = max(length(numG0), length(denG0)) - 1;
phi_F_zf  = zeros(Nparam, nF);
phi_G0_zf = zeros(Nparam, nG0);

% 被控对象/观测/控制状态
G_zf = zeros(nG0, 1);
G0_zf = zeros(nG0, 1);
F_ctrl_zf = zeros(nF, 1);

%% ====== 自适应参数 ======
P = eye(Nparam) * 500;
gamma0 = 1.0;
theta_max_val = 5;
theta = zeros(Nparam, 1);

%% ====== 信号存储 ======
y = zeros(1,Nsim);  u = zeros(1,Nsim);
theta_hist = zeros(Nparam,Nsim);  P_tr_hist = zeros(1,Nsim);
z_hist = zeros(1,Nsim);

fprintf('===== 极速版 2017 JVC CC =====\n');
fprintf('Λ计算: %d次filter/样本 (原210次)\n', Nparam);
fprintf('N=%d, Tsim=%ds, Nsim=%d\n', Nparam, Tsim, Nsim);

%% ====== 仿真循环 ======
tic_sim = tic;
lambda_vals = zeros(Nparam, 1);
phi_reg = zeros(Nparam, 1);

for k = 2:Nsim
    % ---- 1. 被控对象 y(k)=G[u(k-1)]+d(k) ----
    [y_G_k, G_zf] = filter(numG0, denG0, u(k-1), G_zf);
    y(k) = y_G_k + d(k);

    % ---- 2. z(k)=y(k)-G₀[u(k-1)] ----
    [G0_u_k, G0_zf] = filter(numG0, denG0, u(k-1), G0_zf);
    z_k = y(k) - G0_u_k;
    z_hist(k) = z_k;

    % ---- 3. ★ 一遍级联扫描算出全部 Λ_i[z] ★ ----
    % 共享级联链: z → G1 → G1² → ... → G1^N
    % stage s 输出 = G1^s[z] = Λ_{N-s+1}[z]
    y_stage = z_k;
    for stage = 1:Nparam
        [y_stage, Lambda_zf(stage)] = filter(numG1, denG1, y_stage, Lambda_zf(stage));
        lambda_vals(Nparam - stage + 1) = y_stage;
    end

    % F[Λ_i[z]] + G₀[F[Λ_i[z]]] → φ_i (独立状态,无法避免循环)
    for i = 1:Nparam
        [f_out, phi_F_zf(i,:)] = filter(numF, denF, lambda_vals(i), phi_F_zf(i,:)');
        [phi_reg(i), phi_G0_zf(i,:)] = filter(numG0, denG0, f_out, phi_G0_zf(i,:)');
    end

    % ---- 4. 自适应律 ----
    m_s_sq = 1 + gamma0 * (phi_reg.' * phi_reg);
    pred_err = z_k - theta.' * phi_reg;
    epsilon = pred_err / m_s_sq;
    P_phi = P * phi_reg;
    P = P - Ts * (P_phi * P_phi.') / m_s_sq;
    theta_new = theta + Ts * P * epsilon * phi_reg;
    th_norm = norm(theta_new);
    if th_norm > theta_max_val
        theta = theta_new * (theta_max_val / th_norm);
    else
        theta = theta_new;
    end
    theta_hist(:,k) = theta;
    P_tr_hist(k) = trace(P);

    % ---- 5. 控制律 u=-F[θᵀΛ[z]] ----
    Kz = theta.' * lambda_vals;
    [u(k), F_ctrl_zf] = filter(numF, denF, -Kz, F_ctrl_zf);

    if mod(k, round(Nsim/8)) == 0
        fprintf('  t=%.1fs: ||θ||=%.3f, tr(P)=%.2e\n', t(k), th_norm, trace(P));
    end
end
sim_time = toc(tic_sim);

%% ====== 性能评估 ======
final_idx = round(0.5*Nsim):Nsim;
rms_d = sqrt(mean(d(final_idx).^2));
rms_y = sqrt(mean(y(final_idx).^2));
rms_u = sqrt(mean(u(final_idx).^2));
supp_db = 20*log10(rms_y/rms_d);

fprintf('\n========== 性能统计 ==========\n');
fprintf('抑制效果: %.1f dB | 控制RMS: %.0f | ||θ||: %.2f\n', supp_db, rms_u, norm(theta));
fprintf('仿真耗时: %.1f s (%d 样本/s)\n', sim_time, round(Nsim/sim_time));
fprintf('vs filter版: 理论 %.0f× (250→60次filter/样本)\n', 250/60);

%% ====== 绘图 ======
out_dir = 'E:/Code/MATLAB_ANC/ANC-Classic-Control-Simulations/demo3_Robust/demos/';

figure('Position', [100 100 1000 400]);
subplot(1,2,1);
plot(t, d, 'Color', [0.7 0.7 0.7]); hold on;
plot(t, y, 'b-', 'LineWidth', 1);
title(sprintf('输出 vs 扰动 (抑制 %.1f dB)', supp_db)); legend('d','y'); grid on;
subplot(1,2,2);
plot(t, theta_hist(1:5:end,:)');
title(sprintf('θ(t) (耗时%.1fs)', sim_time)); grid on;
saveas(gcf, [out_dir 'VFast_Output.png']);
fprintf('图片已保存。\n');
