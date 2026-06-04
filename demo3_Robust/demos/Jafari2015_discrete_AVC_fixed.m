%% Jafari 2015 TAC 自适应离散 AVC — 修正版
% 对照 About_Book_RobustAdaptive.md (原书 6.6.1节) 完整复现
% 修复: (1)无u限幅 (2)协方差重置 (3)含Δₘ (4)F=100, P₀=20I
clear; close all; clc;

%% ====== 被控对象 (book Eqs) ======
Ts = 1/480;  fs = 1/Ts;  z = tf('z', Ts);
G0 = (-0.00146*(z-0.1438)*(z-1)) / ((z - 0.7096)*(z^2 - 0.04369*z + 0.01392));
[num0, den0] = tfdata(G0, 'v');

% ★ 含未建模动态: G = G₀(1+Δₘ) (book p.179)
Delta_m = -0.0001 / (z + 0.99)^2;
G_true = G0 * (1 + Delta_m);
[num_true, den_true] = tfdata(G_true, 'v');

%% ====== 仿真参数 (book Table 4.1) ======
T_end = 50;  Nsim = T_end*fs;  t = (0:Nsim-1)'*Ts;

% 干扰 (book Eqs)
d1 = 0.7*sin(25*t + pi/3) + 0.5*sin(225*t + pi/4);
d2 = 0.6*sin(85*t - pi/6) + 0.4*sin(125*t + pi/2);
idx_30s = find(t >= 30, 1);
d_total = d1;  d_total(idx_30s:end) = d_total(idx_30s:end) + d2(idx_30s:end);
rng(42);  v = 0.02*randn(Nsim, 1);  v = max(min(v, 0.1), -0.1);
d_total = d_total + v;

% 控制器参数 (book Table 4.1)
N = 50;  F_gain = 100;  gamma0 = 1;
P0 = 20;  P = P0 * eye(N);
theta_max = 50;  theta = zeros(N, 1);

% ★ 协方差重置参数 (book p.164, Eq 6.61)
beta0 = 20;   % 重置值 P = β₀I
beta1 = 0.5;  % 触发阈值 λ_min(P) ≤ β₁
cov_reset_count = 0;

%% ====== 初始化 ======
y = zeros(Nsim,1);  u = zeros(Nsim,1);
zeta_hist = zeros(Nsim,1);  theta_hist = zeros(N,Nsim);
pred_err_hist = zeros(Nsim,1);  P_tr_hist = zeros(Nsim,1);

nG0 = max(length(num0), length(den0))-1;
nG = max(length(num_true), length(den_true))-1;
G_state = zeros(nG,1);  G0_state_zeta = zeros(nG0,1);
G0_state_phi = zeros(N, nG0);
buffer_zeta = zeros(N,1);  u_prev = 0;

fprintf('===== 2015 TAC 修正版 (book 6.6.1) =====\n');
fprintf('N=%d, F=%.0f, γ₀=%.0f, P₀=%.0fI\n', N, F_gain, gamma0, P0);
fprintf('协方差重置: β₀=%.0f, β₁=%.1f\n', beta0, beta1);
fprintf('被控对象: G=G₀(1+Δₘ), 控制器内部: G₀ only\n');

%% ====== 仿真循环 ======
for k = 1:Nsim
    % ---- 1. 物理对象 y(k)=G_true[u(k-1)]+d(k) (含Δₘ!) ----
    [y_no_d, G_state] = filter(num_true, den_true, u_prev, G_state);
    y(k) = y_no_d + d_total(k);

    % ---- 2. 辅佐信号 ζ(k)=y(k)-G₀(z)[u(k-1)] (仅用标称模型!) ----
    [G0_u, G0_state_zeta] = filter(num0, den0, u_prev, G0_state_zeta);
    zeta_k = y(k) - G0_u;
    zeta_hist(k) = zeta_k;

    % ---- 3. 自适应控制 ----
    u_current = 0;
    if t(k) >= 10
        w_vec = flipud(buffer_zeta);

        % φ(k) = G₀(z)·F·α(z)[ζ(k)] (book Eq, F=100)
        phi_vec = zeros(N, 1);
        for i = 1:N
            zd = 0;
            if k > i, zd = zeta_hist(k - i); end
            [go, G0_state_phi(i,:)] = filter(num0, den0, zd, G0_state_phi(i,:)');
            phi_vec(i) = F_gain * go;
        end

        % RLS (book Eq 6.17-6.21)
        pred_err = zeta_k - theta' * phi_vec;
        pred_err_hist(k) = pred_err;
        m2 = 1 + gamma0 * (phi_vec' * phi_vec);
        epsilon = pred_err / m2;

        P_phi = P * phi_vec;
        denom = m2 + phi_vec' * P_phi;
        if denom > 1e-12
            P = P - (P_phi * P_phi') / denom;
        end

        % ★ 协方差重置 (book Eq 6.61, 用trace近似避免eig开销)
        P_tr_hist(k) = trace(P);
        if mod(k, 100) == 0 && P_tr_hist(k) / N < beta1
            P = beta0 * eye(N);
            cov_reset_count = cov_reset_count + 1;
        end

        theta_new = theta + P * epsilon * phi_vec;
        if norm(theta_new) > theta_max
            theta_new = theta_new * (theta_max / norm(theta_new));
        end
        theta = theta_new;

        % 控制律 u(k) = -F·θᵀ·w(k) (book Eq, ★无u限幅!)
        u_current = -F_gain * (theta' * w_vec);
    end

    u(k) = u_current;
    theta_hist(:, k) = theta;
    u_prev = u(k);
    buffer_zeta = [zeta_k; buffer_zeta(1:end-1)];
end

%% ====== 性能评估 ======
idx1 = (t>20 & t<30);  idx2 = (t>40 & t<50);
rms_d1=rms(d_total(idx1)); rms_y1=rms(y(idx1)); att1=20*log10(rms_y1/rms_d1);
rms_d2=rms(d_total(idx2)); rms_y2=rms(y(idx2)); att2=20*log10(rms_y2/rms_d2);

fprintf('\n========== 性能 ==========\n');
fprintf('阶段1 (20-30s, 2频率): 抑制=%.1f dB, u_RMS=%.1f\n', att1, rms(u(idx1)));
fprintf('阶段2 (40-50s, 4频率): 抑制=%.1f dB, u_RMS=%.1f\n', att2, rms(u(idx2)));
fprintf('||θ||=%.2f, tr(P)=%.2e, 协方差重置=%d次\n', norm(theta), trace(P), cov_reset_count);

%% ====== 绘图 ======
figure('Position',[100 100 1200 700]);
subplot(3,2,1); plot(t,d_total,'Color',[.7 .7 .7]); hold on; plot(t,y,'b-');
xline(10,'g--'); xline(30,'r--'); title(sprintf('Output (%.0f/%.0f dB)',att1,att2)); legend('d','y'); grid on; xlim([0 50]);
subplot(3,2,2); plot(t,u,'r-'); xline(10,'g--'); xline(30,'r--'); title('Control u (no saturation)'); grid on;
subplot(3,2,3); plot(t,zeta_hist,'k-'); title('ζ(k)'); grid on;
subplot(3,2,4); plot(t,theta_hist(1:10:end,:)'); title('θ(t)'); grid on;
subplot(3,2,5); win=fs*0.5; rms_s=sqrt(movmean(y.^2,win)); rms_ds=sqrt(movmean(d_total.^2,win));
plot(t,rms_s,'b-'); hold on; plot(t,rms_ds,'Color',[.7 .7 .7]); xline(10,'g--'); xline(30,'r--'); title('Sliding RMS'); legend('y','d'); grid on;
subplot(3,2,6);
[Pyy,fp]=pwelch(y(end-5*fs:end),hann(512),256,512,fs);
[Pdd,~]=pwelch(d_total(end-5*fs:end),hann(512),256,512,fs);
semilogy(fp(fp<=50),Pdd(fp<=50),'Color',[.7 .7 .7]); hold on;
semilogy(fp(fp<=50),Pyy(fp<=50),'b-');
xline(25/(2*pi),'r--'); xline(85/(2*pi),'r--'); xline(125/(2*pi),'r--'); xline(225/(2*pi),'r--');
title('Steady-state spectrum'); legend('d','y'); grid on;

fprintf('Done.\n');
