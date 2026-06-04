%% Jafari 2015 TAC 自适应离散 AVC — 修正版
% 对照 book.md Section 6.6.1 逐项修正
% 修复: (1)F=100静态增益 (2)P(0)=20I (3)regressor含G₀F (4)Δₘ修正
clear; close all; clc;

%% ====== 被控对象 (book.md) ======
Ts = 1/480;  fs = 1/Ts;  z = tf('z', Ts);
G0 = (-0.00146*(z-0.1438)*(z-1)) / ((z - 0.7096)*(z^2 - 0.04369*z + 0.01392));
[num0, den0] = tfdata(G0, 'v');
Delta_m = -0.0001 / (z + 0.99)^2;  % book.md
G_true = G0 * (1 + Delta_m);
[num_true, den_true] = tfdata(G_true, 'v');

%% ====== 仿真参数 (book.md Table 4.1) ======
T_end = 50;  Nsim = T_end*fs;  t = (0:Nsim-1)'*Ts;

d1 = 0.7*sin(25*t + pi/3) + 0.5*sin(225*t + pi/4);
d2 = 0.6*sin(85*t - pi/6) + 0.4*sin(125*t + pi/2);
idx_30s = find(t >= 30, 1);
d_total = d1;  d_total(idx_30s:end) = d_total(idx_30s:end) + d2(idx_30s:end);
rng(42);  v = 0.02*randn(Nsim, 1);  v = max(min(v, 0.1), -0.1);
d_total = d_total + v;

N = 50;  F_gain = 100;  c0 = 1;
P0_val = 20;  P = P0_val * eye(N);
theta_max = 10;  theta = zeros(N, 1);

%% ====== 初始化 ======
y = zeros(Nsim,1);  u = zeros(Nsim,1);
zeta_hist = zeros(Nsim,1);  theta_hist = zeros(N,Nsim);

nG0 = max(length(num0), length(den0))-1;
nG = max(length(num_true), length(den_true))-1;
G_state = zeros(nG,1);  G0_state_zeta = zeros(nG0,1);
G0_state_phi = zeros(N, nG0);
buffer_zeta = zeros(N,1);  u_prev = 0;

fprintf('===== 2015 TAC 修正版 (book.md) =====\n');
fprintf('N=%d, F=%.0f, γ₀=%.0f, P(0)=%.0fI\n', N, F_gain, c0, P0_val);
fprintf('Regressor: φ=F·G₀(z)[α(z)[ζ]], F=100静态增益\n');

%% ====== 仿真循环 ======
for k = 1:Nsim
    [y_no_d, G_state] = filter(num_true, den_true, u_prev, G_state);
    y(k) = y_no_d + d_total(k);
    [G0_u, G0_state_zeta] = filter(num0, den0, u_prev, G0_state_zeta);
    zeta_k = y(k) - G0_u;  zeta_hist(k) = zeta_k;

    u_current = 0;
    if t(k) >= 10
        w_vec = flipud(buffer_zeta);

        % φ(k) = F·G₀(z)[α(z)[ζ(k)]] (book.md Eq.)
        phi_vec = zeros(N, 1);
        for i = 1:N
            zd = 0;
            if k > i, zd = zeta_hist(k - i); end
            [go, G0_state_phi(i,:)] = filter(num0, den0, zd, G0_state_phi(i,:)');
            phi_vec(i) = F_gain * go;
        end

        % RLS (book.md Eq., 标准归一化形式)
        pred_err = zeta_k - theta' * phi_vec;
        m2 = 1 + c0 * (phi_vec' * phi_vec);
        epsilon = pred_err / m2;
        P_phi = P * phi_vec;
        denom = m2 + phi_vec' * P_phi;
        if denom > 1e-12
            P = P - (P_phi * P_phi') / denom;
        end
        theta_new = theta + P * epsilon * phi_vec;
        if norm(theta_new) > theta_max
            theta_new = theta_new * (theta_max / norm(theta_new));
        end
        theta = theta_new;

        % 控制律 u(k) = -F·θᵀ·w(k) (book.md)
        u_current = -F_gain * (theta' * w_vec);
        u_current = max(min(u_current, 100), -100);  % ★ 匹配F=100增益
    end

    u(k) = u_current;  theta_hist(:, k) = theta;
    u_prev = u(k);
    buffer_zeta = [zeta_k; buffer_zeta(1:end-1)];
end

%% ====== 性能 ======
idx1 = (t>20 & t<30);  idx2 = (t>40 & t<50);
rms_d1=rms(d_total(idx1)); rms_y1=rms(y(idx1)); att1=20*log10(rms_y1/rms_d1);
rms_d2=rms(d_total(idx2)); rms_y2=rms(y(idx2)); att2=20*log10(rms_y2/rms_d2);

fprintf('\n========== 性能 ==========\n');
fprintf('阶段1 (20-30s): 抑制=%.1f dB, u_RMS=%.2f\n', att1, rms(u(idx1)));
fprintf('阶段2 (40-50s): 抑制=%.1f dB, u_RMS=%.2f\n', att2, rms(u(idx2)));
fprintf('||θ||=%.3f, tr(P)=%.2e\n', norm(theta), trace(P));

% 频域验证 (修正符号和FIR系数)
w_test = [25, 85, 125, 225] * Ts;
fprintf('\n--- 频域验证 τ₀G₀K ---\n');
for wi = 1:length(w_test)
    % u(z)/ζ(z) = -F·(θ₁z^{-N}+θ₂z^{-(N-1)}+...+θ_Nz^{-1})
    % FIR系数: h_i = -F·θ_{N-i+1}
    h = -F_gain * fliplr(theta');
    Hw = sum(h .* exp(-1j*w_test(wi)*(1:N)));
    G0w = evalfr(G0, exp(1j*w_test(wi)));
    Tc = Hw * G0w;  % G₀(z) * u(z)/ζ(z)
    fprintf('  ω=%.0f: H=%.0f∠%.0f°, τ₀G₀K=%.3f∠%.0f°\n', ...
        w_test(wi)/Ts, abs(Hw), angle(Hw)*180/pi, abs(Tc), angle(Tc)*180/pi);
end

%% ====== 绘图 ======
figure('Position',[100 100 1200 700]);
subplot(3,2,1); plot(t,d_total,'Color',[.7 .7 .7]); hold on; plot(t,y,'b-');
xline(10,'g--'); xline(30,'r--'); title(sprintf('Output (%.0f/%.0f dB)',att1,att2)); legend('d','y'); grid on; xlim([0 50]);
subplot(3,2,2); plot(t,u,'r-'); xline(10,'g--'); xline(30,'r--'); title('Control u'); grid on;
subplot(3,2,3); plot(t,zeta_hist,'k-'); title('ζ(k)'); grid on;
subplot(3,2,4); plot(t,theta_hist(1:10:end,:)'); title('θ(t)'); grid on;
subplot(3,2,5); win=fs*0.5; rms_s=sqrt(movmean(y.^2,win)); rms_ds=sqrt(movmean(d_total.^2,win));
plot(t,rms_s,'b-'); hold on; plot(t,rms_ds,'Color',[.7 .7 .7]); xline(10,'g--'); xline(30,'r--'); title('Sliding RMS'); legend('y','d'); grid on;
subplot(3,2,6);
[Pyy,fp]=pwelch(y(end-fs*5:end),hann(512),256,512,fs);
[Pdd,~]=pwelch(d_total(end-fs*5:end),hann(512),256,512,fs);
semilogy(fp(fp<=50),Pdd(fp<=50),'Color',[.7 .7 .7]); hold on;
semilogy(fp(fp<=50),Pyy(fp<=50),'b-');
xline(25/(2*pi),'r--'); xline(85/(2*pi),'r--'); xline(125/(2*pi),'r--'); xline(225/(2*pi),'r--');
title('Spectrum'); legend('d','y'); grid on;

fprintf('Done.\n');
