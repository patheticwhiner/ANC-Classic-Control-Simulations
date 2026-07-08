function MarinoTomei_ANC_freq_estimator()
%% Marino & Tomei (2016) 自适应扰动抑制 — ANC 次级路径
%
% 严格遵循 docs/MarinoTomei_2016_Derivation.md 的 Eq.(5)-(6).
%
% 控制律:  u = -w0_hat - sum(w_{2i-1}_hat)
% 估计器:  dw/dt = f(w, y, sign{Re[P(j*omega_i)]})
% 自适应:  dtheta_i/dt = eps * sign{Re[P(j*omega_i)]} * w_{2i} * y
%
% 关键: sign{Re[P(j*omega_i)]} 隐含处理被控对象动态.
%       S(z) 在 50/120 Hz 处 Re[S]>0 → 符号均为 +1.

close all; clc;

%% ====== 系统 ======
fs = 4000;  Ts = 1/fs;  Tsim = 8;  N = Tsim * fs;

% 次级路径 S(z) = z^-2 + 1.5*z^-3 - z^-4
S_num = [0, 0, 1, 1.5, -1];  S_len = length(S_num);
S_buf = zeros(S_len, 1);

% 扰动
om1 = 2*pi*50;  om2 = 2*pi*120;  om3 = 2*pi*300;

%% ====== 预计算 sign{Re[S(e^{j*omega*Ts})]} ======
% 来自 Eq.(5)-(6) 的 sign 函数 — 这是 S(z) 知识注入的唯一入口
S_w1 = S_num * exp(-1j * om1*Ts * (0:S_len-1)');  sgn1 = sign(real(S_w1));
S_w2 = S_num * exp(-1j * om2*Ts * (0:S_len-1)');  sgn2 = sign(real(S_w2));
sgn0 = sign(sum(S_num));  % sign[S(1)] = sign[DC gain]
fprintf('=== ANC Marino-Tomei 2016 (sign 预计算) ===\n');
fprintf('sign[S(1)]     = %+d\n', sgn0);
fprintf('sign[Re S(j50)] = %+d  (|S|=%.3f, ∠S=%.1f°)\n', sgn1, abs(S_w1), angle(S_w1)*180/pi);
fprintf('sign[Re S(j120)]= %+d  (|S|=%.3f, ∠S=%.1f°)\n', sgn2, abs(S_w2), angle(S_w2)*180/pi);
fprintf('扰动: 50+120 Hz (+300Hz 未建模)\n');
fprintf('控制律: u = -w0 - w1 - w3  (Eq.5, 无 S 补偿)\n');

%% ====== 估计器参数 ======
k_est = 200;  eps_est = 1.0;   m = 2;
theta = [0.005; 0.03];     % 初始偏离 (真值: (2π·f/fs)² = [0.00617, 0.03553])
w_hat = zeros(2*m+1, 1);   % [w0, w1, w2, w3, w4]

%% ====== 日志 ======
Nlog = 200;  di = max(1, round(N/Nlog));
t_log = zeros(Nlog,1);  f_log = zeros(m,Nlog);
u_log = zeros(Nlog,1);  e_log = zeros(Nlog,1);
li = 1;

for n = 1:N
    tt = (n-1) * Ts;
    d = sin(om1*tt) + sin(om2*tt) + 0.3*sin(om3*tt);
    y_sec = S_num * S_buf;
    e = d + y_sec;

    if n == 1
        t_log(1)=0; f_log(:,1)=[0;0]; u_log(1)=0; e_log(1)=e; li=2;
        continue;
    end

    %% === Eq.(5): 估计器状态更新 (sign{Re[P(j*omega_i)]} 注入) ===
    y = e;
    dw = zeros(2*m+1, 1);
    dw(1) = sgn0 * k_est * y;               % ẇ₀ = sign[P(0)]·k·y
    dw(2) = w_hat(3) + sgn1 * k_est * y;    % ẇ₁ = w₂ + sign[P(jω₁)]·k·y
    dw(3) = -theta(1) * w_hat(2);           % ẇ₂ = -θ₁·w₁
    dw(4) = w_hat(5) + sgn2 * k_est * y;    % ẇ₃ = w₄ + sign[P(jω₂)]·k·y
    dw(5) = -theta(2) * w_hat(4);           % ẇ₄ = -θ₂·w₃

    %% === Eq.(6): 频率自适应 ===
    dth1 =  sgn1 * eps_est * w_hat(3) * y;  % θ̇₁ = sign[P(jω₁)]·ε·w₂·y
    dth2 =  sgn2 * eps_est * w_hat(5) * y;  % θ̇₂ = sign[P(jω₂)]·ε·w₄·y

    w_hat = w_hat + Ts * dw;
    theta(1) = max(theta(1) + Ts * dth1, 1e-8);
    theta(2) = max(theta(2) + Ts * dth2, 1e-8);

    %% === Eq.(5): 控制律 u = -w0 - w1 - w3 ===
    u = -(w_hat(1) + w_hat(2) + w_hat(4));

    %% S(z) 滤波
    S_buf = [u; S_buf(1:end-1)];

    %% 日志
    if mod(n, di) == 0 && li <= Nlog
        t_log(li) = tt;
        f_log(:,li) = sqrt(max(theta,1e-10)) / (2*pi*Ts);
        u_log(li) = u;  e_log(li) = e;
        li = li + 1;
    end
end

%% 输出
fprintf('-----------------------------------------------\n');
fprintf('theta 终值:  [%.6f %.6f]  (真值 [%.4f %.4f])\n', theta, (om1*Ts)^2, (om2*Ts)^2);
fprintf('频率 终值:   [%.1f %.1f] Hz  (真值 [%.0f %.0f])\n', ...
    sqrt(theta(1))*fs/(2*pi), sqrt(theta(2))*fs/(2*pi), om1/(2*pi), om2/(2*pi));
fprintf('控制前 d RMS: %.3f   控制后 e RMS: %.3f\n', sqrt(1), rms(e_log));

%% 绘图
figure('Name','频率估计','Position',[100 80 800 350]);
subplot(1,2,1); plot(t_log,f_log(1,:),'b','LineWidth',1.5); hold on; yline(50,'r--');
ylabel('Hz'); title('f_1 (真值 50 Hz)'); grid on;
subplot(1,2,2); plot(t_log,f_log(2,:),'b','LineWidth',1.5); hold on; yline(120,'r--');
ylabel('Hz'); title('f_2 (真值 120 Hz)'); grid on;

figure('Name','信号','Position',[150 130 800 350]);
subplot(2,1,1); plot(t_log, e_log); title('误差 e(t)'); grid on;
subplot(2,1,2); plot(t_log, u_log); title('控制 u(t)'); grid on;

end
