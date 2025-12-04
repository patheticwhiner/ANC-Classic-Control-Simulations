clear; close all; clc;
addpath('functions');
%% 简化的测试用例
fs = 10000;              % 降低采样频率便于观察
Ts = 1/fs;              
Tsim = 10;               % 缩短仿真时间
Nsim = Tsim*fs;         
t = (0:Nsim-1)*Ts;

%% 简单的测试系统
s = tf('s');
G0 = 0.5*(s-0.2)/(s^2+s+1.25);   
G = G0;                          % 无不确定性
epsi = 0;

%% 简单扰动信号
% 固定频率分量
A = [0.6 0.7];          % 幅值
f = [70 187];           % Hz
phi = [pi/4 pi/2];         % 相位
d = zeros(1,Nsim);
for k=1:length(A)
    d = d + A(k)*sin(2*pi*f(k)*t + phi(k));
end
% 加入小随机噪声
d = d + 0.02*randn(1,Nsim);

%% 简化控制器
% 固定滤波器 F(s)
mu = 300; m = 2; 
k0 = 0.5; alpha = 500;
F = k0 * 2*alpha^2*(s^2+s+1.25)/((s+alpha)^2*(s+0.2));   

Fd = c2d(F, Ts, 'zoh');
[numFd, denFd] = tfdata(Fd, 'v');
% 系统离散化
Gd = c2d(G, Ts, 'zoh');
[numGd, denGd] = tfdata(Gd, 'v');
% G0离散化
G0d = c2d(G0, Ts, 'zoh');
[numG0d, denG0d] = tfdata(G0d, 'v');

%% 简化的Λ(s) - 只用6个参数
Nparam = 6;
lambda_val = 500;

% 这里使用 alpha(z) 作为时移基函数替代 Λ(s)
% alpha(z) = [z^{-N}, z^{-(N-1)}, ..., z^{-1}]^T，脚本中将对每个延迟分量
% 单独通过 F 和 G0 滤波器以得到回归量 φ。

%% 自适应参数 (论文设置)
P = eye(Nparam)*500;
gamma0 = 1.0;

%---------------- 历史缓冲区与变量命名规范化 ----------------%

% 主输出与控制输入
y = zeros(1, Nsim);           % 系统输出
u = zeros(1, Nsim);           % 控制输入
theta = zeros(Nparam,1);      % 自适应参数
phi_reg = zeros(Nparam, 1);   % 回归向量

% 历史缓冲区长度
hist_len = Nparam + 1;        % 历史长度 = 参数数+1

% 滤波器历史缓冲区（统一命名：*_hist_滤波器名）
F_u_hist = zeros(1, hist_len);      % F滤波器输入历史
F_y_hist = zeros(1, hist_len);      % F滤波器输出历史
G_u_hist = zeros(1, hist_len);      % G滤波器输入历史
G_y_hist = zeros(1, hist_len);      % G滤波器输出历史
G0_u_hist = zeros(1, hist_len);     % G0滤波器输入历史
G0_y_hist = zeros(1, hist_len);     % G0滤波器输出历史
% 对于 alpha(z) 的每个延迟分量，保留独立的 F 和 G0 滤波器历史
phi_G0_u_hist = zeros(Nparam, hist_len); % phi链路G0输入历史（每行对应一个延迟分量）
phi_G0_y_hist = zeros(Nparam, hist_len); % phi链路G0输出历史
phi_F_u_hist = zeros(Nparam, hist_len);  % phi链路F输入历史
phi_F_y_hist = zeros(Nparam, hist_len);  % phi链路F输出历史

% 调试与监控变量
plant_output_history = zeros(1, Nsim);   % G系统输出（反噪声）
theta_history = zeros(Nparam, Nsim);     % 参数历史
P_trace_history = zeros(1, Nsim);        % P矩阵迹历史
phi_history = zeros(Nparam, Nsim);       % 回归向量历史

% 关键信号历史（用于监控和绘图）
z_signal_history = zeros(1, Nsim);       % z(t) = y(t) - G0(s)u(t)
epsilon_history = zeros(1, Nsim);        % ε(t) = (z - θ^T φ) / m_s^2
prediction_error_history = zeros(1, Nsim); % 预测误差 z - θ^T φ
m_s_squared_history = zeros(1, Nsim);    % m_s^2 = 1 + γ0 φ^T φ
theta_norm_history = zeros(1, Nsim);     % ||θ(t)||
K_z_history = zeros(1, Nsim);            % K(s,θ(t))z(t) = θ^T(t)φ(t)
phi_norm_history = zeros(1, Nsim);       % ||φ(t)||

fprintf('开始调试仿真...\n');

%% 仿真循环

% --- 新的因果仿真顺序 ---
% 1. 用u(k-1)仿真系统，得到y(k)
% 2. 用y(k)、u(k-1)计算z(k)、phi(k)
% 3. 用z(k)、phi(k)更新theta(k)、P(k)
% 4. 用新theta(k)和观测，计算u(k)
for k = hist_len+1:Nsim
    % 1. 系统仿真：用上一时刻的u(k-1)得到y(k)
    [y_plant_k, G_u_hist, G_y_hist] = IIR_filter(numGd, denGd, u(k-1), G_u_hist, G_y_hist);
    plant_output_history(k) = y_plant_k;
    y(k) = y_plant_k + d(k);

    % 2. 观测信号z(k) = y(k) - G0(s)u(k-1)
    [G0_u_k, G0_u_hist, G0_y_hist] = IIR_filter(numG0d, denG0d, u(k-1), G0_u_hist, G0_y_hist);
    z_k = y(k) - G0_u_k;
    z_signal_history(k) = z_k;

    % 3. phi(k) = G0(z) F(z) alpha(z)[z(k)]，其中 alpha 提供时移基函数
    % alpha(z) = [z^{-Nparam}, z^{-(Nparam-1)}, ..., z^{-1}]^T
    for i = 1:Nparam
        % 获取延迟的 z 值：z(k - i)
        z_delay = 0;
        if k - i >= 1
            z_delay = z_signal_history(k - i);
        end
        % 先通过 F，再通过 G0：为每个延迟分量维护独立历史
        [f_out_i, phi_F_u_hist(i,:), phi_F_y_hist(i,:)] = IIR_filter(numFd, denFd, z_delay, phi_F_u_hist(i,:), phi_F_y_hist(i,:));
        [phi_reg(i), phi_G0_u_hist(i,:), phi_G0_y_hist(i,:)] = IIR_filter(numG0d, denG0d, f_out_i, phi_G0_u_hist(i,:), phi_G0_y_hist(i,:));
    end
    phi_history(:, k) = phi_reg;

    % 4. 参数自适应律更新（离散时间正规化 RLS）
    % 参考文献(6.17-6.21)：采用离散时间的正规化最小二乘递推，而不是
    % 对连续时间方程做 Euler 近似（去掉 Ts 因子）
    % 假设：使用归一化因子 m^2(t)=1+gamma0*phi'phi（同文献），并实现简单的
    % 欧氏投影保证 ||theta|| <= theta_max（文献建议在P-1/2度量下投影，这里为简化采取欧氏投影，若需严格实现可继续改进）。
    m_s_squared = 1 + gamma0 * (phi_reg.' * phi_reg);
    prediction_error = z_k - theta.' * phi_reg;    % ζ(t) - φ^T θ
    epsilon = prediction_error / m_s_squared;       % 归一化误差 ε(t)

    % 规范化 RLS 增益与协方差更新（离散递推，无 Ts 因子）
    denom = m_s_squared + (phi_reg.' * P * phi_reg); % 标量
    if denom <= 0
        denom = 1e-12; % 防止除零/负值（数值保护）
    end
    K = (P * phi_reg) / m_s_squared;                        % 增益向量
    theta_new = theta + K * epsilon;                  % 参数递推（未投影）
    P = P - (P * (phi_reg * phi_reg.') * P) / denom;  % 协方差更新

    % 投影到有界集合 S: ||theta||_2 <= theta_max
    % 这里做简单的欧氏投影；文献中建议按 P^{-1/2} 度量投影（更严格）
    if ~exist('theta_max', 'var')
        theta_max = 1e3; % 合理默认值；可由用户调整
    end
    theta_norm_new = norm(theta_new);
    if theta_norm_new > theta_max
        theta = (theta_new / theta_norm_new) * theta_max;
        % 可选：调整 P 以反映投影（此处不修改 P）
        fprintf('投影触发：||theta||=%.3g > theta_max=%.3g，已缩放。\n', theta_norm_new, theta_max);
    else
        theta = theta_new;
    end

    % 保存历史
    theta_history(:, k) = theta;
    P_trace_history(k) = trace(P);

    % 5. 用新theta(k)和phi(k)计算控制输入u(k)
    K_z_k = theta.' * phi_reg;
    [u(k), F_u_hist, F_y_hist] = IIR_filter(numFd, denFd, -K_z_k, F_u_hist, F_y_hist);

    % 保存关键变量历史数据（用于监控和绘图）
    epsilon_history(k) = epsilon;
    prediction_error_history(k) = prediction_error;
    m_s_squared_history(k) = m_s_squared;
    theta_norm_history(k) = norm(theta);
    K_z_history(k) = K_z_k;
    phi_norm_history(k) = norm(phi_reg);

    % 调试信息
    if k < hist_len+10
        fprintf('k=%d: u(k)=%.6f, y_plant_k=%.6f, ratio=%.2f\n', ...
            k, u(k), y_plant_k, abs(y_plant_k/max(abs(u(k)), 1e-10)));
    end

    % 检查发散
    if any(abs(theta) > 1000) || trace(P) > 1e6
        fprintf('检测到发散！时刻 k=%d\n', k);
        fprintf('theta范数: %.2f\n', norm(theta));
        fprintf('P trace: %.2e\n', trace(P));
        fprintf('z_k: %.6f\n', z_k);
        fprintf('phi_reg范数: %.6f\n', norm(phi_reg));
        fprintf('epsilon: %.6f\n', epsilon);
        break;
    end

    if mod(k, 200) == 0
        fprintf('k=%d: theta_norm=%.3f, P_trace=%.2e, z=%.4f\n', ...
                k, norm(theta), trace(P), z_k);
    end
end


% 还原控制器 K(z) = sum_{i=1..Nparam} theta(i) z^{-i}（使用最后时刻的 theta）
theta_final = theta;
if exist('theta_history','var') && size(theta_history,2) >= Nsim
    theta_final = theta_history(:,end);
end

% 构造 FIR 系数（z^-1 顺序）
b_K = [0; flip(theta_final(:))]; % K(z) = -[0 theta1 theta2 ...] * z^{-0..N}
% b_K 定义为以 z^0, z^-1, ... 的系数向量（第一个系数为0，对应无直流项）
% 绘制频率响应
nfft = 1024;
[Hk, fk] = freqz(b_K, 1, nfft, fs);
figure;
subplot(2,1,1);
semilogx(fk, 20*log10(abs(Hk)), 'k-', 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on;
title('Frequency response of K(z)');
subplot(2,1,2);
semilogx(fk, angle(Hk)*180/pi, 'k-', 'LineWidth', 2);
xlabel('Frequency (Hz)'); ylabel('Phase (deg)'); grid on;
title('Phase of K(z)');


%% ----------- 关键时域对比图 -----------
figure;
subplot(3,2,1);
plot(t, d, 'b-', t, plant_output_history, 'r--', t, y, 'm:','LineWidth',1.2);
legend('扰动 d(t)', 'G输出 y_G(t)', '总输出 y(t)');
title('输出信号对比'); grid on;

subplot(3,2,2);
plot(t, u, 'g-', 'LineWidth',1.2); title('控制信号 u(t)'); grid on;

subplot(3,2,3);
plot(t, z_signal_history, 'k-', 'LineWidth',1.2); title('观测信号 z(t)'); grid on;

subplot(3,2,4);
plot(t, theta_history.'); title('参数演化 θ(t)'); grid on;

subplot(3,2,5);
plot(t, P_trace_history, 'LineWidth',1.2); title('P矩阵迹 trace(P)'); grid on;

subplot(3,2,6);
plot(t, phi_norm_history, 'LineWidth',1.2); title('回归向量范数 ||φ(t)||'); grid on;


%% ----------- 关键收敛与误差对比图 -----------
figure;
subplot(2,2,1);
plot(t, prediction_error_history, 'r-', 'LineWidth',1.2); hold on;
plot(t, epsilon_history, 'b-', 'LineWidth',1.2);
legend('预测误差', '归一化误差ε');
title('误差对比'); grid on;

subplot(2,2,2);
plot(t, m_s_squared_history, 'LineWidth',1.2); title('归一化因子 m_s^2'); grid on;

subplot(2,2,3);
plot(t, theta_norm_history, 'LineWidth',1.2); title('参数范数 ||θ(t)||'); grid on;

subplot(2,2,4);
yyaxis left; plot(t, K_z_history, 'b-', 'LineWidth',1.2); ylabel('K(s,θ)z(t)');
yyaxis right; plot(t, phi_norm_history, 'r-', 'LineWidth',1.2); ylabel('||φ(t)||');
title('K(s,θ)z 与 φ范数'); grid on;


%% ----------- 关键性能指标与频谱 -----------
% 计算抑制效果（最后1秒的RMS值）
final_idx = round(0.8*Nsim):Nsim;
rms_output = sqrt(mean(y(final_idx).^2));
rms_disturbance = sqrt(mean(d(final_idx).^2));
rms_y_G = sqrt(mean(plant_output_history(final_idx).^2));
suppression_ratio = 20*log10(rms_output/rms_disturbance);

figure;
subplot(2,2,1);
bar([rms_disturbance, rms_y_G, rms_output]);
set(gca, 'XTickLabel', {'扰动RMS', 'G输出RMS', '总输出RMS'});
ylabel('RMS Value');
title(sprintf('抑制效果: %.1f dB', suppression_ratio)); grid on;

subplot(2,2,2);
variables = [mean(abs(epsilon_history(final_idx))), ...
            mean(abs(prediction_error_history(final_idx))), ...
            mean(m_s_squared_history(final_idx)), ...
            theta_norm_history(end)];
bar(variables);
set(gca, 'XTickLabel', {'|ε|均值', '|预测误差|均值', 'm_s^2均值', '||θ||终值'});
ylabel('Values');
title('自适应算法关键指标'); grid on;

subplot(2,2,3);
convergence_window = round(0.5*Nsim):Nsim;
plot(t(convergence_window), theta_norm_history(convergence_window), 'b-', 'LineWidth', 2); hold on;
plot(t(convergence_window), P_trace_history(convergence_window)/1000, 'r-', 'LineWidth', 2);
xlabel('Time [s]'); ylabel('Normalized Values');
title('收敛性分析（后半程）');
legend('||θ(t)||', 'trace(P(t))/1000', 'Location', 'best'); grid on;

subplot(2,2,4);
N_fft = 1024;
final_segment = y(end-N_fft+1:end) - mean(y(end-N_fft+1:end));
Y_fft = abs(fft(final_segment));
f_axis = (0:N_fft/2-1) * fs / N_fft;
semilogy(f_axis, Y_fft(1:N_fft/2), 'b-', 'LineWidth', 1.2);
xlabel('频率 [Hz]'); ylabel('幅值');
title('输出信号频谱（稳态段）'); grid on;

% 显示详细性能指标
fprintf('\n=== 自适应控制性能评估 ===\n');
fprintf('扰动RMS: %.4f\n', rms_disturbance);
fprintf('G输出RMS: %.4f\n', rms_y_G);
fprintf('总输出RMS: %.4f\n', rms_output);
fprintf('抑制效果: %.1f dB\n', suppression_ratio);
fprintf('最终参数范数: %.4f\n', norm(theta));
fprintf('最终P矩阵迹: %.4f\n', trace(P));
fprintf('稳态预测误差RMS: %.6f\n', sqrt(mean(prediction_error_history(final_idx).^2)));
fprintf('稳态ε RMS: %.6f\n', sqrt(mean(epsilon_history(final_idx).^2)));
fprintf('稳态φ范数均值: %.4f\n', mean(phi_norm_history(final_idx)));

fprintf('调试完成\n');