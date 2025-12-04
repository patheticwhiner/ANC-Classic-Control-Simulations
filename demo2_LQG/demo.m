%% T12-Qk 非对易性验证实验
clear; clc; close all;

% 基本设置
N = 2000;
fs = 2000;                   % 采样频率
t = (0:N-1)'/fs;

% 固定 LTI 系统 T12(z) = 1 / (1 - 0.5 z^-1)
b_T12 = 1; 
a_T12 = [1 -0.5];
T12_state_A = 0;  % 路径A状态
T12_state_B = 0;  % 路径B状态

% 构造时变参数 Qk(z^-1) = θ1(k) + θ2(k) z^-1
theta1 = 0.5 + 0.5*sin(20*pi*t/0.25);  % 4Hz 变化
theta2 = 0.3 + 0.3*cos(20*pi*t/0.35);  % ~2.8Hz 变化

% 输入信号 r(k)
r = sin(2*pi*100*t) + 0.3*randn(N,1);

% 初始化
yA = zeros(N,1);
yB = zeros(N,1);
r_T12 = zeros(N,1);
q_out = zeros(N,1);

% 路径B先算固定 T12 输出
for k = 2:N
    r_T12(k) = b_T12*r(k-1) - a_T12(2)*r_T12(k-1);
end

% 主循环：计算每个时刻的Qk和输出
for k = 2:N
    % 当前Qk滤波输出
    if k>2
        q_out(k) = theta1(k)*r(k) + theta2(k)*r(k-1);
    else
        q_out(k) = theta1(k)*r(k);
    end
    
    % 路径A：T12[Qk*r]
    yA(k) = b_T12*q_out(k-1) - a_T12(2)*yA(k-1);
    
    % 路径B：Qk[T12*r]
    if k>2
        yB(k) = theta1(k)*r_T12(k) + theta2(k)*r_T12(k-1);
    else
        yB(k) = theta1(k)*r_T12(k);
    end
end

% 差异
diff_y = yA - yB;

%% 绘图
figure;
subplot(3,1,1);
plot(t, yA, 'b', 'DisplayName','T12[Qk r]');
hold on; plot(t, yB, 'r--', 'DisplayName','Qk[T12 r]');
xlabel('时间 (s)'); ylabel('幅度');
title('T_{12}Q_k 与 Q_kT_{12} 的输出对比');
legend; grid on;

subplot(3,1,2);
plot(t, diff_y, 'k', 'DisplayName','差异 (y_A - y_B)');
xlabel('时间 (s)'); ylabel('差值');
title('算子不对易引起的差异');
legend; grid on;

subplot(3,1,3);
plot(t, theta1, 'b', t, theta2, 'r');
xlabel('时间 (s)');
ylabel('\theta_i(k)');
title('时变参数 \theta_1(k), \theta_2(k)');
legend('\theta_1','\theta_2'); grid on;

%% 验证实验：区分非对易性误差 vs 数值误差
fprintf('=== 非对易性 vs 数值误差验证 ===\n');

% 1. 对易性测试：固定参数情况
theta1_fixed = 0.5 * ones(N,1);  % 固定参数
theta2_fixed = 0.3 * ones(N,1);

yA_fixed = zeros(N,1);
yB_fixed = zeros(N,1);
q_out_fixed = zeros(N,1);

% 重新计算固定参数情况
for k = 2:N
    if k>2
        q_out_fixed(k) = theta1_fixed(k)*r(k) + theta2_fixed(k)*r(k-1);
    else
        q_out_fixed(k) = theta1_fixed(k)*r(k);
    end
    
    yA_fixed(k) = b_T12*q_out_fixed(k-1) - a_T12(2)*yA_fixed(k-1);
    
    if k>2
        yB_fixed(k) = theta1_fixed(k)*r_T12(k) + theta2_fixed(k)*r_T12(k-1);
    else
        yB_fixed(k) = theta1_fixed(k)*r_T12(k);
    end
end

diff_y_fixed = yA_fixed - yB_fixed;

% 2. 理论预测：非对易性误差应该与参数变化率相关
% 对于 T12*Qk - Qk*T12，理论上误差与 dθ/dt 成正比
dtheta1_dt = gradient(theta1) * fs;  % 数值微分
dtheta2_dt = gradient(theta2) * fs;
theoretical_error = abs(dtheta1_dt) + abs(dtheta2_dt);  % 简化的理论预测

% 3. 数值精度测试：使用更高精度计算
yA_double = double(yA);
yB_double = double(yB);
diff_y_double = yA_double - yB_double;

%% 统计分析
fprintf('固定参数情况下的差异 RMS: %.2e\n', rms(diff_y_fixed));
fprintf('时变参数情况下的差异 RMS: %.2e\n', rms(diff_y));
fprintf('差异比值 (时变/固定): %.1f\n', rms(diff_y)/rms(diff_y_fixed));

% 计算相关性
corr_theory_actual = corrcoef(theoretical_error(10:end), abs(diff_y(10:end)));
fprintf('理论预测与实际差异的相关系数: %.3f\n', corr_theory_actual(1,2));

%% 增强的绘图
figure('Position', [100 100 1200 800]);

subplot(2,3,1);
plot(t, yA, 'b', 'DisplayName','T12[Qk r]');
hold on; plot(t, yB, 'r--', 'DisplayName','Qk[T12 r]');
xlabel('时间 (s)'); ylabel('幅度');
title('时变参数：T_{12}Q_k 与 Q_kT_{12}');
legend; grid on;

subplot(2,3,2);
plot(t, yA_fixed, 'b', 'DisplayName','T12[Q r] (固定)');
hold on; plot(t, yB_fixed, 'r--', 'DisplayName','Q[T12 r] (固定)');
xlabel('时间 (s)'); ylabel('幅度');
title('固定参数：应该对易');
legend; grid on;

subplot(2,3,3);
semilogy(t, abs(diff_y), 'k', 'DisplayName','时变差异');
hold on; 
semilogy(t, abs(diff_y_fixed), 'g', 'DisplayName','固定差异(数值误差)');
xlabel('时间 (s)'); ylabel('|差异| (对数)');
title('差异幅度对比');
legend; grid on;

subplot(2,3,4);
plot(t, diff_y, 'k', 'DisplayName','实际差异');
hold on;
plot(t, 0.01*theoretical_error, 'r--', 'DisplayName','理论预测×0.01');
xlabel('时间 (s)'); ylabel('差值');
title('实际 vs 理论差异');
legend; grid on;

subplot(2,3,5);
plot(t, dtheta1_dt, 'b', t, dtheta2_dt, 'r');
xlabel('时间 (s)');
ylabel('d\theta/dt');
title('参数变化率');
legend('d\theta_1/dt','d\theta_2/dt'); grid on;

subplot(2,3,6);
scatter(theoretical_error(10:end), abs(diff_y(10:end)), 10, 'filled');
xlabel('理论误差指标');
ylabel('|实际差异|');
title(['相关性: r=' num2str(corr_theory_actual(1,2), '%.3f')]);
grid on;

%% 频域分析验证
fprintf('\n=== 频域验证 ===\n');
% 参数变化的主要频率
[Ptheta, f] = pwelch(theta1-mean(theta1), [], [], [], fs);
[~, max_idx] = max(Ptheta);
main_freq = f(max_idx);

% 差异信号的主要频率
[Pdiff, f_diff] = pwelch(diff_y, [], [], [], fs);
[~, max_idx_diff] = max(Pdiff);
main_freq_diff = f_diff(max_idx_diff);

fprintf('参数变化主频: %.1f Hz\n', main_freq);
fprintf('差异信号主频: %.1f Hz\n', main_freq_diff);

% 如果是真正的非对易性，差异频率应该与参数变化频率相关
