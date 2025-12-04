% 需要检查状态空间表示中是否含有延迟
clear; close all; clc;
%%
load('..\dataset\bltdWhiteNoise_ssmodel.mat');
load('..\dataset\bpf_ssmodel.mat');

%% 耦合系统
% 增广状态矩阵
n = size(Af, 1);    % 原系统状态维度
p = size(Aw, 1);  % 干扰模型状态维度
A = blkdiag(Af, Aw);  % 块对角矩阵 [Af 0; 0 Aw]

% 增广输入矩阵
B = [Bf; zeros(p, size(Bf, 2))];  % 原系统控制输入通道
G = [zeros(n, size(Bw, 2)); Bw];  % 干扰输入通道

% 增广输出矩阵
C = [Cf , Cw];

%% 设计LQR控制器
% 状态权重矩阵 Q（半正定）
Q = C' * C;
% 输入权重矩阵 R（正定）
R = 10e-4;
% 离散 LQR 求解（状态反馈增益）
[K, ~, ~] = dlqr(A, B, Q, R);

% 计算闭环系统的状态矩阵
A_cl = A - B * K;
% 计算闭环系统的特征值
eig_cl = eig(A_cl);
% 检查闭环系统的稳定性
if all(abs(eig_cl) < 1)
    disp('闭环系统是稳定的（所有极点都在单位圆内）。');
else
    disp('闭环系统是不稳定的（存在极点在单位圆外）。');
end

%% 设计卡尔曼滤波器（状态估计器）
% 假设过程噪声协方差矩阵和测量噪声协方差矩阵
Qn = 2000;  % 过程噪声协方差
Rn = 1e-4 * eye(size(C, 1));  % 测量噪声协方差

% 使用离散卡尔曼滤波器求解最优增益矩阵
[L, ~, ~] = dlqe(A, G, C, Qn, Rn);

% 显示卡尔曼滤波器增益
disp('卡尔曼滤波器增益矩阵 L:');
disp(L);

% 计算引入卡尔曼滤波器后的闭环系统状态矩阵
A_cl = A - B * K - L * C;
% 计算闭环系统的特征值
eig_cl = eig(A_cl);
% 显示闭环系统的极点
disp('引入卡尔曼滤波器后的闭环系统极点:');
disp(eig_cl);
% 检查闭环系统的稳定性
if all(abs(eig_cl) < 1)
    disp('引入卡尔曼滤波器后的闭环系统是稳定的（所有极点都在单位圆内）。');
else
    disp('引入卡尔曼滤波器后的闭环系统是不稳定的（存在极点在单位圆外）。');
end

%% 初始化参数
N = 20000;                    % 仿真步数
x = zeros(n + p, 1);    % 增广状态初始化
x_hat = zeros(n + p, 1);    % 状态估计初始化
x_history = zeros(n + p, N);  % 状态轨迹记录
x_hat_history = zeros(n + p, N);  % 状态估计轨迹记录
y_history = zeros(size(Cf, 1), N);  % 输出轨迹记录
u_history = zeros(1, N);  % 控制输入记录

%% 闭环控制仿真(一)
for k = 1:N
    % 生成过程噪声和测量噪声
    e = sqrt(Qn) * randn(size(G, 2), 1);  % 过程噪声
    v = chol(Rn)' * randn(size(C, 1), 1);  % 测量噪声

    % 状态反馈控制律
    u = -K * x_hat;  % 使用估计状态进行反馈控制
    u_history(k) = u;  % 记录控制输入

    % 增广系统动态更新
    x = A * x + B * u + G * e;  % 系统状态更新
    y = C * x + v;  % 输出测量

    % 卡尔曼滤波器更新（状态估计）
    x_hat = A * x_hat + B * u;  % 预测
    x_hat = x_hat + L * (y - C * x_hat);  % 校正

    % 记录状态和输出
    x_history(:, k) = x;
    x_hat_history(:, k) = x_hat;
    y_history(:, k) = y;
end

% 计算干扰信号和反噪声信号
d_history = Cw * x_history(n+1:end, :);  % 干扰信号
anti_history = Cf * x_history(1:end-p, :);  % 反噪声信号

figure;
subplot(2, 1, 1);
plot(1:N, u_history, 'DisplayName', '控制信号'); hold on;
plot(1:N, anti_history, 'DisplayName', '反噪声信号');
xlabel('样本n/samples'); ylabel('幅度');
legend;
grid on; 

subplot(2, 1, 2);
plot(1:N, d_history, 'DisplayName', '干扰信号d'); hold on;
plot(1:N, -anti_history, 'DisplayName', '反噪声信号');
plot(1:N, y_history, 'DisplayName', '残余信号 y');
xlabel('样本n/samples'); ylabel('幅度');
legend;
grid on;

%% 计算声压级及降噪量
p0 = 20e-6;  % 参考声压 (Pa)
% 声压级计算
SPL_d = 20 * log10(rms(d_history(:)) / p0);
SPL_y = 20 * log10(rms(y_history(:)) / p0);
% 降噪量
attenuation = SPL_d - SPL_y;
fprintf('干扰信号声压级: %.2f dB\n', SPL_d);
fprintf('残余信号声压级: %.2f dB\n', SPL_y);
fprintf('降噪量: %.2f dB\n', attenuation);

%% 闭环控制仿真（二）
N = 2000;                    % 仿真步数
xf = zeros(n, 1);          % 原系统状态初始化
xw = zeros(p, 1);         % 干扰模型状态初始化
x_hat = zeros(n + p, 1);      % 增广状态估计初始化
xf_history = zeros(n, N);  % 原系统状态轨迹记录
xw_history = zeros(p, N); % 干扰模型状态轨迹记录
x_hat_history = zeros(n + p, N);  % 增广状态估计轨迹记录
y_history = zeros(size(Cf, 1), N);  % 输出轨迹记录
u_history = zeros(1, N);  % 控制输入记录

for k = 1:N
    % 生成干扰模型的过程噪声
    ew = sqrtm(Qn) * randn(size(Bw, 2), 1);  % 干扰过程噪声
    v = sqrtm(Rn) * randn(size(C, 1), 1);  % 测量噪声

    % 干扰模型动态更新
    xw = Aw * xw + Bw * ew;  % 干扰模型状态更新
    d = Cw * xw;  % 干扰信号

    % 原系统动态更新
    u = -K * x_hat;  % 使用估计状态进行反馈控制
    xf = Af * xf + Bf * u;  % 原系统状态更新
    y = Cf * xf + d + v;  % 系统输出（包含干扰和测量噪声）

    % 卡尔曼滤波器更新（状态估计）
    x_hat = A * x_hat + B * u;  % 预测
    x_hat = x_hat + L * (y - C * x_hat);  % 校正

    % 记录状态和输出
    xf_history(:, k) = xf;
    xw_history(:, k) = xw;
    x_hat_history(:, k) = x_hat;
    y_history(:, k) = y;
    u_history(k) = u;
end

% 计算干扰信号和反噪声信号
d_history = Cw * xw_history;  % 干扰信号
anti_history = Cf * xf_history;  % 反噪声信号

figure;
subplot(2, 1, 1);
plot(1:N, u_history, 'DisplayName', '控制信号'); hold on;
plot(1:N, anti_history, 'DisplayName', '反噪声信号');
xlabel('样本n/samples'); ylabel('幅度');
legend;
grid on; 

subplot(2, 1, 2);
plot(1:N, d_history, 'DisplayName', '干扰信号d'); hold on;
plot(1:N, -anti_history, 'DisplayName', '反噪声信号');
plot(1:N, y_history, 'DisplayName', '残余信号 y');
xlabel('样本n/samples'); ylabel('幅度');
legend;
grid on;

%% 计算声压级及降噪量
p0 = 20e-6;  % 参考声压 (Pa)
% 声压级计算
SPL_d = 20 * log10(rms(d_history(:)) / p0);
SPL_y = 20 * log10(rms(y_history(:)) / p0);
% 降噪量
attenuation = SPL_d - SPL_y;
fprintf('干扰信号声压级: %.2f dB\n', SPL_d);
fprintf('残余信号声压级: %.2f dB\n', SPL_y);
fprintf('降噪量: %.2f dB\n', attenuation);