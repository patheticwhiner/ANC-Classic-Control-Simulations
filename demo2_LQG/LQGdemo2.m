% 使用 LQG 函数改写程序，分开原系统和干扰模型的动态更新
clear; close all; clc;

%% 加载数据
load('..\dataset\bltdWhiteNoise_ssmodel.mat');  % 包含干扰模型 Aw, Bw, Cw
load('..\dataset\bpf_ssmodel.mat');  % 包含原系统 Af, Bf, Cf

%% 耦合系统
% 原系统维度
n = size(Af, 1);    % 原系统状态维度
p = size(Aw, 1);  % 干扰模型状态维度

% 原系统状态空间模型
A_sys = Af;
B_sys = Bf;
C_sys = Cf;

% 干扰模型状态空间模型
A_dist = Aw;
B_dist = Bw;
C_dist = Cw;

% 增广系统状态空间模型
A = blkdiag(A_sys, A_dist);  % 块对角矩阵 [Af 0; 0 Aw]
B = [B_sys; zeros(p, size(B_sys, 2))];  % 控制输入矩阵
G = [zeros(n, size(B_dist, 2)); B_dist];  % 干扰输入矩阵
C = [C_sys, C_dist];  % 增广输出矩阵

%% 设计 LQG 控制器
% 状态权重矩阵 Q（半正定）
Q = C' * C;
% 输入权重矩阵 R（正定）
R = 10e-4;

% 过程噪声协方差矩阵和测量噪声协方差矩阵
Qn = 2 * eye(size(G, 1));  % 过程噪声协方差
Rn = 1e-4 * eye(size(C, 1));  % 测量噪声协方差

% 创建增广系统的状态空间模型
sys = ss(A, [B G], C, 0, 1);  % 离散时间系统，输入包括控制输入和干扰输入

% 使用 LQG 函数设计控制器
lqg_controller = lqg(sys, blkdiag(Q,R), blkdiag(Qn,Rn));  % 设计 LQG 控制器

% 提取 LQG 控制器的增益矩阵
K = lqg_controller.K;  % 状态反馈增益
L = lqg_controller.L;  % 卡尔曼滤波器增益

% 显示 LQG 控制器增益
disp('LQG 控制器状态反馈增益矩阵 K:');
disp(K);
disp('LQG 控制器卡尔曼滤波器增益矩阵 L:');
disp(L);

%% 初始化参数
N = 2000;                    % 仿真步数
x_sys = zeros(n, 1);          % 原系统状态初始化
x_dist = zeros(p, 1);         % 干扰模型状态初始化
x_hat = zeros(n + p, 1);      % 增广状态估计初始化
x_sys_history = zeros(n, N);  % 原系统状态轨迹记录
x_dist_history = zeros(p, N); % 干扰模型状态轨迹记录
x_hat_history = zeros(n + p, N);  % 增广状态估计轨迹记录
y_history = zeros(size(C_sys, 1), N);  % 输出轨迹记录
u_history = zeros(1, N);  % 控制输入记录

%% 闭环控制仿真
for k = 1:N
    % 生成干扰模型的过程噪声
    e_dist = sqrtm(Qn) * randn(size(B_dist, 2), 1);  % 干扰过程噪声
    v = sqrtm(Rn) * randn(size(C, 1), 1);  % 测量噪声

    % 干扰模型动态更新
    x_dist = A_dist * x_dist + B_dist * e_dist;  % 干扰模型状态更新
    d = C_dist * x_dist;  % 干扰信号

    % 原系统动态更新
    u = -K * x_hat;  % 使用估计状态进行反馈控制
    x_sys = A_sys * x_sys + B_sys * u;  % 原系统状态更新
    y = C_sys * x_sys + d + v;  % 系统输出（包含干扰和测量噪声）

    % 卡尔曼滤波器更新（状态估计）
    x_hat = A * x_hat + [B G] * [u; e_dist];  % 预测
    x_hat = x_hat + L * (y - C * x_hat);  % 校正

    % 记录状态和输出
    x_sys_history(:, k) = x_sys;
    x_dist_history(:, k) = x_dist;
    x_hat_history(:, k) = x_hat;
    y_history(:, k) = y;
    u_history(k) = u;
end

% 计算干扰信号和反噪声信号
d_history = C_dist * x_dist_history;  % 干扰信号
anti_history = C_sys * x_sys_history;  % 反噪声信号

%% 绘图
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