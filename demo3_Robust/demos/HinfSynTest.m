%% 1. 定义被控对象和权重函数
clear; clc; close all;

% 定义一个简单的二阶被控对象 G(s)
s = tf('s');
G = 1 / (s^2 + 0.2*s + 1);

% 定义性能权重 W1 (低通滤波器)
% 要求在 0.1 rad/s 以下有-40dB的扰动抑制，在 10 rad/s 以上性能要求放宽
M = 100; % 低频增益 (40 dB)
wb = 1;  % 带宽
A = 1e-4; % 高频增益
W1 = (s/M + wb) / (s + wb*A);
W1.Name = 'W1';

% 定义鲁棒性权重 W3 (高通滤波器)
% 从 100 rad/s 开始，要求T的增益以40dB/dec的速度下降
W3 = tf([1/1.5 10],[1 1000]); % 这是一个示例，实际需要根据模型不确定性来设计
W3.Name = 'W3';

% 在这个例子中，我们不限制控制量，所以 W2 为空
W2 = [];

%% 2. 使用 augw 构建广义对象
P = augw(G, W1, W2, W3);

%% 3. 使用 hinfsyn 综合控制器
NMEAS = 1; % 1个测量输出 y
NCON = 1;  % 1个控制输入 u

[K, CL, gamma] = hinfsyn(P, NMEAS, NCON);

fprintf('H-infinity 综合完成，最终的 gamma = %.4f\n', gamma);

if gamma >= 1
    disp('警告: 设计未满足所有约束，请考虑放宽权重函数。');
else
    disp('成功: 控制器满足所有性能和鲁棒性约束。');
end

%% 4. 验证设计结果
L = G * K; % 开环传递函数
S = feedback(1, L); % 闭环灵敏度函数
T = feedback(L, 1); % 闭环互补灵敏度函数

% 计算未加控制器时的灵敏度函数
S0 = feedback(1, G); % 原始灵敏度
T0 = feedback(G, 1); % 原始互补灵敏度

figure;
% 验证性能: sigma(S) < 1/sigma(W1)
bodemag(S0, S, 1/W1, {1e-2, 1e3});
grid on;
title('灵敏度函数 S(s) vs 性能权重 1/W1(s)');
legend('原始 S0(s)', '闭环 S(s)', '性能边界 1/W1(s)', 'Location', 'best');

figure;
% 验证鲁棒性: sigma(T) < 1/sigma(W3)
bodemag(T0, T, 1/W3, {1e-2, 1e3});
grid on;
title('互补灵敏度函数 T(s) vs 鲁棒性权重 1/W3(s)');
legend('原始 T0(s)', '闭环 T(s)', '鲁棒性边界 1/W3(s)', 'Location', 'best');