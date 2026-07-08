clear; close all; clc;
% 系统参数与初始条件
a1 = 0.5;       % 未知参数，真实值
b2 = 0.5;       % 未知参数，真实值
beta = 0.5;     % 未知参数，真实值
omega1 = 1;     % 扰动频率1
omega2 = 0.5;   % 扰动频率2
omega3 = 4;     % 扰动频率3（仅用于t < 50s）

% 控制器参数
k = 7.1;        % 控制增益
g = 100;        % 参数更新增益
d = [2, 3, 2, 1]; % 矩阵D的参数 [d2, d3, d4, d5]

% 初始条件
x0 = [0.5; 0];          % 初始状态 [x1, x2]
xi1_0 = zeros(4,1);     % 滤波器ξ1初始状态
xi2_0 = zeros(4,1);     % 滤波器ξ2初始状态
chi_hat0 = zeros(4,1);  % 估计器χ初始状态
theta_hat0 = [0.1; 0.1];% 参数估计初始值 [θ1_hat, θ2_hat]

% 仿真时间设置
tspan = [0, 150];       % 仿真时间范围：0到150秒

% 定义赫尔维茨矩阵D（确保特征值在左半平面）
D = [-1, 1, 0, 0;
      0, -d(1), 1, 0;
      0, 0, -d(2), 1;
      0, 0, 0, -d(3)]; % 示例结构，可根据实际调整

% 定义ODE求解选项
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% 运行仿真
[t, X] = ode45(@(t, X) system_dynamics(t, X, a1, b2, beta, omega1, omega2, omega3, k, g, D, d), tspan, ...
              [x0; xi1_0; xi2_0; chi_hat0; theta_hat0], options);

% 提取结果
x1 = X(:,1);            % 系统状态x1
x2 = X(:,2);            % 系统状态x2
e = -X(:,1) - get_w2(t, omega2); % 跟踪误差 e = y - y_r = x1 - (-w2)
u = -k*e - X(:,9) - X(:,5).*X(:,end-1) - X(:,6).*X(:,end); % 控制输入u

% 绘图
figure;
subplot(2,1,1);
plot(t, x1, 'b', t, -get_w2(t, omega2), 'r--');
legend('输出 y=x1', '参考信号 y_r=-w2');
title('输出跟踪性能');
xlabel('时间 (s)');

subplot(2,1,2);
plot(t, e);
title('跟踪误差 e');
xlabel('时间 (s)');

figure;
plot(t, X(:,end-1), 'b', t, X(:,end), 'r');
legend('参数估计 θ1_hat', '参数估计 θ2_hat');
title('参数估计收敛情况');
xlabel('时间 (s)');

% 定义系统动态方程
function dXdt = system_dynamics(t, X, a1, b2, beta, omega1, omega2, omega3, k, g, D, d)
    % 分解状态向量
    x = X(1:2);                 % 系统状态 [x1; x2]
    xi1 = X(3:6);               % 滤波器ξ1状态
    xi2 = X(7:10);              % 滤波器ξ2状态
    chi_hat = X(11:14);         % 估计器χ状态
    theta_hat = X(15:16);       % 参数估计 [θ1_hat; θ2_hat]
    
    % 计算扰动w1和w2
    w1 = get_w1(t, omega1, omega3);
    w2 = get_w2(t, omega2);
    
    % 计算跟踪误差e = y - y_r = x1 - (-w2)
    e = x(1) - (-w2);
    
    % 计算控制输入u
    mu1 = xi1(1);               % μ1 = ξ11
    mu2 = xi2(1);               % μ2 = ξ21
    u = -k*e - chi_hat(1) - mu1*theta_hat(1) - mu2*theta_hat(2);
    
    % 系统动态方程
    dx1dt = x(2) + a1*x(1) + (1/beta)*u + w1;
    dx2dt = (1/beta)*b2*u;
    
    % 滤波器动态方程
    dxi1dt = D*xi1 + [0; u; 0; 0];
    dxi2dt = D*xi2 + [0; 0; 0; u];
    
    % 估计器χ动态方程
    dchi_hatdt = D*chi_hat - [d(1); d(2); d(3); d(4)]*u;
    
    % 参数更新律（投影算子简化实现）
    theta1_hat_dot = g * projection(mu1*e, theta_hat(1), 0, 2); % 假设θ1' ∈ [0,2]
    theta2_hat_dot = g * projection(mu2*e, theta_hat(2), 0, 2); % 假设θ2' ∈ [0,2]
    
    % 整合所有导数
    dXdt = [dx1dt; dx2dt; dxi1dt; dxi2dt; dchi_hatdt; theta1_hat_dot; theta2_hat_dot];
end

% 定义扰动函数w1（分段函数）
function w1 = get_w1(t, omega1, omega3)
    if t < 50
        w1 = sin(omega1*t) + 0.2*sin(omega3*t);
    else
        w1 = sin(omega1*t);
    end
end

% 定义扰动函数w2（分段函数）
function w2 = get_w2(t, omega2)
    if t < 100
        w2 = sin(omega2*t);
    else
        w2 = 0;
    end
end

% 简化投影算子（确保参数在[low, high]范围内）
function theta_dot = projection(update, theta_hat, low, high)
    if theta_hat <= low && update < 0
        theta_dot = 0;
    elseif theta_hat >= high && update > 0
        theta_dot = 0;
    else
        theta_dot = update;
    end
end