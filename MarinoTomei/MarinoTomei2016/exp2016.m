clear; close all; clc;
%%
% 定义系统参数和初始条件
k = 0.5;
epsilon = 0.05;
theta1_0 = 0.3;
theta2_0 = 4.2;

% 定义时间范围
t_span = [0 200];
% 初始条件: [y, w0, w1, w2, w3, w4, theta1, theta2, y_dot]
% 所有初始条件为零，除了theta1和theta2
initial_conditions = [0; 0; 0; 0; 0; 0; theta1_0; theta2_0; 0];

% 使用ode45求解微分方程
[t, y] = ode45(@(t, y) system_dynamics(t, y, k, epsilon), t_span, initial_conditions);

% 提取结果
y_output = y(:,1); % 系统输出y
w0 = y(:,2); % 估计值w0
w1 = y(:,3); % 估计值w1
w2 = y(:,4); % 估计值w2
w3 = y(:,5); % 估计值w3
w4 = y(:,6); % 估计值w4
u_output = -w0 - w1 - w3; % 控制输入u = -w0 - w1 - w3
theta1_est = y(:,7); % 频率估计θ1
theta2_est = y(:,8); % 频率估计θ2

% 计算实际扰动信号用于比较
d_actual = 0.3 + sin(0.5*t) + sin(2*t);

% 绘制结果
figure;
subplot(2,2,1);
plot(t, y_output);
title('y'); axis([0 200 -0.5 1]);
xlabel('时间 (s)');
grid on;

subplot(2,2,2);
plot(t, u_output);
title('u'); axis([0 200 -4 2]);
xlabel('时间 (s)');
grid on;

subplot(2,2,3);
plot(t, theta1_est);
title("$\hat{\theta}_1$", 'Interpreter', 'latex');
xlabel('时间 (s)');
yline(0.5^2, 'r--'); % 添加真实值参考线
legend('估计值', '真实值 (0.25)');
grid on;  axis([0 200 0.2 0.35]);

subplot(2,2,4);
plot(t, theta2_est);
title("$\hat{\theta}_2$", 'Interpreter', 'latex');
xlabel('时间 (s)');
yline(2^2, 'r--'); % 添加真实值参考线
legend('估计值', '真实值 (4)');
grid on;  axis([0 200 3.9 4.2]);

% 添加额外图形以显示扰动和估计
figure;
subplot(2,1,1);
plot(t, d_actual, 'k', t, w0 + w1 + w3, 'r--');
title('扰动信号和估计');
xlabel('时间 (s)');
legend('实际扰动 d(t)', '估计扰动');
grid on;

subplot(2,1,2);
plot(t, y_output + u_output);
title('y+u (剩余误差)');
xlabel('时间 (s)');
grid on;

%% 定义系统动态
function dydt = system_dynamics(t, y, k, epsilon)
    % 状态变量:
    % y(1): 系统输出 y
    % y(2): w0 - 常数项估计
    % y(3): w1 - 第一个正弦项估计的幅值相关量
    % y(4): w2 - 第一个正弦项估计的相位相关量
    % y(5): w3 - 第二个正弦项估计的幅值相关量
    % y(6): w4 - 第二个正弦项估计的相位相关量
    % y(7): theta1 - 第一个频率参数估计 (应该收敛到0.5^2=0.25)
    % y(8): theta2 - 第二个频率参数估计 (应该收敛到2^2=4)
    % y(9): y的导数 (内部状态)
    
    % 扰动信号
    d = 0.3 + sin(0.5*t) + sin(2*t);
    
    % 控制输入
    u = -y(2) - y(3) - y(5);
    
    % 二阶系统动态 y(s) = 1/(s+1)^2 * [u(s) + d(s)]
    % 状态空间表示
    x1 = y(1);             % 系统输出 y
    x1_dot = y(9);         % y的导数 (内部状态)
    x1_ddot = -2*x1_dot - x1 + u + d;  % 二阶系统方程
    
    % 根据控制器方程计算各变量的导数 (公式34)
    dw0 = k * x1;               % 对应 ẇ0 = ky
    dw1 = y(4) + k * x1;        % 对应 ẇ1 = w2 + ky
    dw2 = -y(7) * y(3);         % 对应 ẇ2 = -θ1·w1
    dw3 = y(6) - k * x1;        % 对应 ẇ3 = w4 - ky
    dw4 = -y(8) * y(5);         % 对应 ẇ4 = -θ2·w3
    dtheta1 = epsilon * y(4) * x1; % 对应 θ̇1 = ε·w2·y
    dtheta2 = -epsilon * y(6) * x1; % 对应 θ̇2 = -ε·w4·y
    
    % 组合所有导数
    dydt = [x1_dot; dw0; dw1; dw2; dw3; dw4; dtheta1; dtheta2; x1_ddot];
end