clear; clc; close all;

%% 系统选择 (1-4选择不同示例)
example_choice = 2;  % 修改这个数字来选择不同示例

%% 系统参数配置
switch example_choice
    case 1  % 三阶不稳定最小相位系统 P(s) = 10/[(s+10)(s+1)(s-1)]
        num = 10;
        den = [1 10 -1 -10];
        r = length(den) - length(num);
        y0 = 0.1;
        theta_hat0 = 4.2;
        sys_name = '三阶不稳定最小相位系统';
        
    case 2  % 二阶不稳定最小相位系统 P(s) = 10(s+2)/[(s+10)(s+1)(s-1)]
        num = 10 * [1 2];
        den = [1 10 -1 -10];
        r = length(den) - length(num);
        y0 = -0.1;
        theta_hat0 = 3.8;
        sys_name = '二阶不稳定最小相位系统';
        
    case 3  % 一阶不稳定最小相位系统 P(s) = 10(s+2)²/[(s+10)(s+1)(s-1)]
        num = 10 * [1 4 4];
        den = [1 10 -1 -10];
        r = length(den) - length(num);
        y0 = -0.1;
        theta_hat0 = 4.2;
        sys_name = '一阶不稳定最小相位系统';
        
    case 4  % 含未建模扰动的三阶系统 P(s) = 10/[(s+10)(s+1)(s-1)]
        num = 10;
        den = [1 10 -1 -10];
        r = length(den) - length(num);
        y0 = 0.1;
        theta_hat0 = 4.2;
        sys_name = '含未建模扰动的三阶系统';
end

%% 通用控制器参数
z = [1, 2];              % 控制器极点
zeta = 0.7;              % 阻尼比
omega_n = 1.9;           % 自然频率
epsilon = 0.5;           % 自适应增益
rho_bar = 3;             % 相对阶上界
rho_under = 1;           % 相对阶下界
a = 0.5;                 % 控制器参数a
w_hat0 = [0; 0; 0; 0];   % 初始化估计状态变量
time_span = [0, 100];     % 仿真时间范围

%% 系统分析
[A, B, C, D] = tf2ss(num, den);

% 验证系统可控性和可观测性
if rank(ctrb(A, B)) ~= size(A, 1)
    warning('系统不完全可控!');
end
if rank(obsv(A, C)) ~= size(A, 1)
    warning('系统不完全可观测!');
end

fprintf('=== %s ===\n', sys_name);
fprintf('系统特征值: ');
fprintf('%.3f ', eig(A));
fprintf('\n');

%% 仿真设置和执行
initial_state = [y0; w_hat0; theta_hat0];

% 添加进度显示的输出函数
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);

fprintf('开始仿真...\n');
tic; % 开始计时
[t, x] = ode45(@(t, x) system_dynamics(t, x, A, B, C, z, zeta, omega_n, epsilon, rho_bar, rho_under, a, example_choice), ...
               time_span, initial_state, options);
elapsed_time = toc; % 结束计时
fprintf('\n仿真完成！总耗时: %.2f 秒\n', elapsed_time);

%% 结果提取
y = x(:, 1);
w_hat = x(:, 2:5);
theta_hat = x(:, 6);

% 计算控制输入历史
u_history = zeros(size(t));
for i = 1:length(t)
    % 简化的控制输入计算
    K = rho_bar - rho_under + 1;  % 控制增益
    v = -K * y(i);                % 比例控制
    compensation = a * (w_hat(i,1) + w_hat(i,3));  % 补偿项
    u_history(i) = v + compensation;
end

%% 绘图
figure('Position', [100, 100, 1000, 800]);

subplot(2, 3, 1);
plot(t, y, 'LineWidth', 2);
title('系统输出 y(t)', 'FontSize', 12);
xlabel('时间 (s)'); ylabel('输出 y'); grid on;

subplot(2, 3, 2);
plot(t, theta_hat, 'LineWidth', 2);
title('频率参数估计 \theta(t)', 'FontSize', 12);
xlabel('时间 (s)'); ylabel('\theta 估计值'); grid on;

subplot(2, 3, 3);
plot(t, w_hat(:,1), 'LineWidth', 2);
title('估计状态变量 w_1(t)', 'FontSize', 12);
xlabel('时间 (s)'); ylabel('w_1'); grid on;

subplot(2, 3, 4);
plot(t, w_hat(:,2), 'LineWidth', 2);
title('估计状态变量 w_2(t)', 'FontSize', 12);
xlabel('时间 (s)'); ylabel('w_2'); grid on;

subplot(2, 3, 5);
plot(t, w_hat(:,3), 'LineWidth', 2);
title('估计状态变量 w_3(t)', 'FontSize', 12);
xlabel('时间 (s)'); ylabel('w_3'); grid on;

subplot(2, 3, 6);
plot(t, u_history, 'LineWidth', 2);
title('控制输入 u(t)', 'FontSize', 12);
xlabel('时间 (s)'); ylabel('控制输入 u'); grid on;

sgtitle(sprintf('%s - 自适应控制仿真结果', sys_name), 'FontSize', 14, 'FontWeight', 'bold');

%% 系统动力学函数
function dxdt = system_dynamics(t, x, A, B, C, z, zeta, omega_n, epsilon, rho_bar, rho_under, a, example_choice)
    % 状态变量
    y = x(1);
    w_hat = x(2:5);
    theta_hat = x(6);
    
    % 外部扰动 (仅示例4有扰动)
    d = sin(2*t);
    if example_choice == 4
        d = sin(2*t) + 0.1 * sin(t+0.2);
    end
    
    % 计算控制律
    K = rho_bar - rho_under + 1;  % 控制增益
    v = -K * y;                   % 比例控制
    
    % 计算控制输入
    compensation = a * (w_hat(1) + w_hat(3));  % 基于估计状态的补偿
    u = v + compensation;
    
    % 更新频率参数
    error_signal = w_hat(1) * w_hat(2) + w_hat(3) * w_hat(4);
    theta_dot = -epsilon * error_signal * v;
    
    % 系统状态方程
    disturbance_vector = [0; 0; d];
    if size(A, 1) < 3
        disturbance_vector = disturbance_vector(1:size(A, 1));
    end
    
    system_state = A \ (B * u + disturbance_vector);
    y_dot = C * system_state;
    
    % w_hat 状态更新
    w_hat_dot = [w_hat(2);
                 -theta_hat * w_hat(1) - 2 * zeta * omega_n * w_hat(2) + (omega_n^2) * u;
                 w_hat(4);
                 -theta_hat * w_hat(3) - 2 * zeta * omega_n * w_hat(4) + (omega_n^2) * u];
    
    dxdt = [y_dot; w_hat_dot; theta_dot];
end

