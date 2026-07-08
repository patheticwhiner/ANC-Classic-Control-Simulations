function MarinoTomei_2016_adaptive_freq_est()
%% Marino & Tomei (2016) — 自适应频率估计与扰动抵消
% 基于 Marino, Tomei. "Adaptive disturbance rejection for unknown stable
% linear systems", Transactions of the Institute of Measurement and Control,
% 2016, 38(6): 640-647.
%
% 本脚本合并自:
%   MarinoTomei2016/freq_estimator_core.m  (纯估计器, 无对象动态)
%   MarinoTomei2016/run_freq_estimator.m    (完整闭环, 含二阶对象)
%
% 提供两种运行模式:
%   mode = 'plant'     : 估计器 + 二阶被控对象闭环 (默认, 完整仿真)
%   mode = 'estimator' : 纯估计器, 直接估计扰动信号

close all; clc;

%% ====== 用户设置 ======
run_mode = 'plant';        % 'plant' | 'estimator'

% --- 控制器参数 ---
k       = 0.5;             % 误差反馈增益
epsilon = 0.05;            % 频率自适应增益
theta1_0 = 0.3;            % theta1 初始估计 (真值: omega1^2 = 0.5^2 = 0.25)
theta2_0 = 4.2;            % theta2 初始估计 (真值: omega2^2 = 2^2  = 4.00)

% --- 仿真参数 ---
t_span = [0, 200];         % 仿真时间范围

% --- 扰动定义 ---
% d(t) = A0 + sin(omega1*t) + sin(omega2*t)
A0   = 0.3;                % 直流偏置
om1  = 0.5;                % 频率 1 [rad/s] → theta1_true = om1^2 = 0.25
om2  = 2.0;                % 频率 2 [rad/s] → theta2_true = om2^2 = 4.00

%% ====== 根据模式运行 ======
switch run_mode
    case 'plant'
        run_with_plant();
    case 'estimator'
        run_estimator_only();
    otherwise
        error('Unknown mode: %s. Use ''plant'' or ''estimator''.', run_mode);
end

%% ========================================================================
%  模式 1: 估计器 + 二阶被控对象闭环
%  对象: P(s) = 1 / (s+1)^2
%  y = P(s)*[u + d],  控制器 u = -w0_hat - w1_hat - w3_hat
% ========================================================================
function run_with_plant()
    % --- 初始条件 ---
    % y(1)=y, y(2)=w0, y(3)=w1, y(4)=w2, y(5)=w3, y(6)=w4,
    % y(7)=theta1, y(8)=theta2, y(9)=dy/dt
    y0 = [0; 0; 0; 0; 0; 0; theta1_0; theta2_0; 0];

    % --- 求解 ---
    [t, y] = ode45(@plant_dynamics, t_span, y0);

    % --- 提取结果 ---
    y_out   = y(:,1);       % 系统输出
    w0_est  = y(:,2);       % w0_hat (偏置估计)
    w1_est  = y(:,3);       % w1_hat (第1正弦分量)
    w3_est  = y(:,5);       % w3_hat (第2正弦分量)
    u_ctrl  = -w0_est - w1_est - w3_est;  % 控制输入
    th1_est = y(:,7);       % theta1_hat
    th2_est = y(:,8);       % theta2_hat

    % 实际扰动
    d_actual = A0 + sin(om1*t) + sin(om2*t);

    %% --- 图1: 闭环响应 ---
    figure('Name', '闭环自适应扰动抵消', 'Position', [100, 100, 900, 600]);

    subplot(2,2,1);
    plot(t, y_out, 'b', 'LineWidth', 1.2);
    title('系统输出 y(t)');
    xlabel('时间 (s)'); ylabel('y');
    ylim([-0.5, 1]); grid on;

    subplot(2,2,2);
    plot(t, u_ctrl, 'r', 'LineWidth', 1.2);
    title('控制输入 u(t)');
    xlabel('时间 (s)'); ylabel('u');
    ylim([-4, 2]); grid on;

    subplot(2,2,3);
    plot(t, th1_est, 'b', 'LineWidth', 1.2); hold on;
    yline(om1^2, 'r--', 'LineWidth', 1.2);
    title('\theta_1 估计', 'Interpreter', 'tex');
    xlabel('时间 (s)');
    legend('估计值', sprintf('真值 (%.2f)', om1^2), 'Location', 'best');
    ylim([0.2, 0.35]); grid on;

    subplot(2,2,4);
    plot(t, th2_est, 'b', 'LineWidth', 1.2); hold on;
    yline(om2^2, 'r--', 'LineWidth', 1.2);
    title('\theta_2 估计', 'Interpreter', 'tex');
    xlabel('时间 (s)');
    legend('估计值', sprintf('真值 (%.0f)', om2^2), 'Location', 'best');
    ylim([3.9, 4.2]); grid on;

    %% --- 图2: 扰动估计与残差 ---
    figure('Name', '扰动估计质量', 'Position', [150, 150, 900, 450]);

    d_est = w0_est + w1_est + w3_est;

    subplot(2,1,1);
    plot(t, d_actual, 'k', 'LineWidth', 1.5); hold on;
    plot(t, d_est, 'r--', 'LineWidth', 1.2);
    title('扰动信号与估计');
    xlabel('时间 (s)');
    legend('实际扰动 d(t)', '估计扰动 d_{hat}(t)');
    grid on;

    subplot(2,1,2);
    plot(t, d_actual - d_est, 'b', 'LineWidth', 1.2);
    title('估计误差 d(t) - d_{hat}(t)');
    xlabel('时间 (s)');
    grid on;

    fprintf('=== 闭环仿真完成 ===\n');
    fprintf('theta1 终值: %.4f (真值: %.2f)\n', th1_est(end), om1^2);
    fprintf('theta2 终值: %.4f (真值: %.2f)\n', th2_est(end), om2^2);
    fprintf('输出 y 稳态范围: [%.4f, %.4f]\n', ...
            min(y_out(t > 100)), max(y_out(t > 100)));
end

%% ========================================================================
%  模式 2: 纯估计器 (无被控对象, 直接估计扰动 d(t))
% ========================================================================
function run_estimator_only()
    % --- 初始条件 ---
    % [w0_hat, w1_hat, w2_hat, w3_hat, w4_hat, theta1_hat, theta2_hat]
    y0 = [0; 0; 0; 0; 0; theta1_0; theta2_0];

    % --- 求解 ---
    [t, y] = ode45(@estimator_dynamics, t_span, y0);

    % --- 提取 ---
    w0_hat  = y(:,1);
    w1_hat  = y(:,2);
    w2_hat  = y(:,3);
    w3_hat  = y(:,4);
    w4_hat  = y(:,5);
    th1_hat = y(:,6);
    th2_hat = y(:,7);

    d_actual = A0 + sin(om1*t) + sin(om2*t);
    d_est    = w0_hat + w1_hat + w3_hat;

    %% --- 绘图 ---
    figure('Name', '纯估计器响应', 'Position', [100, 100, 900, 600]);

    subplot(2,2,1);
    plot(t, w0_hat, 'b', 'LineWidth', 1.2);
    title('估计偏置项 w_0');
    xlabel('时间 (s)'); grid on;

    subplot(2,2,2);
    plot(t, w1_hat, 'b', 'LineWidth', 1.2);
    title('估计正弦项 w_1');
    xlabel('时间 (s)'); grid on;

    subplot(2,2,3);
    plot(t, th1_hat, 'b', 'LineWidth', 1.2); hold on;
    yline(om1^2, 'r--', 'LineWidth', 1.2);
    title('\theta_1 估计', 'Interpreter', 'tex');
    xlabel('时间 (s)');
    legend('估计值', sprintf('真值 (%.2f)', om1^2), 'Location', 'best');
    grid on;

    subplot(2,2,4);
    plot(t, th2_hat, 'b', 'LineWidth', 1.2); hold on;
    yline(om2^2, 'r--', 'LineWidth', 1.2);
    title('\theta_2 估计', 'Interpreter', 'tex');
    xlabel('时间 (s)');
    legend('估计值', sprintf('真值 (%.0f)', om2^2), 'Location', 'best');
    grid on;

    %% --- 扰动估计质量 ---
    figure('Name', '估计器: 扰动估计 vs 实际', 'Position', [150, 150, 900, 450]);
    subplot(2,1,1);
    plot(t, d_actual, 'k', 'LineWidth', 1.5); hold on;
    plot(t, d_est, 'r--', 'LineWidth', 1.2);
    title('扰动信号与估计 (纯估计器)');
    xlabel('时间 (s)');
    legend('实际扰动 d(t)', '估计扰动 d_{hat}(t)');
    grid on;

    subplot(2,1,2);
    plot(t, d_actual - d_est, 'b', 'LineWidth', 1.2);
    title('估计误差 d(t) - d_{hat}(t)');
    xlabel('时间 (s)');
    grid on;

    fprintf('=== 纯估计器仿真完成 ===\n');
    fprintf('theta1 终值: %.4f (真值: %.2f)\n', th1_hat(end), om1^2);
    fprintf('theta2 终值: %.4f (真值: %.2f)\n', th2_hat(end), om2^2);
end

%% ========================================================================
%  系统动态: 二阶对象 + 自适应频率估计器 (模式 'plant')
% ========================================================================
function dydt = plant_dynamics(t, y)
    % 状态变量:
    % y(1) = y       系统输出
    % y(2) = w0      偏置估计
    % y(3) = w1      第1正弦分量幅值相关量
    % y(4) = w2      第1正弦分量相位相关量
    % y(5) = w3      第2正弦分量幅值相关量
    % y(6) = w4      第2正弦分量相位相关量
    % y(7) = theta1  频率参数1估计 (→ omega1^2)
    % y(8) = theta2  频率参数2估计 (→ omega2^2)
    % y(9) = y_dot   系统输出导数

    % 扰动信号: d(t) = A0 + sin(omega1*t) + sin(omega2*t)
    d = A0 + sin(om1*t) + sin(om2*t);

    % 控制输入: u = -w0_hat - w1_hat - w3_hat
    u = -y(2) - y(3) - y(5);

    % 二阶对象: P(s) = 1/(s+1)^2, 即 y'' + 2*y' + y = u + d
    x1      = y(1);                    % y
    x1_dot  = y(9);                    % y'
    x1_ddot = -2*x1_dot - x1 + u + d;  % y''

    % 自适应频率估计器动态 (Marino & Tomei, 2016, Eq. 34)
    dw0      =  k * x1;               % w0' = k*y
    dw1      =  y(4) + k * x1;        % w1' = w2 + k*y
    dw2      = -y(7) * y(3);          % w2' = -theta1*w1
    dw3      =  y(6) - k * x1;        % w3' = w4 - k*y
    dw4      = -y(8) * y(5);          % w4' = -theta2*w3
    dtheta1  =  epsilon * y(4) * x1;  % theta1' = epsilon*w2*y
    dtheta2  = -epsilon * y(6) * x1;  % theta2' = -epsilon*w4*y

    dydt = [x1_dot; dw0; dw1; dw2; dw3; dw4; dtheta1; dtheta2; x1_ddot];
end

%% ========================================================================
%  估计器动态: 纯频率估计 (模式 'estimator', 无被控对象)
%  直接以 d(t) 为输入, 无需经过被控对象
% ========================================================================
function dydt = estimator_dynamics(t, y)
    % 状态变量:
    % y(1) = w0_hat    偏置估计
    % y(2) = w1_hat    第1正弦分量幅值相关量
    % y(3) = w2_hat    第1正弦分量相位相关量
    % y(4) = w3_hat    第2正弦分量幅值相关量
    % y(5) = w4_hat    第2正弦分量相位相关量
    % y(6) = theta1_hat    频率参数1估计
    % y(7) = theta2_hat    频率参数2估计

    % 扰动信号 (直接作为输入)
    d_t = A0 + sin(om1*t) + sin(om2*t);

    % 估计输出: d_hat = w0_hat + w1_hat + w3_hat
    d_hat = y(1) + y(2) + y(4);

    % 估计误差
    err = d_t - d_hat;

    % 自适应估计器动态
    dw0     =  k * err;
    dw1     =  y(3) + k * err;
    dw2     = -y(6) * y(2);
    dw3     =  y(5) - k * err;
    dw4     = -y(7) * y(4);
    dtheta1 =  epsilon * y(3) * err;
    dtheta2 = -epsilon * y(5) * err;

    dydt = [dw0; dw1; dw2; dw3; dw4; dtheta1; dtheta2];
end

end
