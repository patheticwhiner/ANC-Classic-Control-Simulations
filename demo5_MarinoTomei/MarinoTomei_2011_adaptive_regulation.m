%% Marino & Tomei (2011) — 自适应调节 (含内部模型参考控制器)
% 合并自 MarinoTomei_euler.m 和 MarinoTomei_ode45.m
% 通过 integrator 开关选择数值求解方法，模型、参数、绘图完全一致。
%
% integrator = 'euler'  : 固定步长欧拉法 (dt=0.01)
% integrator = 'ode45'  : 自适应步长 ode45 (RelTol=1e-6, AbsTol=1e-8)

clear; close all; clc;

%% ====== 用户设置 ======
integrator = 'ode45';      % 'euler' | 'ode45'

%% ====== 系统参数与初始条件 (共用) ======
a1 = 0.5;       % 未知参数，真实值
b2 = 0.5;       % 未知参数，真实值
beta = 0.5;     % 未知参数，真实值
omega1 = 1;     % 扰动频率1
omega2 = 0.5;   % 扰动频率2
omega3 = 4;     % 扰动频率3（仅用于t < 50s）

% 控制器参数
k = 7.100;      % 控制增益
g = 100;        % 参数更新增益
d = [2, 3, 2, 1]; % 矩阵D的参数 [d2, d3, d4, d5]

% 确定m值（基于omega1和omega2）
m = 2;  % 对于两个频率的情况

% 初始条件
x0 = [0.5; 0];          % 初始状态 [x1, x2]
xi1_0 = zeros(4,1);     % 滤波器ξ1初始状态
xi2_0 = zeros(4,1);     % 滤波器ξ2初始状态
chi_hat0 = zeros(4,1);  % 估计器χ初始状态
theta_hat0 = [0.1; 0.1];% 参数估计初始值 [θ1_hat, θ2_hat]

% 定义赫尔维茨矩阵D（确保特征值在左半平面）
D_matrix = [-d(1), 1, 0, 0;
            -d(2), 0, 1, 0;
            -d(3), 0, 0, 1;
            -d(4), 0, 0, 0];

%% ====== 仿真执行 ======
switch integrator
    case 'euler'
        %% --- 固定步长欧拉法 ---
        wc = zeros(2*m+1,1);  % 参考控制器状态
        wc(1) = -x0(1);

        t_start = 0;
        t_end = 150;
        dt = 0.01;
        num_steps = (t_end - t_start) / dt;

        % 状态初始化
        x = x0;
        xi1 = xi1_0;
        xi2 = xi2_0;
        chi_hat = chi_hat0;
        theta_hat = theta_hat0;

        % 初始化结果数组
        t = linspace(t_start, t_end, num_steps+1)';
        x1_history = zeros(size(t));
        e_output_history = zeros(size(t));
        e_input_history = zeros(size(t));
        theta1_hat_history = zeros(size(t));
        theta2_hat_history = zeros(size(t));
        u_history = zeros(size(t));
        u_r_history = zeros(size(t));
        w1_history = zeros(size(t));
        w2_history = zeros(size(t));

        % 投影算子参数
        r_omega = 2;
        epsilon_r = 0.1;

        % 循环仿真
        for i = 1:num_steps+1
            current_t = t(i);

            % 扰动
            w1_val = get_w1(current_t, omega1, omega3);
            w2_val = get_w2(current_t, omega2);

            % 参考信号和调节误差
            y_r_val = -w2_val;
            e_val = x(1) - y_r_val;

            % 控制输入
            mu1 = xi1(1);
            mu2 = xi2(1);
            u_val = -k*e_val - chi_hat(1) - mu1*theta_hat(1) - mu2*theta_hat(2);

            % 参考控制器输入
            u_r_val = wc(1);

            % 存储
            x1_history(i) = x(1);
            e_output_history(i) = e_val;
            e_input_history(i) = u_val - u_r_val;
            theta1_hat_history(i) = theta_hat(1);
            theta2_hat_history(i) = theta_hat(2);
            u_history(i) = u_val;
            u_r_history(i) = u_r_val;
            w1_history(i) = w1_val;
            w2_history(i) = w2_val;

            if i == num_steps+1
                break;
            end

            % 系统状态导数
            dx1dt = x(2) + a1*x(1) + (1/beta)*u_val + w1_val;
            dx2dt = (1/beta)*b2*u_val;

            % 滤波器动态
            dxi1dt = D_matrix*xi1 + [0; u_val; 0; 0];
            dxi2dt = D_matrix*xi2 + [0; 0; 0; u_val];

            % 估计器χ动态
            dchi_hatdt = D_matrix*chi_hat - [d(1); d(2); d(3); d(4)]*u_val;

            % 参数更新律
            phi1 = mu1*e_val;
            phi2 = mu2*e_val;
            theta1_hat_dot = g * projection(phi1, theta_hat(1), r_omega, epsilon_r);
            theta2_hat_dot = g * projection(phi2, theta_hat(2), r_omega, epsilon_r);

            % 参考控制器动态
            theta1_true = omega1^2 + omega2^2;
            theta2_true = omega1^2 * omega2^2;
            Rc = zeros(2*m+1);
            Rc(1,2) = 1;
            Rc(2,1) = -theta1_true;
            Rc(2,3) = 1;
            Rc(3,4) = 1;
            Rc(4,1) = -theta2_true;
            Rc(4,5) = 1;
            dwcdt = Rc * wc;

            % 欧拉更新
            x = x + [dx1dt; dx2dt] * dt;
            xi1 = xi1 + dxi1dt * dt;
            xi2 = xi2 + dxi2dt * dt;
            chi_hat = chi_hat + dchi_hatdt * dt;
            theta_hat(1) = theta_hat(1) + theta1_hat_dot * dt;
            theta_hat(2) = theta_hat(2) + theta2_hat_dot * dt;
            wc = wc + dwcdt * dt;
        end

        % 统一变量名，供后续共用绘图
        e_output = e_output_history;
        e_input  = e_input_history;
        theta1_hat = theta1_hat_history;
        theta2_hat = theta2_hat_history;
        x1 = x1_history;
        u  = u_history;
        w1 = w1_history;
        w2 = w2_history;
        u_r = u_r_history;
        y_r = -w2_history;

    case 'ode45'
        %% --- 自适应步长 ode45 ---
        wc_0 = zeros(2*m+1,1);
        wc_0(1) = -x0(1);

        tspan = [0, 150];
        options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

        [t, X] = ode45(@(t, X) system_dynamics(t, X, a1, b2, beta, ...
            omega1, omega2, omega3, k, g, D_matrix, d), tspan, ...
            [x0; xi1_0; xi2_0; chi_hat0; theta_hat0; wc_0], options);

        % 提取结果
        x1 = X(:,1);
        xi1 = X(:,3:6);
        xi2 = X(:,7:10);
        chi_hat = X(:,11:14);
        theta1_hat = X(:,15);
        theta2_hat = X(:,16);
        wc = X(:,17:end);

        % 计算扰动和导出量
        w1 = arrayfun(@(t_val) get_w1(t_val, omega1, omega3), t);
        w2 = arrayfun(@(t_val) get_w2(t_val, omega2), t);
        y_r = -w2;
        e_output = x1 - y_r;
        u_r = wc(:,1);
        u = -k*e_output - chi_hat(:,1) - xi1(:,1).*theta1_hat - xi2(:,1).*theta2_hat;
        e_input = u - u_r;

    otherwise
        error('Unknown integrator: %s. Use ''euler'' or ''ode45''.', integrator);
end

%% ====== 图1：调节误差与参数估计 (共用) ======
figure;
subplot(2,2,1);
plot(t, e_output, 'b');
title('输出调节误差');
xlabel('时间 (s)'); ylabel('e_{output}');

subplot(2,2,2);
plot(t, e_input, 'r');
title('输入调节误差');
xlabel('时间 (s)'); ylabel('e_{input}');

subplot(2,2,3);
plot(t, theta1_hat, 'b', t, ones(size(t))*(omega1^2 + omega2^2), 'k--');
title('\theta_1 估计');
xlabel('时间 (s)');
legend('估计值', '真实值', 'location', 'best');

subplot(2,2,4);
plot(t, theta2_hat, 'r', t, ones(size(t))*(omega1^2 * omega2^2), 'k--');
title('\theta_2 估计');
xlabel('时间 (s)');
legend('估计值', '真实值', 'location', 'best');

%% ====== 图2：系统输出、扰动与输入参考 (共用) ======
figure;
subplot(2,2,1);
plot(t, x1, 'b', t, y_r, 'r--');
title('系统输出与参考信号');
xlabel('时间 (s)'); ylabel('输出');
legend('y = x1', 'y_r = -w2');

subplot(2,2,2);
plot(t, u, 'k');
title('控制输入');
xlabel('时间 (s)'); ylabel('u');

subplot(2,2,3);
plot(t, w1, 'g', t, w2, 'm');
title('扰动信号');
xlabel('时间 (s)'); ylabel('扰动');
legend('w1', 'w2');

subplot(2,2,4);
plot(t, u_r, 'b');
title('输入参考信号');
xlabel('时间 (s)'); ylabel('参考输入');
legend('u_r');

%% ========================================================================
%  共用辅助函数
%% ========================================================================

function w1 = get_w1(t, omega1, omega3)
    if t < 50
        w1 = sin(omega1*t) + 0.2*sin(omega3*t);
    else
        w1 = sin(omega1*t);
    end
end

function w2 = get_w2(t, omega2)
    if t < 100
        w2 = sin(omega2*t);
    else
        w2 = 0;
    end
end

function phi_proj = projection(phi, theta_hat, r_omega, epsilon_r)
    if isscalar(theta_hat) && isscalar(phi)
        p_r = (theta_hat^2 - r_omega^2)/(epsilon_r^2 + 2*epsilon_r*r_omega);
        grad_p_r = 2*theta_hat/(epsilon_r^2 + 2*epsilon_r*r_omega);
        if p_r <= 0
            phi_proj = phi;
        else
            inner_product = grad_p_r * phi;
            if inner_product <= 0
                phi_proj = phi;
            else
                grad_norm_squared = grad_p_r^2;
                projection_factor = p_r * (grad_p_r^2) / grad_norm_squared;
                phi_proj = (1 - projection_factor) * phi;
            end
        end
    else
        p_r = (norm(theta_hat)^2 - r_omega^2)/(epsilon_r^2 + 2*epsilon_r*r_omega);
        grad_p_r = 2*theta_hat/(epsilon_r^2 + 2*epsilon_r*r_omega);
        if p_r <= 0
            phi_proj = phi;
        else
            inner_product = dot(grad_p_r, phi);
            if inner_product <= 0
                phi_proj = phi;
            else
                grad_norm_squared = norm(grad_p_r)^2;
                projection_matrix = eye(length(theta_hat)) - p_r * (grad_p_r * grad_p_r') / grad_norm_squared;
                phi_proj = projection_matrix * phi;
            end
        end
    end
end

%% ========================================================================
%  ode45 模式专用: 系统动力学函数
%% ========================================================================
function dXdt = system_dynamics(t, X, a1, b2, beta, omega1, omega2, omega3, k, g, D_matrix, d)
    x = X(1:2);
    xi1 = X(3:6);
    xi2 = X(7:10);
    chi_hat = X(11:14);
    theta_hat = X(15:16);
    wc = X(17:end);

    m = 2;

    % 扰动
    w1 = get_w1(t, omega1, omega3);
    w2 = get_w2(t, omega2);

    % 跟踪误差
    e = x(1) - (-w2);

    % 参考控制器动态
    theta1_true = omega1^2 + omega2^2;
    theta2_true = omega1^2 * omega2^2;
    Rc = zeros(2*m+1);
    Rc(1,2) = 1;
    Rc(2,1) = -theta1_true;
    Rc(2,3) = 1;
    Rc(3,4) = 1;
    Rc(4,1) = -theta2_true;
    Rc(4,5) = 1;
    dwcdt = Rc * wc;

    % 控制输入
    mu1 = xi1(1);
    mu2 = xi2(1);
    u = -k*e - chi_hat(1) - mu1*theta_hat(1) - mu2*theta_hat(2);

    % 系统动态
    dx1dt = x(2) + a1*x(1) + (1/beta)*u + w1;
    dx2dt = (1/beta)*b2*u;

    % 滤波器
    dxi1dt = D_matrix*xi1 + [0; u; 0; 0];
    dxi2dt = D_matrix*xi2 + [0; 0; 0; u];

    % 估计器χ
    dchi_hatdt = D_matrix*chi_hat - [d(1); d(2); d(3); d(4)]*u;

    % 参数更新律
    r_omega = 2;
    epsilon_r = 0.1;
    phi1 = mu1*e;
    phi2 = mu2*e;
    theta1_hat_dot = g * projection(phi1, theta_hat(1), r_omega, epsilon_r);
    theta2_hat_dot = g * projection(phi2, theta_hat(2), r_omega, epsilon_r);

    dXdt = [dx1dt; dx2dt; dxi1dt; dxi2dt; dchi_hatdt; theta1_hat_dot; theta2_hat_dot; dwcdt];
end
