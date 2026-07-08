%% Marino & Tomei (2011) — 调节器方程 Benchmark
% 验证自适应控制器在状态空间扰动模型下的性能，含调节器方程求解
% 与理论轨迹对比。
% 合并自 MarinoTomei_bench_rk4.m 和 MarinoTomei_bench_ode45.m
%
% integrator = 'rk4'   : 手写 RK4 (植物/扰动) + 欧拉 (控制状态)
% integrator = 'ode45' : 逐步 ode45 集成所有状态
%
% 注：原 bench_ode45 中 dchi_hatdt 有一处符号错误 (D(:,1) 代替 d)，
%     合并版已修正为与 bench_rk4 一致的写法。

clear; close all; clc;

%% ====== 用户设置 ======
integrator = 'rk4';      % 'rk4' | 'ode45'

%% ====== 系统参数与初始条件 (共用) ======
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
plant_x0 = [0.5; 0];       % 初始状态 [x1, x2]
xi1_0 = zeros(4,1);         % 滤波器ξ1初始状态
xi2_0 = zeros(4,1);         % 滤波器ξ2初始状态
chi_hat0 = zeros(4,1);     % 估计器χ初始状态
theta_hat0 = [0.1; 0.1];   % 参数估计初始值

% 仿真时间设置
tspan = [0, 150];
dt = 0.01;
t = tspan(1):dt:tspan(2);
num_steps = length(t);

%% ====== 扰动的状态空间模型 (共用) ======
% 状态向量 [x1; x2; x3; x4; x5; x6]
% x1,x2: omega1, x3,x4: omega3, x5,x6: omega2

R = zeros(6, 6);
R(1, 2) = 1;  R(2, 1) = -omega1^2;
R(3, 4) = 1;  R(4, 3) = -omega3^2;
R(5, 6) = 1;  R(6, 5) = -omega2^2;

% 输出映射矩阵
P1 = [1, 0, 0.2, 0, 0, 0];  q1 = [0, 0, 0, 0, 1, 0];    % 区间1 [0,50)
P2 = [1, 0, 0,   0, 0, 0];  q2 = [0, 0, 0, 0, 1, 0];    % 区间2 [50,100)
P3 = [1, 0, 0,   0, 0, 0];  q3 = [0, 0, 0, 0, 0, 0];    % 区间3 [100,150]

% 系统矩阵
A = [a1, 1; 0, 0];
b = [1; b2];
C = [1, 0];

%% ====== 求解调节器方程 (共用) ======
fprintf('正在求解调节器方程...\n');
[Pi1, Gamma1] = solveRegulatorEquation(A, b, C, R, P1, q1, beta);
[Pi2, Gamma2] = solveRegulatorEquation(A, b, C, R, P2, q2, beta);
[Pi3, Gamma3] = solveRegulatorEquation(A, b, C, R, P3, q3, beta);
fprintf('区间1 Π1 = '); disp(Pi1);
fprintf('区间2 Π2 = '); disp(Pi2);
fprintf('区间3 Π3 = '); disp(Pi3);

% 扰动系统初始状态
dist_x0 = [0; omega1; 0; omega3; 0; omega2];

% 赫尔维茨矩阵D
D_matrix = [-d(1), 1, 0, 0;
            -d(2), 0, 1, 0;
            -d(3), 0, 0, 1;
            -d(4), 0, 0, 0];

% 投影算子参数
r_omega = 2;
epsilon_r = 0.1;

%% ====== 仿真执行 ======
% 初始化结果数组 (共用)
x1_history     = zeros(size(t));
x2_history     = zeros(size(t));
u_history       = zeros(size(t));
e_history       = zeros(size(t));
w1_history      = zeros(size(t));
w2_history      = zeros(size(t));
theta1_hat_history = zeros(size(t));
theta2_hat_history = zeros(size(t));
xr_history      = zeros(size(t));
ur_history      = zeros(size(t));
e_theory_history = zeros(size(t));

switch integrator
    case 'rk4'
        %% --- 手写 RK4 + 欧拉 ---
        plant_x  = plant_x0;
        xi1      = xi1_0;
        xi2      = xi2_0;
        chi_hat  = chi_hat0;
        theta_hat = theta_hat0;
        dist_x   = dist_x0;

        for i = 1:num_steps
            current_t = t(i);

            % 根据当前时间选择输出矩阵和调节器解
            if current_t < 50
                w1_val = P1 * dist_x;  w2_val = q1 * dist_x;
                Pi = Pi1;  Gamma = Gamma1;
            elseif current_t < 100
                w1_val = P2 * dist_x;  w2_val = q2 * dist_x;
                Pi = Pi2;  Gamma = Gamma2;
            else
                w1_val = P3 * dist_x;  w2_val = q3 * dist_x;
                Pi = Pi3;  Gamma = Gamma3;
            end

            % 理论参考轨迹
            xr = Pi * dist_x;
            ur = Gamma * dist_x;

            % 参考信号和调节误差
            y_r_val = -w2_val;
            e_val = plant_x(1) - y_r_val;

            % 控制输入
            mu1 = xi1(1);
            mu2 = xi2(1);
            u_val = -k*e_val - chi_hat(1) - mu1*theta_hat(1) - mu2*theta_hat(2);

            % 存储
            x1_history(i) = plant_x(1);
            x2_history(i) = plant_x(2);
            u_history(i) = u_val;
            e_history(i) = e_val;
            w1_history(i) = w1_val;
            w2_history(i) = w2_val;
            theta1_hat_history(i) = theta_hat(1);
            theta2_hat_history(i) = theta_hat(2);
            xr_history(i) = xr(1);
            ur_history(i) = ur;
            e_theory_history(i) = xr(1) - y_r_val;

            if i == num_steps, break; end

            % 导数计算
            dx1dt = plant_x(2) + a1*plant_x(1) + (1/beta)*u_val + w1_val;
            dx2dt = (1/beta)*b2*u_val;
            dxi1dt = D_matrix*xi1 + [0; u_val; 0; 0];
            dxi2dt = D_matrix*xi2 + [0; 0; 0; u_val];
            dchi_hatdt = D_matrix*chi_hat - [d(1); d(2); d(3); d(4)]*u_val;
            phi1 = mu1*e_val;
            phi2 = mu2*e_val;
            theta1_dot = g * projection(phi1, theta_hat(1), r_omega, epsilon_r);
            theta2_dot = g * projection(phi2, theta_hat(2), r_omega, epsilon_r);

            h = dt;

            % 植物状态 — 手写 RK4
            k1_x1 = dx1dt;
            k1_x2 = dx2dt;
            k2_x1 = plant_x(2) + h/2*k1_x2 + a1*(plant_x(1) + h/2*k1_x1) + (1/beta)*u_val + w1_val;
            k2_x2 = (1/beta)*b2*u_val;
            k3_x1 = plant_x(2) + h/2*k2_x2 + a1*(plant_x(1) + h/2*k2_x1) + (1/beta)*u_val + w1_val;
            k3_x2 = (1/beta)*b2*u_val;
            k4_x1 = plant_x(2) + h*k3_x2 + a1*(plant_x(1) + h*k3_x1) + (1/beta)*u_val + w1_val;
            k4_x2 = (1/beta)*b2*u_val;

            plant_x(1) = plant_x(1) + (h/6)*(k1_x1 + 2*k2_x1 + 2*k3_x1 + k4_x1);
            plant_x(2) = plant_x(2) + (h/6)*(k1_x2 + 2*k2_x2 + 2*k3_x2 + k4_x2);

            % 控制状态 — 欧拉
            xi1 = xi1 + dxi1dt * dt;
            xi2 = xi2 + dxi2dt * dt;
            chi_hat = chi_hat + dchi_hatdt * dt;
            theta_hat(1) = theta_hat(1) + theta1_dot * dt;
            theta_hat(2) = theta_hat(2) + theta2_dot * dt;

            % 扰动状态 — RK4
            k1 = R * dist_x;
            k2 = R * (dist_x + 0.5*h*k1);
            k3 = R * (dist_x + 0.5*h*k2);
            k4 = R * (dist_x + h*k3);
            dist_x = dist_x + (h/6)*(k1 + 2*k2 + 2*k3 + k4);
        end

        % 统一变量名
        x1 = x1_history;
        u  = u_history;
        e  = e_history;
        w1 = w1_history;
        w2 = w2_history;
        theta1_hat = theta1_hat_history;
        theta2_hat = theta2_hat_history;
        xr = xr_history;
        ur = ur_history;
        e_theory = e_theory_history;

    case 'ode45'
        %% --- 逐步 ode45 集成 ---
        plant_x  = plant_x0;
        xi1      = xi1_0;
        xi2      = xi2_0;
        chi_hat  = chi_hat0;
        theta_hat = theta_hat0;
        dist_x   = dist_x0;

        for i = 1:num_steps
            current_t = t(i);

            % 根据当前时间选择输出矩阵和调节器解
            if current_t < 50
                w1_val = P1 * dist_x;  w2_val = q1 * dist_x;
                Pi = Pi1;  Gamma_val = Gamma1;
            elseif current_t < 100
                w1_val = P2 * dist_x;  w2_val = q2 * dist_x;
                Pi = Pi2;  Gamma_val = Gamma2;
            else
                w1_val = P3 * dist_x;  w2_val = q3 * dist_x;
                Pi = Pi3;  Gamma_val = Gamma3;
            end

            % 理论参考轨迹
            xr_vec = Pi * dist_x;
            ur_val = Gamma_val * dist_x;

            % 参考信号和调节误差
            y_r_val = -w2_val;
            e_val = plant_x(1) - y_r_val;

            % 控制输入
            mu1 = xi1(1);
            mu2 = xi2(1);
            u_val = -k*e_val - chi_hat(1) - mu1*theta_hat(1) - mu2*theta_hat(2);

            % 存储
            x1_history(i) = plant_x(1);
            x2_history(i) = plant_x(2);
            u_history(i) = u_val;
            e_history(i) = e_val;
            w1_history(i) = w1_val;
            w2_history(i) = w2_val;
            theta1_hat_history(i) = theta_hat(1);
            theta2_hat_history(i) = theta_hat(2);
            xr_history(i) = xr_vec(1);
            ur_history(i) = ur_val;
            e_theory_history(i) = xr_vec(1) - y_r_val;

            if i == num_steps, break; end

            % 使用 ode45 求解一个小步长
            tspan_local = [current_t, current_t+dt];
            all_states = [plant_x; xi1; xi2; chi_hat; theta_hat; dist_x];

            [~, states] = ode45(@(t_s, states) bench_dynamics(t_s, states, ...
                a1, b2, beta, u_val, D_matrix, d, g, r_omega, epsilon_r, ...
                R, P1, P2, P3, q1, q2, q3), tspan_local, all_states);

            final_states = states(end, :)';
            plant_x  = final_states(1:2);
            xi1      = final_states(3:6);
            xi2      = final_states(7:10);
            chi_hat  = final_states(11:14);
            theta_hat = final_states(15:16);
            dist_x   = final_states(17:22);
        end

        % 统一变量名
        x1 = x1_history;
        u  = u_history;
        e  = e_history;
        w1 = w1_history;
        w2 = w2_history;
        theta1_hat = theta1_hat_history;
        theta2_hat = theta2_hat_history;
        xr = xr_history;
        ur = ur_history;
        e_theory = e_theory_history;

    otherwise
        error('Unknown integrator: %s. Use ''rk4'' or ''ode45''.', integrator);
end

%% ====== 图1：输出、理论轨迹和误差 (共用) ======
figure;
subplot(2,2,1);
plot(t, x1, 'b-', t, xr, 'r--', t, -w2, 'g-.', 'LineWidth', 1.5);
title('系统输出与理论轨迹');
xlabel('时间 (s)'); ylabel('输出');
legend('实际输出 y=x1', '理论轨迹 xr', '参考信号 -w2', 'Location', 'best');
grid on;

subplot(2,2,2);
plot(t, u, 'b-', t, ur, 'r--', 'LineWidth', 1.5);
title('控制输入对比');
xlabel('时间 (s)'); ylabel('输入');
legend('实际输入 u', '理论输入 ur', 'Location', 'best');
grid on;

subplot(2,2,3);
plot(t, e, 'b-', t, e_theory, 'r--', 'LineWidth', 1.5);
title('跟踪误差对比');
xlabel('时间 (s)'); ylabel('误差');
legend('实际误差 e', '理论误差', 'Location', 'best');
grid on;

subplot(2,2,4);
plot(t, x1 - xr, 'k-', 'LineWidth', 1.5);
title('实际输出与理论轨迹的偏差');
xlabel('时间 (s)'); ylabel('偏差');
grid on;

%% ====== 图2：调节误差与参数估计 (共用) ======
figure;
subplot(2,2,1);
plot(t, e, 'b');
title('输出调节误差');
xlabel('时间 (s)'); ylabel('e_{output}');

subplot(2,2,2);
plot(t, u - ur, 'b');
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

%% ====== 图3：系统输出、扰动与输入 (共用) ======
figure;
subplot(2,2,1);
plot(t, x1, 'b', t, -w2, 'r--');
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
plot(t, ur, 'b');
title('参考输入信号');
xlabel('时间 (s)'); ylabel('u_r(t)');
legend('u_r');

%% ====== 图4：区间切换点分析 (共用) ======
figure;
transition1 = 45:55;
transition2 = 95:105;

subplot(2,2,1);
plot(t(transition1), x1(transition1), 'b-', t(transition1), xr(transition1), 'r--', 'LineWidth', 1.5);
title('区间1-2切换点附近的输出');
xlabel('时间 (s)'); ylabel('输出');
legend('实际输出', '理论轨迹', 'Location', 'best'); grid on;

subplot(2,2,2);
plot(t(transition1), u(transition1), 'b-', t(transition1), ur(transition1), 'r--', 'LineWidth', 1.5);
title('区间1-2切换点附近的控制输入');
xlabel('时间 (s)'); ylabel('输入');
legend('实际输入', '理论输入', 'Location', 'best'); grid on;

subplot(2,2,3);
plot(t(transition2), x1(transition2), 'b-', t(transition2), xr(transition2), 'r--', 'LineWidth', 1.5);
title('区间2-3切换点附近的输出');
xlabel('时间 (s)'); ylabel('输出');
legend('实际输出', '理论轨迹', 'Location', 'best'); grid on;

subplot(2,2,4);
plot(t(transition2), u(transition2), 'b-', t(transition2), ur(transition2), 'r--', 'LineWidth', 1.5);
title('区间2-3切换点附近的控制输入');
xlabel('时间 (s)'); ylabel('输入');
legend('实际输入', '理论输入', 'Location', 'best'); grid on;

%% ========================================================================
%  共用辅助函数
%% ========================================================================

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

function [Pi, Gamma] = solveRegulatorEquation(A, b, C, R, P, q, beta)
    % 数值求解调节器方程: Π*R - A*Π - (1/β)*b*Γ = P,  C*Π = q
    n_state = size(A, 1);
    n_dist  = size(R, 1);
    n_input = 1;

    coef_matrix = zeros(n_state*n_dist + n_dist, n_state*n_dist + n_input*n_dist);

    % Π*R - A*Π
    for i = 1:n_state
        for j = 1:n_dist
            row_idx = (i-1)*n_dist + j;
            for k = 1:n_dist
                col_idx = (i-1)*n_dist + k;
                coef_matrix(row_idx, col_idx) = R(k, j);
            end
            for k = 1:n_state
                for l = 1:n_dist
                    col_idx = (k-1)*n_dist + l;
                    coef_matrix(row_idx, col_idx) = coef_matrix(row_idx, col_idx) - A(i, k) * (l == j);
                end
            end
        end
    end

    % -(1/β)*b*Γ
    for i = 1:n_state
        for j = 1:n_dist
            row_idx = (i-1)*n_dist + j;
            col_idx = n_state*n_dist + j;
            coef_matrix(row_idx, col_idx) = -(1/beta)*b(i);
        end
    end

    % C*Π = q
    for j = 1:n_dist
        row_idx = n_state*n_dist + j;
        for i = 1:n_state
            col_idx = (i-1)*n_dist + j;
            coef_matrix(row_idx, col_idx) = C(i);
        end
    end

    % 右侧向量
    rhs = zeros(n_state*n_dist + n_dist, 1);
    for i = 1:n_state
        for j = 1:n_dist
            row_idx = (i-1)*n_dist + j;
            if i == 1, rhs(row_idx) = P(j); end
        end
    end
    for j = 1:n_dist
        row_idx = n_state*n_dist + j;
        rhs(row_idx) = q(j);
    end

    sol = coef_matrix \ rhs;

    Pi = zeros(n_state, n_dist);
    for i = 1:n_state
        for j = 1:n_dist
            Pi(i, j) = sol((i-1)*n_dist + j);
        end
    end

    Gamma = zeros(1, n_dist);
    for j = 1:n_dist
        Gamma(j) = sol(n_state*n_dist + j);
    end

    % 验证
    residual1 = Pi*R - A*Pi - (1/beta)*b*Gamma - [P; zeros(size(A,1)-1, size(P,2))];
    residual2 = C*Pi - q;
    if norm(residual1, 'fro') > 1e-10 || norm(residual2, 'fro') > 1e-10
        warning('解的残差较大，可能存在数值问题');
    end
end

%% ========================================================================
%  ode45 模式专用: 逐步积分动力学函数
%% ========================================================================
function dstates = bench_dynamics(t_s, states, a1, b2, beta, u_ctrl, D_matrix, d, g, r_omega, epsilon_r, R, P1, P2, P3, q1, q2, q3)
    plant_x  = states(1:2);
    xi1      = states(3:6);
    xi2      = states(7:10);
    chi_hat  = states(11:14);
    theta_hat = states(15:16);
    dist_x   = states(17:22);

    % 扰动输出 (依时间区间)
    if t_s < 50
        w1_val = P1 * dist_x;  w2_val = q1 * dist_x;
    elseif t_s < 100
        w1_val = P2 * dist_x;  w2_val = q2 * dist_x;
    else
        w1_val = P3 * dist_x;  w2_val = q3 * dist_x;
    end

    y_r_val = -w2_val;
    e_val = plant_x(1) - y_r_val;

    % 主系统状态导数
    dx1dt = plant_x(2) + a1*plant_x(1) + (1/beta)*u_ctrl + w1_val;
    dx2dt = (1/beta)*b2*u_ctrl;

    % 滤波器
    dxi1dt = D_matrix*xi1 + [0; u_ctrl; 0; 0];
    dxi2dt = D_matrix*xi2 + [0; 0; 0; u_ctrl];

    % 估计器χ (已修正：与 RK4 版统一使用 d 向量)
    dchi_hatdt = D_matrix*chi_hat - [d(1); d(2); d(3); d(4)]*u_ctrl;

    % 参数更新律
    mu1 = xi1(1);
    mu2 = xi2(1);
    phi1 = mu1*e_val;
    phi2 = mu2*e_val;
    theta1_dot = g * projection(phi1, theta_hat(1), r_omega, epsilon_r);
    theta2_dot = g * projection(phi2, theta_hat(2), r_omega, epsilon_r);

    % 扰动系统状态导数
    ddist_x = R * dist_x;

    dstates = [dx1dt; dx2dt; dxi1dt; dxi2dt; dchi_hatdt; theta1_dot; theta2_dot; ddist_x];
end
