clear; close all; clc;

%% 系统参数与初始条件
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
plant_x = [0.5; 0];       % 初始状态 [x1, x2]
xi1 = zeros(4,1);         % 滤波器ξ1初始状态
xi2 = zeros(4,1);         % 滤波器ξ2初始状态
chi_hat = zeros(4,1);     % 估计器χ初始状态
theta_hat = [0.1; 0.1];   % 参数估计初始值 [θ1_hat, θ2_hat]

% 仿真时间设置
tspan = [0, 150];         % 仿真时间范围：0到150秒
dt = 0.01;                % 时间步长
t = tspan(1):dt:tspan(2);
num_steps = length(t);

%% 构建扰动的状态空间模型
% 状态向量 [x1; x2; x3; x4; x5; x6]
% x1,x2: 对应omega1的正弦分量
% x3,x4: 对应omega3的正弦分量
% x5,x6: 对应omega2的正弦分量

% 构建完整的状态空间矩阵R（6x6）
R = zeros(6, 6);
% omega1对应的状态方程
R(1, 2) = 1;
R(2, 1) = -omega1^2;
% omega3对应的状态方程
R(3, 4) = 1;
R(4, 3) = -omega3^2;
% omega2对应的状态方程
R(5, 6) = 1;
R(6, 5) = -omega2^2;

% 输出映射矩阵（区间1 [0, 50)）
P1 = [1, 0, 0.2, 0, 0, 0];  % w1 = x1 + 0.2*x3
q1 = [0, 0, 0, 0, 1, 0];    % w2 = x5

% 输出映射矩阵（区间2 [50, 100)）
P2 = [1, 0, 0, 0, 0, 0];    % w1 = x1
q2 = [0, 0, 0, 0, 1, 0];    % w2 = x5

% 输出映射矩阵（区间3 [100, 150]）
P3 = [1, 0, 0, 0, 0, 0];    % w1 = x1
q3 = [0, 0, 0, 0, 0, 0];    % w2 = 0

% 定义系统矩阵用于求解调节器方程
A = [a1, 1; 0, 0];          % 系统状态矩阵
b = [1; b2];                % 输入矩阵
C = [1, 0];                 % 输出矩阵

%% 求解三个区间的调节器方程: ΠR = AΠ + 1/β*b*Γ + P, q = CΠ
fprintf('正在求解调节器方程...\n');
% 区间1求解
[Pi1, Gamma1] = solveRegulatorEquation(A, b, C, R, P1, q1, beta);
fprintf('区间1 [0, 50) 调节器方程解:\n');
fprintf('Π1 = \n'); disp(Pi1);
fprintf('Γ1 = \n'); disp(Gamma1);

% 区间2求解
[Pi2, Gamma2] = solveRegulatorEquation(A, b, C, R, P2, q2, beta);
fprintf('区间2 [50, 100) 调节器方程解:\n');
fprintf('Π2 = \n'); disp(Pi2);
fprintf('Γ2 = \n'); disp(Gamma2);

% 区间3求解
[Pi3, Gamma3] = solveRegulatorEquation(A, b, C, R, P3, q3, beta);
fprintf('区间3 [100, 150] 调节器方程解:\n');
fprintf('Π3 = \n'); disp(Pi3);
fprintf('Γ3 = \n'); disp(Gamma3);

% 扰动系统初始状态
dist_x = [0; omega1; 0; omega3; 0; omega2]; % 初始值使得所有信号在t=0时为0

%% 定义赫尔维茨矩阵D（确保特征值在左半平面）
D = [-d(1), 1, 0, 0;
     -d(2), 0, 1, 0;
     -d(3), 0, 0, 1;
     -d(4), 0, 0, 0]; % 示例结构，可根据实际调整

% 投影算子参数
r_omega = 2;      % 参数约束半径
epsilon_r = 0.1;  % 投影算子平滑参数

% 初始化结果数组
x1_history = zeros(size(t));
x2_history = zeros(size(t));
xi1_history = zeros(length(t), 4);
xi2_history = zeros(length(t), 4);
chi_hat_history = zeros(length(t), 4);
theta1_hat_history = zeros(size(t));
theta2_hat_history = zeros(size(t));
w1_history = zeros(size(t));
w2_history = zeros(size(t));
u_history = zeros(size(t));
e_history = zeros(size(t));
dist_x_history = zeros(length(t), 6);
% 初始化理论参考轨迹和控制输入
xr_history = zeros(size(t));
ur_history = zeros(size(t));
e_theory_history = zeros(size(t)); % 理论误差

% 循环仿真
for i = 1:num_steps
    current_t = t(i);
    
    % 存储当前状态
    dist_x_history(i,:) = dist_x';
    
    % 根据当前时间确定使用哪个输出矩阵和调节器方程解
    if current_t < 50
        % 区间1: [0, 50)
        w1 = P1 * dist_x;
        w2 = q1 * dist_x;
        Pi = Pi1;
        Gamma = Gamma1;
    elseif current_t < 100
        % 区间2: [50, 100)
        w1 = P2 * dist_x;
        w2 = q2 * dist_x;
        Pi = Pi2;
        Gamma = Gamma2;
    else
        % 区间3: [100, 150]
        w1 = P3 * dist_x;
        w2 = q3 * dist_x;
        Pi = Pi3;
        Gamma = Gamma3;
    end
    
    % 计算理论参考轨迹和控制输入
    xr = Pi * dist_x;         % 理论状态轨迹
    ur = Gamma * dist_x;      % 理论控制输入
    
    % 存储理论轨迹
    xr_history(i) = xr(1);    % 只记录x1分量
    ur_history(i) = ur;
    
    % 计算参考信号和调节误差
    y_r = -w2;
    e = plant_x(1) - y_r;
    e_theory = xr(1) - y_r;   % 理论误差
    
    % 计算控制输入
    mu1 = xi1(1);
    mu2 = xi2(1);
    u = -k*e - chi_hat(1) - mu1*theta_hat(1) - mu2*theta_hat(2);
    
    % 存储状态和结果
    x1_history(i) = plant_x(1);
    x2_history(i) = plant_x(2);
    xi1_history(i,:) = xi1';
    xi2_history(i,:) = xi2';
    chi_hat_history(i,:) = chi_hat';
    theta1_hat_history(i) = theta_hat(1);
    theta2_hat_history(i) = theta_hat(2);
    w1_history(i) = w1;
    w2_history(i) = w2;
    u_history(i) = u;
    e_history(i) = e;
    e_theory_history(i) = e_theory;
    
    % 如果是最后一步，跳出循环
    if i == num_steps
        break;
    end
    
    % 计算系统状态导数
    dx1dt = plant_x(2) + a1*plant_x(1) + (1/beta)*u + w1;
    dx2dt = (1/beta)*b2*u;
    
    % 滤波器动态方程
    dxi1dt = D*xi1 + [0; u; 0; 0];
    dxi2dt = D*xi2 + [0; 0; 0; u];
    
    % 估计器χ动态方程
    dchi_hatdt = D*chi_hat - [d(1); d(2); d(3); d(4)]*u;
    
    % 参数更新律
    phi1 = mu1*e;
    phi2 = mu2*e;
    theta1_hat_dot = g * projection(phi1, theta_hat(1), r_omega, epsilon_r);
    theta2_hat_dot = g * projection(phi2, theta_hat(2), r_omega, epsilon_r);
    
    % 扰动系统状态导数
    ddist_x = R * dist_x;
    
    % 使用四阶龙格-库塔法更新状态
    h = dt;
    
    % 主系统状态更新
    k1_x1 = dx1dt;
    k1_x2 = dx2dt;
    k2_x1 = plant_x(2) + h/2*k1_x2 + a1*(plant_x(1) + h/2*k1_x1) + (1/beta)*u + w1;
    k2_x2 = (1/beta)*b2*u;
    k3_x1 = plant_x(2) + h/2*k2_x2 + a1*(plant_x(1) + h/2*k2_x1) + (1/beta)*u + w1;
    k3_x2 = (1/beta)*b2*u;
    k4_x1 = plant_x(2) + h*k3_x2 + a1*(plant_x(1) + h*k3_x1) + (1/beta)*u + w1;
    k4_x2 = (1/beta)*b2*u;
    
    plant_x(1) = plant_x(1) + (h/6)*(k1_x1 + 2*k2_x1 + 2*k3_x1 + k4_x1);
    plant_x(2) = plant_x(2) + (h/6)*(k1_x2 + 2*k2_x2 + 2*k3_x2 + k4_x2);
    
    % 控制系统状态更新（简化为欧拉法）
    xi1 = xi1 + dxi1dt * dt;
    xi2 = xi2 + dxi2dt * dt;
    chi_hat = chi_hat + dchi_hatdt * dt;
    theta_hat(1) = theta_hat(1) + theta1_hat_dot * dt;
    theta_hat(2) = theta_hat(2) + theta2_hat_dot * dt;
    
    % 扰动系统状态更新（使用RK4）
    k1 = R * dist_x;
    k2 = R * (dist_x + 0.5 * h * k1);
    k3 = R * (dist_x + 0.5 * h * k2);
    k4 = R * (dist_x + h * k3);
    
    dist_x = dist_x + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end

%% 绘制图1：输出、理论轨迹和误差
figure;
subplot(2,2,1);
plot(t, x1_history, 'b-', t, xr_history, 'r--', t, -w2_history, 'g-.', 'LineWidth', 1.5);
title('系统输出与理论轨迹');
xlabel('时间 (s)');
ylabel('输出');
legend('实际输出 y=x1', '理论轨迹 xr', '参考信号 -w2', 'Location', 'best');
grid on;

subplot(2,2,2);
plot(t, u_history, 'b-', t, ur_history, 'r--', 'LineWidth', 1.5);
title('控制输入对比');
xlabel('时间 (s)');
ylabel('输入');
legend('实际输入 u', '理论输入 ur', 'Location', 'best');
grid on;

subplot(2,2,3);
plot(t, e_history, 'b-', t, e_theory_history, 'r--', 'LineWidth', 1.5);
title('跟踪误差对比');
xlabel('时间 (s)');
ylabel('误差');
legend('实际误差 e', '理论误差', 'Location', 'best');
grid on;

subplot(2,2,4);
plot(t, x1_history - xr_history, 'k-', 'LineWidth', 1.5);
title('实际输出与理论轨迹的偏差');
xlabel('时间 (s)');
ylabel('偏差');
grid on;

%% 绘制图2：调节误差与参数估计
figure;
subplot(2,2,1);
plot(t, e_history, 'b');
title('输出调节误差');
xlabel('时间 (s)');
ylabel('e_{output}');

subplot(2,2,2);
plot(t, u_history-ur_history, 'b');
title('输出调节误差');
xlabel('时间 (s)');
ylabel('e_{input}');

subplot(2,2,3);
plot(t, theta1_hat_history, 'b', t, ones(size(t))*(omega1^2 + omega2^2), 'k--');
title('θ_1估计');
xlabel('时间 (s)');
legend('估计值', '真实值','location','best');

subplot(2,2,4);
plot(t, theta2_hat_history, 'r', t, ones(size(t))*(omega1^2 * omega2^2), 'k--');
title('θ_2估计');
xlabel('时间 (s)');
legend('估计值', '真实值','location','best');

%% 绘制图3：系统输出、扰动与输入
figure;
subplot(2,2,1);
plot(t, x1_history, 'b', t, -w2_history, 'r--');
title('系统输出与参考信号');
xlabel('时间 (s)');
ylabel('输出');
legend('y = x1', 'y_r = -w2');

subplot(2,2,2);
plot(t, u_history, 'k');
title('控制输入');
xlabel('时间 (s)');
ylabel('u');

subplot(2,2,3);
plot(t, w1_history, 'g', t, w2_history, 'm');
title('扰动信号');
xlabel('时间 (s)');
ylabel('扰动');
legend('w1', 'w2');

subplot(2,2,4);
plot(t, ur_history, 'b');
title('参考输入信号');
xlabel('时间 (s)');
ylabel('u_r(t)');
legend('u_r');

% 分析区间切换点
figure;
transition1 = 45:55;  % 第一个区间切换点附近
transition2 = 95:105; % 第二个区间切换点附近

subplot(2,2,1);
plot(t(transition1), x1_history(transition1), 'b-', t(transition1), xr_history(transition1), 'r--', 'LineWidth', 1.5);
title('区间1-2切换点附近的输出');
xlabel('时间 (s)');
ylabel('输出');
legend('实际输出', '理论轨迹', 'Location', 'best');
grid on;

subplot(2,2,2);
plot(t(transition1), u_history(transition1), 'b-', t(transition1), ur_history(transition1), 'r--', 'LineWidth', 1.5);
title('区间1-2切换点附近的控制输入');
xlabel('时间 (s)');
ylabel('输入');
legend('实际输入', '理论输入', 'Location', 'best');
grid on;

subplot(2,2,3);
plot(t(transition2), x1_history(transition2), 'b-', t(transition2), xr_history(transition2), 'r--', 'LineWidth', 1.5);
title('区间2-3切换点附近的输出');
xlabel('时间 (s)');
ylabel('输出');
legend('实际输出', '理论轨迹', 'Location', 'best');
grid on;

subplot(2,2,4);
plot(t(transition2), u_history(transition2), 'b-', t(transition2), ur_history(transition2), 'r--', 'LineWidth', 1.5);
title('区间2-3切换点附近的控制输入');
xlabel('时间 (s)');
ylabel('输入');
legend('实际输入', '理论输入', 'Location', 'best');
grid on;

%% 投影算子的实现
function phi_proj = projection(phi, theta_hat, r_omega, epsilon_r)
    % 检查是否为标量情况
    if isscalar(theta_hat) && isscalar(phi)
        % 标量情况的处理
        % 计算p_r(theta_hat)
        p_r = (theta_hat^2 - r_omega^2)/(epsilon_r^2 + 2*epsilon_r*r_omega);
        
        % 计算梯度 grad p_r(theta_hat)
        grad_p_r = 2*theta_hat/(epsilon_r^2 + 2*epsilon_r*r_omega);
        
        % 根据条件判断使用哪种投影方式
        if p_r <= 0
            % 情况1：p_r(theta_hat) <= 0
            phi_proj = phi;
        else
            % 计算内积 <grad p_r(theta_hat), phi>
            inner_product = grad_p_r * phi;
            
            if inner_product <= 0
                % 情况2：p_r(theta_hat) > 0 且 <grad p_r(theta_hat), phi> <= 0
                phi_proj = phi;
            else
                % 情况3：p_r(theta_hat) > 0 且 <grad p_r(theta_hat), phi> > 0
                % 计算投影矩阵并应用到phi
                grad_norm_squared = grad_p_r^2;
                projection_factor = p_r * (grad_p_r^2) / grad_norm_squared;
                phi_proj = (1 - projection_factor) * phi;
            end
        end
    else
        % 向量情况的处理
        % 计算p_r(theta_hat)
        p_r = (norm(theta_hat)^2 - r_omega^2)/(epsilon_r^2 + 2*epsilon_r*r_omega);
        
        % 计算梯度 grad p_r(theta_hat)
        grad_p_r = 2*theta_hat/(epsilon_r^2 + 2*epsilon_r*r_omega);
        
        % 根据条件判断使用哪种投影方式
        if p_r <= 0
            % 情况1：p_r(theta_hat) <= 0
            phi_proj = phi;
        else
            % 计算内积 <grad p_r(theta_hat), phi>
            inner_product = dot(grad_p_r, phi);
            
            if inner_product <= 0
                % 情况2：p_r(theta_hat) > 0 且 <grad p_r(theta_hat), phi> <= 0
                phi_proj = phi;
            else
                % 情况3：p_r(theta_hat) > 0 且 <grad p_r(theta_hat), phi> > 0
                % 计算投影矩阵并应用到phi
                grad_norm_squared = norm(grad_p_r)^2;
                projection_matrix = eye(length(theta_hat)) - p_r * (grad_p_r * grad_p_r') / grad_norm_squared;
                phi_proj = projection_matrix * phi;
            end
        end
    end
end

%% 数值方法求解调节器方程
function [Pi, Gamma] = solveRegulatorEquation(A, b, C, R, P, q, beta)
    % 使用数值方法求解调节器方程
    % 方程: Π*R - A*Π - (1/β)*b*Γ = P, C*Π = q
    
    n_state = size(A, 1);      % 系统状态维数
    n_dist = size(R, 1);       % 扰动状态维数
    n_input = 1;               % 输入维数
    
    % 创建系数矩阵
    coef_matrix = zeros(n_state*n_dist + n_dist, n_state*n_dist + n_input*n_dist);
    
    % 构造Sylvester方程部分
    for i = 1:n_state
        for j = 1:n_dist
            row_idx = (i-1)*n_dist + j;
            
            % Π*R项
            for k = 1:n_dist
                col_idx = (i-1)*n_dist + k;
                coef_matrix(row_idx, col_idx) = R(k, j);
            end
            
            % -A*Π项
            for k = 1:n_state
                for l = 1:n_dist
                    col_idx = (k-1)*n_dist + l;
                    coef_matrix(row_idx, col_idx) = coef_matrix(row_idx, col_idx) - A(i, k) * (l == j);
                end
            end
        end
    end
    
    % 填充输入项 -(1/β)*b*Γ
    for i = 1:n_state
        for j = 1:n_dist
            row_idx = (i-1)*n_dist + j;
            col_idx = n_state*n_dist + j;
            coef_matrix(row_idx, col_idx) = -(1/beta)*b(i);
        end
    end
    
    % 填充输出约束 C*Π = q
    for j = 1:n_dist
        row_idx = n_state*n_dist + j;
        for i = 1:n_state
            col_idx = (i-1)*n_dist + j;
            coef_matrix(row_idx, col_idx) = C(i);
        end
    end
    
    % 构造右侧向量
    rhs = zeros(n_state*n_dist + n_dist, 1);
    
    % 填充Sylvester方程右侧 (P)
    for i = 1:n_state
        for j = 1:n_dist
            row_idx = (i-1)*n_dist + j;
            if i == 1
                rhs(row_idx) = P(j);  % P行对应的值
            else
                rhs(row_idx) = 0;     % 其他行为0
            end
        end
    end
    
    % 填充输出约束右侧 (q)
    for j = 1:n_dist
        row_idx = n_state*n_dist + j;
        rhs(row_idx) = q(j);
    end
    
    % 求解线性方程组
    sol = coef_matrix \ rhs;
    
    % 提取解：Π和Γ
    Pi = zeros(n_state, n_dist);
    for i = 1:n_state
        for j = 1:n_dist
            idx = (i-1)*n_dist + j;
            Pi(i, j) = sol(idx);
        end
    end
    
    Gamma = zeros(1, n_dist);
    for j = 1:n_dist
        idx = n_state*n_dist + j;
        Gamma(j) = sol(idx);
    end
    
    % 验证解的正确性
    residual1 = Pi*R - A*Pi - (1/beta)*b*Gamma - [P; zeros(size(A,1)-1, size(P,2))];
    residual2 = C*Pi - q;
    
    if norm(residual1, 'fro') > 1e-10 || norm(residual2, 'fro') > 1e-10
        warning('解的残差较大，可能存在数值问题');
    end
end
