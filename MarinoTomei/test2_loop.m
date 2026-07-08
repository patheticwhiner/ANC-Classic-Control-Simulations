clear; close all; clc;
% 系统参数与初始条件
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
x = [0.5; 0];          % 初始状态 [x1, x2]
xi1 = zeros(4,1);     % 滤波器ξ1初始状态
xi2 = zeros(4,1);     % 滤波器ξ2初始状态
chi_hat = zeros(4,1);  % 估计器χ初始状态
theta_hat = [0.1; 0.1];% 参数估计初始值 [θ1_hat, θ2_hat]
wc = zeros(2*m+1,1);  % 参考控制器状态初始条件 (设置第一个元素为初始e以匹配初始控制)
wc(1) = -x(1); % 设置初始值使得初始u_r匹配原方程

% 仿真时间设置
t_start = 0;
t_end = 150;
dt = 0.01;        % 时间步长
num_steps = (t_end - t_start) / dt;

% 定义赫尔维茨矩阵D（确保特征值在左半平面）
D = [-d(1), 1, 0, 0;
     -d(2), 0, 1, 0;
     -d(3), 0, 0, 1;
     -d(4), 0, 0, 0]; % 示例结构，可根据实际调整

% 初始化结果数组
t = linspace(t_start, t_end, num_steps+1)';
x1_history = zeros(size(t));
x2_history = zeros(size(t));
xi1_history = zeros(length(t), 4);
xi2_history = zeros(length(t), 4);
chi_hat_history = zeros(length(t), 4);
theta1_hat_history = zeros(size(t));
theta2_hat_history = zeros(size(t));
wc_history = zeros(length(t), 2*m+1);
u_history = zeros(size(t));
u_r_history = zeros(size(t));
e_output_history = zeros(size(t));
e_input_history = zeros(size(t));
w1_history = zeros(size(t));
w2_history = zeros(size(t));

% 投影算子参数
r_omega = 2;      % 参数约束半径
epsilon_r = 0.1;  % 投影算子平滑参数

% 循环仿真
for i = 1:num_steps+1
    current_t = t(i);
    
    % 计算当前时刻的扰动
    w1 = get_w1(current_t, omega1, omega3);
    w2 = get_w2(current_t, omega2);
    
    % 计算参考信号和调节误差
    y_r = -w2;
    e = x(1) - y_r;
    
    % 计算控制输入
    mu1 = xi1(1);
    mu2 = xi2(1);
    u = -k*e - chi_hat(1) - mu1*theta_hat(1) - mu2*theta_hat(2);
    
    % 计算参考控制器输入
    u_r = wc(1);
    
    % 存储当前状态和计算结果
    x1_history(i) = x(1);
    x2_history(i) = x(2);
    xi1_history(i,:) = xi1';
    xi2_history(i,:) = xi2';
    chi_hat_history(i,:) = chi_hat';
    theta1_hat_history(i) = theta_hat(1);
    theta2_hat_history(i) = theta_hat(2);
    wc_history(i,:) = wc';
    u_history(i) = u;
    u_r_history(i) = u_r;
    e_output_history(i) = e;
    e_input_history(i) = u - u_r;
    w1_history(i) = w1;
    w2_history(i) = w2;
    
    % 如果是最后一步，跳出循环
    if i == num_steps+1
        break;
    end
    
    % 计算系统状态导数
    dx1dt = x(2) + a1*x(1) + (1/beta)*u + w1;
    dx2dt = (1/beta)*b2*u;
    
    % 滤波器动态方程
    dxi1dt = D*xi1 + [0; u; 0; 0];
    dxi2dt = D*xi2 + [0; 0; 0; u];
    
    % 估计器χ动态方程
    dchi_hatdt = D*chi_hat - [d(1); d(2); d(3); d(4)]*u;
    
    % 参数更新律（使用投影算子）
    phi1 = mu1*e;
    phi2 = mu2*e;
    theta1_hat_dot = g * projection(phi1, theta_hat(1), r_omega, epsilon_r);
    theta2_hat_dot = g * projection(phi2, theta_hat(2), r_omega, epsilon_r);
    
    % 构造R_c矩阵并计算wc的导数
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
    
    % 使用欧拉法更新状态
    x(1) = x(1) + dx1dt * dt;
    x(2) = x(2) + dx2dt * dt;
    xi1 = xi1 + dxi1dt * dt;
    xi2 = xi2 + dxi2dt * dt;
    chi_hat = chi_hat + dchi_hatdt * dt;
    theta_hat(1) = theta_hat(1) + theta1_hat_dot * dt;
    theta_hat(2) = theta_hat(2) + theta2_hat_dot * dt;
    wc = wc + dwcdt * dt;
end

% 绘制图1：调节误差与参数估计
figure;
subplot(2,2,1);
plot(t, e_output_history, 'b');
title('输出调节误差');
xlabel('时间 (s)');
ylabel('e_{output}');

subplot(2,2,2);
plot(t, e_input_history, 'r');
title('输入调节误差');
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

% 绘制图2：系统输出、扰动与输入参考
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
plot(t, u_r_history, 'b');
title('输入参考信号');
xlabel('时间 (s)');
ylabel('参考输入');
legend('u_r');

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

% 投影算子的实现
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
