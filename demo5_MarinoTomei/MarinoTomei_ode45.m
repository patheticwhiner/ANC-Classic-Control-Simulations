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
x0 = [0.5; 0];          % 初始状态 [x1, x2]
xi1_0 = zeros(4,1);     % 滤波器ξ1初始状态
xi2_0 = zeros(4,1);     % 滤波器ξ2初始状态
chi_hat0 = zeros(4,1);  % 估计器χ初始状态
theta_hat0 = [0.1; 0.1];% 参数估计初始值 [θ1_hat, θ2_hat]
wc_0 = zeros(2*m+1,1);  % 参考控制器状态初始条件 (设置第一个元素为初始e以匹配初始控制)
wc_0(1) = -x0(1); % 设置初始值使得初始u_r匹配原方程

% 仿真时间设置
tspan = [0, 150];       % 仿真时间范围：0到150秒

% 定义赫尔维茨矩阵D（确保特征值在左半平面）
D = [-d(1), 1, 0, 0;
     -d(2), 0, 1, 0;
     -d(3), 0, 0, 1;
     -d(4), 0, 0, 0]; % 示例结构，可根据实际调整

% 定义ODE求解选项
options = odeset('RelTol', 1e-6, 'AbsTol', 1e-8);

% 运行仿真
[t, X] = ode45(@(t, X) system_dynamics(t, X, a1, b2, beta, omega1, omega2, omega3, k, g, D, d), tspan, ...
              [x0; xi1_0; xi2_0; chi_hat0; theta_hat0; wc_0], options);

% 提取结果
x1 = X(:,1);            % 系统状态x1
x2 = X(:,2);            % 系统状态x2
xi1 = X(:,3:6);         % 滤波器ξ1状态
xi2 = X(:,7:10);        % 滤波器ξ2状态
chi_hat = X(:,11:14);   % 估计器χ状态
theta1_hat = X(:,15);   % θ1估计
theta2_hat = X(:,16);   % θ2估计
wc = X(:,17:end);       % 参考控制器状态

% 提取扰动信号
w1 = arrayfun(@(t) get_w1(t, omega1, omega3), t);
w2 = arrayfun(@(t) get_w2(t, omega2), t);
y_r = -w2;              % 参考信号 y_r = -w2

% 导出结果
e_output = x1 - y_r;    % 输出调节误差 e = y - y_r = x1 - (-w2)
gamma_c = [1, zeros(1, 2*m)]; % 定义gamma_c = [1 0 ... 0]
% u_r = wc * gamma_c'; % 参考输入 u_r = gamma_c * w_c
u_r = wc(:,1);
u = -k*e_output - chi_hat(:,1) - xi1(:,1).*theta1_hat - xi2(:,1).*theta2_hat; % 控制输入u
e_input = u - u_r;      % 输入调节误差（假设理想输入为 -k*e）

% 绘制图1：调节误差与参数估计
figure;
subplot(2,2,1);
plot(t, e_output, 'b');
title('输出调节误差');
xlabel('时间 (s)');
ylabel('e_{output}');

subplot(2,2,2);
plot(t, e_input, 'r');
title('输入调节误差');
xlabel('时间 (s)');
ylabel('e_{input}');

subplot(2,2,3);
plot(t, theta1_hat, 'b', t, ones(size(t))*(omega1^2 + omega2^2), 'k--');
title('θ_1估计');
xlabel('时间 (s)');
legend('估计值', '真实值','location','best');

subplot(2,2,4);
plot(t, theta2_hat, 'r', t, ones(size(t))*(omega1^2 * omega2^2), 'k--');
title('θ_2估计');
xlabel('时间 (s)');
legend('估计值', '真实值','location','best');

% 绘制图2：系统输出、扰动与输入参考
figure;
subplot(2,2,1);
plot(t, x1, 'b', t, y_r, 'r--');
title('系统输出与参考信号');
xlabel('时间 (s)');
ylabel('输出');
legend('y = x1', 'y_r = -w2');

subplot(2,2,2);
plot(t, u, 'k');
title('控制输入');
xlabel('时间 (s)');
ylabel('u');

subplot(2,2,3);
plot(t, w1, 'g', t, w2, 'm');
title('扰动信号');
xlabel('时间 (s)');
ylabel('扰动');
legend('w1', 'w2');

subplot(2,2,4);
plot(t, u_r, 'b');
title('输入参考信号');
xlabel('时间 (s)');
ylabel('参考输入');
legend('u_r');

% 定义系统动态方程
function dXdt = system_dynamics(t, X, a1, b2, beta, omega1, omega2, omega3, k, g, D, d)
    % 分解状态向量
    x = X(1:2);                 % 系统状态 [x1; x2]
    xi1 = X(3:6);               % 滤波器ξ1状态
    xi2 = X(7:10);              % 滤波器ξ2状态
    chi_hat = X(11:14);         % 估计器χ状态
    theta_hat = X(15:16);       % 参数估计 [θ1_hat; θ2_hat]
    
    % 计算m值（基于omega1和omega2）
    m = 2;  % 对于两个频率的情况
    
    % 提取参考控制器状态
    wc = X(17:end);             % 参考控制器状态 w_c
    
    % 计算扰动w1和w2
    w1 = get_w1(t, omega1, omega3);
    w2 = get_w2(t, omega2);
    
    % 计算跟踪误差e = y - y_r = x1 - (-w2)
    e = x(1) - (-w2);
    
    % 构造R_c矩阵
    % 计算theta_1和theta_2的真实值
    theta1_true = omega1^2 + omega2^2;
    theta2_true = omega1^2 * omega2^2;
    Rc = zeros(2*m+1);
    Rc(1,2) = 1;
    Rc(2,1) = -theta1_true;
    Rc(2,3) = 1;
    Rc(3,4) = 1;
    Rc(4,1) = -theta2_true;
    Rc(4,5) = 1;
    
    % 计算w_c的导数
    dwcdt = Rc * wc;
    
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
    
    % 参数更新律（使用投影算子）
    r_omega = 2;      % 参数约束半径
    epsilon_r = 0.1;  % 投影算子平滑参数
    
    % 对每个参数单独应用投影算子
    phi1 = mu1*e;
    phi2 = mu2*e;
    theta1_hat_dot = g * projection(phi1, theta_hat(1), r_omega, epsilon_r);
    theta2_hat_dot = g * projection(phi2, theta_hat(2), r_omega, epsilon_r);
    
    % 整合所有导数
    dXdt = [dx1dt; dx2dt; dxi1dt; dxi2dt; dchi_hatdt; theta1_hat_dot; theta2_hat_dot; dwcdt];
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

