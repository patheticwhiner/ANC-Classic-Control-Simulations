clear; close all; clc;

% 系统参数定义
a1 = 0.5;       % 系统参数
b2 = 0.5;       % 系统参数
beta = 0.5;     % 系统参数
omega1 = 1;     % 扰动频率1
omega2 = 0.5;   % 扰动频率2
omega3 = 4;     % 扰动频率3（仅用于区间1）

% 定义系统矩阵
A = [a1, 1; 0, 0];  % 系统状态矩阵
b = [1; b2];        % 输入矩阵
C = [1, 0];         % 输出矩阵

fprintf('=== 系统矩阵 ===\n');
fprintf('A = \n');
disp(A);
fprintf('b = \n');
disp(b);
fprintf('C = \n');
disp(C);
fprintf('beta = %.2f\n\n', beta);

% 构建扰动的状态空间矩阵R（6x6）
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

% 解调节器方程: ΠR = AΠ + 1/β*b*Γ + P, q = CΠ
% 对三个不同区间分别求解

% ==== 区间1: [0, 50) - w1包含omega1和omega3, w2包含omega2 ====
fprintf('=== 区间1: [0, 50) ===\n');
fprintf('扰动: w1 = sin(%.2f*t) + 0.2*sin(%.2f*t), w2 = sin(%.2f*t)\n', omega1, omega3, omega2);

% 输出映射矩阵
P1 = [1, 0, 0.2, 0, 0, 0];  % w1 = x1 + 0.2*x3
q1 = [0, 0, 0, 0, 1, 0];    % w2 = x5

fprintf('扰动系统状态空间矩阵 R:\n');
disp(R);
fprintf('扰动输出矩阵 P1 (对应w1):\n');
disp(P1);
fprintf('扰动输出矩阵 q1 (对应w2):\n');
disp(q1);

% 直接使用函数求解调节器方程
[Pi1, Gamma1] = solveRegulatorEquation(A, b, C, R, P1, q1, beta);

% 显示解
fprintf('\n区间1解:\n');
fprintf('Π1 = \n');
disp(Pi1);
fprintf('Γ1 = \n');
disp(Gamma1);

% 验证解的正确性
verifyRegulatorSolution(A, b, C, R, P1, q1, beta, Pi1, Gamma1);

% ==== 区间2: [50, 100) - w1只包含omega1, w2包含omega2 ====
fprintf('\n=== 区间2: [50, 100) ===\n');
fprintf('扰动: w1 = sin(%.2f*t), w2 = sin(%.2f*t)\n', omega1, omega2);

% 输出映射矩阵（区间2）
P2 = [1, 0, 0, 0, 0, 0];    % w1 = x1
q2 = [0, 0, 0, 0, 1, 0];    % w2 = x5

fprintf('扰动输出矩阵 P2 (对应w1):\n');
disp(P2);
fprintf('扰动输出矩阵 q2 (对应w2):\n');
disp(q2);

% 直接使用函数求解调节器方程
[Pi2, Gamma2] = solveRegulatorEquation(A, b, C, R, P2, q2, beta);

% 显示解
fprintf('\n区间2解:\n');
fprintf('Π2 = \n');
disp(Pi2);
fprintf('Γ2 = \n');
disp(Gamma2);

% 验证解的正确性
verifyRegulatorSolution(A, b, C, R, P2, q2, beta, Pi2, Gamma2);

% ==== 区间3: [100, 150] - w1只包含omega1, w2为0 ====
fprintf('\n=== 区间3: [100, 150] ===\n');
fprintf('扰动: w1 = sin(%.2f*t), w2 = 0\n', omega1);

% 输出映射矩阵（区间3）
P3 = [1, 0, 0, 0, 0, 0];    % w1 = x1
q3 = [0, 0, 0, 0, 0, 0];    % w2 = 0

fprintf('扰动输出矩阵 P3 (对应w1):\n');
disp(P3);
fprintf('扰动输出矩阵 q3 (对应w2):\n');
disp(q3);

% 直接使用函数求解调节器方程
[Pi3, Gamma3] = solveRegulatorEquation(A, b, C, R, P3, q3, beta);

% 显示解
fprintf('\n区间3解:\n');
fprintf('Π3 = \n');
disp(Pi3);
fprintf('Γ3 = \n');
disp(Gamma3);

% 验证解的正确性
verifyRegulatorSolution(A, b, C, R, P3, q3, beta, Pi3, Gamma3);

% 对比三种情况的解
fprintf('\n=== 三种情况的Π矩阵对比 ===\n');
fprintf('Π1(区间1) - Π2(区间2) 的差异范数: %.6e\n', norm(Pi1 - Pi2, 'fro'));
fprintf('Π2(区间2) - Π3(区间3) 的差异范数: %.6e\n', norm(Pi2 - Pi3, 'fro'));
fprintf('Π1(区间1) - Π3(区间3) 的差异范数: %.6e\n\n', norm(Pi1 - Pi3, 'fro'));

fprintf('=== 三种情况的Γ向量对比 ===\n');
fprintf('Γ1(区间1) - Γ2(区间2) 的差异范数: %.6e\n', norm(Gamma1 - Gamma2, 'fro'));
fprintf('Γ2(区间2) - Γ3(区间3) 的差异范数: %.6e\n', norm(Gamma2 - Gamma3, 'fro'));
fprintf('Γ1(区间1) - Γ3(区间3) 的差异范数: %.6e\n', norm(Gamma1 - Gamma3, 'fro'));

% 绘制理论轨迹图
t = 0:0.01:150;
w = zeros(6, length(t));
w(:,1) = [0; omega1; 0; omega3; 0; omega2]; % 初始状态

% 使用RK4求解扰动系统
for i = 1:length(t)-1
    dt = t(i+1) - t(i);
    k1 = R * w(:,i);
    k2 = R * (w(:,i) + 0.5*dt*k1);
    k3 = R * (w(:,i) + 0.5*dt*k2);
    k4 = R * (w(:,i) + dt*k3);
    w(:,i+1) = w(:,i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

% 计算各区间的理论轨迹
x_ideal = zeros(2, length(t));
u_ideal = zeros(1, length(t));
w1_out = zeros(1, length(t));
w2_out = zeros(1, length(t));

for i = 1:length(t)
    if t(i) < 50
        % 区间1
        x_ideal(:,i) = Pi1 * w(:,i);
        u_ideal(i) = Gamma1 * w(:,i);
        w1_out(i) = P1 * w(:,i);
        w2_out(i) = q1 * w(:,i);
    elseif t(i) < 100
        % 区间2
        x_ideal(:,i) = Pi2 * w(:,i);
        u_ideal(i) = Gamma2 * w(:,i);
        w1_out(i) = P2 * w(:,i);
        w2_out(i) = q2 * w(:,i);
    else
        % 区间3
        x_ideal(:,i) = Pi3 * w(:,i);
        u_ideal(i) = Gamma3 * w(:,i);
        w1_out(i) = P3 * w(:,i);
        w2_out(i) = q3 * w(:,i);
    end
end

% 绘制理论轨迹
figure;
subplot(2,2,1);
plot(t, x_ideal(1,:), 'b-', t, -w2_out, 'r--', 'LineWidth', 1.5);
title('理论输出轨迹');
xlabel('时间 (s)');
ylabel('输出');
legend('系统输出 x_1', '参考信号 -w_2');
grid on;

subplot(2,2,2);
plot(t, u_ideal, 'k-', 'LineWidth', 1.5);
title('理论控制输入');
xlabel('时间 (s)');
ylabel('u');
grid on;

subplot(2,2,3);
plot(t, w1_out, 'g-', t, w2_out, 'm-', 'LineWidth', 1.5);
title('扰动信号');
xlabel('时间 (s)');
ylabel('扰动');
legend('w_1', 'w_2');
grid on;

subplot(2,2,4);
plot(t, x_ideal(1,:) - (-w2_out), 'r-', 'LineWidth', 1.5);
title('理论跟踪误差 e = x_1 - (-w_2)');
xlabel('时间 (s)');
ylabel('误差');
grid on;

% 打印扰动系统的特征值
fprintf('\n=== 扰动系统特征值 ===\n');
eig_R = eig(R);
fprintf('R的特征值:\n');
disp(eig_R);

% 选取几个特殊时间点观察稳态解
fprintf('\n=== 特殊时间点的系统状态 ===\n');
special_times = [25, 75, 125]; % 三个区间中间点

for i = 1:length(special_times)
    t_idx = find(t >= special_times(i), 1);
    fprintf('t = %.0f:\n', special_times(i));
    fprintf('w(t) = [%.4f; %.4f; %.4f; %.4f; %.4f; %.4f]\n', w(1,t_idx), w(2,t_idx), w(3,t_idx), w(4,t_idx), w(5,t_idx), w(6,t_idx));
    fprintf('w1(t) = %.4f, w2(t) = %.4f\n', w1_out(t_idx), w2_out(t_idx));
    fprintf('x(t) = [%.4f; %.4f]\n', x_ideal(1,t_idx), x_ideal(2,t_idx));
    fprintf('u(t) = %.4f\n\n', u_ideal(t_idx));
end

% 使用MATLAB的solve函数尝试求解符号解
fprintf('\n=== 使用符号解法求解简化的调节器方程 ===\n');
try
    % 简化问题为2维扰动系统
    syms s1 s2 s3 s4 g1 g2 real
    Pi_sym = [s1 s2; s3 s4];
    Gamma_sym = [g1 g2];
    
    % 简化的R矩阵（只考虑omega1）
    R_simple = [0 1; -omega1^2 0];
    P_simple = [1 0];
    
    % 构建方程
    eq1 = Pi_sym*R_simple - A*Pi_sym - (1/beta)*b*Gamma_sym == [P_simple; 0 0];
    eq2 = C*Pi_sym == [1 0];
    
    % 联立求解
    solutions = solve([eq1(:); eq2(:)], [s1 s2 s3 s4 g1 g2]);
    
    % 显示符号解
    Pi_sol = [solutions.s1 solutions.s2; solutions.s3 solutions.s4];
    Gamma_sol = [solutions.g1 solutions.g2];
    
    fprintf('符号解 (只考虑omega1):\n');
    fprintf('Π = \n');
    disp(Pi_sol);
    fprintf('Γ = \n');
    disp(Gamma_sol);
catch
    fprintf('符号解法求解失败，可能是方程过于复杂。\n');
end

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
end

function verifyRegulatorSolution(A, b, C, R, P, q, beta, Pi, Gamma)
    % 验证解是否满足调节器方程
    residual1 = Pi*R - A*Pi - (1/beta)*b*Gamma;
    residual2 = C*Pi - q;
    
    % 计算残差
    fprintf('验证方程1 (Π*R - A*Π - 1/β*b*Γ = P)，残差范数: %.6e\n', norm(residual1 - [P; zeros(1,size(P,2))], 'fro'));
    fprintf('验证方程2 (C*Π = q)，残差范数: %.6e\n', norm(residual2 - q, 'fro'));
    
    % 检查残差是否在可接受范围内
    if norm(residual1 - [P; zeros(1,size(P,2))], 'fro') < 1e-10 && norm(residual2 - q, 'fro') < 1e-10
        fprintf('解验证成功！残差在数值精度范围内。\n');
    else
        fprintf('警告：解的残差超出预期，可能存在数值问题。\n');
    end
end
