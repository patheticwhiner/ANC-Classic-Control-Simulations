clear; close all; clc;

% 系统参数定义
omega1 = 1;     % 扰动频率1
omega2 = 0.5;   % 扰动频率2
omega3 = 4;     % 扰动频率3（仅用于t < 50s）

% 时间设置
tspan = [0 150];
dt = 0.01;
t = tspan(1):dt:tspan(2);

% 分段定义三个时间区间
% 区间1: [0, 50) - w1包含omega1和omega3, w2包含omega2
% 区间2: [50, 100) - w1只包含omega1, w2包含omega2
% 区间3: [100, 150] - w1只包含omega1, w2为0

% ==== 正确的状态空间建模方法 ====
% 对于正弦信号 sin(omega*t)，其状态空间模型为：
% x1' = x2
% x2' = -omega^2 * x1
% y = x1
% 初始条件为 x1(0) = 0, x2(0) = omega

% 创建合并的状态空间模型
% 状态向量 [x1; x2; x3; x4; x5; x6]
% x1,x2: 对应omega1的正弦分量
% x3,x4: 对应omega3的正弦分量
% x5,x6: 对应omega2的正弦分量

% 构建完整的状态空间矩阵R（6x6）
% 初始条件下，各频率对应的正弦信号均为sin(omega*t)
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

% 初始条件 - 假设所有信号在t=0时为0，导数为对应的角频率
x_init = [0; omega1; 0; omega3; 0; omega2];

% 创建结果存储数组
x = zeros(6, length(t));
w1 = zeros(size(t));
w2 = zeros(size(t));

% 赋初值
x(:, 1) = x_init;
w1(1) = P1 * x_init;  % 区间1的w1输出
w2(1) = q1 * x_init;  % 区间1的w2输出

% 使用四阶龙格-库塔法求解状态空间模型，提高精度
for i = 1:length(t)-1
    current_t = t(i);
    h = dt;
    
    % 标准RK4步骤
    k1 = R * x(:, i);
    k2 = R * (x(:, i) + 0.5 * h * k1);
    k3 = R * (x(:, i) + 0.5 * h * k2);
    k4 = R * (x(:, i) + h * k3);
    
    x(:, i+1) = x(:, i) + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
    
    % 根据当前时间确定使用哪个输出矩阵
    if current_t + dt < 50
        % 区间1: [0, 50)
        w1(i+1) = P1 * x(:, i+1);
        w2(i+1) = q1 * x(:, i+1);
    elseif current_t + dt < 100
        % 区间2: [50, 100)
        w1(i+1) = P2 * x(:, i+1);
        w2(i+1) = q2 * x(:, i+1);
    else
        % 区间3: [100, 150]
        w1(i+1) = P3 * x(:, i+1);
        w2(i+1) = q3 * x(:, i+1);
    end
end

% 生成直接计算的参考信号用于比较
w1_direct = zeros(size(t));
w2_direct = zeros(size(t));

for i = 1:length(t)
    current_t = t(i);
    
    % 直接使用分段函数计算w1
    if current_t < 50
        w1_direct(i) = sin(omega1*current_t) + 0.2*sin(omega3*current_t);
    else
        w1_direct(i) = sin(omega1*current_t);
    end
    
    % 直接使用分段函数计算w2
    if current_t < 100
        w2_direct(i) = sin(omega2*current_t);
    else
        w2_direct(i) = 0;
    end
end

% 计算误差
error_w1 = w1 - w1_direct;
error_w2 = w2 - w2_direct;
max_error_w1 = max(abs(error_w1));
max_error_w2 = max(abs(error_w2));
rms_error_w1 = sqrt(mean(error_w1.^2));
rms_error_w2 = sqrt(mean(error_w2.^2));

% 绘制比较图
figure;
subplot(2,2,1);
plot(t, w1, 'b-', t, w1_direct, 'r--', 'LineWidth', 1.5);
title('w1(t) = Pw(t) 比较');
xlabel('时间 (s)');
ylabel('幅值');
legend('状态空间模型', '直接计算');
grid on;

subplot(2,2,2);
plot(t, w2, 'b-', t, w2_direct, 'r--', 'LineWidth', 1.5);
title('w2(t) = -qw 比较');
xlabel('时间 (s)');
ylabel('幅值');
legend('状态空间模型', '直接计算');
grid on;

subplot(2,2,3);
plot(t, error_w1, 'k-', 'LineWidth', 1.5);
title(['w1(t) 误差 (最大误差: ', num2str(max_error_w1), ', RMS: ', num2str(rms_error_w1), ')']);
xlabel('时间 (s)');
ylabel('误差');
grid on;

subplot(2,2,4);
plot(t, error_w2, 'k-', 'LineWidth', 1.5);
title(['w2(t) 误差 (最大误差: ', num2str(max_error_w2), ', RMS: ', num2str(rms_error_w2), ')']);
xlabel('时间 (s)');
ylabel('误差');
grid on;

% 分别绘制状态轨迹
figure;
subplot(3,1,1);
plot(t, x(1,:), 'b-', t, x(2,:), 'r--', 'LineWidth', 1.5);
title('状态轨迹: omega1相关状态');
xlabel('时间 (s)');
ylabel('状态值');
legend('x1 (位置)', 'x2 (速度)');
grid on;

subplot(3,1,2);
plot(t, x(3,:), 'b-', t, x(4,:), 'r--', 'LineWidth', 1.5);
title('状态轨迹: omega3相关状态');
xlabel('时间 (s)');
ylabel('状态值');
legend('x3 (位置)', 'x4 (速度)');
grid on;

subplot(3,1,3);
plot(t, x(5,:), 'b-', t, x(6,:), 'r--', 'LineWidth', 1.5);
title('状态轨迹: omega2相关状态');
xlabel('时间 (s)');
ylabel('状态值');
legend('x5 (位置)', 'x6 (速度)');
grid on;

% 输出误差分析结果到控制台
fprintf('==== 误差分析 ====\n');
fprintf('w1(t) 最大误差: %e\n', max_error_w1);
fprintf('w1(t) RMS误差: %e\n', rms_error_w1);
fprintf('w2(t) 最大误差: %e\n', max_error_w2);
fprintf('w2(t) RMS误差: %e\n', rms_error_w2);

% 验证RMS误差的波形
figure;
t_zoom = t(t>=45 & t<=55);  % 区间1到区间2的过渡
w1_zoom = w1(t>=45 & t<=55);
w1_direct_zoom = w1_direct(t>=45 & t<=55);
error_w1_zoom = error_w1(t>=45 & t<=55);

subplot(2,1,1);
plot(t_zoom, w1_zoom, 'b-', t_zoom, w1_direct_zoom, 'r--', 'LineWidth', 1.5);
title('区间过渡处w1(t)比较 (t=45~55s)');
xlabel('时间 (s)');
ylabel('幅值');
legend('状态空间模型', '直接计算');
grid on;

subplot(2,1,2);
plot(t_zoom, error_w1_zoom, 'k-', 'LineWidth', 1.5);
title('区间过渡处w1(t)误差');
xlabel('时间 (s)');
ylabel('误差');
grid on;

% 合并的状态空间模型分析
fprintf('\n==== 合并状态空间模型分析 ====\n');
fprintf('状态空间矩阵R:\n');
disp(R);
fprintf('区间1 [0, 50) 输出矩阵P1和q1:\n');
disp(P1);
disp(q1);
fprintf('区间2 [50, 100) 输出矩阵P2和q2:\n');
disp(P2);
disp(q2);
fprintf('区间3 [100, 150] 输出矩阵P3和q3:\n');
disp(P3);
disp(q3);