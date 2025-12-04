% 需要检查状态空间表示中是否含有延迟
clear; close all; clc;
%%
load('..\dataset\bltdWhiteNoise_ssmodel.mat');
load('..\dataset\ARMAX_SYSID_30303022.mat');
% 将ARMAX模型转换为状态空间模型
% ARMAX模型是一个idpoly对象，包含A, B, C多项式系数。
% 我们使用ssdata函数将其转换为状态空间表示。
% Af, Bf, Cf 分别是次能馈通路的系统矩阵，输入矩阵和输出矩阵
% Gf 是噪声模型矩阵
[Af, Bf, Cf, ~, Gf] = ssdata(ARMAXmodel.model);
fs = ARMAXmodel.fs;
Ts = 1/fs;

%% 耦合系统
% 增广状态矩阵
n = size(Af, 1);    % 原系统状态维度
p = size(Aw, 1);  % 干扰模型状态维度
A = blkdiag(Af, Aw);  % 块对角矩阵 [Af 0; 0 Aw]

% 增广输入矩阵
B = [Bf; zeros(p, size(Bf, 2))];  % 原系统控制输入通道
G = [zeros(n, size(Bw, 2)); Bw];  % 干扰输入通道

% 增广输出矩阵
C = [Cf , Cw];

%% 设计LQR控制器
% 状态权重矩阵 Q（半正定）
Q = C' * C;
% 输入权重矩阵 R（正定）
R = 1e-8;
% 离散 LQR 求解（状态反馈增益）
[K, ~, ~] = dlqr(A, B, Q, R);

% 计算闭环系统的状态矩阵
A_cl = A - B * K;
% 计算闭环系统的特征值
eig_cl = eig(A_cl);
% 检查闭环系统的稳定性
if all(abs(eig_cl) < 1)
    disp('闭环系统是稳定的（所有极点都在单位圆内）。');
else
    disp('闭环系统是不稳定的（存在极点在单位圆外）。');
end

%% 设计卡尔曼滤波器（状态估计器）
% 假设过程噪声协方差矩阵和测量噪声协方差矩阵
Qn = 2;  % 过程噪声协方差
Rn = 1e-10 * eye(size(C, 1));  % 测量噪声协方差

% 使用离散卡尔曼滤波器求解最优增益矩阵
[L, ~, ~] = dlqe(A, G, C, Qn, Rn);

% 显示卡尔曼滤波器增益
disp('卡尔曼滤波器增益矩阵 L:');
disp(L);

% 计算引入卡尔曼滤波器后的闭环系统状态矩阵
A_cl = A - B * K - L * C;
% 计算闭环系统的特征值
eig_cl = eig(A_cl);
% 显示闭环系统的极点
disp('引入卡尔曼滤波器后的闭环系统极点:');
disp(eig_cl);
% 检查闭环系统的稳定性
if all(abs(eig_cl) < 1)
    disp('引入卡尔曼滤波器后的闭环系统是稳定的（所有极点都在单位圆内）。');
else
    disp('引入卡尔曼滤波器后的闭环系统是不稳定的（存在极点在单位圆外）。');
end

%% 初始化参数
N = 20000;                    % 仿真步数
x = zeros(n + p, 1);    % 增广状态初始化
x_hat = zeros(n + p, 1);    % 状态估计初始化
x_history = zeros(n + p, N);  % 状态轨迹记录
x_hat_history = zeros(n + p, N);  % 状态估计轨迹记录
y_history = zeros(size(Cf, 1), N);  % 输出轨迹记录
u_history = zeros(1, N);  % 控制输入记录

%% 闭环控制仿真(一)
for k = 1:N
    % 生成过程噪声和测量噪声
    e = sqrt(Qn) * randn(size(G, 2), 1);  % 过程噪声
    v = chol(Rn)' * randn(size(C, 1), 1);  % 测量噪声

    % 状态反馈控制律
    u = -K * x_hat;  % 使用估计状态进行反馈控制
    u_history(k) = u;  % 记录控制输入

    % 增广系统动态更新
    x = A * x + B * u + G * e;  % 系统状态更新
    y = C * x + v;  % 输出测量

    % 卡尔曼滤波器更新（状态估计）
    x_hat = A * x_hat + B * u;  % 预测
    x_hat = x_hat + L * (y - C * x_hat);  % 校正

    % 记录状态和输出
    x_history(:, k) = x;
    x_hat_history(:, k) = x_hat;
    y_history(:, k) = y;
end

% 计算干扰信号和反噪声信号
d_history = Cw * x_history(n+1:end, :);  % 干扰信号
anti_history = Cf * x_history(1:end-p, :);  % 反噪声信号

figure;
subplot(2, 1, 1);
plot(1:N, u_history, 'DisplayName', '控制信号'); hold on;
plot(1:N, anti_history, 'DisplayName', '反噪声信号');
xlabel('样本n/samples'); ylabel('幅度');
legend;
grid on; 

subplot(2, 1, 2);
plot(1:N, d_history, 'DisplayName', '干扰信号d'); hold on;
plot(1:N, -anti_history, 'DisplayName', '反噪声信号');
plot(1:N, y_history, 'DisplayName', '残余信号 y');
xlabel('样本n/samples'); ylabel('幅度');
legend;
grid on;

% 计算声压级及降噪量
p0 = 20e-6;  % 参考声压 (Pa)
% 声压级计算
SPL_d = 20 * log10(rms(d_history(:)) / p0);
SPL_y = 20 * log10(rms(y_history(:)) / p0);
% 降噪量
attenuation = SPL_d - SPL_y;
fprintf('干扰信号声压级: %.2f dB\n', SPL_d);
fprintf('残余信号声压级: %.2f dB\n', SPL_y);
fprintf('降噪量: %.2f dB\n', attenuation);

%% 闭环控制仿真（二）
N = 20000;                    % 仿真步数
xf = zeros(n, 1);          % 原系统状态初始化
xw = zeros(p, 1);         % 干扰模型状态初始化
x_hat = zeros(n + p, 1);      % 增广状态估计初始化
xf_history = zeros(n, N);  % 原系统状态轨迹记录
xw_history = zeros(p, N); % 干扰模型状态轨迹记录
x_hat_history = zeros(n + p, N);  % 增广状态估计轨迹记录
y_history = zeros(size(Cf, 1), N);  % 输出轨迹记录
u_history = zeros(1, N);  % 控制输入记录

for k = 1:N
    % 生成干扰模型的过程噪声
    ew = sqrtm(Qn) * randn(size(Bw, 2), 1);  % 干扰过程噪声
    v = sqrtm(Rn) * randn(size(C, 1), 1);  % 测量噪声

    % 干扰模型动态更新
    xw = Aw * xw + Bw * ew;  % 干扰模型状态更新
    d = Cw * xw;  % 干扰信号

    % 原系统动态更新
    u = -K * x_hat;  % 使用估计状态进行反馈控制
    xf = Af * xf + Bf * u;  % 原系统状态更新
    y = Cf * xf + d + v;  % 系统输出（包含干扰和测量噪声）

    % 卡尔曼滤波器更新（状态估计）
    x_hat = A * x_hat + B * u;  % 预测
    x_hat = x_hat + L * (y - C * x_hat);  % 校正

    % 记录状态和输出
    xf_history(:, k) = xf;
    xw_history(:, k) = xw;
    x_hat_history(:, k) = x_hat;
    y_history(:, k) = y;
    u_history(k) = u;
end

% 计算干扰信号和反噪声信号
d_history = Cw * xw_history;  % 干扰信号
anti_history = Cf * xf_history;  % 反噪声信号

%%
figure;
subplot(2, 1, 1);
plot(1:N, u_history, 'DisplayName', '控制信号'); hold on;
plot(1:N, anti_history, 'DisplayName', '反噪声信号');
xlabel('样本n/samples'); ylabel('幅度');
legend;
grid on; 

subplot(2, 1, 2);
plot(1:N, d_history, 'DisplayName', '干扰信号d'); hold on;
plot(1:N, -anti_history, 'DisplayName', '反噪声信号');
plot(1:N, y_history, 'DisplayName', '残余信号 y');
xlabel('样本n/samples'); ylabel('幅度');
legend;
grid on;

% 计算声压级及降噪量
p0 = 20e-6;  % 参考声压 (Pa)
% 声压级计算
SPL_d = 20 * log10(rms(d_history(:)) / p0);
SPL_y = 20 * log10(rms(y_history(:)) / p0);
% 降噪量
attenuation = SPL_d - SPL_y;
fprintf('干扰信号声压级: %.2f dB\n', SPL_d);
fprintf('残余信号声压级: %.2f dB\n', SPL_y);
fprintf('降噪量: %.2f dB\n', attenuation);

%% 闭环控制仿真（三）Q 参数化自适应（未完成）
% 基于"闭环控制仿真（二）"的时序结构，加入 Q 参数化自适应项 u_Q(k)
% u(k) = -K*x_hat + u_Q(k)
% u_Q(k) = theta' * phi(k)
% phi(k) 来自带通滤波后的残差 r(k)=y(k)-y_hat(k)
% theta 由 RLS 更新
% Closed-loop A matrix for the plant when u = -K*x + u_q
% Note: the correct sign is A_cl_pl = A - B*K
A_cl_pl = A - B * K;
B_t12 = B;                         % input channel for the Q-module (u_q)
C_t12 = C;                         % output channel
D_t12 = zeros(size(C_t12,1), size(B_t12,2));

% Build state-space for T12 (transfer from u_q to output)
T12_ss = ss(A_cl_pl, B_t12, C_t12, D_t12, Ts);

% Convert to transfer-function form. For MIMO systems tf returns a matrix of SISO TFs.
T12_tf = tf(T12_ss);

% Robust extraction: for MIMO models extract each SISO element separately.
% Preallocate cell arrays
[ny, nu] = size(T12_tf);
numT12 = cell(ny, nu);
denT12 = cell(ny, nu);
for i = 1:ny
    for j = 1:nu
        % tfdata on a SISO element always returns vectors with 'v' option
        [ni, di] = tfdata(T12_tf(i,j), 'v');
        % store as row vectors
        numT12{i,j} = ni(:).';
        denT12{i,j} = di(:).';
    end
end

% If overall system is SISO, expose numeric vectors for convenience
if ny == 1 && nu == 1
    num_t12 = numT12{1,1};
    den_t12 = denT12{1,1};
end

N = 20000;                       % 仿真步数
xf = zeros(n, 1);                % 原系统状态
xw = zeros(p, 1);                % 干扰模型状态
x_hat = zeros(n + p, 1);         % 增广状态估计

% ---------- Q 模块参数 ----------
nq = 2;                          % Q 模块阶数
theta = zeros(nq, 1);            % 参数初值
P = 1e3 * eye(nq);               % 协方差初值
lambda = 0.98;                   % 遗忘因子

% 带通滤波器 H(z)：限制自适应频带
bp_low = 50; bp_high = 300;
Wn = [bp_low, bp_high] / (fs/2);
[b_h, a_h] = butter(4, Wn, 'bandpass');

% ---------- 数据缓存 ----------
y_history = zeros(size(Cf, 1), N);
u_history = zeros(1, N);
theta_history = zeros(nq, N);
d_history = zeros(1, N);
anti_history = zeros(1, N);
psi_buf = zeros(nq, 1);  % 自适应回归向量缓冲区

% ---------- T12 预测缓冲（用于补偿 u_q 对输出的影响） ----------
uq = 0;
num_t12 = num_t12(:).';
den_t12 = den_t12(:).';
% 构造合成滤波器（将 T12 与带通 H 串联），用于对残差 r 进行带限处理
b_comb = conv(num_t12, b_h);
a_comb = conv(den_t12, a_h);
% 初始化 IIR 滤波器历史缓存
u_comb_hist = zeros(1, max(1, length(b_comb)-1));     % 合成滤波器输入历史
y_comb_hist = zeros(1, max(1, length(a_comb)-1));     % 合成滤波器输出历史

for k = 1:N
    % 生成干扰模型的过程噪声
    ew = sqrtm(Qn) * randn(size(Bw, 2), 1);  % 干扰过程噪声
    v = sqrtm(Rn) * randn(size(C, 1), 1);

    % --- 1. 干扰模型动态更新 ---
    xw = Aw * xw + Bw * ew;                  % 干扰状态更新
    d = Cw * xw;                             % 干扰信号

    
    % --- 3. 残差信号与带通滤波 ---
    y_hat = C * x_hat;                       % 估计输出
    r = y - y_hat;                           % 残差
    
    % 使用自定义 IIR_filter 替代 filter 函数
    [psi_k, u_comb_hist, y_comb_hist] = IIR_filter(b_comb, a_comb, r, u_comb_hist, y_comb_hist);
    psi_buf = [psi_k; psi_buf(1:end-1)];
    % 根据公式，phi(k) = -[psi(k); psi(k-1); ...]
    phi = -psi_buf;
        
    % --- 自适应 Q 输出 ---
    u_q = theta' * phi;
    % u_q = 0;

    % --- 控制输入 ---
    u = -K * x_hat + u_q;
    % --- 原系统状态更新 ---
    xf = Af * xf + Bf * u;
    % --- 输出测量（包含干扰 + 测量噪声）---
    y = Cf * xf + d + v;                     % 实际输出
    
    % --- 7. 增广估计更新（预测 + 校正）---
    x_hat = A * x_hat + B * u; % 预测
    x_hat = x_hat + L * (y - C * x_hat); % 校正
    
    % --- 8. RLS 参数更新 ---
    % 利用公式 epsilon = ( y - (T12*Q - Q*T12)*r ) / denom
    denom = 1 + phi' * P * phi;
    epsilon = (y + 0) / denom; 
    theta = theta + P * phi * epsilon;
    P = (1/lambda) * (P - (P * (phi * phi') * P) / denom);
    
    % --- 9. 记录 ---
    u_history(k) = u;
    y_history(:, k) = y;
    d_history(:, k) = d;
    anti_history(:, k) = Cf * xf;
    theta_history(:, k) = theta;
end

% ---------- 绘图 ----------
figure;
subplot(3,1,1);
plot(d_history, 'DisplayName', '干扰信号 d'); hold on;
plot(-anti_history, '--', 'DisplayName', '反噪声');
plot(y_history, 'DisplayName', '系统输出 y');
xlabel('样本'); ylabel('幅值'); legend; grid on; title('输出与干扰信号');

subplot(3,1,2);
plot(u_history, 'DisplayName', '控制输入 u');
xlabel('样本'); ylabel('幅度'); legend; grid on; title('控制信号（含自适应项）');

subplot(3,1,3);
plot(theta_history');
xlabel('样本'); ylabel('\theta_i'); grid on;
title('Q 参数 \theta 收敛轨迹');

% ---------- 计算声压级及降噪量 ----------
p0 = 20e-6;
SPL_d = 20 * log10(rms(d_history(:)) / p0);
SPL_y = 20 * log10(rms(y_history(:)) / p0);
attenuation = SPL_d - SPL_y;
fprintf('干扰信号声压级: %.2f dB\n', SPL_d);
fprintf('残余信号声压级: %.2f dB\n', SPL_y);
fprintf('降噪量: %.2f dB\n', attenuation);

%%
function y = FIR(filter, X)
    % Robust FIR: ensure column vectors, pad X with zeros if shorter than filter
    b = filter(:);
    x = X(:);
    Lb = length(b);
    if length(x) < Lb
        x = [x; zeros(Lb - length(x), 1)];
    end
    % return scalar dot-product
    y = b.' * x(1:Lb);
end

function [y_out, u_hist_new, y_hist_new] = IIR_filter(num, den, u_in, u_hist, y_hist)
    % 实现IIR滤波器: H(z) = num(z)/den(z)
    % 差分方程: den(z)Y(z) = num(z)U(z)
    % 即: a0*y[n] + a1*y[n-1] + ... = b0*u[n] + b1*u[n-1] + ...
    
    % 更新输入历史（新样本加入队首，长度与输入一致）
    u_hist_new = [u_in, u_hist(1:end-1)];
    % FIR计算时临时补零或截断
    if length(u_hist_new) < length(num)
        u_hist_pad = [u_hist_new, zeros(1, length(num)-length(u_hist_new))];
    else
        u_hist_pad = u_hist_new(1:length(num));
    end
    if length(y_hist) < length(den)-1
        y_hist_pad = [y_hist, zeros(1, length(den)-1-length(y_hist))];
    else
        y_hist_pad = y_hist(1:length(den)-1);
    end
    % 计算输出：y[n] = (num部分 - den[1:end]部分) / den[0]
    num_part = FIR(num, u_hist_pad);
    den_part = FIR(den(2:end), y_hist_pad);  % 不包括den(1)，因为那是当前输出项
    y_out = (num_part - den_part) / den(1);
    % 更新输出历史（长度与输入一致）
    y_hist_new = [y_out, y_hist(1:end-1)];
end