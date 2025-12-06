% 需要检查状态空间表示中是否含有延迟
clear; close all; clc;

%% 1. 加载数据与模型构建
load('..\dataset\bltdWhiteNoise_ssmodel.mat');
load('..\dataset\bpf_ssmodel.mat');

% --- 耦合系统 (Plant) ---
n = size(Af, 1);    % 原系统状态维度
p = size(Aw, 1);    % 干扰模型状态维度
A = blkdiag(Af, Aw); 

% 增广输出矩阵 C (测量通道)
C = [Cf , Cw];
ny = size(C, 1);

% 增广输入矩阵构建
% B_aug: 控制输入通道 (u)
B_aug = [Bf; zeros(p, size(Bf, 2))];  
% G_aug: 干扰/噪声输入通道 (w)
G_aug = [zeros(n, size(Bw, 2)); Bw];  

nu = size(B_aug, 2); % 控制输入维度
nw = size(G_aug, 2); % 噪声输入维度
nx = size(A, 1);     % 总状态维度

% D 矩阵
D_aug = zeros(ny, nu + nw);

%% 2. 核心修正：构建带明确分组的系统对象
% 只有在 ss() 构造函数中直接指定 InputGroup，才能最稳妥地被 lqg 识别
inputs = [B_aug, G_aug];

sys_plant = ss(A, inputs, C, D_aug, -1, ...
    'InputGroup', struct('controls', 1:nu, 'noise', (nu+1):(nu+nw)), ...
    'OutputGroup', struct('measured', 1:ny));

fprintf('系统构建完成：\n 状态 nx=%d\n 控制输入 nu=%d (Index: %s)\n 噪声输入 nw=%d (Index: %s)\n', ...
    nx, nu, mat2str(1:nu), nw, mat2str((nu+1):(nu+nw)));

%% 3. 权重矩阵设计
% --- LQR 权重 (QXU) ---
% 维度应当是 (nx + nu) x (nx + nu) = 17 x 17
Q_lqr = C' * C;           % 状态权重 Q
R_lqr = 10e-4 * eye(nu);  % 输入权重 R
QXU = blkdiag(Q_lqr, R_lqr);

% --- Kalman 滤波器噪声协方差 (QWV) ---
% 维度应当是 (nw + ny) x (nw + ny) = 2 x 2
Qn = 2000 * eye(nw);      % 过程噪声协方差 W
Rn = 1e-4 * eye(ny);      % 测量噪声协方差 V
QWV = blkdiag(Qn, Rn);

%% 4. 控制器求解 (提供两种方式)
METHOD = 'LQG_FUNC'; % 可选 'LQG_FUNC' 或 'MANUAL'

if strcmp(METHOD, 'LQG_FUNC')
    disp('--- 方式 A: 使用 lqg 函数求解 ---');
    try
        % lqg 会根据 sys_plant.InputGroup 自动拆分 QXU 和 QWV
        Reg = lqg(sys_plant, QXU, QWV);
        disp('lqg 函数求解成功！');
    catch ME
        warning('lqg 函数调用失败，原因: %s', ME.message);
        disp('尝试切换到方式 B (手动设计)...');
        METHOD = 'MANUAL';
    end
end

if strcmp(METHOD, 'MANUAL')
    disp('--- 方式 B: 使用 lqr + kalman + lqgreg 分步求解 ---');
    % 1. 设计 LQR 增益 K
    % 提取 B (控制) 和 G (噪声)
    B_ctrl = sys_plant.B(:, 1:nu);
    G_noise = sys_plant.B(:, (nu+1):end);
    
    K = dlqr(sys_plant.A, B_ctrl, Q_lqr, R_lqr);
    
    % 2. 设计 Kalman 滤波器
    % 构造噪声模型: x[k+1] = Ax + Bu + Gw, y = Cx + v
    % kalman 函数需要区分 B (已知输入) 和 G (噪声输入)
    % 这里我们创建一个仅用于滤波设计的 plant
    sys_kf = ss(sys_plant.A, [B_ctrl, G_noise], sys_plant.C, [zeros(ny,nu), zeros(ny,nw)], -1);
    
    % 指定哪些是已知输入(Controls)，哪些是过程噪声(Noise)
    % kalman(sys, Qn, Rn) 默认最后面的输入是噪声，但最好显式指定
    sys_kf.InputGroup.Known = 1:nu;
    sys_kf.InputGroup.Noise = (nu+1):(nu+nw);
    
    [kest, L, P] = kalman(sys_kf, Qn, Rn);
    
    % 3. 组合成 LQG 调节器
    Reg = lqgreg(kest, K);
    disp('手动设计完成。');
end

% 提取控制器参数用于仿真
[A_reg, B_reg, C_reg, D_reg] = ssdata(Reg);

%% 5. 闭环控制仿真
N = 20000;
x = zeros(nx, 1);                 % 物理系统状态
x_reg = zeros(size(A_reg, 1), 1); % 控制器状态

u_history = zeros(nu, N);
y_history = zeros(ny, N);
d_history_rec = zeros(1, N);      % 记录干扰(用于对比)

fprintf('开始 %d 步仿真...\n', N);
for k = 1:N
    % 1. 环境噪声生成
    w = sqrt(Qn(1,1)) * randn(nw, 1);
    v = sqrt(Rn(1,1)) * randn(ny, 1);

    % 2. 物理系统输出 (y = Cx + v)
    y = C * x + v;

    % 3. 控制器计算 (u = C_reg*x_reg + D_reg*y)
    u = C_reg * x_reg + D_reg * y;

    % 4. 状态更新
    % 控制器: x_reg[k+1]
    x_reg_next = A_reg * x_reg + B_reg * y;
    
    % 物理系统: x[k+1] = Ax + Bu + Gw
    x_next = A * x + B_aug * u + G_aug * w;
    
    % 更新
    x_reg = x_reg_next;
    x = x_next;

    % 5. 数据记录
    u_history(:, k) = u;
    y_history(:, k) = y;
    
    % 提取纯干扰分量 (d = Cw * xw)
    xw = x(n+1:end); 
    d_history_rec(k) = Cw * xw;
end

%% 6. 结果分析与绘图
p0 = 20e-6; 
SPL_d = 20 * log10(rms(d_history_rec) / p0);
SPL_y = 20 * log10(rms(y_history(:)) / p0);
attenuation = SPL_d - SPL_y;

fprintf('\n--- 最终结果 ---\n');
fprintf('干扰信号声压级: %.2f dB\n', SPL_d);
fprintf('残差信号声压级: %.2f dB\n', SPL_y);
fprintf('降噪量: %.2f dB\n', attenuation);

figure('Name', 'LQG Control Performance', 'Color', 'w');
subplot(2, 1, 1);
plot(u_history(1,:), 'b', 'LineWidth', 0.5); 
ylabel('幅度'); title('控制输出 u[k]'); grid on;
xlim([0 N]);

subplot(2, 1, 2);
plot(d_history_rec, 'k', 'DisplayName', '开环干扰 d[k]'); hold on;
plot(y_history(1,:), 'r', 'DisplayName', '闭环残差 y[k]');
ylabel('幅度'); xlabel('样本 k');
title(sprintf('降噪效果 (Att: %.2f dB)', attenuation));
legend('Location', 'best'); grid on;
xlim([0 N]);