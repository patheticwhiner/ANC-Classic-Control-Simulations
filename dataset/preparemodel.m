%%
clear; close all; clc;

%% 导入模型
addpath('functions');
load('dataset\ARMAX_SYSID_30303022.mat');
% 定义对象的标称模型
d = ARMAXmodel.orders(4);
A_poly = ARMAXmodel.model.A;
B_poly = ARMAXmodel.model.B;
B_poly = B_poly(d+1:end);
nA = ARMAXmodel.orders(1);
nB = ARMAXmodel.orders(2);
fs = ARMAXmodel.fs;
ts = 1/fs;

zplane(B_poly,A_poly);

%% 转换为能观标准型 (Observable Canonical Form)
a_coeffs = A_poly(2:end)'; % 取出 a_1, ..., a_nA 并转为列向量
b_coeffs = [B_poly, zeros(1, nA - nB)]'; % 转为列向量并补零

% 构建能观标准型的状态空间矩阵
Ao = zeros(nA);
if nA > 0
    Ao(:, 1) = -a_coeffs; % 第一列
    if nA > 1
        Ao(1:nA-1, 2:nA) = eye(nA-1); % 右上角的单位矩阵
    end
end

Bo = b_coeffs;
Co = [1, zeros(1, nA-1)];
Do = 0; % 假设为严格真系统

% 创建状态空间对象以便于分析
sys_obs = ss(Ao, Bo, Co, Do, ts);

fprintf('能观标准型状态空间矩阵构建完成。\n');

%% 验证 A_poly 和 B_poly 的互质性（考虑数值误差）
[is_coprime, gcd_poly, common_roots, analysis] = check_coprimality(A_poly, B_poly);

%% 检查系统的强可镇定性
[isStrong, unstable_zeros, multiplicities, Bminus, Bplus] = checkStrongStabilizability(B_poly, A_poly);

%% 绘制灵敏度函数
% 1. 准备系统 (Plant)
% 将 ARMAX 多项式转换为 LTI 传递函数对象
% 注意: armax 模型系数是 z^-1 的多项式，tf需要z的多项式，因此需要翻转系数
plant_tf = tf(fliplr(B_poly), fliplr(A_poly), ts);
% 2. 准备控制器 (Controller)
controller_tf = tf(1, 1, ts);
% 3. 准备频率向量
w_rad_s = logspace(-1, log10(pi/ts), fs/2); 
% 4. 调用函数绘制灵敏度图
fprintf('正在调用 plotSensitivity 函数...\n');
plotSensitivity(plant_tf, controller_tf, w_rad_s);
sgtitle('灵敏度分析 (单位反馈 K=1)', 'FontWeight', 'bold');

%% 分析系统零极点
% 分析零点和极点
[poles, zers, pole_data, zero_data] = analyze_zeros_poles(a_coeffs, b_coeffs, fs);
% 打印结果表格
print_analysis_table(poles, zers, pole_data, zero_data, fs);

%% 观测器设计
fprintf('\n=== 观测器设计 ===\n');

% 能观性检查
obsv_matrix = obsv(Ao, Co);
rank_obsv = rank(obsv_matrix);
cond_obsv = cond(obsv_matrix); % 计算条件数
fprintf('能观性矩阵的秩: %d (系统阶数: %d)\n', rank_obsv, nA);
fprintf('能观性矩阵的条件数: %.2e\n', cond_obsv);
if cond_obsv > 1e10
    fprintf('警告: 系统几乎不能观 (条件数非常大)，极点配置法很可能失败。\n');
end

if rank_obsv == nA && is_coprime
    fprintf('✓ 系统理论上完全能观，开始设计观测器。\n');
    
    % --- 方法1: LQR/Kalman滤波器设计 (更稳健，推荐) ---
    fprintf('\n--- 使用LQR方法设计 ---\n');
    Q_obs = eye(nA) * 10;  % 状态噪声权重 (可调)
    R_obs = 1;           % 测量噪声权重 (可调)
    try
        L_lqr = lqr(Ao', Co', Q_obs, R_obs)';
        lqr_obs_poles = eig(Ao - L_lqr*Co);
        fprintf('LQR观测器设计成功。\n');
        % ... (显示极点)
    catch ME
        fprintf('✗ LQR观测器设计失败: %s\n', ME.message);
        L_lqr = [];
    end
    
else
    fprintf('✗ 系统不完全能观，无法设计全阶观测器。\n');
end

%% 状态反馈控制器设计
fprintf('\n=== 状态反馈控制器设计 ===\n');

% 检查能控性
ctrb_matrix = ctrb(Ao, Bo);
rank_ctrb = rank(ctrb_matrix);
fprintf('能控性矩阵的秩: %d (系统阶数: %d)\n', rank_ctrb, nA);

if rank_ctrb == nA
    fprintf('✓ 系统完全能控，可以设计状态反馈控制器\n');
    
    % LQR方法设计控制器
    if exist('lqr', 'file') == 2
        Q_ctrl = eye(nA)*1e-6;  % 状态权重矩阵
        R_ctrl = 1000;        % 控制输入权重
        [K_lqr, ~, lqr_cl_poles] = lqr(Ao, Bo, Q_ctrl, R_ctrl);
        K = K_lqr; % 增加此行
        actual_cl_poles = lqr_cl_poles; % 增加此行
        
        fprintf('LQR闭环极点: ');
        for i = 1:length(lqr_cl_poles)
            if isreal(lqr_cl_poles(i))
                fprintf('%.4f ', lqr_cl_poles(i));
            else
                fprintf('%.4f±%.4fj ', real(lqr_cl_poles(i)), abs(imag(lqr_cl_poles(i)))); 
            end
        end
        fprintf('\n');
    end
    
    % 计算前馈增益（用于跟踪参考信号）
    Nbar = 1 / (Co * inv(eye(nA) - Ao + Bo*K) * Bo); % 用K_lqr或K
    fprintf('前馈增益 Nbar = %.4f\n', Nbar);
    
else
    fprintf('✗ 系统不完全能控，无法设计状态反馈控制器\n');
    K = [];
    K_lqr = [];
    Nbar = [];
end

%% 输出反馈控制器设计
fprintf('\n=== 输出反馈控制器设计 ===\n');

if rank_obsv == nA && rank_ctrb == nA
    fprintf('✓ 系统能观能控，可以设计输出反馈控制器\n');
    
    % 构建增广系统矩阵用于分离原理
    A_aug = [Ao - Bo*K, Bo*K; 
             zeros(nA, nA), Ao - L*Co];
    
    % 验证分离原理：增广系统的极点应该是控制器极点和观测器极点的并集
    aug_poles = eig(A_aug);
    expected_poles = [actual_cl_poles; actual_obs_poles];
    
    fprintf('增广系统极点: ');
    for i = 1:length(aug_poles)
        if isreal(aug_poles(i))
            fprintf('%.4f ', aug_poles(i));
        else
            fprintf('%.4f±%.4fj ', real(aug_poles(i)), abs(imag(aug_poles(i))));
        end
    end
    fprintf('\n');
    
    fprintf('✓ 输出反馈控制器设计完成\n');
    fprintf('  控制律: u = -K*x_hat + Nbar*r\n');
    fprintf('  观测器: x_hat_dot = (Ao - L*Co)*x_hat + Bo*u + L*y\n');
    
else
    fprintf('✗ 无法设计输出反馈控制器\n');
end

%% 保存设计结果
design_results = struct();
design_results.Ao = Ao;
design_results.Bo = Bo;
design_results.Co = Co;
design_results.Do = Do;
design_results.sys_obs = sys_obs;

if exist('L', 'var') && ~isempty(L)
    design_results.L = L;
    design_results.observer_poles = actual_obs_poles;
end

if exist('L_lqr', 'var') && ~isempty(L_lqr)
    design_results.L_lqr = L_lqr;
end

if exist('K', 'var') && ~isempty(K)
    design_results.K = K;
    design_results.closed_loop_poles = actual_cl_poles;
    design_results.Nbar = Nbar;
end

if exist('K_lqr', 'var') && ~isempty(K_lqr)
    design_results.K_lqr = K_lqr;
end

design_results.is_observable = (rank_obsv == nA);
design_results.is_controllable = (rank_ctrb == nA);
design_results.is_coprime = is_coprime;

fprintf('\n=== 设计结果已保存到 design_results 结构体 ===\n');

%% 开环系统分析
fprintf('\n=== 开环系统分析 ===\n');
figure('Name', 'Open-Loop System Analysis');

% 1. 零极点图
subplot(2,2,1);
zplane(Bo, Ao); % 注意：对于状态空间，zplane需要不同的输入，但我们用多项式更直观
pzmap(sys_obs);
title('Pole-Zero Map');
grid on;

% 2. 阶跃响应
subplot(2,2,2);
step(sys_obs);
title('Step Response');
grid on;

% 3. 脉冲响应
subplot(2,2,3);
impulse(sys_obs);
title('Impulse Response');
grid on;

% 4. 波特图
subplot(2,2,4);
bode(sys_obs);
title('Bode Plot');
grid on;

sgtitle('Open-Loop System Analysis', 'FontSize', 14, 'FontWeight', 'bold');

%% 互质性检验函数
function [is_coprime, gcd_poly, common_roots, analysis] = check_coprimality(A, B, varargin)
    % CHECK_COPRIMALITY 检验两个多项式是否互质（考虑数值误差）
    %
    % 输入:
    %   A, B - 多项式系数向量（降幂排列）
    %   varargin - 可选参数:
    %     'tolerance' - 数值容差（默认: 1e-8）
    %     'method' - 检验方法 ('gcd', 'roots', 'both') （默认: 'both'）
    %
    % 输出:
    %   is_coprime - 逻辑值，是否互质
    %   gcd_poly - 最大公约式多项式
    %   common_roots - 共同根
    %   analysis - 分析结果结构体
    
    % 解析输入参数
    p = inputParser;
    addParameter(p, 'tolerance', 1e-8, @isnumeric);
    addParameter(p, 'method', 'both', @ischar);
    parse(p, varargin{:});
    
    tol = p.Results.tolerance;
    method = p.Results.method;
    
    % 初始化输出
    common_roots = [];
    analysis = struct();
    analysis.tolerance = tol;
    
    % 方法1: 基于多项式GCD
    if ismember(method, {'gcd', 'both'})
        try
            % 使用符号工具箱（如果可用）
            if exist('gcd', 'file') == 2
                % 转换为符号多项式
                syms z;
                A_sym = poly2sym(A, z);
                B_sym = poly2sym(B, z);
                gcd_sym = gcd(A_sym, B_sym);
                gcd_poly = sym2poly(gcd_sym);
            else
                % 使用数值方法近似GCD
                gcd_poly = numerical_gcd(A, B, tol);
            end
        catch
            % 备用数值方法
            gcd_poly = numerical_gcd(A, B, tol);
        end
        
        % 判断GCD是否为常数
        gcd_order = length(gcd_poly) - 1;
        is_coprime_gcd = (gcd_order == 0) || (abs(gcd_poly(1)) < tol && gcd_order == 1);
        
        analysis.gcd_order = gcd_order;
        analysis.gcd_poly = gcd_poly;
    else
        is_coprime_gcd = true;
        gcd_poly = 1;
        analysis.gcd_order = 0;
    end
    
    % 方法2: 基于根的比较
    if ismember(method, {'roots', 'both'})
        % 计算多项式的根
        roots_A = roots(A);
        roots_B = roots(B);
        
        % 去除数值噪声很大的根
        roots_A = roots_A(abs(roots_A) < 1e6);
        roots_B = roots_B(abs(roots_B) < 1e6);
        
        % 寻找共同根
        common_roots = find_common_roots(roots_A, roots_B, tol);
        
        % 计算最小根距离
        min_dist = inf;
        for i = 1:length(roots_A)
            for j = 1:length(roots_B)
                dist = abs(roots_A(i) - roots_B(j));
                if dist < min_dist
                    min_dist = dist;
                end
            end
        end
        
        is_coprime_roots = isempty(common_roots);
        analysis.min_root_distance = min_dist;
        analysis.common_roots_count = length(common_roots);
    else
        is_coprime_roots = true;
        analysis.min_root_distance = inf;
    end
    
    % 综合判断
    if strcmp(method, 'both')
        is_coprime = is_coprime_gcd && is_coprime_roots;
        analysis.method_agreement = (is_coprime_gcd == is_coprime_roots);
    elseif strcmp(method, 'gcd')
        is_coprime = is_coprime_gcd;
    else
        is_coprime = is_coprime_roots;
    end
    
    analysis.is_coprime_gcd = is_coprime_gcd;
    analysis.is_coprime_roots = is_coprime_roots;
    fprintf('\n=== 互质性分析结果 ===\n');
    if is_coprime
        fprintf('✓ A(z) 和 B(z) 近似互质，系统能观\n');
    else
        fprintf('✗ A(z) 和 B(z) 不互质，系统不能观\n');
        fprintf('共同根的数量: %d\n', length(common_roots));
        fprintf('共同根: ');
        for i = 1:length(common_roots)
            if isreal(common_roots(i))
                fprintf('%.4f ', common_roots(i));
            else
                fprintf('%.4f±%.4fj ', real(common_roots(i)), abs(imag(common_roots(i))));
            end
        end
        fprintf('\n');
    end
    fprintf('结果分析:\n');
    fprintf('- GCD 多项式阶数: %d\n', analysis.gcd_order);
    fprintf('- 最小根距离: %.2e\n', analysis.min_root_distance);
    fprintf('- 数值容差: %.2e\n', analysis.tolerance);
end

function common_roots = find_common_roots(roots_A, roots_B, tol)
    % 寻找两组根中的共同根
    common_roots = [];
    used_B = false(size(roots_B));
    
    for i = 1:length(roots_A)
        for j = 1:length(roots_B)
            if ~used_B(j) && abs(roots_A(i) - roots_B(j)) < tol
                common_roots(end+1) = roots_A(i);
                used_B(j) = true;
                break;
            end
        end
    end
end

function gcd_poly = numerical_gcd(A, B, tol)
    % 数值方法计算多项式GCD（简化版本）
    % 这是一个基础实现，实际应用中可能需要更复杂的算法
    
    % 移除前导零
    A = A(find(abs(A) > tol, 1):end);
    B = B(find(abs(B) > tol, 1):end);
    
    if isempty(A), A = 0; end
    if isempty(B), B = 0; end
    
    % 简单情况处理
    if length(A) == 1 || length(B) == 1
        gcd_poly = 1;
        return;
    end
    
    % 使用Sylvester矩阵方法的简化版本
    % 计算两个多项式根的最小距离
    try
        roots_A = roots(A);
        roots_B = roots(B);
        
        min_dist = inf;
        for i = 1:length(roots_A)
            for j = 1:length(roots_B)
                dist = abs(roots_A(i) - roots_B(j));
                if dist < min_dist
                    min_dist = dist;
                end
            end
        end
        
        if min_dist < tol
            % 存在共同根，不互质
            gcd_poly = [1, 0]; % 表示一阶多项式
        else
            gcd_poly = 1; % 互质
        end
    catch
        % 如果出错，保守地假设不互质
        gcd_poly = [1, 0];
    end
end