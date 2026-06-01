% RUN_RST_EMOPSO  基于 ε-MOPSO 的 RST 控制器参数整定主仿真脚本
%
% 工作流程:
%   1. 加载/定义被控对象模型 (A, B, d)
%   2. 提取不稳定零点 β_i 和不稳定极点 α_k
%   3. 设定 P_D, H_R, H_S, n_X, n_Y
%   4. 配置 ε-MOPSO 算法参数
%   5. 构建 sys_info 结构体和包装目标函数
%   6. 调用 eMOPSO_core 进行多目标优化
%   7. 可视化 Pareto 前沿
%   8. 从 Pareto 前沿选择拐点解
%   9. 保存优化结果到 .mat 文件
%
% 参考文献:
%   [1] Sierra, M.R. and Coello Coello, C.A., "Improving PSO-Based
%       Multi-objective Optimization Using Crowding, Mutation and
%       ε-Dominance", EMO 2005, LNCS 3410, pp. 505-519, 2005.
%   [2] Zames, G. and Francis, B.A., "Feedback, Minimax Sensitivity, and
%       Optimal Robustness", IEEE Trans. AC, 28(5), pp. 585-601, 1983.

clear; close all; clc;

% 确保路径可达
addpath('../functions');
addpath('.');

%% ====================================================================
% 第 1 节：被控对象定义
% ====================================================================

% ---- 方式 a) 手动定义（用于测试） ----
% 示例：二阶系统 G(z) = z^{-d} * B(z^{-1}) / A(z^{-1})
%   B(z^{-1}) = b0 + b1*z^{-1} + ...
%   A(z^{-1}) = 1  + a1*z^{-1} + a2*z^{-2} + ...
B = [0.2, 0.15];      % 分子系数 (降幂: b0, b1, ...)
A = [1, -1.2, 0.45];  % 分母系数 (降幂: 1, a1, a2, ...)
d = 1;                 % 纯延迟步数
Ts = 1;                % 采样时间

% ---- 方式 b) 从文件加载（需要时取消注释） ----
% load('dataset/ARMAX_SYSID_30303022.mat');
% % 假设加载的结构体包含字段 A, B, d, Ts
% % 若结构体名称不同，请相应调整
% A  = model.A;
% B  = model.B;
% d  = model.d;
% Ts = model.Ts;

fprintf('========== 被控对象模型 ==========\n');
fprintf('A(z^{-1}) = ');
fprintf('%+.4f z^{-%d} ', [A; (0:length(A)-1)]);
fprintf('\nB(z^{-1}) = ');
fprintf('%+.4f z^{-%d} ', [B; (0:length(B)-1)]);
fprintf('\n延迟步数 d = %d, 采样时间 Ts = %.4f\n', d, Ts);

%% ====================================================================
% 第 2 节：系统分析 — 提取不稳定零极点
% ====================================================================

% 计算系统零极点
sys_zeros = roots(B);
sys_poles = roots(A);

% 提取不稳定零点 (|z| > 1)
beta = sys_zeros(abs(sys_zeros) > 1);
% 提取不稳定极点 (|z| > 1)
alpha = sys_poles(abs(sys_poles) > 1);

% 如果无情不稳定零点，添加一个虚拟零点 β_virtual = ∞
% (在实际计算中用一个足够大的值表示，或在目标函数中跳过)
if isempty(beta)
    fprintf('\n注意: 系统无不稳定零点 (最小相位系统)。\n');
    fprintf('将使用默认虚拟零点 β_virtual = 1e6 以保证目标函数维度一致。\n');
    beta = 1e6;
end

% 如果无不确定极点，添加虚拟极点
if isempty(alpha)
    fprintf('注意: 系统无不稳定极点 (开环稳定系统)。\n');
    alpha = 1.0001;  % 虚拟极点，略大于1
end

fprintf('\n========== 系统零极点分析 ==========\n');
fprintf('不稳定零点 (|z|>1): %d 个\n', length(beta));
for i = 1:length(beta)
    fprintf('  β_%d = %.4f + %.4fi, |β_%d| = %.4f\n', ...
        i, real(beta(i)), imag(beta(i)), i, abs(beta(i)));
end
fprintf('不稳定极点 (|z|>1): %d 个\n', length(alpha));
for i = 1:length(alpha)
    fprintf('  α_%d = %.4f + %.4fi, |α_%d| = %.4f\n', ...
        i, real(alpha(i)), imag(alpha(i)), i, abs(alpha(i)));
end

%% ====================================================================
% 第 3 节：RST 控制器预设参数
% ====================================================================

% 主导极点多项式 P_D (示例: 二阶 Butterworth 型, 截止频率 ω_c ≈ 0.3π)
% P_D(z) = 1 - 2*exp(-ζ*ω_n*Ts)*cos(ω_n*Ts*sqrt(1-ζ^2))*z^{-1} + exp(-2*ζ*ω_n*Ts)*z^{-2}
zeta_P = 0.8;
omega_P = 0.25 * pi;  % 主导极点自然频率
PD1 = -2 * exp(-zeta_P * omega_P) * cos(omega_P * sqrt(1 - zeta_P^2));
PD2 = exp(-2 * zeta_P * omega_P);
P_D = [1, PD1, PD2];

fprintf('\n========== 主导极点配置 ==========\n');
fprintf('ζ = %.2f, ω_n = %.2f rad/sample\n', zeta_P, omega_P);
fprintf('P_D(z^{-1}) = 1 + %.4f z^{-1} + %.4f z^{-2}\n', PD1, PD2);
fprintf('P_D 极点: %.4f ± %.4fi\n', ...
    real(exp(-zeta_P*omega_P)*exp(1j*omega_P*sqrt(1-zeta_P^2))), ...
    abs(imag(exp(-zeta_P*omega_P)*exp(1j*omega_P*sqrt(1-zeta_P^2)))));

% 控制器固定因子
H_R = [1, 1];    % 分子: 1 + z^{-1} (用于嵌入积分作用或高频衰减)
H_S = [1, -1];   % 分母: 1 - z^{-1} (用于嵌入积分作用消除稳态误差)

% X 和 Y 多项式阶数 (根据 Bezout 方程解存在条件)
% n_S = n_B + n_HR - 1, n_R = n_A + n_HS - 1
% 此处 nX, nY 是 X 和 Y 多项式的阶数（即根的数量）
nX = 2;  % Bezout可解: deg(P_D*X) ≤ deg(A*HS)+deg(B*HR)-1 = 4
nY = 3;  % 适中的复杂度

fprintf('\n========== 控制器预设参数 ==========\n');
fprintf('H_R(z^{-1}) = ');
fprintf('%+.4f z^{-%d} ', [H_R; (0:length(H_R)-1)]);
fprintf('\nH_S(z^{-1}) = ');
fprintf('%+.4f z^{-%d} ', [H_S; (0:length(H_S)-1)]);
fprintf('\nX 多项式阶数 nX = %d\n', nX);
fprintf('Y 多项式阶数 nY = %d\n', nY);

%% ====================================================================
% 第 4 节：eMOPSO 算法参数配置
% ====================================================================

% 决策变量维度
n_var = nX + nY;

% 目标函数数量: m_f = (不稳定零点数) + 1(延迟) + 1(Bezout残差)
n_obj = length(beta) + 2;

% 搜索空间边界: 收紧到 ±0.95 防止根过于靠近单位圆
lb = -0.95 * ones(1, n_var);   % 收紧搜索空间
ub =  0.95 * ones(1, n_var);   % 防止根过于靠近单位圆

% ε-MOPSO 参数 (增大种群和迭代以改善收敛性)
options = struct();
options.n_pop       = 40;     % 增大种群保持多样性
options.k_max       = 80;     % 增加迭代改善收敛
options.epsilon     = 0.02;   % 更细的ε值
options.c1          = 1.5;    % 认知加速系数
options.c2          = 1.7;    % 社会加速系数
options.w_max       = 0.9;    % 最大惯性权重
options.w_min       = 0.3;    % 降低最小惯性权重加速收敛
options.archive_max = 80;     % 增大Archive容量
options.verbose     = true;   % 开启详细输出

fprintf('\n========== ε-MOPSO 参数配置 ==========\n');
fprintf('决策变量维度: %d\n', n_var);
fprintf('目标函数数量: %d\n', n_obj);
fprintf('种群大小: %d\n', options.n_pop);
fprintf('最大迭代: %d\n', options.k_max);
fprintf('ε 值: %.4f\n', options.epsilon);
fprintf('搜索空间: [%.2f, %.2f]\n', lb(1), ub(1));

%% ====================================================================
% 第 5 节：构建系统信息与目标函数包装器
% ====================================================================

% 构建 sys_info 结构体
sys_info = struct();
sys_info.A      = A;
sys_info.B      = B;
sys_info.d      = d;
sys_info.P_D    = P_D;
sys_info.H_R    = H_R;
sys_info.H_S    = H_S;
sys_info.nX     = nX;
sys_info.nY     = nY;
sys_info.beta   = beta(:);
sys_info.alpha  = alpha(:);
sys_info.Ts     = Ts;
sys_info.nFreq  = 500;   % 适中的频率分辨率

% 包装目标函数 (含约束惩罚)
% 重要: eMOPSO_core 调用 obj_func(x) 返回列向量
obj_func = @(theta) wrapped_objective(theta, sys_info);

fprintf('\n========== 开始优化 ==========\n');

%% ====================================================================
% 第 6 节：运行优化
% ====================================================================

tic;
[archive, archive_f, history] = eMOPSO_core(obj_func, n_var, n_obj, lb, ub, options);
elapsed_time = toc;

fprintf('\n========== 优化完成 ==========\n');
fprintf('耗时: %.2f 秒\n', elapsed_time);
fprintf('Pareto 前沿解数量: %d\n', size(archive, 1));

%% ====================================================================
% 第 7 节：Pareto 前沿可视化
% ====================================================================

figure('Name', 'Pareto Front - RST eMOPSO', 'Position', [100, 100, 800, 600]);

if n_obj == 2
    % 2目标情况 — 二维散点图
    scatter(archive_f(:,1), archive_f(:,2), 40, 'b', 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('f_1 (Zero matching error)');
    ylabel('f_2 (Delay integral error)');
    title(sprintf('Pareto Front (ε=%.3f, %d solutions)', options.epsilon, size(archive_f,1)));
    grid on;

elseif n_obj >= 3
    % 3目标情况 — 三维散点图
    scatter3(archive_f(:,1), archive_f(:,2), archive_f(:,3), 40, archive_f(:,3), 'filled', ...
        'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
    xlabel('f_1 (Zero error)');
    ylabel('f_2 (Zero error)');
    zlabel('f_3 (Delay integral error)');
    title(sprintf('Pareto Front (ε=%.3f, %d solutions)', options.epsilon, size(archive_f,1)));
    grid on;
    view(45, 30);
    colorbar;
end

% 收敛曲线
figure('Name', 'Convergence', 'Position', [150, 150, 500, 400]);
plot(history, 'b-', 'LineWidth', 1.5);
xlabel('Iteration');
ylabel('Archive Size');
title('Archive Size vs Iteration');
grid on;

%% ====================================================================
% 第 8 节：选择 Pareto 最优解
% ====================================================================

% 选择拐点解 (knee point) — 用归一化加权和的最小值
f_range = max(archive_f) - min(archive_f);
f_range(f_range == 0) = 1;  % 防止除零
f_norm = (archive_f - min(archive_f)) ./ f_range;
[~, knee_idx] = min(sum(f_norm, 2));
theta_opt = archive(knee_idx, :)';
f_opt = archive_f(knee_idx, :)';

fprintf('\n========== 选定最优解 ==========\n');
fprintf('解索引: %d (共 %d 个 Pareto 解)\n', knee_idx, size(archive, 1));
fprintf('目标函数值: ');
fprintf('%.6f ', f_opt);
fprintf('\n');
fprintf('决策变量 θ (前%d个 X根, 后%d个 Y根):\n', nX, nY);
fprintf('  X roots = [');
fprintf('%.4f ', theta_opt(1:nX));
fprintf(']\n');
fprintf('  Y roots = [');
fprintf('%.4f ', theta_opt(nX+1:end));
fprintf(']\n');

%% ====================================================================
% 第 9 节：保存结果
% ====================================================================

% 确保 output 目录存在
if ~exist('output', 'dir'), mkdir('output'); end
save('output/RST_eMOPSO_results.mat', ...
    'archive', 'archive_f', 'history', ...
    'theta_opt', 'f_opt', ...
    'sys_info', 'options', 'elapsed_time');
fprintf('\n结果已保存到: output/RST_eMOPSO_results.mat\n');

% =========================================================================
% 本地包装函数：目标函数 + 约束惩罚
% =========================================================================

function phi = wrapped_objective(theta, sys_info)
% WRAPPED_OBJECTIVE  包装目标函数并施加约束惩罚
%   调用 RST_objective 计算原始目标值，再调用 RST_constraints 计算约束违例，
%   最后使用 penalty_function 施加惩罚。
%
%   输入:
%     theta    - 决策变量向量 (列向量)
%     sys_info - 系统信息结构体
%
%   输出:
%     phi      - 惩罚后的增广目标函数值 (列向量)

    % 确保 theta 为列向量
    theta = theta(:);

    % 计算原始目标函数，放大1000倍以匹配 epsilon=0.05 的量级
    f = RST_objective(theta, sys_info) * 1000;

    % 计算约束
    [g1, g2] = RST_constraints(theta, sys_info);

    % 施加惩罚 (Lambda=0.5，等比放大以匹配缩放后的目标函数)
    phi = penalty_function(f, g1, g2, 0.5, 1);
end
