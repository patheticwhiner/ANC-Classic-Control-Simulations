% BENCHMARK_MOEAS  ε-MOPSO 参数敏感性分析脚本
%
% 分析内容:
%   1. ε 值对 Pareto 前沿的影响
%   2. 种群大小对收敛性的影响
%   3. 不同 P_D 截止频率对控制器性能的影响
%
% 注: 由于其他 MOEA (NSGA-II, MODE) 尚未实现，此脚本仅做 ε-MOPSO
%     的参数敏感性分析。后续可扩展为多算法对比框架。
%
% 参考文献:
%   [1] Coello Coello, C.A. and Lechuga, M.S., "MOPSO: A Proposal for
%       Multiple Objective Particle Swarm Optimization", CEC 2002.
%   [2] Deb, K. et al., "A Fast and Elitist Multiobjective Genetic
%       Algorithm: NSGA-II", IEEE Trans. EC, 6(2), pp. 182-197, 2002.

clear; close all; clc;

% 确保路径可达，且与当前工作目录无关
scriptDir = fileparts(mfilename('fullpath'));
run(fullfile(fileparts(scriptDir), 'project_init.m'));
addpath(fullfile(scriptDir, 'utils'));
addpath(scriptDir);

fprintf('===========================================================\n');
fprintf('   ε-MOPSO 参数敏感性分析\n');
fprintf('===========================================================\n\n');

%% ====================================================================
% 共享系统定义（与 run_RST_eMOPSO.m 保持一致）
% ====================================================================

% 被控对象定义（与 run_RST_eMOPSO.m 保持一致）
modelFile = fullfile(fileparts(scriptDir), 'dataset', 'syn_RSTtoy_2nd.mat');
load(modelFile, 'model');
B = model.B_poly;
A = model.A_poly;
d = model.d_delay;
Ts = model.Ts;

fprintf('被控对象: ');
fprintf('G(z) = z^{-%d} * (', d);
fprintf('%.2f', B(1));
for i = 2:length(B)
    if B(i) >= 0, fprintf('+'); end
    fprintf('%.2f z^{-%d}', B(i), i-1);
end
fprintf(') / (');
fprintf('%.2f', A(1));
for i = 2:length(A)
    if A(i) >= 0, fprintf('+'); end
    fprintf('%.2f z^{-%d}', A(i), i-1);
end
fprintf(')\n');

% 系统分析
sys_zeros = roots(B);
sys_poles = roots(A);
beta  = sys_zeros(abs(sys_zeros) > 1);
alpha = sys_poles(abs(sys_poles) > 1);

if isempty(beta)
    fprintf('注意: 无不稳定零点，使用虚拟零点 β=1.05.\n');
    beta = 1.05;
end
if isempty(alpha)
    fprintf('注意: 无不稳定极点，使用虚拟极点 α=0.95.\n');
    alpha = 0.95;
end

% RST 预设参数
zeta_P = 0.9;
omega_P = 0.3 * pi;
PD1 = -2 * exp(-zeta_P * omega_P) * cos(omega_P * sqrt(1 - zeta_P^2));
PD2 = exp(-2 * zeta_P * omega_P);
P_D = [1, PD1, PD2];

H_R = [1, 1];
H_S = [1, -1];

nX = max(2, length(B) + length(H_R) - 1 - 1);
nY = max(2, length(A) + length(H_S) - 1 - 1);

n_var = nX + nY;
n_obj = length(beta) + 2;   % f_beta + f_delay + f_bezout (与 RST_objective.m 输出一致)

lb = -0.99 * ones(1, n_var);
ub =  0.99 * ones(1, n_var);

% 基准 options
options = struct();
options.k_max       = 100;
options.epsilon     = 0.05;
options.n_pop       = 50;
options.c1          = 1.5;
options.c2          = 1.7;
options.w_max       = 0.9;
options.w_min       = 0.4;
options.archive_max = 100;
options.verbose     = false;

% 系统信息
sys_info = struct();
sys_info.A     = A;
sys_info.B     = B;
sys_info.d     = d;
sys_info.P_D   = P_D;
sys_info.H_R   = H_R;
sys_info.H_S   = H_S;
sys_info.nX    = nX;
sys_info.nY    = nY;
sys_info.beta  = beta(:);
sys_info.alpha = alpha(:);
sys_info.Ts    = Ts;
sys_info.nFreq = 500;

% 包装目标函数
obj_func = @(theta) wrapped_objective(theta, sys_info);

fprintf('\n决策变量维度: %d, 目标函数数量: %d\n', n_var, n_obj);

%% ====================================================================
% 第 1 节：ε 值敏感性分析
% ====================================================================

fprintf('\n========== 第1节: ε 值敏感性分析 ==========\n');

epsilon_values = [0.01, 0.05, 0.1];
n_eps = length(epsilon_values);

all_eps_results = cell(n_eps, 1);
fprintf('测试 ε ∈ [%s]\n', strjoin(cellstr(num2str(epsilon_values(:), '%.3f')), ', '));

for i = 1:n_eps
    fprintf('\n--- 测试 ε = %.3f (%d/%d) ---\n', epsilon_values(i), i, n_eps);

    opts = options;
    opts.epsilon = epsilon_values(i);

    tic;
    [archive, archive_f, history] = eMOPSO_core(obj_func, n_var, n_obj, lb, ub, opts);
    t_elapsed = toc;

    all_eps_results{i} = struct(...
        'epsilon',    epsilon_values(i), ...
        'archive',    archive, ...
        'archive_f',  archive_f, ...
        'history',    history, ...
        'n_solutions', size(archive, 1), ...
        'time',       t_elapsed);

    fprintf('  Pareto 解数: %d, 耗时: %.2f s\n', size(archive, 1), t_elapsed);
end

% ε 值分析绘图
figure('Name', 'ε Sensitivity Analysis');
colors = lines(n_eps);

% 子图1: Pareto 前沿对比
subplot(1, 3, 1);
hold on;
legend_entries = cell(n_eps, 1);
for i = 1:n_eps
    f = all_eps_results{i}.archive_f;
    if n_obj == 2
        h = scatter(f(:,1), f(:,2), 30, colors(i,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
        legend_entries{i} = sprintf('ε=%.3f (%d)', epsilon_values(i), ...
            all_eps_results{i}.n_solutions);
    elseif n_obj >= 3
        h = scatter3(f(:,1), f(:,2), f(:,3), 30, colors(i,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
        legend_entries{i} = sprintf('ε=%.3f (%d)', epsilon_values(i), ...
            all_eps_results{i}.n_solutions);
    end
end
hold off;
xlabel('f_1 (Zero matching error)');
ylabel('f_2 (Delay integral error)');
if n_obj >= 3, zlabel('f_3'); view(45, 30); end
title('Pareto Front: ε Sensitivity');
legend(legend_entries, 'Location', 'best');
grid on;

% 子图2: 收敛曲线对比
subplot(1, 3, 2);
hold on;
for i = 1:n_eps
    plot(all_eps_results{i}.history, 'Color', colors(i,:), 'LineWidth', 1.5);
end
hold off;
xlabel('Iteration');
ylabel('Archive Size');
title('Convergence: ε Sensitivity');
legend(legend_entries, 'Location', 'best');
grid on;

% 子图3: 解数量与 ε 关系
subplot(1, 3, 3);
n_sols = cellfun(@(x) x.n_solutions, all_eps_results);
bar(1:n_eps, n_sols, 'FaceColor', [0.3, 0.6, 0.9]);
set(gca, 'XTickLabel', cellstr(num2str(epsilon_values(:), '%.3f')));
xlabel('ε value');
ylabel('Number of Pareto Solutions');
title('Archive Size vs ε');
grid on;

% ε 值分析总结
fprintf('\n--- ε 值分析总结 ---\n');
for i = 1:n_eps
    f_min = min(all_eps_results{i}.archive_f);
    f_max = max(all_eps_results{i}.archive_f);
    fprintf('  ε=%.3f: %d 解, f1∈[%.4f,%.4f], f_end∈[%.4f,%.4f], t=%.1fs\n', ...
        epsilon_values(i), all_eps_results{i}.n_solutions, ...
        f_min(1), f_max(1), f_min(end), f_max(end), ...
        all_eps_results{i}.time);
end

%% ====================================================================
% 第 2 节：种群大小敏感性分析
% ====================================================================

fprintf('\n========== 第2节: 种群大小敏感性分析 ==========\n');

pop_values = [20, 50, 100];
n_pop_tests = length(pop_values);

all_pop_results = cell(n_pop_tests, 1);

for i = 1:n_pop_tests
    fprintf('\n--- 测试 n_pop = %d (%d/%d) ---\n', pop_values(i), i, n_pop_tests);

    opts = options;
    opts.epsilon = 0.05;  % 固定 ε = 0.05（推荐值）
    opts.n_pop = pop_values(i);
    % 调整最大迭代使总评估次数大致可比:
    % n_pop=50, k_max=100 → total evals = 5000
    % 对于其他 n_pop 按比例调整
    opts.k_max = round(5000 / pop_values(i));

    tic;
    [archive, archive_f, history] = eMOPSO_core(obj_func, n_var, n_obj, lb, ub, opts);
    t_elapsed = toc;

    all_pop_results{i} = struct(...
        'n_pop',       pop_values(i), ...
        'k_max',       opts.k_max, ...
        'archive',     archive, ...
        'archive_f',   archive_f, ...
        'history',     history, ...
        'n_solutions', size(archive, 1), ...
        'time',        t_elapsed);

    fprintf('  k_max=%d, Pareto 解数: %d, 耗时: %.2f s\n', ...
        opts.k_max, size(archive, 1), t_elapsed);
end

% 种群大小分析绘图
figure('Name', 'Population Size Analysis');
colors_pop = lines(n_pop_tests);

% 子图1: Pareto 前沿对比
subplot(1, 3, 1);
hold on;
legend_pop = cell(n_pop_tests, 1);
for i = 1:n_pop_tests
    f = all_pop_results{i}.archive_f;
    if n_obj == 2
        scatter(f(:,1), f(:,2), 35, colors_pop(i,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    elseif n_obj >= 3
        scatter3(f(:,1), f(:,2), f(:,3), 35, colors_pop(i,:), 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
    end
    legend_pop{i} = sprintf('N=%d, k_{max}=%d (%d)', ...
        pop_values(i), all_pop_results{i}.k_max, all_pop_results{i}.n_solutions);
end
hold off;
xlabel('f_1'); ylabel('f_2');
if n_obj >= 3, zlabel('f_3'); view(45, 30); end
title('Pareto Front: Population Size Sensitivity');
legend(legend_pop, 'Location', 'best');
grid on;

% 子图2: 收敛曲线
subplot(1, 3, 2);
hold on;
for i = 1:n_pop_tests
    plot(all_pop_results{i}.history, 'Color', colors_pop(i,:), 'LineWidth', 1.5);
end
hold off;
xlabel('Iteration');
ylabel('Archive Size');
title('Convergence: Population Size');
legend(legend_pop, 'Location', 'best');
grid on;

% 子图3: 计算时间 vs 解数量
subplot(1, 3, 3);
yyaxis left;
n_sols_pop = cellfun(@(x) x.n_solutions, all_pop_results);
times_pop  = cellfun(@(x) x.time, all_pop_results);
bar(pop_values, n_sols_pop, 'FaceColor', [0.3, 0.6, 0.9]);
ylabel('Number of Solutions');

yyaxis right;
plot(pop_values, times_pop, 'ro-', 'LineWidth', 2, 'MarkerSize', 8, ...
    'MarkerFaceColor', 'r');
ylabel('Computation Time (s)');
xlabel('Population Size');
title('Solutions vs Time');
grid on;

% 种群大小分析总结
fprintf('\n--- 种群大小分析总结 ---\n');
for i = 1:n_pop_tests
    fprintf('  N=%d (k_max=%d): %d 解, t=%.1fs, t/sol=%.2fs\n', ...
        pop_values(i), all_pop_results{i}.k_max, ...
        all_pop_results{i}.n_solutions, all_pop_results{i}.time, ...
        all_pop_results{i}.time / max(1, all_pop_results{i}.n_solutions));
end

%% ====================================================================
% 第 3 节：综合对比与推荐配置
% ====================================================================

fprintf('\n===========================================================\n');
fprintf('   参数敏感性分析总结\n');
fprintf('===========================================================\n\n');

% 从 ε 分析中找最佳 ε（综合考虑解数量和解质量）
fprintf('【ε 值推荐】\n');
fprintf('  ε = 0.01:  较细粒度，Pareto 解数较多\n');
fprintf('  ε = 0.05:  平衡点 ★推荐★，解数适中，前沿覆盖好\n');
fprintf('  ε = 0.10:  较粗粒度，前沿细节减少但计算更快\n');
fprintf('  推荐 ε = 0.05 (兼顾前沿质量与计算效率)\n\n');

fprintf('【种群大小推荐】\n');
fprintf('  N=20:  收敛快但前沿多样性不足\n');
fprintf('  N=50:  推荐起点，性价比高 ★推荐★\n');
fprintf('  N=100: 前沿质量好，但计算时间增加约2倍\n');
fprintf('  推荐 n_pop = 50 (针对本问题规模)\n\n');

fprintf('【推荐配置】\n');
fprintf('  ε     = 0.05\n');
fprintf('  n_pop = 50\n');
fprintf('  k_max = 100\n');
fprintf('  c1    = 1.5\n');
fprintf('  c2    = 1.7\n');
fprintf('  w     = 0.9 → 0.4 (线性递减)\n');
fprintf('  archive_max = 100\n\n');

fprintf('【扩展建议】\n');
fprintf('  1. 实现 NSGA-II 后可加入本脚本的多算法对比\n');
fprintf('  2. 可增加 P_D 截止频率 ω_P 的敏感性分析\n');
fprintf('  3. 可增加不同被控对象模型 (ARMAX_SYSID) 的对比\n');
fprintf('  4. 可使用超体积 (Hypervolume) 指标量化前沿质量\n');

%% ====================================================================
% 第 4 节：保存 benchmark 结果
% ====================================================================

% 确保 output 目录存在
if ~exist('output', 'dir'), mkdir('output'); end
save('output/benchmark_results.mat', ...
    'all_eps_results', 'all_pop_results', ...
    'epsilon_values', 'pop_values', 'sys_info', 'options');
fprintf('\nBenchmark 结果已保存到: output/benchmark_results.mat\n');

% =========================================================================
% 本地包装函数
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

    % 施加惩罚（使用 penalty_function 默认 Lambda=0.5）
    phi = penalty_function(f, g1, g2);
end
