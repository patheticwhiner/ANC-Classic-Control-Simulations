function metrics = evaluateFrontier(archive_f, func_num, options)
% EVALUATEFRONTIER 多目标优化 Pareto 前沿量化评估
%
%   对 ε-MOPSO 求得的 Pareto 前沿进行多维度的定量评估，输出结构化指标
%   以便客观比较不同参数配置下的算法性能。
%
%   输入:
%     archive_f  - 目标函数值矩阵 (n_solutions × n_obj)
%     func_num   - 测试函数编号 (1–4)
%     options    - 评估选项结构体 (可选)
%       .n_reference  - 参考前沿采样点数 (默认: 500)
%       .verbose      - 是否打印详细输出 (默认: true)
%
%   输出:
%     metrics     - 评估指标结构体，包含以下字段:
%       .n_solutions   - 前沿解数量
%       .hv            - 超体积 (Hypervolume)，越大越好
%       .gd            - 世代距离 (GD)，越小越好
%       .igd           - 反世代距离 (IGD)，越小越好
%       .spread        - 分布均匀度 (Spread Δ)，越接近0越好
%       .coverage      - 目标空间覆盖率，越大越好
%       .f_range       - 各目标的范围 [min; max]
%       .converged      - 是否达到可行收敛 (逻辑值)
%
%   参考文献:
%     Zitzler, E. and Thiele, L., "Multiobjective Optimization Using
%     Evolutionary Algorithms — A Comparative Case Study", PPSN 1998.
%     Deb, K. et al., "A Fast and Elitist Multiobjective Genetic Algorithm:
%     NSGA-II", IEEE TEC, 2002.

    if nargin < 3
        options = struct();
    end
    if ~isfield(options, 'n_reference'), options.n_reference = 500; end
    if ~isfield(options, 'verbose'),     options.verbose = true;     end

    n_sol = size(archive_f, 1);
    n_obj = size(archive_f, 2);

    % =====================================================================
    % 生成理论 Pareto 参考前沿
    % =====================================================================
    [ref_front, ref_desc] = generateReferenceFront(func_num, options.n_reference);

    % =====================================================================
    % 基本统计
    % =====================================================================
    metrics = struct();
    metrics.n_solutions = n_sol;
    metrics.f_range = [min(archive_f,[],1); max(archive_f,[],1)];

    % =====================================================================
    % Hypervolume (HV) — 以 nadir 点为参考
    % =====================================================================
    nadir = max(archive_f, [], 1) * 1.1 + 0.1;  % 稍大于最大值的参考点
    metrics.hv = computeHypervolume(archive_f, nadir);

    % =====================================================================
    % Generational Distance (GD) — 平均距离到参考前沿
    % =====================================================================
    metrics.gd = computeGD(archive_f, ref_front);

    % =====================================================================
    % Inverted Generational Distance (IGD)
    % =====================================================================
    metrics.igd = computeIGD(archive_f, ref_front);

    % =====================================================================
    % Spread (Δ) — 分布均匀度
    % =====================================================================
    metrics.spread = computeSpread(archive_f);

    % =====================================================================
    % 覆盖率: 用 convex hull 面积/体积近似
    % =====================================================================
    metrics.coverage = computeCoverage(archive_f, ref_front);

    % =====================================================================
    % 收敛判定
    % =====================================================================
    metrics.converged = (n_sol > 0 && metrics.gd < 1.0);

    % =====================================================================
    % 打印评估报告
    % =====================================================================
    if options.verbose
        printReport(metrics, func_num, ref_desc);
    end

end

% =========================================================================
% 理论参考前沿生成
% =========================================================================
function [ref_front, desc] = generateReferenceFront(func_num, n_ref)
% GENERATEREFERENCEFRONT  生成各测试函数的理论 Pareto 前沿离散采样

    switch func_num
        case 1  % F1: f2 = 1 - sqrt(f1), f1 ∈ [0,1]
            f1 = linspace(0, 1, n_ref)';
            f2 = 1 - sqrt(f1);
            ref_front = [f1, f2];
            desc = 'f2 = 1 − √f1 (凸曲线)';

        case 2  % F2: f2 = 2 - √f1 - f1·sin(10πf1), f1 ∈ [0,1]  (Appendix A.2)
            f1 = linspace(0, 1, n_ref)';
            f2 = 2 - sqrt(f1) - f1 .* sin(10 * pi * f1);
            ref_front = [f1, f2];
            desc = 'f2 = 2 − √f1 − f1·sin(10πf1) (振荡凹曲线)';

        case 3  % F3: 单位球面第一象限, f1²+f2²+f3² = 1 (r=1固定)
            n_per_dim = ceil(sqrt(n_ref));
            theta = linspace(0, pi/2, n_per_dim)';
            phi   = linspace(0, pi/2, n_per_dim)';
            [T, P] = meshgrid(theta, phi);
            ref_front = [cos(T(:)).*cos(P(:)), cos(T(:)).*sin(P(:)), sin(T(:))];
            % 均匀采样回 n_ref 点
            if size(ref_front,1) > n_ref
                idx = round(linspace(1, size(ref_front,1), n_ref));
                ref_front = ref_front(idx, :);
            end
            desc = '单位球面第一象限, f1²+f2²+f3² = 1';

        case 4  % F4: f3 = 3.5 - Σ2xi·sin(3πxi), x3→0  (Appendix A.4)
            %   f1=x1, f2=x2, f3=3.5-2x1·sin(3πx1)-2x2·sin(3πx2)-2x3·sin(3πx3)
            %   Pareto最优: x3→0, f3=3.5-2f1·sin(3πf1)-2f2·sin(3πf2)
            n_per_dim = ceil(sqrt(n_ref));
            f1 = linspace(0, 1, n_per_dim)';
            f2 = linspace(0, 1, n_per_dim)';
            [F1, F2] = meshgrid(f1, f2);
            F3_min = 3.5 - 2*F1(:).*sin(3*pi*F1(:)) - 2*F2(:).*sin(3*pi*F2(:));
            ref_front = [F1(:), F2(:), F3_min];
            desc = 'f3 = 3.5 − Σ2xi·sin(3πxi), x3→0 (振荡曲面)';

        otherwise
            error('无效的函数编号 %d', func_num);
    end
end

% =========================================================================
% Hypervolume 计算 (Lebesgue 测度近似, 蒙特卡洛法)
% =========================================================================
function hv = computeHypervolume(archive_f, ref_point)
% COMPUTEHYPERVOLUME  蒙特卡洛法计算超体积
%   使用随机采样点在目标空间中的比例近似 Lebesgue 测度。

    n_obj = size(archive_f, 2);
    f_min = min(archive_f, [], 1);

    % 蒙特卡洛采样 (维度自适应)
    n_samples = min(100000, 10000 * 2^n_obj);
    samples = rand(n_samples, n_obj) .* (ref_point - f_min) + f_min;

    count = 0;
    for i = 1:size(samples, 1)
        % 检查是否有任何解支配该采样点
        dominated = false;
        for j = 1:size(archive_f, 1)
            if all(archive_f(j,:) <= samples(i,:)) && any(archive_f(j,:) < samples(i,:))
                dominated = true;
                break;
            end
        end
        if dominated
            count = count + 1;
        end
    end

    volume = prod(ref_point - f_min);
    hv = volume * count / n_samples;
end

% =========================================================================
% Generational Distance (GD)
% =========================================================================
function gd = computeGD(archive_f, ref_front)
% COMPUTEGD  世代距离: 每个解到最近参考点的平均欧氏距离

    if isempty(archive_f) || isempty(ref_front)
        gd = Inf;
        return;
    end

    n = size(archive_f, 1);
    dists = zeros(n, 1);
    for i = 1:n
        d = sqrt(sum((ref_front - archive_f(i,:)).^2, 2));
        dists(i) = min(d);
    end
    gd = sqrt(sum(dists.^2)) / n;
end

% =========================================================================
% Inverted Generational Distance (IGD)
% =========================================================================
function igd = computeIGD(archive_f, ref_front)
% COMPUTEIGD  反世代距离: 每个参考点到最近解的平均距离

    if isempty(archive_f) || isempty(ref_front)
        igd = Inf;
        return;
    end

    n = size(ref_front, 1);
    dists = zeros(n, 1);
    for i = 1:n
        d = sqrt(sum((archive_f - ref_front(i,:)).^2, 2));
        dists(i) = min(d);
    end
    igd = sum(dists) / n;
end

% =========================================================================
% Spread (Δ) — 分布均匀度
% =========================================================================
function delta = computeSpread(archive_f)
% COMPUTESPREAD  NSGA-II 风格的分布度量
%   Δ = 0 表示完全均匀分布, Δ > 0 表示非均匀

    if size(archive_f, 1) < 2
        delta = Inf;
        return;
    end

    n_sol = size(archive_f, 1);
    n_obj = size(archive_f, 2);

    % 计算所有解之间的成对欧氏距离
    D = zeros(n_sol);
    for i = 1:n_sol
        for j = i+1:n_sol
            D(i,j) = norm(archive_f(i,:) - archive_f(j,:));
            D(j,i) = D(i,j);
        end
    end

    % 平均距离
    d_mean = sum(D(:)) / (n_sol * (n_sol - 1));

    % 距离标准差 / 平均距离 (归一化)
    if d_mean > 0
        delta = std(nonzeros(D(:))) / d_mean;
    else
        delta = Inf;
    end
end

% =========================================================================
% 覆盖率 (目标空间填充比例)
% =========================================================================
function coverage = computeCoverage(archive_f, ref_front)
% COMPUTECOVERAGE  计算前沿在参考前沿上的覆盖率
%   通过将目标空间分割为网格，计算被覆盖的网格比例。

    n_obj = size(archive_f, 2);

    % 确定包围盒
    f_min = min([archive_f; ref_front], [], 1);
    f_max = max([archive_f; ref_front], [], 1);
    f_range = f_max - f_min;
    f_range(f_range < 1e-6) = 1e-6;

    % 网格分辨率 (自适应)
    grid_size = min(20, ceil(1000^(1/n_obj)));
    edges = cell(1, n_obj);
    for m = 1:n_obj
        edges{m} = linspace(f_min(m), f_max(m), grid_size + 1);
    end

    % 网格占用检查
    if n_obj == 2
        [E1, E2] = meshgrid(1:grid_size, 1:grid_size);
        grid_cells = [E1(:), E2(:)];
    elseif n_obj == 3
        [E1, E2, E3] = meshgrid(1:grid_size, 1:grid_size, 1:grid_size);
        grid_cells = [E1(:), E2(:), E3(:)];
    else
        coverage = NaN;
        return;
    end

    n_occupied = 0;
    for c = 1:size(grid_cells, 1)
        cell_center = zeros(1, n_obj);
        for m = 1:n_obj
            cell_center(m) = (edges{m}(grid_cells(c,m)) + edges{m}(grid_cells(c,m)+1)) / 2;
        end
        % 检查此网格中心是否被任何解支配
        for j = 1:size(archive_f, 1)
            if all(archive_f(j,:) <= cell_center)
                n_occupied = n_occupied + 1;
                break;
            end
        end
    end

    coverage = n_occupied / size(grid_cells, 1);
end

% =========================================================================
% 评估报告输出
% =========================================================================
function printReport(metrics, func_num, ref_desc)
% PRINTREPORT  格式化输出前沿质量评估报告

    fprintf('\n');
    fprintf('============================================================\n');
    fprintf('  F%d Pareto 前沿质量评估\n', func_num);
    fprintf('============================================================\n');
    fprintf('  参考前沿: %s\n', ref_desc);
    fprintf('\n');

    % 评级判定
    if metrics.converged
        gd_grade = gradeGD(metrics.gd);
        igd_grade = gradeIGD(metrics.igd);
        spread_grade = gradeSpread(metrics.spread);

        fprintf('  ┌─────────────────────────────────────────────────────┐\n');
        fprintf('  │  指标           │  数值        │  评级              │\n');
        fprintf('  ├─────────────────────────────────────────────────────┤\n');
        fprintf('  │  解数量         │  %4d         │  —                 │\n', metrics.n_solutions);
        fprintf('  │  GD (世代距离)  │  %10.5f  │  %s  │\n', metrics.gd, gd_grade);
        fprintf('  │  IGD (反世代)   │  %10.5f  │  %s  │\n', metrics.igd, igd_grade);
        fprintf('  │  Spread (Δ)     │  %10.4f  │  %s  │\n', metrics.spread, spread_grade);
        fprintf('  │  HV (超体积)    │  %10.5f  │  —                 │\n', metrics.hv);
        fprintf('  │  覆盖率         │  %7.1f%%   │  —                 │\n', metrics.coverage * 100);
        fprintf('  └─────────────────────────────────────────────────────┘\n');

        % 综合判定
        overall = (metrics.gd < 0.1) && (metrics.igd < 0.2) && (metrics.spread < 1.0);
        if overall
            fprintf('\n  ✅ 综合评估: 前沿质量良好\n');
        else
            fprintf('\n  ⚠️  综合评估: 前沿质量需改进\n');
        end
    else
        fprintf('  ❌ 算法未收敛或解集为空\n');
        fprintf('     解数量: %d, GD: %.4f\n', metrics.n_solutions, metrics.gd);
    end

    % 各目标范围
    fprintf('\n  目标空间范围:\n');
    for m = 1:size(metrics.f_range, 2)
        fprintf('    f_%d ∈ [%.4f, %.4f]\n', m, metrics.f_range(1,m), metrics.f_range(2,m));
    end

    fprintf('============================================================\n\n');

end

function grade = gradeGD(gd)
    if gd < 0.01,      grade = '★★★★★ 极优';
    elseif gd < 0.05,  grade = '★★★★  优秀';
    elseif gd < 0.1,   grade = '★★★   良好';
    elseif gd < 0.5,   grade = '★★    可接受';
    else               grade = '★     差';
    end
end

function grade = gradeIGD(igd)
    if igd < 0.02,     grade = '★★★★★ 极优';
    elseif igd < 0.1,  grade = '★★★★  优秀';
    elseif igd < 0.2,  grade = '★★★   良好';
    elseif igd < 0.5,  grade = '★★    可接受';
    else               grade = '★     差';
    end
end

function grade = gradeSpread(spread)
    if spread < 0.3,   grade = '★★★★★ 均匀';
    elseif spread < 0.6, grade = '★★★★  较均匀';
    elseif spread < 1.0, grade = '★★★   一般';
    elseif spread < 2.0, grade = '★★    聚集';
    else                grade = '★     严重聚集';
    end
end
