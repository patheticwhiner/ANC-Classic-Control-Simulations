function [archive, archive_f, history] = eMOPSO_core(obj_func, n_var, n_obj, lb, ub, options)
% EMOPSO_CORE ε-多目标粒子群优化算法核心
%   实现了基于 ε-支配 的 MOPSO 算法，用于求解多目标优化问题的
%   Pareto 最优解集。
%
%   输入:
%     obj_func  - 目标函数句柄，接口: f = obj_func(x)
%                 其中 x 为列向量 (n_var × 1), f 为列向量 (n_obj × 1)
%     n_var     - 决策变量维度
%     n_obj     - 目标函数数量
%     lb        - 决策变量下界 (行向量或列向量, 长度 n_var)
%     ub        - 决策变量上界 (行向量或列向量, 长度 n_var)
%     options   - 算法参数结构体 (可选), 包含以下字段:
%       .n_pop       - 种群大小 (默认: 50)
%       .k_max       - 最大迭代次数 (默认: 200)
%       .epsilon     - ε支配阈值 (默认: 0.01)
%       .c1          - 认知加速系数 (默认: 1.5)
%       .c2          - 社会加速系数 (默认: 1.7)
%       .w_max       - 最大惯性权重 (默认: 0.9)
%       .w_min       - 最小惯性权重 (默认: 0.4)
%       .archive_max - Archive最大容量 (默认: 100)
%       .verbose     - 是否显示进度 (默认: true)
%
%   输出:
%     archive    - Archive中的决策变量矩阵 (n_archive × n_var)
%     archive_f  - Archive中的目标函数值矩阵 (n_archive × n_obj)
%     history    - 每代Archive大小记录 (k_max × 1)
%
%   参考文献:
%     [1] Sierra, M.R. and Coello Coello, C.A., "Improving PSO-Based
%         Multi-objective Optimization Using Crowding, Mutation and
%         ε-Dominance", EMO 2005, LNCS 3410, pp. 505-519, 2005.
%     [2] Mostaghim, S. and Teich, J., "Strategies for Finding Good
%         Local Guides in Multi-objective Particle Swarm Optimization
%         (MOPSO)", SIS 2003, pp. 26-33, 2003.

    % =====================================================================
    % 输入参数验证与默认值设置
    % =====================================================================
    if nargin < 6, options = struct(); end
    if nargin < 5
        error('eMOPSO_core:NotEnoughInputs', ...
              '至少需要5个输入参数。');
    end

    % 确保 lb, ub 为行向量
    lb = lb(:)';
    ub = ub(:)';

    % 检查维度一致性
    if length(lb) ~= n_var || length(ub) ~= n_var
        error('eMOPSO_core:DimensionMismatch', ...
              'lb和ub的长度必须等于 n_var (%d)。', n_var);
    end

    % 设置默认参数值
    if ~isfield(options, 'n_pop'),       options.n_pop       = 50;   end
    if ~isfield(options, 'k_max'),       options.k_max       = 200;  end
    if ~isfield(options, 'epsilon'),     options.epsilon     = 0.01; end
    if ~isfield(options, 'c1'),          options.c1          = 1.5;  end
    if ~isfield(options, 'c2'),          options.c2          = 1.7;  end
    if ~isfield(options, 'w_max'),       options.w_max       = 0.9;  end
    if ~isfield(options, 'w_min'),       options.w_min       = 0.4;  end
    if ~isfield(options, 'archive_max'), options.archive_max = 100;  end
    if ~isfield(options, 'verbose'),     options.verbose     = true; end

    % 提取参数
    n_pop       = options.n_pop;
    k_max       = options.k_max;
    epsilon     = options.epsilon;
    c1          = options.c1;
    c2          = options.c2;
    w_max       = options.w_max;
    w_min       = options.w_min;
    archive_max = options.archive_max;
    verbose     = options.verbose;

    % =====================================================================
    % 步骤1: 初始化种群
    % =====================================================================
    rng_val = rng;  % 保存随机数状态以便复现

    % 在上下界范围内均匀随机初始化
    pop = lb + rand(n_pop, n_var) .* (ub - lb);         % (n_pop × n_var)

    % 初始化速度为0
    vel = zeros(n_pop, n_var);                           % (n_pop × n_var)

    % =====================================================================
    % 步骤2: 评估初始种群
    % =====================================================================
    pop_f = zeros(n_pop, n_obj);
    for i = 1:n_pop
        pop_f(i, :) = obj_func(pop(i, :)')';            % 转置为列向量再转回行
    end

    % =====================================================================
    % 步骤3: 初始化个体最优 (pbest)
    % =====================================================================
    pbest    = pop;                                      % (n_pop × n_var)
    pbest_f  = pop_f;                                    % (n_pop × n_obj)

    % =====================================================================
    % 步骤4: 初始化Archive (外部存档)
    % =====================================================================
    [archive, archive_f] = initialize_archive(pop, pop_f, epsilon, archive_max);

    % =====================================================================
    % 步骤5: 主循环
    % =====================================================================
    history = zeros(k_max, 1);

    if verbose
        fprintf('=====================================================\n');
        fprintf('  ε-MOPSO 优化开始\n');
        fprintf('  种群大小: %d, 最大迭代: %d, ε: %.4f\n', ...
                n_pop, k_max, epsilon);
        fprintf('=====================================================\n');
    end

    for k = 1:k_max
        % ---------------------------------------------------------
        % 5.1 计算当前惯性权重 (线性递减)
        % ---------------------------------------------------------
        w = w_max - (w_max - w_min) * k / k_max;

        % ---------------------------------------------------------
        % 5.2 对每个粒子进行更新
        % ---------------------------------------------------------
        for i = 1:n_pop
            % --- 5.2.1 选择全局引导 (gbest) ---
            gbest = select_gbest(archive, archive_f, ...
                                 pbest_f(i, :), epsilon);

            % --- 5.2.2 更新速度和位置 ---
            [new_x, new_v] = update_particle(pop(i, :), vel(i, :), ...
                                             pbest(i, :), gbest, ...
                                             w, c1, c2, lb, ub);

            % --- 5.2.3 评估新位置 ---
            new_f = obj_func(new_x')';                   % 转为行向量

            % --- 5.2.4 更新个体最优 (ε-支配判断) ---
            if epsilon_dominates(new_f', pbest_f(i, :)', epsilon)
                pbest(i, :)   = new_x;
                pbest_f(i, :) = new_f;
            elseif ~epsilon_dominates(pbest_f(i, :)', new_f', epsilon)
                % 互不支配时，随机选择
                if rand() < 0.5
                    pbest(i, :)   = new_x;
                    pbest_f(i, :) = new_f;
                end
            end

            % --- 5.2.5 用新解更新Archive ---
            [archive, archive_f] = update_archive(archive, archive_f, ...
                                                   new_x, new_f, ...
                                                   epsilon, archive_max);

            % 更新种群
            pop(i, :) = new_x;
            vel(i, :) = new_v;
        end

        % ---------------------------------------------------------
        % 5.3 记录当前Archive大小
        % ---------------------------------------------------------
        history(k) = size(archive, 1);

        % ---------------------------------------------------------
        % 5.4 显示进度
        % ---------------------------------------------------------
        if verbose && (mod(k, 10) == 0 || k == 1 || k == k_max)
            fprintf('  迭代 %4d/%d, Archive大小: %3d\n', ...
                    k, k_max, history(k));
        end
    end

    if verbose
        fprintf('=====================================================\n');
        fprintf('  ε-MOPSO 优化完成\n');
        fprintf('  最终Archive大小: %d\n', size(archive, 1));
        fprintf('=====================================================\n');
    end

end

% =========================================================================
% 子函数: ε-支配判断
% =========================================================================
function result = epsilon_dominates(f1, f2, epsilon)
% EPSILON_DOMINATES 判断 f1 是否 ε-支配 f2
%   f1 ε-支配 f2 ⇔ ∀m: f1(m)/(1+ε) ≤ f2(m) 且 ∃m: f1(m)/(1+ε) < f2(m)
%
%   输入:
%     f1, f2  - 目标函数值列向量
%     epsilon - ε值
%   输出:
%     result  - true 如果 f1 ε-支配 f2

    % 确保为列向量
    f1 = f1(:);
    f2 = f2(:);

    % 缩放因子
    scale = 1 + epsilon;

    % 检查条件1: ∀m, f1(m)/scale ≤ f2(m)
    cond1 = all(f1 / scale <= f2);

    % 检查条件2: ∃m, f1(m)/scale < f2(m)
    cond2 = any(f1 / scale < f2);

    result = cond1 && cond2;
end

% =========================================================================
% 子函数: Box索引计算
% =========================================================================
function B = box_index(f, epsilon)
% BOX_INDEX 计算目标向量在ε-网格中的Box索引
%   B_m = floor(log(f_m) / log(1 + epsilon))
%
%   输入:
%     f       - 目标函数值向量 (行向量或列向量)
%     epsilon - ε值
%   输出:
%     B       - Box索引行向量

    % 防止 log(0)
    f = max(f, 1e-12);

    B = floor(log(f(:)') / log(1 + epsilon));
end

% =========================================================================
% 子函数: 初始化Archive
% =========================================================================
function [archive, archive_f] = initialize_archive(pop, pop_f, epsilon, archive_max)
% INITIALIZE_ARCHIVE 从初始种群构建Archive（使用ε-支配规则）

    [n_pop, n_var] = size(pop);
    n_obj = size(pop_f, 2);

    % 使用临时列表构建
    archive    = zeros(0, n_var);
    archive_f  = zeros(0, n_obj);

    for i = 1:n_pop
        [archive, archive_f] = update_archive(archive, archive_f, ...
                                               pop(i, :), pop_f(i, :), ...
                                               epsilon, archive_max);
    end
end

% =========================================================================
% 子函数: Archive更新
% =========================================================================
function [archive, archive_f] = update_archive(archive, archive_f, ...
                                                new_x, new_f, epsilon, archive_max)
% UPDATE_ARCHIVE 用新解更新Archive
%   更新规则:
%     1. 如果new_x与Archive中某解在同一Box，保留目标值更小的
%     2. 如果new_x ε-支配某Archive解的Box，替换之
%     3. 如果超过archive_max，随机删除拥挤Box中的解

    n_archive = size(archive, 1);

    % 确保 new_x, new_f 为行向量
    new_x = new_x(:)';
    new_f = new_f(:)';

    % 计算新解的Box索引
    new_box = box_index(new_f, epsilon);

    % 标记需要删除的Archive索引
    to_remove = false(n_archive, 1);

    for j = 1:n_archive
        arch_box = box_index(archive_f(j, :), epsilon);

        % 规则1: 同一Box — 保留目标值更小的
        if isequal(new_box, arch_box)
            if sum(new_f) < sum(archive_f(j, :))
                to_remove(j) = true;
            else
                % 新解不优于Archive解，不添加
                return;
            end
        end

        % 规则2: 新解 ε-支配 Archive解
        if epsilon_dominates(new_f', archive_f(j, :)', epsilon)
            to_remove(j) = true;
        end

        % 如果Archive解 ε-支配新解，则不添加新解
        if epsilon_dominates(archive_f(j, :)', new_f', epsilon)
            return;
        end
    end

    % 删除被支配的解
    archive(to_remove, :)   = [];
    archive_f(to_remove, :) = [];

    % 添加新解
    archive   = [archive; new_x];
    archive_f = [archive_f; new_f];

    % 规则3: 容量控制 — 随机删除拥挤Box中的解
    if size(archive, 1) > archive_max
        % 计算每个解的Box拥挤度（同Box中的解数量）
        n_arch = size(archive, 1);
        boxes = zeros(n_arch, 1);
        box_map = containers.Map('KeyType', 'char', 'ValueType', 'double');

        for j = 1:n_arch
            b_key = mat2str(box_index(archive_f(j, :), epsilon));
            if isKey(box_map, b_key)
                box_map(b_key) = box_map(b_key) + 1;
            else
                box_map(b_key) = 1;
            end
        end

        % 找出最拥挤的Box
        keys_list = box_map.keys;
        max_count = 0;
        max_key = '';
        for jj = 1:length(keys_list)
            if box_map(keys_list{jj}) > max_count
                max_count = box_map(keys_list{jj});
                max_key = keys_list{jj};
            end
        end

        % 从最拥挤的Box中随机删除一个解
        candidates = [];
        for j = 1:n_arch
            if isequal(mat2str(box_index(archive_f(j, :), epsilon)), max_key)
                candidates = [candidates; j]; %#ok<AGROW>
            end
        end

        if ~isempty(candidates)
            remove_idx = candidates(randi(length(candidates)));
            archive(remove_idx, :)   = [];
            archive_f(remove_idx, :) = [];
        else
            % 回退: 随机删除
            remove_idx = randi(size(archive, 1));
            archive(remove_idx, :)   = [];
            archive_f(remove_idx, :) = [];
        end
    end
end

% =========================================================================
% 子函数: gbest选择 (Sigma方法)
% =========================================================================
function gbest = select_gbest(archive, archive_f, particle_f, epsilon) %#ok<INUSD>
% SELECT_GBEST 使用Sigma方法选择全局引导
%   对粒子计算其sigma向量:
%     σ_i = (f1² - f2², f1² - f3², ..., f_{M-1}² - f_M²)
%   在Archive中选择sigma向量欧氏距离最近的粒子作为gbest

    n_archive = size(archive, 1);

    if n_archive == 0
        error('eMOPSO_core:EmptyArchive', 'Archive为空，无法选择gbest。');
    end

    if n_archive == 1
        gbest = archive(1, :);
        return;
    end

    n_obj = size(archive_f, 2);

    % 计算粒子的sigma向量
    sigma_particle = compute_sigma(particle_f(:)');

    % 计算Archive中每个解的sigma向量
    sigma_archive = zeros(n_archive, n_obj - 1);
    for j = 1:n_archive
        sigma_archive(j, :) = compute_sigma(archive_f(j, :));
    end

    % 计算欧氏距离
    diff = sigma_archive - repmat(sigma_particle, n_archive, 1);
    distances = sqrt(sum(diff.^2, 2));

    % 选择距离最近的
    [~, idx] = min(distances);
    gbest = archive(idx, :);
end

function sigma = compute_sigma(f)
% COMPUTE_SIGMA 计算目标向量f的sigma值
%   σ = (f1² - f2², f1² - f3², ..., f_{M-1}² - f_M²)
    M = length(f);
    if M < 2
        sigma = 0;
        return;
    end
    sigma = zeros(1, M - 1);
    for m = 1:(M - 1)
        sigma(m) = f(1)^2 - f(m + 1)^2;
    end
end

% =========================================================================
% 子函数: 单个粒子更新
% =========================================================================
function [new_x, new_v] = update_particle(x, v, pbest, gbest, w, c1, c2, lb, ub)
% UPDATE_PARTICLE 更新单个粒子的速度和位置
%   v_new = w*v + c1*r1*(pbest - x) + c2*r2*(gbest - x)
%   x_new = x + v_new
%   边界处理: 反射策略

    r1 = rand(size(x));
    r2 = rand(size(x));

    % 速度更新
    new_v = w * v + c1 * r1 .* (pbest - x) + c2 * r2 .* (gbest - x);

    % 位置更新
    new_x = x + new_v;

    % 边界反射处理
    for d = 1:length(x)
        while new_x(d) < lb(d) || new_x(d) > ub(d)
            if new_x(d) < lb(d)
                % 下界反射
                new_x(d) = lb(d) + (lb(d) - new_x(d));
                new_v(d) = -new_v(d);  % 速度反向
            elseif new_x(d) > ub(d)
                % 上界反射
                new_x(d) = ub(d) - (new_x(d) - ub(d));
                new_v(d) = -new_v(d);  % 速度反向
            end
        end
    end
end
