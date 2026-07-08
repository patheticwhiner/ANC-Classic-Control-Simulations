function run_eMOPSO(case_key)
% RUN_EMOPSO  统一 ε-MOPSO 多目标优化入口脚本
%
%   用法: run_eMOPSO(case_key)
%
%   支持的 case_key:
%     'RST_toy'     - RST 控制器优化 (二阶教学模型)
%     'RST_armax'   - RST 控制器优化 (ARMAX 实测声学管道模型, 完整30阶)
%     'F1'          - 标准测试函数 F1 (2目标, Rastrigin型)
%     'F2'          - 标准测试函数 F2 (2目标, ZDT型, 30维)
%     'F3'          - 标准测试函数 F3 (3目标, 球坐标参数化)
%     'F4'          - 标准测试函数 F4 (3目标, 正弦振荡非连续前沿)
%
%   省略 case_key 或输入无效值时打印帮助信息。
%
%   参考:
%     About_RST_eMOPSO_spec.md — RST 控制器整定数学规范
%     testFunctions.m          — 标准测试函数定义 (附录A)

    % =====================================================================
    % 路径初始化
    % =====================================================================
    % 获取脚本所在目录，确保路径解析独立于当前工作目录
    script_path = mfilename('fullpath');
    script_dir = fileparts(script_path);
    if isempty(script_dir), script_dir = pwd; end
    addpath(fullfile(script_dir, '..', 'functions'));
    addpath(fullfile(script_dir, 'utils'));
    addpath(fullfile(script_dir, '..', 'dataset'));
    addpath(script_dir);
    output_dir = fullfile(script_dir, 'output');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    % 无参数调用：打印帮助
    if nargin < 1 || isempty(case_key)
        print_usage();
        return;
    end

    % =====================================================================
    % 统一分派
    % =====================================================================
    ck = lower(case_key);

    switch ck
        case {'rst_toy', 'rst_armax'}
            RST_engine(ck);

        case {'f1', 'f2', 'f3', 'f4'}
            func_num = str2double(ck(2));
            TestFunction_engine(func_num);

        otherwise
            fprintf('错误: 未知的 case_key ''%s''。\n\n', case_key);
            print_usage();
    end
end

% =========================================================================
% 引擎 1: RST 控制器优化引擎
% =========================================================================
function RST_engine(case_key)
% RST_ENGINE  RST 控制器 ε-MOPSO 优化引擎
%   处理被控对象定义、约束惩罚、优化求解、控制器合成与验证全流程。

    % 与主函数一致的 output 目录
    script_path = mfilename('fullpath');
    script_dir = fileparts(script_path);
    if isempty(script_dir), script_dir = pwd; end
    output_dir = fullfile(script_dir, 'output');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    fprintf('===========================================================\n');
    fprintf('   ε-MOPSO RST 控制器优化 — %s\n', upper(case_key));
    fprintf('===========================================================\n\n');

    %% === 第1节: 被控对象定义 ===

    if strcmp(case_key, 'rst_toy')
        % 二阶教学模型
        modelFile = fullfile(script_dir, '..', 'dataset', 'syn_RSTtoy_2nd.mat');
        load(modelFile, 'model');
        B  = model.B_poly;
        A  = model.A_poly;
        d  = model.d_delay;
        Ts = model.Ts;
        fprintf('模型来源: 二阶教学模型 (syn_RSTtoy_2nd.mat)\n');
    else  % rst_armax
        % 加载完整 ARMAX(30,30,30,22) 实测模型（不降阶）
        [B, A, d, Ts] = load_armax_full();
    end

    fprintf('A(z^{-1}) = ');
    fprintf('%+.4f z^{-%d} ', [A; (0:length(A)-1)]);
    fprintf('\nB(z^{-1}) = ');
    fprintf('%+.4f z^{-%d} ', [B; (0:length(B)-1)]);
    fprintf('\n延迟 d=%d, Ts=%.4f\n', d, Ts);

    %% === 第2节: 系统分析 — 提取不稳定零极点 ===

    sys_zeros = roots(B);
    sys_poles = roots(A);

    beta  = sys_zeros(abs(sys_zeros) > 1);
    alpha = sys_poles(abs(sys_poles) > 1);

    if isempty(beta)
        fprintf('\n注意: 无不稳定零点，使用虚拟零点 β=1.05。\n');
        beta = 1.05;
    end
    if isempty(alpha)
        fprintf('注意: 无不稳定极点，使用虚拟极点 α=0.95。\n');
        alpha = 0.95;
    end

    fprintf('\n不稳定零点: %d 个\n', length(beta));
    for i = 1:length(beta)
        fprintf('  β_%d = %.4f+%.4fi, |β|=%.4f\n', ...
            i, real(beta(i)), imag(beta(i)), abs(beta(i)));
    end
    fprintf('不稳定极点: %d 个\n', length(alpha));

    %% === 第3节: RST 预设参数 ===

    % 主导极点 P_D (二阶 Butterworth 型)
    zeta_P  = 0.8;
    omega_P = 0.25 * pi;
    PD1 = -2 * exp(-zeta_P * omega_P) * cos(omega_P * sqrt(1 - zeta_P^2));
    PD2 = exp(-2 * zeta_P * omega_P);
    P_D = [1, PD1, PD2];

    % 控制器固定因子
    H_R = [1, 1];
    H_S = [1, -1];

    % X/Y 多项式阶数选择
    %   设计原理 (基于实验验证):
    %     RST_toy  (二阶教学): nX=2, nY=2, 总变量 4, P_DX 4阶 → 1个零约束 ✅
    %     RST_armax (30阶实测): nextra = nah+nbh-1 - (deg(P_D)+nX)
    %       nX=5  → nextra=75 (90%), max|R|=5.1×10⁵ ❌
    %       nX=25 → nextra=55 (66%), max|R|=6.7×10⁴ ↓87%
    %       nX=50 → nextra=30 (36%), max|R|预计~10³  目标
    %       物理: 每减少1个零约束 = 少一次Toeplitz递归放大
    if strcmp(case_key, 'rst_toy')
        nX = 2;
        nY = 2;
    else  % rst_armax
        nX = 50;
        nY = 12;
    end

    fprintf('\nP_D: 1 %+.4f z^{-1} %+.4f z^{-2}\n', PD1, PD2);
    fprintf('H_R = [1, 1], H_S = [1, -1]\n');
    fprintf('nX=%d, nY=%d (整形多项式阶数)\n', nX, nY);

    %% === 第4节: ε-MOPSO 参数配置 ===

    n_var = nX + nY;
    n_obj = length(beta) + 2;  % β匹配 + 延迟积分 + Bezout残差

    lb = -0.95 * ones(1, n_var);
    ub =  0.95 * ones(1, n_var);

    % 根据问题规模自适应调整参数
    if strcmp(case_key, 'rst_toy')
        n_pop_val = 40;
        k_max_val = 80;
        eps_val   = 0.02;
    else  % rst_armax — 62维(nX=50+nY=12) + 9个β + 2个附加目标 = 11目标
        n_pop_val = 200;
        k_max_val = 400;
        eps_val   = 0.02;
    end

    options = struct();
    options.n_pop       = n_pop_val;
    options.k_max       = k_max_val;
    options.epsilon     = eps_val;
    options.c1          = 1.5;
    options.c2          = 1.7;
    options.w_max       = 0.9;
    options.w_min       = 0.3;
    options.archive_max = 80;
    options.verbose     = true;

    fprintf('\nε-MOPSO: n_var=%d, n_obj=%d, n_pop=%d, k_max=%d, ε=%.3f\n', ...
        n_var, n_obj, n_pop_val, k_max_val, eps_val);

    %% === 第5节: 构建系统信息与目标函数 ===

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

    % 包装目标函数（含约束惩罚）
    obj_func = @(theta) wrapped_objective(theta, sys_info);

    fprintf('\n========== 开始优化 ==========\n');

    %% === 第6节: 运行优化 ===

    tic;
    [archive, archive_f, history] = eMOPSO_core(obj_func, n_var, n_obj, lb, ub, options);
    elapsed_time = toc;

    fprintf('\n========== 优化完成 ==========\n');
    fprintf('耗时: %.2f s, Pareto 解数: %d\n', elapsed_time, size(archive, 1));

    %% === 第7节: Pareto 前沿可视化 ===

    figure('Name', sprintf('Pareto Front — %s', case_key));

    if n_obj == 2
        scatter(archive_f(:,1), archive_f(:,2), 40, 'b', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        xlabel('f_1'); ylabel('f_2');
    elseif n_obj >= 3
        scatter3(archive_f(:,1), archive_f(:,2), archive_f(:,3), 40, ...
            archive_f(:,3), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
        xlabel('f_1'); ylabel('f_2'); zlabel('f_3');
        view(45, 30); colorbar;
    end
    title(sprintf('Pareto Front — %s (ε=%.3f, %d solutions)', ...
        case_key, options.epsilon, size(archive_f,1)));
    grid on;

    % 收敛曲线
    figure('Name', sprintf('Convergence — %s', case_key));
    plot(history, 'b-', 'LineWidth', 1.5);
    xlabel('Iteration'); ylabel('Archive Size');
    title(sprintf('Convergence — %s', case_key));
    grid on;

    %% === 第8节: 选择 Pareto 最优解 ===

    f_range = max(archive_f) - min(archive_f);
    f_range(f_range == 0) = 1;
    f_norm = (archive_f - min(archive_f)) ./ f_range;
    [~, knee_idx] = min(sum(f_norm, 2));
    theta_opt = archive(knee_idx, :)';
    f_opt = archive_f(knee_idx, :)';

    fprintf('\n拐点解索引: %d/%d\n', knee_idx, size(archive,1));
    fprintf('目标值: '); fprintf('%.6f ', f_opt); fprintf('\n');
    fprintf('X roots: '); fprintf('%.4f ', theta_opt(1:nX)); fprintf('\n');
    fprintf('Y roots: '); fprintf('%.4f ', theta_opt(nX+1:end)); fprintf('\n');

    %% === 第9节: 保存结果 ===

    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    switch case_key
        case 'rst_toy'
            out_file = fullfile(output_dir, 'RST_eMOPSO_results.mat');
        case 'rst_armax'
            out_file = fullfile(output_dir, 'RST_eMOPSO_armax_results.mat');
    end

    save(out_file, 'archive', 'archive_f', 'history', ...
        'theta_opt', 'f_opt', 'sys_info', 'options', 'elapsed_time');
    fprintf('\n结果已保存到: %s\n', out_file);

    %% === 第10节: 控制器合成 ===

    fprintf('\n========== 控制器合成 ==========\n');
    [R, S, T, ~] = postprocess_RST(theta_opt, sys_info, 'Plot', true);

    %% === 第11节: 导出为 controller_validator 兼容格式 ===

    export_for_validator(out_file, 'output');

    %% === 第12节: 自动验证控制器性能 ===

    validate_results = validate_controller(out_file, 'Plot', true);
end

% =========================================================================
% ARMAX 模型加载（完整30阶，不降阶）
% =========================================================================
function [B, A, d, Ts] = load_armax_full()
% LOAD_ARMAX_FULL  加载 ARMAX(30,30,30,22) 实测模型，直接使用完整30阶
%
%   输出:
%     B, A  - 分子/分母多项式系数 (降幂), 完整30阶
%     d     - 延迟步数
%     Ts    - 采样时间

    fprintf('正在加载 ARMAX(30,30,30,22) 实测模型 (完整30阶, 不降阶)...\n');
    info = DataManager('armax_30303022');

    if ~strcmp(info.type, 'armax_model')
        error('数据源类型错误: 期望 armax_model, 实际 %s', info.type);
    end

    % 提取多项式
    A = info.model.A;
    B_full = info.model.B;
    d  = info.orders(4);
    Ts = 1 / info.fs;

    fprintf('  原始模型: ARMAX(%d,%d,%d,%d), fs=%d Hz\n', ...
        info.orders(1), info.orders(2), info.orders(3), d, info.fs);

    % 去除延迟前导零 (B_full 前 d 项为零, 对应 z^0 ~ z^{-(d-1)})
    B = B_full(d+1:end);

    % 确保 A 首项为 1
    if abs(A(1) - 1) > 1e-10
        B = B / A(1);
        A = A / A(1);
    end

    fprintf('  A 阶数: %d, B 阶数: %d (含延迟 d=%d)\n', length(A)-1, length(B)-1, d);
end

% =========================================================================
% 引擎 2: 标准测试函数引擎
% =========================================================================
function TestFunction_engine(func_num)
% TESTFUNCTION_ENGINE  标准多目标测试函数 ε-MOPSO 验证引擎
%   对 F1-F4 标准测试函数运行 ε-MOPSO 算法。
%   这些函数没有约束，直接传递原始目标值给 eMOPSO_core。

    % 与主函数一致的 output 目录
    script_path = mfilename('fullpath');
    script_dir = fileparts(script_path);
    if isempty(script_dir), script_dir = pwd; end
    output_dir = fullfile(script_dir, 'output');
    if ~exist(output_dir, 'dir'), mkdir(output_dir); end

    fprintf('===========================================================\n');
    fprintf('   ε-MOPSO 标准测试 — F%d\n', func_num);
    fprintf('===========================================================\n\n');

    %% === 问题维度配置 ===

    switch func_num
        case 1  % F1: 2目标, n=2, x₂∈[-30,30] (Appendix A.1)
            n_var = 2;
            n_obj = 2;
            fprintf('F1: 2目标, %d维, Rastrigin型 g 函数 (x2∈[-30,30])\n', n_var);
            fprintf('Pareto前沿理论: f2 = 1 - sqrt(f1), f1∈[0,1]\n');
        case 2  % F2: 2目标, n=30, ZDT型+振荡h (Appendix A.2)
            n_var = 30;
            n_obj = 2;
            fprintf('F2: 2目标, %d维, ZDT型 g + 振荡 h 函数\n', n_var);
            fprintf('Pareto前沿理论: f2 = 2 - sqrt(f1) - f1*sin(10πf1)\n');
        case 3  % F3: 3目标, n=3, r=1+x₃∈[1,2] (Appendix A.3)
            n_var = 3;
            n_obj = 3;
            fprintf('F3: 3目标, %d维, 球坐标 (r=1+x₃∈[1,2])\n', n_var);
            fprintf('Pareto前沿理论: 单位球面第一象限 f1^2+f2^2+f3^2=1\n');
        case 4  % F4: 3目标, n=3, f3=3.5-Σ2xi·sin(3πxi) (Appendix A.4)
            n_var = 3;
            n_obj = 3;
            fprintf('F4: 3目标, %d维, 振荡曲面\n', n_var);
            fprintf('Pareto前沿理论: f3 = 3.5 - 2f1·sin(3πf1) - 2f2·sin(3πf2)\n');
    end

    %% === 搜索空间 ===

    lb = zeros(1, n_var);
    ub = ones(1, n_var);

    %% === ε-MOPSO 自适应参数 ===

    options = struct();
    options.c1          = 1.5;
    options.c2          = 1.7;
    options.w_max       = 0.9;
    options.w_min       = 0.4;
    options.verbose     = false;
    options.p_mut       = 0;      % 默认不启用多项式变异
    options.eta_m       = 20;

    switch func_num
        case 1  % F1: 默认参数即可
            options.n_pop       = 100;
            options.k_max       = 200;
            options.epsilon     = 0.01;
            options.archive_max = 100;

        case 2  % F2: 30维高维 — 增大种群+迭代+变异
            options.n_pop       = 200;
            options.k_max       = 500;
            options.epsilon     = 0.01;
            options.archive_max = 100;
            options.p_mut       = 0.1;   % 多项式变异维持多样性

        case 3  % F3: 3目标球面 — 增大Archive
            options.n_pop       = 100;
            options.k_max       = 200;
            options.epsilon     = 0.01;
            options.archive_max = 150;

        case 4  % F4: 3目标非连续 — 减小ε+增大Archive
            options.n_pop       = 100;
            options.k_max       = 200;
            options.epsilon     = 0.005;  % 更细的ε网格
            options.archive_max = 200;    % 更大容量防止饱和
    end

    fprintf('\nε-MOPSO: n_pop=%d, k_max=%d, ε=%.4f, p_mut=%.2f\n', ...
        options.n_pop, options.k_max, options.epsilon, options.p_mut);

    %% === 构建目标函数 (直接使用 testFunctions, 无约束) ===

    obj_func = @(x) testFunctions(x, func_num);

    fprintf('\n========== 开始优化 ==========\n');

    %% === 运行优化 ===

    tic;
    [archive, archive_f, history] = eMOPSO_core(obj_func, n_var, n_obj, lb, ub, options);
    elapsed_time = toc;

    fprintf('========== 优化完成 ==========\n');
    fprintf('耗时: %.2f s, Pareto 解数: %d\n', elapsed_time, size(archive, 1));

    %% === Pareto 前沿可视化 ===

    figure('Name', sprintf('Pareto Front — F%d', func_num));

    if n_obj == 2
        scatter(archive_f(:,1), archive_f(:,2), 30, 'r', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
        xlabel(sprintf('f_1')); ylabel(sprintf('f_2'));
    elseif n_obj == 3
        scatter3(archive_f(:,1), archive_f(:,2), archive_f(:,3), 30, ...
            archive_f(:,3), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
        xlabel('f_1'); ylabel('f_2'); zlabel('f_3');
        view(45, 30); colorbar;
    end
    title(sprintf('Pareto Front — F%d (ε=%.3f, %d solutions)', ...
        func_num, options.epsilon, size(archive_f,1)));
    grid on;

    % 收敛曲线
    figure('Name', sprintf('Convergence — F%d', func_num));
    plot(history, 'r-', 'LineWidth', 1.5);
    xlabel('Iteration'); ylabel('Archive Size');
    title(sprintf('Convergence — F%d', func_num));
    grid on;

    %% === 保存结果 ===

    out_file = fullfile(output_dir, sprintf('test_F%d_results.mat', func_num));
    save(out_file, 'archive', 'archive_f', 'history', ...
        'options', 'elapsed_time', 'func_num');
    fprintf('\n结果已保存到: %s\n', out_file);

    %% === 前沿质量评估 (定量指标) ===

    fprintf('\n--- F%d 前沿统计 ---\n', func_num);
    f_min = min(archive_f, [], 1);
    f_max = max(archive_f, [], 1);
    for m = 1:n_obj
        fprintf('  f_%d ∈ [%.4f, %.4f]\n', m, f_min(m), f_max(m));
    end

    % 调用量化评估工具
    metrics = evaluateFrontier(archive_f, func_num);

    % 叠加参考前沿到 Pareto 图（F1/F2 为2目标，F3/F4 为3目标）
    figure('Name', sprintf('Frontier Quality — F%d', func_num));
    if n_obj == 2
        % 实际前沿
        scatter(archive_f(:,1), archive_f(:,2), 30, 'r', 'filled', ...
            'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
        hold on;
        % 理论参考前沿
        f1_ref = linspace(f_min(1), f_max(1), 500);
        switch func_num
            case 1, f2_ref = 1 - sqrt(f1_ref);
            case 2, f2_ref = 2 - sqrt(f1_ref) - f1_ref .* sin(10 * pi * f1_ref);
        end
        plot(f1_ref, f2_ref, 'b-', 'LineWidth', 2);
        legend('ε-MOPSO', '理论Pareto前沿', 'Location', 'best');
        xlabel('f_1'); ylabel('f_2');
        title(sprintf('F%d | GD=%.4f | IGD=%.4f | Spread=%.3f | %d解', ...
            func_num, metrics.gd, metrics.igd, metrics.spread, metrics.n_solutions));
    elseif n_obj == 3
        % 实际前沿
        scatter3(archive_f(:,1), archive_f(:,2), archive_f(:,3), 30, ...
            archive_f(:,3), 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.3);
        hold on;
        % 理论参考前沿 (内联生成，3目标)
        n_per_dim = ceil(sqrt(200));
        switch func_num
            case 3
                % F3: 单位球面第一象限 (r=1, x₃=0)
                theta = linspace(0, pi/2, n_per_dim)';
                phi   = linspace(0, pi/2, n_per_dim)';
                [T, P] = meshgrid(theta, phi);
                rf = [cos(T(:)).*cos(P(:)), cos(T(:)).*sin(P(:)), sin(T(:))];
                if size(rf,1) > 200
                    idx = randperm(size(rf,1), 200);
                    rf = rf(idx, :);
                end
            case 4
                % F4: f3 = 3.5 - 2f1·sin(3πf1) - 2f2·sin(3πf2)  (x₃=0)
                f1_r = linspace(0, 1, n_per_dim)';
                f2_r = linspace(0, 1, n_per_dim)';
                [F1, F2] = meshgrid(f1_r, f2_r);
                F3_min = 3.5 - 2*F1(:).*sin(3*pi*F1(:)) - 2*F2(:).*sin(3*pi*F2(:));
                rf = [F1(:), F2(:), F3_min];
                if size(rf,1) > 200
                    idx = randperm(size(rf,1), 200);
                    rf = rf(idx, :);
                end
        end
        scatter3(rf(:,1), rf(:,2), rf(:,3), 15, ...
            [0 0.4 0.8], 'o', 'MarkerFaceAlpha', 0.3, 'MarkerEdgeColor', 'none');
        legend('ε-MOPSO', '理论Pareto前沿', 'Location', 'best');
        xlabel('f_1'); ylabel('f_2'); zlabel('f_3');
        view(45, 30); colorbar;
        title(sprintf('F%d | GD=%.4f | IGD=%.4f | Spread=%.3f | %d解', ...
            func_num, metrics.gd, metrics.igd, metrics.spread, metrics.n_solutions));
    end
    grid on;
end

% =========================================================================
% 本地包装函数: 目标函数 + 约束惩罚 (RST 引擎专用)
% =========================================================================
function phi = wrapped_objective(theta, sys_info)
% WRAPPED_OBJECTIVE  包装 RST 目标函数并施加约束惩罚
%
%   先放大原始目标 1000 倍以匹配 epsilon 量级，
%   再使用乘性惩罚 (penalty_function) 处理约束违例。

    theta = theta(:);

    % 原始目标函数 × 1000
    f = RST_objective(theta, sys_info) * 1000;

    % 约束计算
    [g1, g2] = RST_constraints(theta, sys_info);

    % 施加乘性惩罚 (默认 Lambda=0.5)
    phi = penalty_function(f, g1, g2);
end

% =========================================================================
% 帮助信息
% =========================================================================
function print_usage()
    fprintf('===========================================================\n');
    fprintf('  run_eMOPSO — 统一 ε-MOPSO 多目标优化入口\n');
    fprintf('===========================================================\n\n');
    fprintf('用法: run_eMOPSO(''<case_key>'')\n\n');
    fprintf('可用 case_key:\n');
    fprintf('  ''RST_toy''     — RST 控制器优化 (二阶教学模型)\n');
    fprintf('  ''RST_armax''   — RST 控制器优化 (ARMAX 实测模型, 完整30阶)\n');
    fprintf('  ''F1''          — 标准测试函数 F1 (2目标, Rastrigin型)\n');
    fprintf('  ''F2''          — 标准测试函数 F2 (2目标, ZDT型, 30维)\n');
    fprintf('  ''F3''          — 标准测试函数 F3 (3目标, 球坐标参数化)\n');
    fprintf('  ''F4''          — 标准测试函数 F4 (3目标, 正弦振荡非连续前沿)\n\n');
    fprintf('示例:\n');
    fprintf('  run_eMOPSO(''RST_toy'')\n');
    fprintf('  run_eMOPSO(''F1'')\n');
    fprintf('  run_eMOPSO(''RST_armax'')\n\n');
    fprintf('输出文件位于 output/ 目录:\n');
    fprintf('  RST_eMOPSO_results.mat         — RST_toy 优化结果\n');
    fprintf('  RST_eMOPSO_armax_results.mat   — RST_armax 优化结果\n');
    fprintf('  test_F1_results.mat ~ F4       — 标准测试函数优化结果\n');
end
