function results = analyze_armax_model(modelSource)
% analyze_armax_model - ARMAX模型分析：能观标准型转换、互质性检验、零极点分析
%
% 用法:
%   results = analyze_armax_model('armax_30303022')
%
% 输入:
%   modelSource - DataManager 可识别的 ARMAX 模型名称 (如 'armax_30303022')
%
% 输出:
%   results - struct，包含:
%     .Ao, .Bo, .Co, .Do       - 能观标准型状态空间矩阵
%     .sys_obs                  - 能观标准型 ss 对象
%     .is_coprime               - A(z), B(z) 是否互质
%     .is_strongly_stabilizable - 是否强可镇定
%     .poles, .zeros            - 零极点位置
%     .obsv_rank, .obsv_cond    - 能观性矩阵秩与条件数
%     .plant_tf                 - 被控对象传递函数

    if nargin < 1
        modelSource = 'armax_30303022';
    end

    %% 1. 加载模型
    fprintf('=== ARMAX 模型分析: %s ===\n', modelSource);
    info = DataManager(modelSource);

    if ~strcmp(info.type, 'armax_model')
        error('数据源 ''%s'' 不是 ARMAX 模型 (type=%s)', modelSource, info.type);
    end

    A_poly = info.model.A;
    B_poly = info.model.B;
    d      = info.orders(4);
    B_poly = B_poly(d+1:end);
    nA     = info.orders(1);
    nB     = info.orders(2);
    fs     = info.fs;
    ts     = 1/fs;

    fprintf('模型: ARMAX(%d,%d,%d,%d), fs=%d Hz\n', ...
        info.orders(1), info.orders(2), info.orders(3), d, fs);

    %% 2. 转换为能观标准型 (Observable Canonical Form)
    a_coeffs = A_poly(2:end)';
    b_coeffs = [B_poly, zeros(1, nA - nB)]';

    Ao = zeros(nA);
    if nA > 0
        Ao(:, 1) = -a_coeffs;
        if nA > 1
            Ao(1:nA-1, 2:nA) = eye(nA-1);
        end
    end

    Bo = b_coeffs;
    Co = [1, zeros(1, nA-1)];
    Do = 0;

    sys_obs = ss(Ao, Bo, Co, Do, ts);
    fprintf('能观标准型状态空间构建完成。\n');

    %% 3. 互质性检验
    [is_coprime, gcd_poly, common_roots, cp_analysis] = ...
        check_coprimality_local(A_poly, B_poly);
    fprintf('互质性: %s\n', ternary(is_coprime, '✓ 互质 → 系统能观', '✗ 不互质 → 系统不能观'));

    %% 4. 强可镇定检查
    addpath(fullfile(fileparts(mfilename('fullpath'))));
    [isStrong, unstable_zeros] = checkStrongStabilizability(B_poly, A_poly);
    fprintf('强可镇定: %s\n', ternary(isStrong, '✓ 是', '✗ 否 (存在不稳定零点)'));

    %% 5. 灵敏度函数
    plant_tf = tf(fliplr(B_poly), fliplr(A_poly), ts);
    controller_tf = tf(1, 1, ts);
    w_rad_s = logspace(-1, log10(pi/ts), floor(fs/2));
    plotSensitivity(plant_tf, controller_tf, w_rad_s);
    sgtitle(sprintf('灵敏度分析 (模型: %s, 单位反馈 K=1)', modelSource), 'FontWeight', 'bold');

    %% 6. 零极点分析
    [poles, zers, pole_data, zero_data] = analyze_zeros_poles(a_coeffs, b_coeffs, fs);
    print_analysis_table(poles, zers, pole_data, zero_data, fs);

    %% 7. 能观性矩阵分析
    obsv_matrix = obsv(Ao, Co);
    rank_obsv = rank(obsv_matrix);
    cond_obsv = cond(obsv_matrix);
    fprintf('\n能观性矩阵: 秩=%d (系统阶数=%d), 条件数=%.2e\n', rank_obsv, nA, cond_obsv);
    if cond_obsv > 1e10
        fprintf('警告: 系统几乎不能观 (条件数极大)。\n');
    end

    %% 8. 开环系统图形分析
    figure('Name', sprintf('开环系统分析: %s', modelSource));
    subplot(2,2,1);
    pzmap(sys_obs);
    title('零极点图'); grid on;

    subplot(2,2,2);
    step(sys_obs);
    title('阶跃响应'); grid on;

    subplot(2,2,3);
    impulse(sys_obs);
    title('脉冲响应'); grid on;

    subplot(2,2,4);
    bode(sys_obs);
    title('波特图'); grid on;
    sgtitle(sprintf('开环系统分析: %s', modelSource), 'FontSize', 14, 'FontWeight', 'bold');

    %% 整理输出
    results = struct();
    results.Ao = Ao;  results.Bo = Bo;  results.Co = Co;  results.Do = Do;
    results.sys_obs = sys_obs;
    results.is_coprime = is_coprime;
    results.is_strongly_stabilizable = isStrong;
    results.poles = poles;
    results.zeros = zers;
    results.obsv_rank = rank_obsv;
    results.obsv_cond = cond_obsv;
    results.plant_tf = plant_tf;
    results.coprime_analysis = cp_analysis;

    fprintf('\n=== 分析完成 ===\n');
end

%% 本地工具函数

function s = ternary(cond, t, f)
    if cond, s = t; else, s = f; end
end

function [is_coprime, gcd_poly, common_roots, analysis] = check_coprimality_local(A, B, varargin)
    % 检验两个多项式是否互质（考虑数值误差）
    p = inputParser;
    addParameter(p, 'tolerance', 1e-8, @isnumeric);
    addParameter(p, 'method', 'both', @ischar);
    parse(p, varargin{:});

    tol = p.Results.tolerance;
    method = p.Results.method;
    common_roots = [];
    analysis = struct();
    analysis.tolerance = tol;

    % 方法1: GCD
    if ismember(method, {'gcd', 'both'})
        try
            syms z;
            A_sym = poly2sym(A, z);
            B_sym = poly2sym(B, z);
            gcd_sym = gcd(A_sym, B_sym);
            gcd_poly = sym2poly(gcd_sym);
        catch
            gcd_poly = numerical_gcd_local(A, B, tol);
        end
        gcd_order = length(gcd_poly) - 1;
        is_coprime_gcd = (gcd_order == 0) || (abs(gcd_poly(1)) < tol && gcd_order == 1);
        analysis.gcd_order = gcd_order;
        analysis.gcd_poly = gcd_poly;
    else
        is_coprime_gcd = true;
        gcd_poly = 1;
        analysis.gcd_order = 0;
    end

    % 方法2: 根比较
    if ismember(method, {'roots', 'both'})
        roots_A = roots(A);
        roots_B = roots(B);
        roots_A = roots_A(abs(roots_A) < 1e6);
        roots_B = roots_B(abs(roots_B) < 1e6);
        common_roots = find_common_roots_local(roots_A, roots_B, tol);

        min_dist = inf;
        for i = 1:length(roots_A)
            for j = 1:length(roots_B)
                dist = abs(roots_A(i) - roots_B(j));
                if dist < min_dist, min_dist = dist; end
            end
        end
        is_coprime_roots = isempty(common_roots);
        analysis.min_root_distance = min_dist;
        analysis.common_roots_count = length(common_roots);
    else
        is_coprime_roots = true;
        analysis.min_root_distance = inf;
    end

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

    fprintf('\n=== 互质性分析 ===\n');
    if ~is_coprime
        fprintf('共同根数量: %d\n', length(common_roots));
    end
    fprintf('GCD多项式阶数: %d, 最小根距离: %.2e, 容差: %.2e\n', ...
        analysis.gcd_order, analysis.min_root_distance, analysis.tolerance);
end

function common_roots = find_common_roots_local(roots_A, roots_B, tol)
    common_roots = [];
    used_B = false(size(roots_B));
    for i = 1:length(roots_A)
        for j = 1:length(roots_B)
            if ~used_B(j) && abs(roots_A(i) - roots_B(j)) < tol
                common_roots(end+1) = roots_A(i); %#ok<AGROW>
                used_B(j) = true;
                break;
            end
        end
    end
end

function gcd_poly = numerical_gcd_local(A, B, tol)
    A = A(find(abs(A) > tol, 1):end);
    B = B(find(abs(B) > tol, 1):end);
    if isempty(A), A = 0; end
    if isempty(B), B = 0; end
    if length(A) == 1 || length(B) == 1
        gcd_poly = 1;
        return;
    end
    try
        roots_A = roots(A);
        roots_B = roots(B);
        min_dist = inf;
        for i = 1:length(roots_A)
            for j = 1:length(roots_B)
                dist = abs(roots_A(i) - roots_B(j));
                if dist < min_dist, min_dist = dist; end
            end
        end
        if min_dist < tol
            gcd_poly = [1, 0];
        else
            gcd_poly = 1;
        end
    catch
        gcd_poly = [1, 0];
    end
end
