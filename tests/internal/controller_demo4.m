function result = controller_demo4(signals, test_name, variant, params)
%CONTROLLER_DEMO4 Demo4: epsilon-MOPSO sensitivity-shaping RST (+ Youla-Q RLS).
%
%   用法:
%     result = controller_demo4(signals, test_name, variant, params)
%
%   signals 由 load_cylinder1dm_signals 提供，本函数只做计算。
%   设计流程: eMOPSO 搜索 X/Y 多项式根 (Zames-Francis 灵敏度整形目标 +
%   模值/延迟裕度约束惩罚) → Pareto 存档确定性选解 → bezoutd_reg 反解
%   R/S → 全阶模型仿真。优化与约束复用 demo4_Robust/utils 的
%   eMOPSO_core / RST_objective / RST_constraints / penalty_function /
%   bezoutd_reg (路径由 tests/startup.m 提供)。
%
%   variant:
%     'fixed'    — 仅冻结中央控制器 (离线整形设计)
%     'adaptive' — 在 eMOPSO 中央控制器上叠加 Youla-Q RLS 外环，
%                  自适应律与 demo1 完全相同 (Landau 2017 结构)，
%                  用于隔离比较"中央控制器设计方法"的影响
%
%   params 专用字段 (默认值):
%     .nX, .nY        — X/Y 整形多项式阶数 (4, 2)
%                       阶数约束: deg(P_D)+nX <= deg(A·H_S)+deg(B~·H_R)-1，
%                       cylinder1dm 模型上无陷波时 nX<=9，带二阶陷波时 nX<=11
%     .f_notch        — H_S 内模陷波频率 Hz，[] 表示无陷波 ([])
%     .zeta_notch     — 陷波阻尼，>0 保证控制器严格 Schur (0.01)
%     .H_R, .H_S      — 控制器固定因子 ([1], [1])
%                       注: run_eMOPSO 默认 H_S=[1,-1] 会把控制器极点放在
%                       单位圆上 (积分器)，调节问题不需要，默认去掉。
%     .P_D_zeta, .P_D_omega — 主导极点 Butterworth 参数 (0.8, 0.25*pi)
%     .seed           — eMOPSO 随机种子，固定保证可复现 (42)
%     .n_pop, .k_max, .epsilon, .archive_max — eMOPSO 预算 (40/80/0.01/50)
%     .score_band     — 无陷波时选解评分的目标频带 Hz ([100 500])
%     .nQ, .lambda1, .lambda2, .F_diag — Youla-Q RLS 参数，仅 adaptive
%                       variant 使用，语义与 demo1 相同 (4, 0.98, 1, 1e-3)
%     .actuator_limit — 归一化执行器硬限幅 (5)
%     .control_scale — 控制输出缩放，默认 1

    %% ====================================================================
    % §0. 参数解析
    %% ====================================================================

    if nargin < 2 || isempty(test_name), test_name = 'T1'; end
    if nargin < 3 || isempty(variant), variant = 'fixed'; end
    if ~ismember(variant, {'fixed', 'adaptive'})
        error('Demo4:UnknownVariant', ...
            'Demo4 仅支持 fixed/adaptive variant，收到 ''%s''。', variant);
    end

    if nargin < 4, params = struct(); end
    if ~isfield(params, 'nX'),          params.nX = 4; end
    if ~isfield(params, 'nY'),          params.nY = 2; end
    if ~isfield(params, 'f_notch'),     params.f_notch = []; end
    if ~isfield(params, 'zeta_notch'),  params.zeta_notch = 0.01; end
    if ~isfield(params, 'H_R'),         params.H_R = 1; end
    if ~isfield(params, 'H_S'),         params.H_S = 1; end
    if ~isfield(params, 'P_D_zeta'),    params.P_D_zeta = 0.8; end
    if ~isfield(params, 'P_D_omega'),   params.P_D_omega = 0.25 * pi; end
    if ~isfield(params, 'seed'),        params.seed = 42; end
    if ~isfield(params, 'n_pop'),       params.n_pop = 40; end
    if ~isfield(params, 'k_max'),       params.k_max = 80; end
    if ~isfield(params, 'epsilon'),     params.epsilon = 0.01; end
    if ~isfield(params, 'archive_max'), params.archive_max = 50; end
    if ~isfield(params, 'score_band'),  params.score_band = [100 500]; end
    if ~isfield(params, 'nQ'),          params.nQ = 4; end
    if ~isfield(params, 'lambda1'),     params.lambda1 = 0.98; end
    if ~isfield(params, 'lambda2'),     params.lambda2 = 1; end
    if ~isfield(params, 'F_diag'),      params.F_diag = 1e-3; end
    if ~isfield(params, 'actuator_limit'), params.actuator_limit = 5; end
    if ~isfield(params, 'control_scale'), params.control_scale = 1; end

    %% ====================================================================
    % §1. 模型加载 (同 demo1)
    %% ====================================================================
    if isfield(signals, 'model_data')
        modelData = signals.model_data;
    else
        info = DataManager(signals.model);
        modelData = struct('A', info.model.A, 'B', info.model.B, ...
            'orders', info.orders, 'fs', info.fs);
    end
    A = modelData.A(:).';          % 行向量
    B_full = modelData.B(:).';     % 含 d 步延迟前导零
    d = modelData.orders(4);
    fs = modelData.fs;
    Ts = 1/fs;
    B_coeffs = B_full(d+1:end);    % sys_info.B 约定: 不含延迟前导零

    %% ====================================================================
    % §2. 构建 sys_info (与 run_eMOPSO.m RST_engine 相同的接线)
    %% ====================================================================

    % H_S 内模陷波 (阻尼同 demo1 的 second_order_factor)
    H_S_aug = params.H_S(:).';
    if ~isempty(params.f_notch)
        H_S_aug = conv(H_S_aug, ...
            second_order_factor(params.f_notch, params.zeta_notch, Ts));
    end

    % 主导极点 P_D (二阶 Butterworth 型)
    zeta_P = params.P_D_zeta;
    omega_P = params.P_D_omega;
    PD1 = -2 * exp(-zeta_P * omega_P) * cos(omega_P * sqrt(1 - zeta_P^2));
    PD2 = exp(-2 * zeta_P * omega_P);

    % 不稳定零极点 (无则用虚拟值，同 run_eMOPSO.m:100-114)
    sys_zeros = roots(B_coeffs);
    sys_poles = roots(A);
    beta = sys_zeros(abs(sys_zeros) > 1);
    alpha = sys_poles(abs(sys_poles) > 1);
    if isempty(beta), beta = 1.05; end
    if isempty(alpha), alpha = 0.95; end

    sys_info = struct();
    sys_info.A     = A;
    sys_info.B     = B_coeffs;
    sys_info.d     = d;
    sys_info.P_D   = [1, PD1, PD2];
    sys_info.H_R   = params.H_R(:).';
    sys_info.H_S   = H_S_aug;
    sys_info.nX    = params.nX;
    sys_info.nY    = params.nY;
    sys_info.beta  = beta(:);
    sys_info.alpha = alpha(:);
    sys_info.Ts    = Ts;
    sys_info.nFreq = 500;

    %% ====================================================================
    % §3. 设计缓存: eMOPSO 只按 design_key 跑一次
    %% ====================================================================
    persistent ctrl_cache;
    if isempty(ctrl_cache)
        ctrl_cache = struct('key', '', 'design', [], 'archive_f', [], ...
            'selection', []);
    end

    fNotchKey = 0;
    if ~isempty(params.f_notch), fNotchKey = params.f_notch; end
    design_key = sprintf('%s_nX%d_nY%d_fn%g_z%g_HR%s_HS%s_s%d_p%d_k%d_e%g', ...
        signals.model, params.nX, params.nY, fNotchKey, params.zeta_notch, ...
        mat2str(params.H_R, 4), mat2str(params.H_S, 4), params.seed, ...
        params.n_pop, params.k_max, params.epsilon);

    if strcmp(ctrl_cache.key, design_key)
        design = ctrl_cache.design;
        archive_f = ctrl_cache.archive_f;
        selection = ctrl_cache.selection;
    else
        %% ================================================================
        % §4. eMOPSO 多目标优化 (固定种子，可复现)
        %% ================================================================
        fprintf('--- Demo4: eMOPSO 灵敏度整形 RST 设计 ---\n');
        fprintf('  nX=%d, nY=%d, notch=%s, n_pop=%d, k_max=%d, seed=%d\n', ...
            params.nX, params.nY, mat2str(fNotchKey), ...
            params.n_pop, params.k_max, params.seed);

        n_var = params.nX + params.nY;
        n_obj = length(sys_info.beta) + 2;   % β 匹配 + 延迟积分 + Bezout 残差
        lb = -0.95 * ones(1, n_var);
        ub =  0.95 * ones(1, n_var);

        options = struct('n_pop', params.n_pop, 'k_max', params.k_max, ...
            'epsilon', params.epsilon, 'c1', 1.5, 'c2', 1.7, ...
            'w_max', 0.9, 'w_min', 0.3, ...
            'archive_max', params.archive_max, 'verbose', false);

        obj_func = @(theta) wrapped_objective(theta, sys_info);

        rng_state = rng;
        rng(params.seed);
        % penalty_function 在约束违例时逐次告警，优化中会产生数千条
        % 带回溯的输出；这里临时静音，结束后恢复。
        warn_state = warning('off', 'all');
        [archive, archive_f] = eMOPSO_core(obj_func, ...
            n_var, n_obj, lb, ub, options);
        warning(warn_state);
        rng(rng_state);

        if isempty(archive)
            error('Demo4:EmptyArchive', 'eMOPSO 未返回任何 Pareto 解。');
        end

        %% ================================================================
        % §5. Pareto 存档确定性选解
        %% ================================================================
        scoreFreq = params.f_notch;   % 有陷波: 按陷波频率抑制评分
        [design, selection] = select_pareto_solution(archive, archive_f, ...
            sys_info, scoreFreq, params.score_band, fs);

        fprintf(['  选解 #%d/%d: score=%.2f, Ms=%.2f dB, ' ...
            '设计抑制=%.2f dB, 闭环半径=%.5f, 控制器半径=%.5f\n'], ...
            selection.index, size(archive, 1), selection.score, ...
            selection.Ms_db, selection.suppression_design_db, ...
            design.closed_loop_radius, design.controller_radius);

        ctrl_cache.key = design_key;
        ctrl_cache.design = design;
        ctrl_cache.archive_f = archive_f;
        ctrl_cache.selection = selection;
    end

    %% ====================================================================
    % §6. 选择测试信号并仿真 (全阶模型 + 执行器限幅)
    %% ====================================================================
    test_sig = signals.(test_name);
    d_sig  = test_sig.d(:)';
    y_open = test_sig.y_open(:)';

    switch variant
        case 'fixed'
            [y_closed, u, u_demand, clipped] = sim_fixed_rst( ...
                d_sig, A, B_full, design.R0, design.S0, ...
                params.actuator_limit, params.control_scale);

        case 'adaptive'
            % Youla-Q RLS 外环 (与 demo1 相同的 Landau 2017 结构)。
            % 中央控制器换成 eMOPSO 设计的 R0/S0，P0 用其实际闭环特征
            % 多项式；设计中的内模/固定因子已并入 R0/S0，因此这里
            % HR=HS=1 (纯 YK 增量形式)。
            fprintf('  Q-RLS: nQ=%d, λ₁=%.3f, λ₂=%.1f, F₀=%.0e\n', ...
                params.nQ, params.lambda1, params.lambda2, params.F_diag);
            B_star = B_full(2:end);
            [y_closed, u, Q_final, u_demand, clipped, qNormHistory, ...
                adaptationError] = sim_adaptive_landau(d_sig, A, B_full, ...
                B_star, design.R0, design.S0, 1, 1, design.P_cl, params);
            fprintf('  Q_final: [%s], ||Q||=%.3f\n', ...
                num2str(Q_final, '%.3f '), norm(Q_final));
    end

    %% ====================================================================
    % §7. 标准化指标 + demo4 专用字段
    %% ====================================================================
    meta = struct('demo', 'demo4', 'variant', variant, 'test', test_name);
    result = compute_metrics(y_open, y_closed, u, test_sig, meta);

    result.extra.R0 = design.R0;
    result.extra.S0 = design.S0;
    result.extra.P0 = design.P_cl;
    result.extra.A = A;
    result.extra.B = B_full;
    result.extra.deg_R0 = length(design.R0) - 1;
    result.extra.deg_S0 = length(design.S0) - 1;
    result.extra.design = struct( ...
        'closed_loop_radius', design.closed_loop_radius, ...
        'controller_radius', design.controller_radius, ...
        'Ms_db', selection.Ms_db, ...
        'suppression_design_db', selection.suppression_design_db, ...
        'score', selection.score, ...
        'selected_index', selection.index, ...
        'n_pareto_solutions', selection.n_solutions, ...
        'stable_solutions', selection.n_stable, ...
        'nX', params.nX, 'nY', params.nY, ...
        'n_pop', params.n_pop, 'k_max', params.k_max, ...
        'seed', params.seed, ...
        'f_notch', params.f_notch, ...
        'zeta_notch', params.zeta_notch);
    result.extra.archive_f = archive_f;
    result.extra.syp_omega = design.omega;
    result.extra.syp_mag = design.Syp_mag;
    result.extra.u_demand = u_demand;
    result.extra.u_demand_max = max(abs(u_demand));
    result.extra.saturation_count = sum(clipped);
    result.extra.actuator_limit = params.actuator_limit;
    result.extra.t = test_sig.t;
    result.extra.y_open = y_open;
    result.extra.y_closed = y_closed;
    result.extra.u = u;

    if strcmp(variant, 'adaptive')
        result.extra.Q_final  = Q_final;
        result.extra.nQ       = params.nQ;
        result.extra.lambda1  = params.lambda1;
        result.extra.lambda2  = params.lambda2;
        result.extra.q_norm_history = qNormHistory;
        result.extra.adaptation_error = adaptationError;
    end

    fprintf('  Demo4 %-9s | %s: %.1f dB | t_conv=%.2fs | max|u|=%.3f\n', ...
        variant, test_name, result.supp_db, result.t_conv_s, result.u_max);
end

%% ========================================================================
% 辅助函数
%% ========================================================================

function phi = wrapped_objective(theta, sys_info)
% 包装 RST 目标函数并施加约束惩罚 (基于 run_eMOPSO.m 本地函数)
%   先放大原始目标 1000 倍以匹配 epsilon 量级，
%   再使用乘性惩罚 (penalty_function) 处理约束违例。
%
%   在原实现之上追加控制器内部稳定性惩罚: Zames-Francis 目标和
%   模值/延迟裕度约束都不约束 Bezout 解 S' 的根，实测该对象上
%   整个 Pareto 存档的 S' 根半径都在 1.01-1.02，导致无可实现解。
%   这里对 max|roots(S')| >= 1 的解按超出量乘性放大全部目标，
%   把搜索推向内部稳定区域 (与 demo1 拒绝非 Schur S 的要求一致)。
    theta = theta(:);
    f = RST_objective(theta, sys_info) * 1000;
    [g1, g2] = RST_constraints(theta, sys_info);
    phi = penalty_function(f, g1, g2);

    try
        X = real(poly(theta(1:sys_info.nX)));
        B_tilde = [zeros(1, sys_info.d), sys_info.B];
        [~, S_prime] = bezoutd(sys_info.A, B_tilde, ...
            sys_info.H_S, sys_info.H_R, conv(sys_info.P_D, X));
        sRadius = max(abs(roots(S_prime)));
        if sRadius >= 1
            phi = phi * (2 + 10 * (sRadius - 1));
        end
    catch
        phi = phi * 10;   % Bezout 失败视为严重违例
    end
end

function design = synthesize_rst(theta, sys_info)
% 由决策变量反解 RST 控制器并计算闭环特征多项式与输出灵敏度。
%   与 postprocess_RST 的合成数学一致 (P_c = A·H_S·S' + z^{-d}B·H_R·R'
%   = P_D·X)，但静默、无绘图，且直接返回含 H 因子的完整控制器
%   R0 = H_R·R'，S0 = H_S·S'。
    theta = theta(:);
    X = real(poly(theta(1:sys_info.nX)));
    P_DX = conv(sys_info.P_D, X);
    B_tilde = [zeros(1, sys_info.d), sys_info.B];

    % evalc 吞掉 bezoutd_reg 的逐次诊断打印 (选解循环会调用几十次)
    [~, R_prime, S_prime] = evalc( ...
        'bezoutd_reg(sys_info.A, B_tilde, sys_info.H_S, sys_info.H_R, P_DX)');
    R_prime = trimPolynomial(R_prime(:).', "left");
    S_prime = trimPolynomial(S_prime(:).', "left");

    R0 = conv(sys_info.H_R, R_prime);
    S0 = conv(sys_info.H_S, S_prime);
    P_cl = addPolynomials(conv(sys_info.A, S0), conv(B_tilde, R0));

    omega = linspace(0, pi, 1024);
    Syp_mag = abs(freqz(conv(sys_info.A, S0), P_cl, omega));
    Sup_mag = abs(freqz(conv(sys_info.A, R0), P_cl, omega));   % 输入灵敏度

    design = struct('R0', R0, 'S0', S0, 'P_cl', P_cl, ...
        'closed_loop_radius', root_radius(P_cl), ...
        'controller_radius', root_radius(S0), ...
        'omega', omega, 'Syp_mag', Syp_mag(:).', 'Sup_mag', Sup_mag(:).');
end

function [best, selection] = select_pareto_solution(archive, archive_f, ...
        sys_info, scoreFreq, scoreBand, fs)
% 从 Pareto 存档中确定性选解:
%   1. 过滤闭环与控制器均严格 Schur 的解
%   2. 评分 = 抑制指标 - 2·max(0, Ms_dB - 6) - max(0, Sup峰值_dB - 40)
%      有陷波: 抑制取陷波频率处 -20log10|Syp|
%      无陷波: 抑制取 score_band 频带内平均 -20log10|Syp|
%      输入灵敏度峰值惩罚与 demo1 设计评分一致，约束控制代价
%   3. 全部不稳定时回退 knee-point (归一化目标向量行和最小) 并告警
    nSolutions = size(archive, 1);
    best = [];
    selection = struct('index', 0, 'score', -inf, 'Ms_db', NaN, ...
        'suppression_design_db', NaN, 'n_solutions', nSolutions, ...
        'n_stable', 0);

    for index = 1:nSolutions
        try
            candidate = synthesize_rst(archive(index, :), sys_info);
        catch
            continue;   % Bezout 失败的解直接跳过
        end
        if candidate.closed_loop_radius >= 1 - 1e-7 ...
                || candidate.controller_radius >= 1 - 1e-7
            continue;
        end
        selection.n_stable = selection.n_stable + 1;

        Ms_db = max(20*log10(candidate.Syp_mag + eps));
        Sup_peak_db = max(20*log10(candidate.Sup_mag + eps));
        if ~isempty(scoreFreq)
            [~, idx] = min(abs(candidate.omega - 2*pi*scoreFreq/fs));
            supp_db = -20*log10(max(candidate.Syp_mag(idx), eps));
        else
            bandMask = candidate.omega >= 2*pi*scoreBand(1)/fs ...
                & candidate.omega <= 2*pi*scoreBand(2)/fs;
            supp_db = mean(-20*log10(max(candidate.Syp_mag(bandMask), eps)));
        end
        score = supp_db - 2*max(0, Ms_db - 6) - max(0, Sup_peak_db - 40);

        if score > selection.score
            best = candidate;
            selection.index = index;
            selection.score = score;
            selection.Ms_db = Ms_db;
            selection.suppression_design_db = supp_db;
        end
    end

    if isempty(best)
        warning('Demo4:NoStablePareto', ...
            'Pareto 存档中无稳定解，回退 knee-point 选择。');
        f_min = min(archive_f, [], 1);
        f_range = max(archive_f, [], 1) - f_min;
        f_range(f_range == 0) = 1;
        f_norm = (archive_f - f_min) ./ f_range;
        [~, kneeIndex] = min(sum(f_norm, 2));
        best = synthesize_rst(archive(kneeIndex, :), sys_info);
        selection.index = kneeIndex;
        selection.Ms_db = max(20*log10(best.Syp_mag + eps));
        selection.score = -inf;
    end
end

function factor = second_order_factor(frequency, zeta, Ts)
% 阻尼二阶内模因子 (同 demo1 设计)
wn = 2*pi*frequency;
radius = exp(-zeta*wn*Ts);
wd = wn*sqrt(max(0, 1-zeta^2));
factor = [1, -2*radius*cos(wd*Ts), radius^2];
end

function [y, u, u_demand, clipped] = sim_fixed_rst( ...
        d, A, B, R, S, actuatorLimit, controlScale)
% 固定 RST 控制器仿真 (同 demo1)
%   被控对象: y(k) = d(k) + (B/A)*u(k)
%   控制律:   S(z⁻¹)·u(k) = -R(z⁻¹)·y(k)
    N = length(d);
    y = zeros(1, N);
    u = zeros(1, N);
    u_demand = zeros(1, N);
    clipped = false(1, N);

    nA = length(A);  nB = length(B);
    nR = length(R);  nS = length(S);
    buf_len = max([nA, nB, nR, nS]) + 1;

    u_buf    = zeros(1, buf_len);
    anti_buf = zeros(1, buf_len);
    y_buf    = zeros(1, buf_len);

    for k = 1:N
        anti_new = FIR(B(2:end), u_buf) - FIR(A(2:end), anti_buf);
        y_new = d(k) + anti_new;

        u_request = -(FIR(R, [y_new, y_buf(1:end-1)]) + ...
                  FIR(S(2:end), u_buf)) / S(1);
        u_request = controlScale * u_request;
        [u_new, clipped(k)] = apply_actuator_limit(u_request, actuatorLimit);

        anti_buf = [anti_new, anti_buf(1:end-1)];
        u_buf    = [u_new,    u_buf(1:end-1)];
        y_buf    = [y_new,    y_buf(1:end-1)];

        y(k) = y_new;
        u(k) = u_new;
        u_demand(k) = u_request;
    end
end

function [u_actual, clipped] = apply_actuator_limit(u_demand, limit)
% 归一化执行器硬限幅
    u_actual = min(max(u_demand, -limit), limit);
    clipped = abs(u_actual - u_demand) > 1e-10;
end

function [y, u, Q_final, u_demand, clipped, qNormHistory, ...
        adaptationError] = sim_adaptive_landau(...
    disturbance, A, B, B_star, Ro, So, HR, HS, Po, params)
% R₀/S₀ + Youla-Q RLS 自适应 (与 demo1 完全相同的 Landau 2017 结构)
%
%   关键: 抗噪路径用 B(2:end)，扰动观测用 B_star = B(2:end) → 延迟一致
%         P₀ 为中央控制器的实际闭环特征多项式
%         RLS 协方差: F_new = λ₁·inv(F) + λ₂·Φ'Φ; F = inv(F_new)

    nQ = params.nQ;
    lambda1 = params.lambda1;
    lambda2 = params.lambda2;
    F_diag = params.F_diag;

    N = length(disturbance);

    % 多项式长度 (B(2:end) 用于抗噪和观测)
    B_landau = B(2:end);
    nA = length(A);  nB = length(B_landau);
    nRo = length(Ro);  nSo = length(So);
    nPo = length(Po);

    % 滚动缓冲区长度
    buf_len = max([nA, nB, nRo, nSo, nPo, nQ+3]) + 5;

    % 信号缓冲区 (Landau 滚动风格)
    U  = zeros(1, buf_len);
    Y  = zeros(1, buf_len);
    Yf = zeros(1, buf_len);
    W  = zeros(1, buf_len);
    W1 = zeros(1, buf_len);
    W2 = zeros(1, buf_len);

    % RLS 状态
    Q_vec = zeros(1, nQ + 1);
    Phi   = zeros(1, nQ + 1);
    F_mat = F_diag * eye(nQ + 1);

    % 预计算多项式卷积
    HR_HS = conv(HR, HS);
    HR_HS_Bstar = conv(HR_HS, B_star);

    % 输出记录
    y = zeros(1, N);
    u = zeros(1, N);
    u_demand = zeros(1, N);
    clipped = false(1, N);
    qNormHistory = zeros(1, N);
    adaptationError = zeros(1, N);

    for k = 1:N
        % ---- 1. 被控对象: y = d + (B/A)*u ----
        antinoise = FIR(B_landau, U) - FIR(A(2:end), Yf);
        Yf = [antinoise, Yf(1:end-1)];

        y_new = antinoise + disturbance(k);
        Y = [y_new, Y(1:end-1)];

        % ---- 2. 扰动观测: w = A*y - B_star*u ----
        w_new = FIR(A, Y) - FIR(B_star, U);
        W = [w_new, W(1:end-1)];

        % ---- 3. 滤波扰动: w1 = (S₀/P₀)*w ----
        w1_new = (FIR(So, W) - FIR(Po(2:end), W1)) / Po(1);
        W1 = [w1_new, W1(1:end-1)];

        % ---- 4. RLS 先验/后验误差 ----
        e0 = w1_new - FIR(Q_vec, Phi);
        e_post = e0 / (1 + Phi * F_mat * Phi');

        % ---- 5. RLS 更新 ----
        Q_vec = Q_vec + (F_mat * Phi' * e_post)';

        F_new = lambda1 * inv(F_mat) + lambda2 * (Phi' * Phi);
        F_new = F_new + 1e-8 * eye(nQ + 1);
        F_mat = inv(F_new);

        if any(isnan(F_mat(:))) || any(isinf(F_mat(:))) || rcond(F_mat) < 1e-12
            Q_vec = zeros(1, nQ + 1);
            F_mat = F_diag * eye(nQ + 1);
        end

        % ---- 6. 控制律 ----
        HR_HS_Q = conv(HR_HS, Q_vec);
        u_request = -FIR(Ro, Y) - FIR(HR_HS_Q, W) - FIR(So(2:end), U);
        u_request = u_request / So(1);
        u_request = params.control_scale * u_request;
        [u_new, clipped(k)] = apply_actuator_limit( ...
            u_request, params.actuator_limit);
        U = [u_new, U(1:end-1)];

        % ---- 7. 回归向量: w2 = (H_R*H_S*B_star/P₀)*w ----
        w2_new = (FIR(HR_HS_Bstar, W) - FIR(Po(2:end), W2)) / Po(1);
        W2 = [w2_new, W2(1:end-1)];
        Phi = [w2_new, Phi(1:end-1)];

        % 记录
        y(k) = y_new;
        u(k) = u_new;
        u_demand(k) = u_request;
        qNormHistory(k) = norm(Q_vec);
        adaptationError(k) = e_post;
    end

    Q_final = Q_vec;
end

function y = FIR(b, x)
% FIR 点积滤波: y = Σ b(i)·x(i)
    L = min(length(b), length(x));
    y = b(1:L) * x(1:L)';
end

function radius = root_radius(polynomial)
if numel(polynomial) <= 1
    radius = 0;
else
    radius = max(abs(roots(polynomial)));
end
end
