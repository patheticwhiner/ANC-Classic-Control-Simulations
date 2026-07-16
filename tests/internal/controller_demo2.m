function result = controller_demo2(signals, test_name, variant, params)
%CONTROLLER_DEMO2 Demo2: augmented LQG plus feedback adaptation.
%
%   增广模式 (params.augmented=true, 默认):
%     构建扰动状态空间模型 → 增广 ARMAX 状态 → LQG 在增广系统上设计
%     → Kalman 观测器可分离扰动与对象状态 (Qian.md / 钱梵梵 2022 框架)
%
%   裸模式 (params.augmented=false):
%     仅在 ARMAX 状态空间上设计 LQG (无扰动模型)
%
%   params:
%     .augmented — 是否使用增广扰动模型 (默认 true)
%     .Q_lqr, .R_lqr — LQR 权重；空参数时使用当前场景冻结值
%     .Qn_plant  — 对象状态过程噪声 (默认 1e-4)
%     .Qn_dist   — 扰动状态过程噪声 (默认 1)
%     .Rn        — 测量噪声 (默认 0.01)
%     .nQ, .lambda, .F_init, .bp_band — filtered-x Q-RLS 参数
%     .adaptation_method — 'rls' (默认) 或 'lms'
%     .lms_epsilon, .lms_leakage — normalized-LMS 稳定化参数
%     .adaptive_warmup_seconds — 自适应启动前固定 LQG 预热时间
%     .adaptation_gain — RLS 参数更新松弛系数
%     .theta_norm_limit — Q 参数向量范数上限
%     .actuator_limit — 归一化执行器硬限幅，默认 5
%     .control_scale, .ramp_seconds — 控制增益缩放与启动渐入

    if nargin < 2 || isempty(test_name), test_name = 'T1'; end
    if nargin < 3 || isempty(variant), variant = 'fixed'; end
    if nargin < 4, params = struct(); end
    if isempty(fieldnames(params))
        params = frozen_demo2_defaults(test_name, variant);
    end
    rnProvided = isfield(params, 'Rn');
    if ~isfield(params, 'augmented'), params.augmented = true; end
    if ~isfield(params, 'Q_lqr'),     params.Q_lqr     = 1; end
    if ~isfield(params, 'R_lqr'),     params.R_lqr     = 1; end
    if ~isfield(params, 'Qn_plant'),  params.Qn_plant  = 1e-4; end
    if ~isfield(params, 'Qn_dist'),   params.Qn_dist   = 1; end
    if ~isfield(params, 'Rn'),        params.Rn        = 0.01; end
    if ~isfield(params, 'nQ'),        params.nQ        = 6; end
    if ~isfield(params, 'lambda'),    params.lambda    = 0.98; end
    if ~isfield(params, 'F_init'),    params.F_init    = 1; end
    if ~isfield(params, 'bp_band'),   params.bp_band   = [200, 600]; end
    if ~isfield(params, 'adaptive_warmup_seconds'), params.adaptive_warmup_seconds = 0.5; end
    if ~isfield(params, 'adaptation_gain'), params.adaptation_gain = 1; end
    if ~isfield(params, 'adaptation_method'), params.adaptation_method = 'rls'; end
    if ~isfield(params, 'lms_epsilon'), params.lms_epsilon = 1e-6; end
    if ~isfield(params, 'lms_leakage'), params.lms_leakage = 0; end
    if ~isfield(params, 'theta_norm_limit'), params.theta_norm_limit = 10; end
    if ~isfield(params, 'rls_regularization'), params.rls_regularization = 1e-9; end
    if ~isfield(params, 'actuator_limit'), params.actuator_limit = 5; end
    if ~isfield(params, 'control_scale'), params.control_scale = 1; end
    if ~isfield(params, 'ramp_seconds'), params.ramp_seconds = 0.1; end
    params.adaptation_method = lower(char(params.adaptation_method));
    if ~ismember(params.adaptation_method, {'rls', 'lms'})
        error('Demo2:UnknownAdaptationMethod', ...
            'Unknown adaptation method: %s', params.adaptation_method);
    end

    %% ====================================================================
    % §1. 模型加载
    %% ====================================================================
    if isfield(signals, 'model_data')
        modelData = signals.model_data;
    else
        info = DataManager(signals.model);
        modelData = struct('A', info.model.A, 'B', info.model.B, ...
            'orders', info.orders, 'fs', info.fs);
    end
    fs = modelData.fs;
    Ts = 1/fs;

    [Af, Bf, Cf] = deterministic_state_space(modelData, Ts);
    n_plant = size(Af, 1);

    %% ====================================================================
    % §2. 构建增广扰动模型
    %% ====================================================================
    test_sig = signals.(test_name);
    if ~rnProvided && strcmp(test_sig.type, 'bandlimited_noise')
        params.Rn = 0.1;
    end

    if params.augmented
        [A_dist, C_dist, n_dist, G_dist_init] = build_disturbance_model(test_sig, fs, Ts);

        % 增广状态空间: x̃ = [x_plant; x_dist]
        A_aug = blkdiag(Af, A_dist);
        B_aug = [Bf; zeros(n_dist, 1)];
        C_aug = [Cf, C_dist];
        % 过程噪声: 对象噪声 + 扰动激励 (窄带噪声注入扰动状态)
        G_aug = blkdiag(eye(n_plant), G_dist_init);

        A_design = A_aug;
        B_design = B_aug;
        C_design = C_aug;
        G_design = G_aug;
        n_total   = n_plant + n_dist;

        Qn = blkdiag(params.Qn_plant * eye(n_plant), ...
                     params.Qn_dist * eye(size(G_dist_init,2)));
        Rn = params.Rn * eye(size(C_design, 1));
    else
        % 裸 LQG: 对确定性通道重新设计状态估计器。
        A_design = Af;  B_design = Bf;  C_design = Cf;
        G_design = eye(n_plant);
        n_total   = n_plant;  n_dist = 0;

        Qn = params.Qn_plant * eye(n_plant);
        Rn = params.Rn * eye(size(Cf, 1));
    end

    %% ====================================================================
    % §3. LQG 设计
    %% ====================================================================
    persistent lqg_cache;
    if isempty(lqg_cache)
        lqg_cache = struct('key', '', 'K', [], 'L', [], 'A_cl', [], ...
            'A_design', [], 'B_design', [], 'C_design', [], ...
            'n_plant', [], 'n_dist', []);
    end

    cache_key = sprintf(['%s_%s_aug%d_Q%.12g_R%.12g_Qp%.12g_' ...
        'Qd%.12g_Rn%.12g_scale%.12g'], ...
        signals.model, test_name, params.augmented, ...
        params.Q_lqr, params.R_lqr, params.Qn_plant, params.Qn_dist, ...
        params.Rn, params.control_scale);

    if ~strcmp(lqg_cache.key, cache_key)
        fprintf('--- Demo2: LQG (n_plant=%d, n_dist=%d) ---\n', n_plant, n_dist);

        % LQR
        Q_lqr = params.Q_lqr * (C_design' * C_design);
        R_lqr = params.R_lqr;
        [K, ~, ~] = dlqr(A_design, B_design, Q_lqr, R_lqr);
        K = params.control_scale * K;

        [L, ~, ~] = dlqe(A_design, G_design, C_design, Qn, Rn);

        A_cl = A_design - B_design * K - L * C_design;
        if max(abs(eig(A_cl))) >= 1
            warning('Demo2: A_cl 不稳定 (max|eig|=%.4f)', max(abs(eig(A_cl))));
        end

        lqg_cache.key       = cache_key;
        lqg_cache.K         = K;
        lqg_cache.L         = L;
        lqg_cache.A_cl      = A_cl;
        lqg_cache.A_design  = A_design;
        lqg_cache.B_design  = B_design;
        lqg_cache.C_design  = C_design;
        lqg_cache.n_plant   = n_plant;
        lqg_cache.n_dist    = n_dist;

        fprintf('  n=%d, max|eig(A_cl)|=%.4f\n', n_total, max(abs(eig(A_cl))));
    else
        K         = lqg_cache.K;
        L         = lqg_cache.L;
        A_cl      = lqg_cache.A_cl;
        A_design  = lqg_cache.A_design;
        B_design  = lqg_cache.B_design;
        C_design  = lqg_cache.C_design;
        n_plant   = lqg_cache.n_plant;
        n_dist    = lqg_cache.n_dist;
    end

    %% ====================================================================
    % §4. 仿真
    %% ====================================================================
    d_sig  = test_sig.d(:)';
    y_open = test_sig.y_open(:)';
    [A_t12, B_t12, C_t12] = build_lqg_t12( ...
        Af, Bf, Cf, A_design, B_design, C_design, K, L);
    t12_radius = max(abs(eig(A_t12)));

    % 使用完整对象动力学 (仅植物部分 + 扰动叠加)
    switch variant
        case 'fixed'
            [y_closed, u, u_demand, clipped] = sim_lqg_augmented(...
                Af, Bf, Cf, n_plant, n_dist, ...
                A_design, B_design, C_design, K, L, d_sig, ...
                params.actuator_limit, round(params.ramp_seconds * fs));

        case 'adaptive'
            [y_closed, u, adaptive, u_demand, clipped] = ...
                sim_lqg_augmented_adaptive(...
                Af, Bf, Cf, n_plant, n_dist, ...
                A_design, B_design, C_design, K, L, ...
                A_t12, B_t12, C_t12, d_sig, fs, params);

        otherwise
            error('Demo2: 未知 variant ''%s''', variant);
    end

    %% ====================================================================
    % §5. 指标
    %% ====================================================================
    meta = struct('demo', 'demo2', 'variant', variant, 'test', test_name);
    result = compute_metrics(y_open, y_closed, u, test_sig, meta);

    result.extra.augmented  = params.augmented;
    result.extra.n_dist     = n_dist;
    result.extra.Q_lqr      = params.Q_lqr;
    result.extra.R_lqr      = params.R_lqr;
    result.extra.control_scale = params.control_scale;
    result.extra.ramp_seconds = params.ramp_seconds;
    result.extra.adaptation_method = params.adaptation_method;
    result.extra.max_eig_Acl = max(abs(eig(A_cl)));
    result.extra.max_eig_T12 = t12_radius;
    result.extra.K = K;
    result.extra.L = L;
    result.extra.A_cl = A_cl;
    % --- dSPACE/Simulink 移植所需的设计矩阵（纯输出附加，不影响仿真数值） ---
    result.extra.Af = Af;
    result.extra.Bf = Bf;
    result.extra.Cf = Cf;
    result.extra.A_design = A_design;
    result.extra.B_design = B_design;
    result.extra.C_design = C_design;
    result.extra.n_plant = n_plant;
    result.extra.A_t12 = A_t12;
    result.extra.B_t12 = B_t12;
    result.extra.C_t12 = C_t12;
    result.extra.ramp_samples = round(params.ramp_seconds * fs);
    result.extra.u_demand = u_demand;
    result.extra.u_demand_max = max(abs(u_demand));
    result.extra.saturation_count = sum(clipped);
    result.extra.actuator_limit = params.actuator_limit;
    result.extra.t = test_sig.t;
    result.extra.y_open = y_open;
    result.extra.y_closed = y_closed;
    result.extra.u = u;

    if strcmp(variant, 'adaptive')
        result.extra.nQ         = params.nQ;
        result.extra.lambda     = params.lambda;
        result.extra.F_init     = params.F_init;
        result.extra.bp_band    = params.bp_band;
        result.extra.adaptation_gain = params.adaptation_gain;
        % --- 移植所需：与 sim_lqg_augmented_adaptive 相同方式重算带通系数 ---
        Wn_extra = params.bp_band / (fs/2);
        Wn_extra = max(1e-3, min(Wn_extra, 0.999));
        [result.extra.bp_b, result.extra.bp_a] = butter(4, Wn_extra, 'bandpass');
        result.extra.warmup_samples = round(params.adaptive_warmup_seconds * fs);
        result.extra.theta_norm_limit = params.theta_norm_limit;
        result.extra.lms_epsilon = params.lms_epsilon;
        result.extra.lms_leakage = params.lms_leakage;
        result.extra.rls_regularization = params.rls_regularization;
        result.extra.norm_theta = norm(adaptive.theta_final);
        result.extra.theta_final = adaptive.theta_final;
        result.extra.theta_norm_history = adaptive.theta_norm_history;
        result.extra.adaptation_error = adaptive.adaptation_error;
        result.extra.filtered_residual = adaptive.filtered_residual;
        result.extra.q_output = adaptive.q_output;
        result.extra.filtered_regressor_norm = adaptive.filtered_regressor_norm;
        result.extra.adaptation_update_count = adaptive.update_count;
    end

    fprintf('  Demo2 %-9s | %s: %.1f dB | t_conv=%.2fs | max|u|=%.3f\n', ...
        variant, test_name, result.supp_db, result.t_conv_s, result.u_max);
end

function params = frozen_demo2_defaults(testName, variant)
%FROZEN_DEMO2_DEFAULTS Match the parameters used by the stage report.
params = struct();
if strcmp(testName, 'T1') && strcmp(variant, 'fixed')
    params = struct( ...
        'Q_lqr', 1, 'R_lqr', 1e-4, ...
        'Qn_plant', 1e-4, 'Qn_dist', 0.1, 'Rn', 0.1, ...
        'control_scale', 1, 'ramp_seconds', 0.1);
elseif strcmp(testName, 'T1') && strcmp(variant, 'adaptive')
    params = struct( ...
        'Q_lqr', 1, 'R_lqr', 1e-4, ...
        'Qn_plant', 1e-4, 'Qn_dist', 0.1, 'Rn', 0.1, ...
        'control_scale', 1, 'ramp_seconds', 0.1, ...
        'nQ', 4, 'lambda', 0.995, 'F_init', 0.1, 'adaptation_gain', 0.1, ...
        'adaptation_method', 'rls', 'theta_norm_limit', 8);
elseif strcmp(testName, 'T2') && strcmp(variant, 'adaptive')
    params = struct( ...
        'Q_lqr', 1, 'R_lqr', 1e-2, ...
        'Qn_plant', 1e-4, 'Qn_dist', 0.1, 'Rn', 0.1, ...
        'control_scale', 1, 'ramp_seconds', 0.1, ...
        'bp_band', [280 440], 'adaptive_warmup_seconds', 0.5, ...
        'theta_norm_limit', 8, 'nQ', 32, ...
        'lambda', 0.98, 'F_init', 0.1, 'adaptation_gain', 0.15, ...
        'adaptation_method', 'rls');
elseif strcmp(testName, 'T3') && strcmp(variant, 'fixed')
    params = struct( ...
        'Q_lqr', 1, 'R_lqr', 1e-2, ...
        'Qn_plant', 1e-4, 'Qn_dist', 0.1, 'Rn', 0.1, ...
        'control_scale', 1, 'ramp_seconds', 0.1);
end
end

function [A_t12, B_t12, C_t12] = build_lqg_t12( ...
        Af_plant, Bf_plant, Cf_plant, ...
        A_design, B_design, C_design, K, L)
% Exact closed-loop map from additive Q output u_q to measured output y.
    n_plant = size(Af_plant, 1);
    n_design = size(A_design, 1);

    A_t12 = [Af_plant, -Bf_plant * K; ...
             L * Cf_plant, A_design - L * C_design - B_design * K];
    B_t12 = [Bf_plant; B_design];
    C_t12 = [Cf_plant, zeros(size(Cf_plant, 1), n_design)];

    if size(A_t12, 1) ~= n_plant + n_design
        error('Demo2:T12DimensionMismatch', 'T12 联合闭环状态维度不一致。');
    end
end

function [A, B, C] = deterministic_state_space(model, Ts)
    persistent cachedModel cachedTs cachedA cachedB cachedC;
    if ~isempty(cachedModel) && isequal(model, cachedModel) && Ts == cachedTs
        A = cachedA;
        B = cachedB;
        C = cachedC;
        return;
    end

    numerator = model.B;
    denominator = model.A;

    plant = ss(tf(numerator, denominator, Ts, 'Variable', 'z^-1'));
    [A, B, C] = ssdata(plant);
    cachedModel = model;
    cachedTs = Ts;
    cachedA = A;
    cachedB = B;
    cachedC = C;
end

%% ====================================================================
% 扰动模型构建
%% ====================================================================
function [A_dist, C_dist, n_dist, G_dist] = build_disturbance_model(test_sig, fs, Ts)
    switch test_sig.type
        case 'fixed_sine'
            f0 = test_sig.f_hz;
            omega = 2*pi*f0*Ts;
            r_d = 0.999;
            A_dist = r_d * [cos(omega), -sin(omega); sin(omega), cos(omega)];
            C_dist = [1, 0];
            G_dist = eye(2);  % 过程噪声驱动两个扰动状态
            n_dist = 2;

        case 'linear_chirp'
            f_center = mean(test_sig.f_range);
            omega = 2*pi*f_center*Ts;
            zeta_d = 0.05;
            r = exp(-zeta_d*omega);
            A_dist = r * [cos(omega), -sin(omega); sin(omega), cos(omega)];
            C_dist = [1, 0];
            G_dist = eye(2);
            n_dist = 2;

        case 'bandlimited_noise'
            f_range = test_sig.f_range;
            Wn = f_range / (fs/2);
            Wn = max(1e-3, min(Wn, 0.999));
            [b_bp, a_bp] = butter(2, Wn, 'bandpass');
            [A_dist, ~, C_dist, ~] = tf2ss(b_bp, a_bp);
            n_dist = size(A_dist, 1);
            G_dist = eye(n_dist);

        otherwise
            omega = 2*pi*334*Ts;
            r_d = 0.999;
            A_dist = r_d * [cos(omega), -sin(omega); sin(omega), cos(omega)];
            C_dist = [1, 0];
            G_dist = eye(2);
            n_dist = 2;
    end
end

%% ====================================================================
% 仿真引擎
%% ====================================================================
function [y, u, u_demand, clipped] = sim_lqg_augmented(...
        Af_plant, Bf_plant, Cf_plant, n_plant, n_dist, ...
        A_design, B_design, C_design, K, L, d, actuatorLimit, rampSamples)
% 仿真: 真实扰动 d(k) 直接加到输出 (观测器内部用增广模型估计)
    N = length(d);
    x_plant = zeros(n_plant, 1);
    x_hat   = zeros(n_plant + n_dist, 1);
    y = zeros(1, N);
    u = zeros(1, N);
    u_demand = zeros(1, N);
    clipped = false(1, N);

    for k = 1:N
        % 测量: y = plant_output + REAL disturbance
        y_k = Cf_plant * x_plant + d(k);

        % 控制律 (作用于增广状态估计)
        ramp = min(1, max(0, (k - 1) / max(1, rampSamples)));
        u_request = -ramp * K * x_hat;
        [u_k, clipped(k)] = apply_actuator_limit(u_request, actuatorLimit);

        % Kalman 更新 (增广观测器学习分离扰动与对象状态)
        x_hat = (A_design - L * C_design) * x_hat ...
            + B_design * u_k + L * y_k;

        % 对象状态更新
        x_plant = Af_plant * x_plant + Bf_plant * u_k;

        y(k) = y_k;
        u(k) = u_k;
        u_demand(k) = u_request;
    end
end

function [y, u, adaptive, u_demand, clipped] = sim_lqg_augmented_adaptive(...
        Af_plant, Bf_plant, Cf_plant, n_plant, n_dist, ...
        A_design, B_design, C_design, K, L, ...
        A_t12, B_t12, C_t12, d, fs, params)
% 增广 LQG + filtered-x Youla-Q RLS。
% H(z)r 生成 Q 输出基；T12(z)H(z)r 仅用于参数梯度。

    n_total = n_plant + n_dist;
    N = length(d);
    x_plant = zeros(n_plant, 1);
    x_hat   = zeros(n_total, 1);
    y = zeros(1, N);
    u = zeros(1, N);
    u_demand = zeros(1, N);
    clipped = false(1, N);
    q_output = zeros(1, N);
    filtered_residual = zeros(1, N);
    adaptation_error = zeros(1, N);
    theta_norm_history = zeros(1, N);
    filtered_regressor_norm = zeros(1, N);

    bp = params.bp_band;
    Wn = bp / (fs/2);
    Wn = max(1e-3, min(Wn, 0.999));
    [b_h, a_h] = butter(4, Wn, 'bandpass');
    filter_state = zeros(max(length(a_h), length(b_h)) - 1, 1);

    nQ = params.nQ;
    lambda = params.lambda;
    theta = zeros(nQ, 1);
    P_mat = params.F_init * eye(nQ);
    q_basis = zeros(nQ, 1);
    t12_states = zeros(size(A_t12, 1), nQ);
    warmup_samples = round(params.adaptive_warmup_seconds * fs);
    ramp_samples = round(params.ramp_seconds * fs);
    update_count = 0;

    for k = 1:N
        y_k = Cf_plant * x_plant + d(k);

        r_k = y_k - C_design * x_hat;
        [r_filt, filter_state] = filter(b_h, a_h, r_k, filter_state);
        q_basis = [r_filt; q_basis(1:end-1)];

        filtered_x = (C_t12 * t12_states).';
        y_Q_k = theta' * q_basis;
        ramp = min(1, max(0, (k - 1) / max(1, ramp_samples)));
        u_request = ramp * (-K * x_hat + y_Q_k);
        [u_k, clipped(k)] = apply_actuator_limit( ...
            u_request, params.actuator_limit);

        if k >= warmup_samples && ~clipped(k)
            if strcmp(params.adaptation_method, 'lms')
                phi_energy = filtered_x' * filtered_x;
                lms_step = params.adaptation_gain / ...
                    max(params.lms_epsilon + phi_energy, 1e-12);
                theta = (1 - params.lms_leakage) * theta ...
                    - lms_step * filtered_x * y_k;
            else
                denom = lambda + filtered_x' * P_mat * filtered_x;
                gain = P_mat * filtered_x / max(denom, 1e-10);
                theta = theta - params.adaptation_gain * gain * y_k;
                P_mat = (P_mat - gain * filtered_x' * P_mat) / lambda;
            end
            update_count = update_count + 1;

            theta_norm = norm(theta);
            if theta_norm > params.theta_norm_limit
                theta = theta * (params.theta_norm_limit / theta_norm);
            end
        end
        if strcmp(params.adaptation_method, 'rls')
            P_mat = 0.5 * (P_mat + P_mat') ...
                + params.rls_regularization * eye(nQ);
        end
        if any(isnan(theta)) || any(isinf(theta))
            theta = zeros(nQ, 1);  P_mat = params.F_init * eye(nQ);
        end

        x_hat = (A_design - L * C_design) * x_hat ...
            + B_design * u_k + L * y_k;
        x_plant = Af_plant * x_plant + Bf_plant * u_k;
        t12_states = A_t12 * t12_states + B_t12 * q_basis.';

        y(k) = y_k;
        u(k) = u_k;
        u_demand(k) = u_request;
        q_output(k) = y_Q_k;
        filtered_residual(k) = r_filt;
        adaptation_error(k) = y_k;
        theta_norm_history(k) = norm(theta);
        filtered_regressor_norm(k) = norm(filtered_x);
    end

    adaptive = struct( ...
        'theta_final', theta, ...
        'theta_norm_history', theta_norm_history, ...
        'adaptation_error', adaptation_error, ...
        'filtered_residual', filtered_residual, ...
        'q_output', q_output, ...
        'filtered_regressor_norm', filtered_regressor_norm, ...
        'update_count', update_count);
end

function [u_actual, clipped] = apply_actuator_limit(u_demand, limit)
% 归一化执行器硬限幅
    u_actual = min(max(u_demand, -limit), limit);
    clipped = abs(u_actual - u_demand) > 1e-10;
end
