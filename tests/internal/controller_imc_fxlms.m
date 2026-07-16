function result = controller_imc_fxlms(signals, test_name, variant, params)
%CONTROLLER_IMC_FXLMS  Feedback IMC-FxLMS baseline controller.
%
%   result = controller_imc_fxlms(signals, test_name, variant, params)
%
%   纯反馈 ANC 基线: 通过 IMC 结构从误差信号重构内部参考, 使用
%   filtered-x LMS 自适应 FIR 控制滤波器。无参考麦克风, 信息结构
%   与 Marino-Tomei 频估器一致。
%
%   variant:
%     'fixed'    — 冻结 FIR 系数 (从零初始化, 不更新)
%     'adaptive' — online filtered-x LMS 更新
%
%   params 专用字段:
%     .N_fir          — 控制 FIR 阶数，默认 64
%     .mu             — LMS 步长，默认 0.005
%     .delta          — NLMS 正则化项，默认 1e-4
%     .actuator_limit — 归一化执行器硬限幅，默认 5
%     .ramp_seconds   — 控制输出启动渐入时间，默认 0.1
%     .theta_init     — FIR 初始系数向量，默认 zeros(N_fir,1)
%     .norm_limit     — 系数范数上限，默认 Inf

    %% ====================================================================
    % §0. 参数解析与模型加载
    %% ====================================================================

    if nargin < 2 || isempty(test_name), test_name = 'T1'; end
    if nargin < 3 || isempty(variant),   variant   = 'adaptive'; end
    if nargin < 4,                       params    = struct(); end

    if ~isfield(params, 'N_fir'),          params.N_fir          = 64; end
    if ~isfield(params, 'mu'),             params.mu             = 0.005; end
    if ~isfield(params, 'delta'),          params.delta          = 1e-4; end
    if ~isfield(params, 'actuator_limit'), params.actuator_limit = 5; end
    if ~isfield(params, 'ramp_seconds'),   params.ramp_seconds   = 0.1; end
    if ~isfield(params, 'norm_limit'),     params.norm_limit     = Inf; end
    if ~isfield(params, 'theta_init'),     params.theta_init     = []; end

    % --- 加载次级路径模型 ---
    if isfield(signals, 'model_data')
        modelData = signals.model_data;
    else
        info = DataManager(signals.model);
        modelData = struct('A', info.model.A, 'B', info.model.B, ...
            'orders', info.orders, 'fs', info.fs);
    end
    A_sec = modelData.A(:).';
    B_sec = modelData.B(:).';
    d_delay = modelData.orders(4);
    fs = modelData.fs;
    Ts = 1 / fs;

    B_poly = B_sec;
    nA = length(A_sec) - 1;
    nB = length(B_poly) - 1;
    N_fir = params.N_fir;

    % --- FIR 初始系数 ---
    if isempty(params.theta_init)
        theta_init = zeros(N_fir, 1);
    else
        theta_init = params.theta_init(:);
        if length(theta_init) > N_fir
            theta_init = theta_init(1:N_fir);
        elseif length(theta_init) < N_fir
            theta_init = [theta_init; zeros(N_fir - length(theta_init), 1)];
        end
    end

    %% ====================================================================
    % §1. 加载测试信号
    %% ====================================================================

    test_sig = signals.(test_name);
    d_sig = test_sig.d(:)';
    Nsim = length(d_sig);
    Tsim = test_sig.Tsim;
    y_open = test_sig.y_open(:)';

    %% ====================================================================
    % §2. 自适应增益
    %% ====================================================================

    switch variant
        case 'fixed'
            mu_use = 0;
        case 'adaptive'
            mu_use = params.mu;
        otherwise
            error('IMC-FxLMS: 未知 variant ''%s''', variant);
    end

    %% ====================================================================
    % §3. IMC-FxLMS 仿真
    %% ====================================================================

    [y_closed, u, u_demand, clipped, theta_final, coeff_time, coeff_norm] = ...
        sim_imc_fxlms_feedback(d_sig, Tsim, fs, N_fir, ...
        mu_use, params.delta, theta_init, ...
        A_sec, B_poly, nA, nB, ...
        params.actuator_limit, round(params.ramp_seconds * fs), ...
        params.norm_limit);

    %% ====================================================================
    % §4. 计算标准化指标
    %% ====================================================================

    meta = struct('demo', 'imc_fxlms', 'variant', variant, 'test', test_name);
    result = compute_metrics(y_open, y_closed, u, test_sig, meta);

    %% ====================================================================
    % §5. 追加专用字段
    %% ====================================================================

    result.extra.N_fir           = N_fir;
    result.extra.mu              = mu_use;
    result.extra.delta           = params.delta;
    result.extra.actuator_limit  = params.actuator_limit;
    result.extra.ramp_seconds    = params.ramp_seconds;
    result.extra.u_demand        = u_demand;
    result.extra.u_demand_max    = max(abs(u_demand));
    result.extra.saturation_count = sum(clipped);
    result.extra.t               = test_sig.t;
    result.extra.y_open          = y_open;
    result.extra.y_closed        = y_closed;
    result.extra.u               = u;
    result.extra.theta_final     = theta_final;
    result.extra.norm_theta      = norm(theta_final);
    result.extra.coefficient_time = coeff_time;
    result.extra.coefficient_norm = coeff_norm;
    result.extra.numeric_failure_count = sum(~isfinite(u_demand));
    if isfield(modelData, 'path_family')
        result.extra.path_family = modelData.path_family;
    else
        result.extra.path_family = 'unknown';
    end

    fprintf('  IMC-FxLMS %-7s | %s: %.1f dB | t_conv=%.2fs | max|u|=%.3f | ||w||=%.2f\n', ...
        variant, test_name, result.supp_db, result.t_conv_s, ...
        result.u_max, norm(theta_final));
end

%% ====================================================================
%  仿真引擎: 反馈 IMC-FxLMS
%% ====================================================================

function [y, u, u_demand, clipped, theta_final, coeff_time, coeff_norm] = ...
        sim_imc_fxlms_feedback(d_sig, Tsim, fs, N_fir, ...
        mu, delta, theta_init, ...
        A_sec, B_poly, nA, nB, ...
        actuator_limit, ramp_samples, norm_limit)
%SIM_IMC_FXLMS_FEEDBACK  Pure-feedback IMC-FxLMS simulation engine.
%
%   Information structure: error microphone only, no reference sensor.
%   Internal reference: d_hat(n) = e(n) - g_hat(n) where g_hat = S_hat(z)*u.
%   With perfect secondary-path model (S_hat = S), d_hat = d_p.
%
%   b_plant = B_poly(2:end) includes (d_delay-1) leading zeros for the
%   pure delay — identical to the demo3 secondary-path pattern.

    Nsim = round(Tsim * fs);

    % ---- 状态 ----
    y = zeros(1, Nsim);
    u = zeros(1, Nsim);
    u_demand = zeros(1, Nsim);
    clipped = false(1, Nsim);
    g = zeros(1, Nsim);        % 真实 (也是建模) 次级路径输出: S(z)*u = S_hat(z)*u
    d_hat = zeros(1, Nsim);    % 估计的主噪声 (内部参考)
    x_f = zeros(1, Nsim);      % filtered-x 信号: S_hat(z)*d_hat

    theta = theta_init(:);

    % ---- b_plant 模式 (与 demo3 一致) ----
    b_plant = B_poly(2:end);
    Lb = length(b_plant);
    a_plant = A_sec;

    % ---- 帧记录 ----
    frame_len = max(32, round(0.05 * fs));
    n_frames = ceil(Nsim / frame_len);
    coeff_time = zeros(n_frames, 1);
    coeff_norm = zeros(n_frames, 1);
    fi = 1;

    k_start = max([N_fir, nA + Lb + 2, ramp_samples + 10]);

    for k = k_start:Nsim
        % ---- 1. 次级路径: g = filter(B_poly, A_sec, u) ----
        g(k) = -a_plant(2:end) * g(k-1:-1:max(1,k-nA))' ...
               + b_plant * u(k-1:-1:max(1,k-Lb))';

        % ---- 2. 误差麦克风 ----
        e = d_sig(k) + g(k);
        y(k) = e;

        % ---- 3. 内部参考重构: d_hat = e - g (S_hat = S) ----
        d_hat(k) = e - g(k);

        % ---- 4. Filtered-x: x_f = S_hat(z)*d_hat = filter(B_poly, A_sec, d_hat) ----
        x_f(k) = -a_plant(2:end) * x_f(k-1:-1:max(1,k-nA))' ...
                 + b_plant * d_hat(k-1:-1:max(1,k-Lb))';

        % ---- 5. NLMS 更新 ----
        if k > N_fir && ~clipped(max(1, k-1))
            x_f_vec = x_f(k-(1:N_fir))';
            d_hat_vec = d_hat(k-(1:N_fir))';

            if mu > 0
                nrm = x_f_vec' * x_f_vec + delta;
                theta = theta + mu * e * x_f_vec / nrm;

                theta_norm = norm(theta);
                if isfinite(norm_limit) && theta_norm > norm_limit
                    theta = theta * (norm_limit / theta_norm);
                end
            end

            u_request = -(theta' * d_hat_vec);
        else
            u_request = 0;
        end

        % ---- 6. 启动渐入 + 限幅 ----
        ramp = min(1, max(0, (k - k_start) / max(1, ramp_samples)));
        u_request = ramp * u_request;

        if ~isfinite(u_request)
            u_request = 0;
        end

        [u(k), clipped(k)] = apply_limit(u_request, actuator_limit);
        u_demand(k) = u_request;

        % ---- 帧记录 ----
        if mod(k, frame_len) == 0 && fi <= n_frames
            coeff_time(fi) = k / fs;
            coeff_norm(fi) = norm(theta);
            fi = fi + 1;
        end
    end

    theta_final = theta;
end

%% ====================================================================
%  执行器限幅
%% ====================================================================

function [u_actual, was_clipped] = apply_limit(u_demand, limit)
    u_actual = min(max(u_demand, -limit), limit);
    was_clipped = abs(u_actual - u_demand) > 1e-10;
end
