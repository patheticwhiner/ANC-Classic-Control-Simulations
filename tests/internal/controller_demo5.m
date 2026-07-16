function result = controller_demo5(signals, test_name, variant, params)
%CONTROLLER_DEMO5 Marino-Tomei adaptive frequency estimator for ANC.
%
%   result = controller_demo5(signals, test_name, variant, params)
%
%   基于 Marino & Tomei (2016) 自适应内模的窄带 ANC 控制器。
%   纯反馈结构——仅使用误差麦克风，无需参考麦克风。
%   在线估计窄带噪声频率并通过自适应内模生成抵消信号。
%
%   variant:
%     'fixed'    — 冻结初始频率猜测 (epsilon=0)，展示无自适应的后果
%     'adaptive' — 在线频率估计与自适应抵消
%
%   params 专用字段:
%     .q              — 建模频率分量数，默认 1
%     .k              — 比例误差反馈增益，默认 0.08
%     .epsilon        — 频率自适应增益，默认 1e-5
%     .freq_init_hz   — 初始频率估计 [Hz]，默认 T1.f_hz * 0.9
%     .freq_min_hz    — 频率下界 [Hz]，默认 20
%     .freq_max_hz    — 频率上界 [Hz]，默认 500
%     .method         — 'euler' 或 'exact' (默认) 离散化
%     .actuator_limit — 归一化执行器硬限幅，默认 5
%     .ramp_seconds   — 控制输出启动渐入时间，默认 0.1
%     .dc_cancel      — 是否包含 DC 抵消项 w0，默认 false
%     .case_a_min_ratio — Case A 最小 |Re(S)|/|S|，默认 0.1
%     .feedback_sign  — 'paper' 使用离散有效 sign(Re{S_eff})，'benchmark' 使用其相反号
%     .output_timing  — 'updated' (默认) 输出更新后状态，与参考 Demo 一致；
%                       'previous' 仅用于离散相位边界诊断

    %% ====================================================================
    % §0. 参数解析与模型加载
    %% ====================================================================

    if nargin < 2 || isempty(test_name), test_name = 'T1'; end
    if nargin < 3 || isempty(variant),   variant   = 'adaptive'; end
    if nargin < 4,                       params    = struct(); end

    if ~isfield(params, 'q'),              params.q              = 1; end
    if ~isfield(params, 'actuator_limit'), params.actuator_limit = 5; end
    if ~isfield(params, 'ramp_seconds'),   params.ramp_seconds   = 0.1; end
    if ~isfield(params, 'method'),         params.method         = 'exact'; end
    if ~isfield(params, 'dc_cancel'),      params.dc_cancel      = false; end
    if ~isfield(params, 'case_a_min_ratio'), params.case_a_min_ratio = 0.1; end
    if ~isfield(params, 'feedback_sign'),   params.feedback_sign = 'paper'; end
    if ~isfield(params, 'output_timing'),   params.output_timing = 'updated'; end

    % --- 加载次级路径模型 S(z) = B(z)/A(z) * z^(-d) ---
    if isfield(signals, 'model_data')
        modelData = signals.model_data;
    else
        info = DataManager(signals.model);
        modelData = struct('A', info.model.A, 'B', info.model.B, ...
            'orders', info.orders, 'fs', info.fs);
    end
    A_sec = modelData.A(:).';
    B_sec = modelData.B(:).';
    fs = modelData.fs;
    Ts = 1 / fs;

    % B_poly 保留延迟前导零，与 demo3 的 b_plant 模式一致
    B_poly = B_sec;
    nA = length(A_sec) - 1;
    nB = length(B_poly) - 1;

    % --- 设计频率与参数默认值 ---
    if ~isfield(params, 'freq_init_hz') || isempty(params.freq_init_hz)
        if isfield(signals.(test_name), 'f_hz')
            f_nominal = signals.(test_name).f_hz;
        else
            f_nominal = 357;
        end
        % 初始估计偏离真值 10%，展示自适应修正能力
        params.freq_init_hz = f_nominal * 0.9;
    end
    if ~isfield(params, 'freq_min_hz')
        params.freq_min_hz = 20;
    end
    if ~isfield(params, 'freq_max_hz')
        params.freq_max_hz = min(500, 0.45 * fs);
    end
    if ~isfield(params, 'k')
        params.k = 0.08;
    end
    if ~isfield(params, 'epsilon')
        params.epsilon = 1e-5;
    end

    freq_init = params.freq_init_hz(:);
    if isscalar(freq_init)
        freq_init = repmat(freq_init, params.q, 1);
    end
    q = params.q;
    if numel(freq_init) ~= q
        error('controller_demo5:FrequencyCountMismatch', ...
            'freq_init_hz must be scalar or contain q=%d entries.', q);
    end

    % --- 数字频率参数: theta = (omega*Ts)^2 = (2*pi*f/fs)^2 ---
    theta_init = (2 * pi * freq_init / fs).^2;
    theta_min  = (2 * pi * params.freq_min_hz / fs)^2;
    theta_max  = (2 * pi * params.freq_max_hz / fs)^2;

    % --- feedback polarity ---
    % Keep the paper convention and the benchmark-convention probe explicit;
    % the sign cannot be inferred from a pole-radius check alone because the
    % ANC output direction also determines whether the residual is cancelled.
    params.feedback_sign = validatestring(params.feedback_sign, ...
        {'paper', 'benchmark'});
    params.output_timing = validatestring(params.output_timing, ...
        {'updated', 'previous'});
    sgn = zeros(q, 1);
    rawSgn = zeros(q, 1);
    for i = 1:q
        omega = 2 * pi * freq_init(i);
        z_i = exp(1j * omega * Ts);
        S_w = (B_poly * (z_i.^(-(0:nB))).') / ...
              (A_sec * (z_i.^(-(0:nA))).');
        rawSgn(i) = sign(real(S_w));
        if rawSgn(i) == 0
            rawSgn(i) = 1;
        end
        S_eff = effective_case_a_response(S_w, omega, fs, ...
            params.output_timing, params.method);
        effectiveSign = sign(real(S_eff));
        if effectiveSign == 0
            effectiveSign = 1;
        end
        if strcmp(params.feedback_sign, 'paper')
            sgn(i) = effectiveSign;
        else
            sgn(i) = -effectiveSign;
        end
    end

    % DC 增益符号: sign(S(1)) = sign(sum(B) / sum(A))
    rawSgn0 = sign(sum(B_poly) / sum(A_sec));
    if rawSgn0 == 0, rawSgn0 = 1; end
    if strcmp(params.feedback_sign, 'paper')
        sgn0 = rawSgn0;
    else
        sgn0 = -rawSgn0;
    end

    % Case A requires a reliable, constant sign of the discrete effective
    % real response near every frequency visited by the disturbance. A chirp
    % that crosses the effective zero is outside this implementation's
    % theoretical information structure.
    [caseARatio, caseASign, caseAFrequencies, caseARawRatio, caseARawSign] = ...
        case_a_diagnostics(test_name, signals.(test_name), A_sec, B_poly, fs, ...
        params.output_timing, params.method);
    if strcmp(params.feedback_sign, 'paper')
        initialSignMatches = all(ismember(caseASign, sgn));
    else
        initialSignMatches = all(ismember(-caseASign, sgn));
    end
    caseAApplicable = all(caseARatio >= params.case_a_min_ratio) ...
        && numel(unique(caseASign)) == 1 && initialSignMatches;

    % --- Euler 稳定性检查 ---
    if strcmp(params.method, 'euler')
        omegaTs_max = 2 * pi * max(freq_init) / fs;
        if omegaTs_max > 0.3
            warning(['Demo5: Euler method unstable when omega*Ts=%.2f > 0.3. ' ...
                'Use method=''exact'' instead.'], omegaTs_max);
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
    % §2. 自适应增益 (fixed variant 强制 epsilon=0)
    %% ====================================================================

    switch variant
        case 'fixed'
            epsilon = 0;
        control_gain = params.k;
        case 'adaptive'
            epsilon = params.epsilon;
        control_gain = params.k;
        otherwise
            error('Demo5: 未知 variant ''%s''', variant);
    end

    %% ====================================================================
    % §3. Marino-Tomei 频率估计器仿真
    %% ====================================================================

    [y_closed, u, u_demand, clipped, f_est_log, t_log, theta_final] = ...
        sim_marino_tomei(d_sig, Tsim, fs, Ts, q, ...
        theta_init, theta_min, theta_max, sgn, sgn0, ...
        control_gain, epsilon, params.method, ...
        A_sec, B_poly, nA, nB, ...
        params.actuator_limit, round(params.ramp_seconds * fs), ...
        params.dc_cancel, params.output_timing);

    %% ====================================================================
    % §4. 计算标准化指标
    %% ====================================================================

    meta = struct('demo', 'demo5', 'variant', variant, 'test', test_name);
    result = compute_metrics(y_open, y_closed, u, test_sig, meta);

    %% ====================================================================
    % §5. 追加 demo5 专用字段
    %% ====================================================================

    result.extra.actuator_limit  = params.actuator_limit;
    result.extra.ramp_seconds    = params.ramp_seconds;
    result.extra.u_demand        = u_demand;
    result.extra.u_demand_max    = max(abs(u_demand));
    result.extra.saturation_count = sum(clipped);
    result.extra.t               = test_sig.t;
    result.extra.y_open          = y_open;
    result.extra.y_closed        = y_closed;
    result.extra.u               = u;
    result.extra.q               = q;
    result.extra.k               = control_gain;
    result.extra.epsilon         = epsilon;
    result.extra.freq_init_hz    = freq_init;
    % --- dSPACE/Simulink 移植所需的内模配置（纯输出附加，不影响仿真数值） ---
    result.extra.theta_init      = theta_init;
    result.extra.theta_min       = theta_min;
    result.extra.theta_max       = theta_max;
    result.extra.freq_min_hz     = params.freq_min_hz;
    result.extra.freq_max_hz     = params.freq_max_hz;
    result.extra.dc_cancel       = params.dc_cancel;
    result.extra.ramp_samples    = round(params.ramp_seconds * fs);
    result.extra.k_start         = max([nA + (length(B_poly) - 1) + 2, ...
        round(params.ramp_seconds * fs) + 10]);
    result.extra.A_sec           = A_sec;
    result.extra.B_poly          = B_poly;
    result.extra.freq_final_hz   = sqrt(theta_final) * fs / (2 * pi);
    result.extra.theta_final     = theta_final;
    result.extra.method           = params.method;
    result.extra.feedback_sign    = params.feedback_sign;
    result.extra.output_timing    = params.output_timing;
    result.extra.f_est_log       = f_est_log;
    result.extra.t_log           = t_log;
    result.extra.sign_S          = sgn;
    result.extra.sign_S0         = sgn0;
    result.extra.raw_sign_S      = rawSgn;
    result.extra.raw_sign_S0     = rawSgn0;
    result.extra.case_a_ratio    = caseARatio;
    result.extra.case_a_sign     = caseASign;
    result.extra.case_a_raw_ratio = caseARawRatio;
    result.extra.case_a_raw_sign  = caseARawSign;
    result.extra.case_a_response = 'discrete_effective';
    result.extra.case_a_phase_rule = 'updated:+alpha/2, previous:-alpha/2';
    if isfield(modelData, 'path_family')
        result.extra.path_family = modelData.path_family;
    else
        result.extra.path_family = 'unknown';
    end
    result.extra.case_a_frequencies_hz = caseAFrequencies;
    result.extra.case_a_min_ratio = min(caseARatio);
    result.extra.case_a_threshold = params.case_a_min_ratio;
    result.extra.case_a_sign_consistent = numel(unique(caseASign)) == 1;
    result.extra.case_a_initial_sign_matches = initialSignMatches;
    result.extra.case_a_applicable = caseAApplicable;
    result.extra.numeric_failure_count = sum(~isfinite(u_demand));
    result.extra.frequency_bound_hit = any( ...
        result.extra.freq_final_hz <= params.freq_min_hz * (1 + 1e-9) | ...
        result.extra.freq_final_hz >= params.freq_max_hz * (1 - 1e-9));

    if strcmp(test_name, 'T1')
        result.extra.frequency_error_final_hz = ...
            result.extra.freq_final_hz(1) - test_sig.f_hz;
    elseif strcmp(test_name, 'T2') && isfield(test_sig, 'f_inst')
        result.extra.f_inst = test_sig.f_inst(:);
        fTruth = interp1(test_sig.t(:), test_sig.f_inst(:), t_log(:), ...
            'linear', 'extrap');
        fEstimate = f_est_log(1, :).';
        steady = t_log(:) >= 3;
        result.extra.frequency_tracking_error_rms_hz = ...
            sqrt(mean((fEstimate(steady) - fTruth(steady)).^2));
    end

    fprintf(['  Demo5 %-9s | %s: %.1f dB | t_conv=%.2fs | max|u|=%.3f ' ...
        '| f_final=%.2f Hz | Case A=%d\n'], ...
        variant, test_name, result.supp_db, result.t_conv_s, ...
        result.u_max, result.extra.freq_final_hz(1), caseAApplicable);
end

function [ratio, signs, frequencies, rawRatio, rawSigns] = case_a_diagnostics( ...
        testName, testSignal, A, B, fs, outputTiming, method)
if strcmp(testName, 'T1') && isfield(testSignal, 'f_hz')
    frequencies = testSignal.f_hz;
elseif isfield(testSignal, 'f_inst')
    frequencies = linspace(min(testSignal.f_inst), max(testSignal.f_inst), 121);
elseif isfield(testSignal, 'f_range')
    frequencies = linspace(min(testSignal.f_range), max(testSignal.f_range), 121);
else
    frequencies = linspace(20, min(500, 0.45 * fs), 121);
end

response = zeros(size(frequencies));
effectiveResponse = zeros(size(frequencies));
for index = 1:numel(frequencies)
    z = exp(1j * 2 * pi * frequencies(index) / fs);
    response(index) = (B * (z.^(-(0:numel(B)-1))).') / ...
        (A * (z.^(-(0:numel(A)-1))).');
    effectiveResponse(index) = effective_case_a_response( ...
        response(index), 2 * pi * frequencies(index), fs, ...
        outputTiming, method);
end
rawRatio = abs(real(response)) ./ max(abs(response), eps);
rawSigns = sign(real(response));
ratio = abs(real(effectiveResponse)) ./ max(abs(effectiveResponse), eps);
signs = sign(real(effectiveResponse));
rawSigns(rawSigns == 0) = 1;
signs(signs == 0) = 1;
end

function response = effective_case_a_response(response, omega, fs, outputTiming, method)
% Exact inner-model updates have a half-sample phase relative to e(k).
% Keep the raw criterion for Euler, which has no equivalent exact phase.
if strcmp(method, 'exact')
    alpha = omega / fs;
    if strcmp(outputTiming, 'updated')
        response = response * exp(1j * alpha / 2);
    else
        response = response * exp(-1j * alpha / 2);
    end
end
end

%% ====================================================================
%  仿真引擎: Marino-Tomei 自适应频率估计 ANC
%% ====================================================================

function [y, u, u_demand, clipped, f_est_log, t_log, theta_final] = ...
        sim_marino_tomei(d_sig, Tsim, fs, Ts, q, ...
        theta_init, theta_min, theta_max, sgn, sgn0, ...
        control_gain, epsilon, method, ...
        A_sec, B_poly, nA, nB, ...
        actuator_limit, ramp_samples, dc_cancel, output_timing)
%SIM_MARINO_TOMEI  Marino-Tomei 自适应频率估计 ANC 仿真引擎.
%
%   次级路径 g = filter(B_poly, A_sec, u)，其中 B_poly 已包含延迟前导零。
%   与 demo3 一致：b_plant = B_poly(2:end), Lb = nB。
%
%   b_plant 的前 (d_delay-1) 个元素为零，自动实现 z^{-d} 延迟。

    Nsim = round(Tsim * fs);

    % ---- 状态初始化 ----
    y = zeros(1, Nsim);
    u = zeros(1, Nsim);
    u_demand = zeros(1, Nsim);
    clipped = false(1, Nsim);
    g = zeros(1, Nsim);

    % 内模状态
    w0 = 0;
    w = zeros(2, q);
    theta = theta_init(:);

    % ---- b_plant 模式 (与 demo3 一致) ----
    b_plant = B_poly(2:end);
    Lb = length(b_plant);
    a_plant = A_sec;

    % ---- 日志 ----
    log_interval = max(1, round(Nsim / 800));
    n_log = floor(Nsim / log_interval) + 1;
    f_est_log = zeros(q, n_log);
    t_log = zeros(n_log, 1);
    li = 1;

    % ---- 控制启动 ----
    k_start = max([nA + Lb + 2, ramp_samples + 10]);

    for sampleIndex = 1:Nsim
        % ====== 1. 次级路径: g = filter(B_poly, A_sec, u) ======
        if sampleIndex > k_start
            g(sampleIndex) = -a_plant(2:end) * ...
                g(sampleIndex-1:-1:max(1,sampleIndex-nA))' ...
                + b_plant * u(sampleIndex-1:-1:max(1,sampleIndex-Lb))';
        end

        % ====== 2. 误差麦克风 ======
        e = d_sig(sampleIndex) + g(sampleIndex);
        y(sampleIndex) = e;

        % ====== 3. Marino-Tomei 控制器 ======
        if sampleIndex >= k_start
            w_old = w;
            w0_old = w0;
            theta_old = theta;

            if dc_cancel
                w0 = w0 + sgn0 * control_gain * e;
            end

            for i = 1:q
                switch method
                    case 'euler'
                        w(1,i) = w_old(1,i) + Ts * ...
                            (w_old(2,i) + sgn(i) * control_gain * e);
                        w(2,i) = w_old(2,i) + Ts * (-theta_old(i) * w_old(1,i));

                    case 'exact'
                        alpha = sqrt(max(theta_old(i), 1e-12));
                        ca = cos(alpha);
                        sa = sin(alpha);

                        w(1,i) = ca * w_old(1,i) + (sa / alpha) * w_old(2,i) ...
                               + sgn(i) * control_gain * (sa / alpha) * e;
                        w(2,i) = -alpha * sa * w_old(1,i) + ca * w_old(2,i) ...
                               + sgn(i) * control_gain * (ca - 1) * e;

                    otherwise
                        error('Unknown method: %s.', method);
                end

                if epsilon > 0
                    theta(i) = theta_old(i) + epsilon * sgn(i) * w_old(2,i) * e;
                    theta(i) = min(max(theta(i), theta_min), theta_max);
                end
            end

            if strcmp(output_timing, 'previous')
                u_request = -(w0_old + sum(w_old(1,:)));
            else
                u_request = -(w0 + sum(w(1,:)));
            end

            ramp = min(1, max(0, ...
                (sampleIndex - k_start) / max(1, ramp_samples)));
            u_request = ramp * u_request;
        else
            u_request = 0;
        end

        u_demand(sampleIndex) = u_request;
        [u(sampleIndex), clipped(sampleIndex)] = ...
            apply_limit(u_request, actuator_limit);

        % ====== 日志 ======
        if mod(sampleIndex, log_interval) == 0 || sampleIndex == 1
            t_log(li) = (sampleIndex - 1) * Ts;
            f_est_log(:, li) = sqrt(theta) * fs / (2 * pi);
            li = li + 1;
        end
    end

    t_log = t_log(1:li-1);
    f_est_log = f_est_log(:, 1:li-1);
    theta_final = theta;
end

%% ====================================================================
%  执行器限幅
%% ====================================================================

function [u_actual, was_clipped] = apply_limit(u_demand, limit)
    if isnan(u_demand)
        u_actual = 0;
        was_clipped = true;
    else
        u_actual = min(max(u_demand, -limit), limit);
        was_clipped = ~isfinite(u_demand) || ...
            abs(u_actual - u_demand) > 1e-10;
    end
end
