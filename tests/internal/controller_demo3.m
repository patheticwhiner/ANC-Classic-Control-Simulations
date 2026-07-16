function result = controller_demo3(signals, test_name, variant, params)
%CONTROLLER_DEMO3 Demo3: F(z)+FIR + Jafari NLMS 自适应.
%
%   用法:
%     result = controller_demo3(signals, test_name, variant, params)
%
%   signals 由 load_cylinder1dm_signals 提供，本函数只做计算。
%
%   variant:
%     'fixed'    — 固定 FIR K(z,θ) 从 H∞ 插值设计
%     'adaptive' — NLMS 在线更新 FIR，从 fixed θ 初始化，
%                  使用 H(z) 滤波误差作 regressor (Jafari 方法)
%
%   params 专用字段:
%     .N_fir    — 固定 FIR 阶数，默认按采样率缩放 (120@48kHz)
%     .N_nlms   — NLMS 自适应 FIR 阶数，默认同 N_fir
%     .mu_nlms  — NLMS 步长，默认 0.05
%     .control_scale — 控制输出缩放，默认 1
%     .theta_norm_limit — 自适应 FIR 范数上限，默认 Inf
%     .adaptive_structure — 'jafari' (default) or 'imc_fxnlms'
%     .f_design — 固定 FIR 设计频率 (Hz)，默认取 T1.f_hz
%     .actuator_limit — 归一化执行器硬限幅，默认 5
%     .ramp_seconds — 控制输出启动渐入时间，默认 0.1 s

    %% ====================================================================
    % §0. 参数解析与模型加载
    %% ====================================================================

    if nargin < 2 || isempty(test_name), test_name = 'T1'; end
    if nargin < 3 || isempty(variant), variant = 'fixed'; end
    if nargin < 4, params = struct(); end
    if ~isfield(params, 'N_fir'),   params.N_fir   = []; end
    if ~isfield(params, 'N_nlms'),  params.N_nlms  = []; end
    if ~isfield(params, 'mu_nlms'), params.mu_nlms = 0.05; end
    if ~isfield(params, 'control_scale'), params.control_scale = 1; end
    if ~isfield(params, 'theta_norm_limit'), params.theta_norm_limit = Inf; end
    if ~isfield(params, 'adaptive_structure'), params.adaptive_structure = 'jafari'; end
    if ~isfield(params, 'f_design'), params.f_design = []; end
    if ~isfield(params, 'actuator_limit'), params.actuator_limit = 5; end
    if ~isfield(params, 'ramp_seconds'), params.ramp_seconds = 0.1; end

    if isfield(signals, 'model_data')
        modelData = signals.model_data;
    else
        info = DataManager(signals.model);
        modelData = struct('A', info.model.A, 'B', info.model.B, ...
            'orders', info.orders, 'fs', info.fs);
    end
    A = modelData.A;
    B_poly = modelData.B;
    d = modelData.orders(4);
    fs = modelData.fs;
    Ts = 1/fs;
    if isempty(params.N_fir)
        params.N_fir = max(16, round(120 * fs / 48000));
    end
    if isempty(params.N_nlms)
        params.N_nlms = max(16, round(64 * fs / 48000));
    end
    if isempty(params.f_design)
        params.f_design = signals.T1.f_hz;
    end

    B = B_poly(d+1:end);
    nA = length(A) - 1;
    nB = length(B) - 1;

    %% ====================================================================
    % §1. F(z) 频谱展平 + 固定 FIR 控制器设计
    %% ====================================================================

    persistent ctrl_cache;
    if isempty(ctrl_cache)
        ctrl_cache = struct('model', '', 'theta_fixed', [], ...
            'F_num', [], 'F_den', [], 'H_num', [], 'H_den', [], ...
            'N_fir', [], 'gamma', []);
    end

    design_key = sprintf('%s_N%d_f%.3f', ...
        signals.model, params.N_fir, params.f_design);

    if ~strcmp(ctrl_cache.model, design_key)
        fprintf('--- Demo3: F(z) 频谱展平 + FIR 设计 ---\n');

        % NMP 零点 → 单位圆镜像反射
        z_orig = roots(B);
        z_nmp  = z_orig(abs(z_orig) > 1);
        z_mp   = z_orig(abs(z_orig) <= 1);
        B_tilde = real(poly([z_mp; 1./conj(z_nmp)]));
        B_tilde = B_tilde * (B(1) / B_tilde(1));

        % 一阶低通 L(z)
        flattenCutoff = min(2000, 0.4 * fs);
        rho = exp(-2*pi*flattenCutoff*Ts);
        L_num = 1 - rho;
        L_den = [1, -rho];

        % F(z) = A / (L * B_tilde)
        F_num = A;
        F_den = conv(L_den, B_tilde);
        F_den = F_den / F_den(1);

        % 展平后次级通道: H(z) = (z^{-d}*B*F)/A
        H_num = conv([zeros(1,d), B], F_num);
        H_den = conv(A, F_den);

        % 固定 FIR: 插值约束 |G_eff(ω₀)·K(ω₀)| = 1
        f0 = params.f_design;
        omega0 = 2*pi*f0 / fs;

        N_fir = params.N_fir;
        zi0 = exp(1j*omega0);
        G_eff_0 = polyval_freqz(H_num, H_den, omega0);
        firDelays = 1:N_fir;
        Aeq = [real(G_eff_0 * zi0.^(-firDelays));
               imag(G_eff_0 * zi0.^(-firDelays))];
        theta_fixed = Aeq \ [1; 0];

        % γ 值
        Nfreq = 2000;
        w_grid = linspace(0, pi, Nfreq);
        G_eff_grid = freqz(H_num, H_den, w_grid).';
        gamma = max(abs(G_eff_grid(:) .* ...
            (exp(-1j*firDelays'*w_grid)' * theta_fixed)));

        % 缓存
        ctrl_cache.model       = design_key;
        ctrl_cache.theta_fixed = theta_fixed;
        ctrl_cache.F_num       = F_num;
        ctrl_cache.F_den       = F_den;
        ctrl_cache.H_num       = H_num;
        ctrl_cache.H_den       = H_den;
        ctrl_cache.N_fir       = N_fir;
        ctrl_cache.gamma       = gamma;

        fprintf('  f₀=%.1f Hz, N_fir=%d, ||θ||=%.2f, γ=%.4f\n', ...
            f0, N_fir, norm(theta_fixed), gamma);
    else
        theta_fixed = ctrl_cache.theta_fixed;
        F_num       = ctrl_cache.F_num;
        F_den       = ctrl_cache.F_den;
        H_num       = ctrl_cache.H_num;
        H_den       = ctrl_cache.H_den;
        N_fir       = ctrl_cache.N_fir;
        gamma       = ctrl_cache.gamma;
    end

    %% ====================================================================
    % §2. 被控对象参数
    %% ====================================================================
    a_plant = A;
    b_plant = B_poly(2:end);
    Lb = length(b_plant);
    nF_num = length(F_num);
    nF_den = length(F_den) - 1;
    k_start = max(N_fir, nA + Lb) + 1;

    %% ====================================================================
    % §3. 选择测试信号 + 仿真
    %% ====================================================================
    test_sig = signals.(test_name);
    d_sig    = test_sig.d(:)';
    Nsim     = length(d_sig);
    Tsim     = test_sig.Tsim;
    y_open   = test_sig.y_open(:)';

    switch variant
        case 'fixed'
            [y_closed, u, u_demand, clipped] = sim_fixed_jafari( ...
                d_sig, Tsim, fs, Ts, ...
                N_fir, theta_fixed, a_plant, b_plant, nA, Lb, ...
                nF_num, nF_den, F_num, F_den, k_start, ...
                params.actuator_limit, round(params.ramp_seconds * fs), params.control_scale, params.theta_norm_limit);

        case 'adaptive'
            N_adapt = params.N_nlms;
            if strcmp(params.adaptive_structure, 'imc_fxnlms')
                theta_init = zeros(N_adapt, 1);
                [y_closed, u, theta_final, u_demand, clipped, ...
                    coefficientTime, coefficientNorm, internalError] = ...
                    sim_imc_fxnlms(d_sig, Tsim, fs, N_adapt, ...
                    params.mu_nlms, theta_init, a_plant, b_plant, nA, Lb, ...
                    params.actuator_limit, round(params.ramp_seconds * fs), ...
                    params.control_scale, params.theta_norm_limit);
            else
                % 从固定 FIR 初始化 Jafari NLMS (补零或截断)
                if N_adapt > N_fir
                    theta_init = [theta_fixed; zeros(N_adapt - N_fir, 1)];
                else
                    theta_init = theta_fixed(1:N_adapt);
                end
                [y_closed, u, theta_final, u_demand, clipped, ...
                    coefficientTime, coefficientNorm, internalError] = ...
                    sim_nlms_jafari(d_sig, Tsim, fs, Ts, ...
                    N_adapt, params.mu_nlms, theta_init, ...
                    a_plant, b_plant, nA, Lb, ...
                    nF_num, nF_den, F_num, F_den, ...
                    H_num, H_den, k_start, ...
                    params.actuator_limit, round(params.ramp_seconds * fs), ...
                    params.control_scale, params.theta_norm_limit);
            end

        otherwise
            error('Demo3: 未知 variant ''%s''', variant);
    end

    %% ====================================================================
    % §4. 计算标准化指标
    %% ====================================================================
    meta = struct('demo', 'demo3', 'variant', variant, 'test', test_name);
    result = compute_metrics(y_open, y_closed, u, test_sig, meta);

    %% ====================================================================
    % §5. 追加 demo3 专用字段
    %% ====================================================================
    result.extra.N_fir = N_fir;
    result.extra.gamma = gamma;
    result.extra.f_design = params.f_design;
    result.extra.ramp_seconds = params.ramp_seconds;
    result.extra.theta_fixed = theta_fixed;
    result.extra.F_num = F_num;
    result.extra.F_den = F_den;
    result.extra.H_num = H_num;
    result.extra.H_den = H_den;
    result.extra.u_demand = u_demand;
    result.extra.u_demand_max = max(abs(u_demand));
    result.extra.saturation_count = sum(clipped);
    result.extra.actuator_limit = params.actuator_limit;
    result.extra.t = test_sig.t;
    result.extra.y_open = y_open;
    result.extra.y_closed = y_closed;
    result.extra.u = u;

    if strcmp(variant, 'adaptive')
        result.extra.N_adapt  = params.N_nlms;
        result.extra.mu_nlms = params.mu_nlms;
        result.extra.adaptive_structure = params.adaptive_structure;
        result.extra.norm_theta = norm(theta_final);
        result.extra.coefficient_time = coefficientTime;
        result.extra.coefficient_norm = coefficientNorm;
        result.extra.internal_error = internalError;
        result.extra.theta_init = theta_init;
    end

    fprintf('  Demo3 %-9s | %s: %.1f dB | t_conv=%.2fs | max|u|=%.3f\n', ...
        variant, test_name, result.supp_db, result.t_conv_s, result.u_max);
end

function [y, u, thetaFinal, uDemand, clipped, coefficientTime, ...
        coefficientNorm, internalError] = sim_imc_fxnlms( ...
        d_sig, Tsim, fs, N_adapt, mu, theta_init, ...
        a_plant, b_plant, nA, Lb, actuatorLimit, rampSamples, ...
        controlScale, thetaNormLimit)
%SIM_IMC_FXNLMS Standard source-reference filtered-x NLMS baseline.

    Nsim = round(Tsim * fs);
    y = zeros(1, Nsim);
    u = zeros(1, Nsim);
    uDemand = zeros(1, Nsim);
    clipped = false(1, Nsim);
    g = zeros(1, Nsim);
    filteredReference = zeros(1, Nsim);
    theta = theta_init(:);

    frameLength = max(32, round(0.05 * fs));
    frameCount = ceil(Nsim / frameLength);
    coefficientTime = zeros(frameCount, 1);
    coefficientNorm = zeros(frameCount, 1);
    frameIndex = 1;
    kStart = max([N_adapt, nA, Lb]) + 1;

    for k = kStart:Nsim
        g(k) = -a_plant(2:end) * g(k-1:-1:k-nA)' ...
               + b_plant * u(k-1:-1:k-Lb)';
        y(k) = g(k) + d_sig(k);

        filteredReference(k) = ...
            -a_plant(2:end) * filteredReference(k-1:-1:k-nA)' ...
            + b_plant * d_sig(k-1:-1:k-Lb)';

        referenceVector = d_sig(k-(1:N_adapt))';
        filteredVector = filteredReference(k-(1:N_adapt))';
        if ~clipped(k-1)
            theta = theta + mu * y(k) * filteredVector ...
                / (filteredVector' * filteredVector + 0.01);
            thetaNorm = norm(theta);
            if isfinite(thetaNormLimit) && thetaNorm > thetaNormLimit
                theta = theta * (thetaNormLimit / thetaNorm);
            end
        end

        ramp = min(1, max(0, (k - kStart) / max(1, rampSamples)));
        uRequest = -controlScale * ramp * (theta' * referenceVector);
        if ~isfinite(uRequest), uRequest = 0; end
        [u(k), clipped(k)] = apply_actuator_limit(uRequest, actuatorLimit);
        uDemand(k) = uRequest;

        if mod(k, frameLength) == 0 && frameIndex <= frameCount
            coefficientTime(frameIndex) = k / fs;
            coefficientNorm(frameIndex) = norm(theta);
            frameIndex = frameIndex + 1;
        end
    end

    internalError = y;
    thetaFinal = theta;
end

%% ====================================================================
% 仿真引擎 (提取自 JafariANC_RealAcoustic.m)
%% ====================================================================

function [y, u, u_demand, clipped] = sim_fixed_jafari( ...
        d_sig, Tsim, fs, Ts, N_fir, theta_fixed, ...
        a_plant, b_plant, nA, Lb, nF_num, nF_den, F_num, F_den, ...
        k_start, actuatorLimit, rampSamples, controlScale, thetaNormLimit)
    Nsim = round(Tsim * fs);
    y = zeros(1, Nsim);
    u = zeros(1, Nsim);
    u_demand = zeros(1, Nsim);
    clipped = false(1, Nsim);
    g = zeros(1, Nsim);
    z = zeros(1, Nsim);
    xF = zeros(nF_num, 1);
    yF = zeros(max(1, nF_den), 1);

    for k = k_start:Nsim
        % 次级路径: g = (B/A)*u (抗噪声)
        g(k) = -a_plant(2:end) * g(k-1:-1:max(1,k-nA))' ...
               + b_plant * u(k-1:-1:max(1,k-Lb))';
        % 残差 = 主噪声 + 抗噪声
        y(k) = g(k) + d_sig(k);
        z(k) = y(k) - g(k);

        if k > N_fir
            v = theta_fixed' * z(k-(1:N_fir))';
            xF = [v; xF(1:nF_num-1)];
            uk = F_num * xF - F_den(2:end) * yF(1:nF_den);
            yF = [uk; yF(1:end-1)];
            ramp = min(1, max(0, (k - k_start) / max(1, rampSamples)));
            u_request = -controlScale * ramp * uk;
            if ~isfinite(u_request)
                u_request = 0;
            end
            [u(k), clipped(k)] = apply_actuator_limit(u_request, actuatorLimit);
            u_demand(k) = u_request;
        end
    end
end

function [y, u, thetaFinal, uDemand, clipped, coefficientTime, ...
        coefficientNorm, internalError] = sim_nlms_jafari( ...
        d_sig, Tsim, fs, Ts, N_adapt, mu, theta_init, ...
        a_plant, b_plant, nA, Lb, nF_num, nF_den, F_num, F_den, ...
        H_num, H_den, k_start, actuatorLimit, rampSamples, controlScale, thetaNormLimit)
%SIM_NLMS_JAFARI  Jafari NLMS-on-FIR 自适应仿真引擎.
%
%   固定 FIR theta_fixed 作为 NLMS 初始值，在线更新 theta。
%   H(z) = B·F(z)/A 用于生成滤波误差 zf(k) 作为 NLMS regressor，
%   FIR 输出 v = theta' * z(k-N_adapt+1:k) 经 F(z) 产生控制量。
%
%   参考: Jafari & Ioannou, "Adaptive Noise Cancellation,"
%         Rejection of Sinusoidal and Narrowband Disturbances.

    Nsim = round(Tsim * fs);
    y = zeros(1, Nsim);
    u = zeros(1, Nsim);
    uDemand = zeros(1, Nsim);
    clipped = false(1, Nsim);
    g = zeros(1, Nsim);
    z = zeros(1, Nsim);
    zf = zeros(1, Nsim);

    % NLMS 状态
    theta = theta_init(:);
    th = zeros(N_adapt, Nsim);

    % F(z) 滤波器状态
    xF = zeros(nF_num, 1);
    yF = zeros(max(1, nF_den), 1);

    % H(z) 滤波器状态
    xH = zeros(length(H_num) - 1, 1);
    yH = zeros(length(H_den) - 1, 1);

    % 记录收敛曲线 (每 0.05 s 取一帧)
    frameLength = max(32, round(0.05 * fs));
    frameCount = ceil(Nsim / frameLength);
    coefficientTime = zeros(frameCount, 1);
    coefficientNorm = zeros(frameCount, 1);
    internalError = zeros(1, Nsim);
    frameIndex = 1;

    regulatorOnset = 0;

    for k = k_start:Nsim
        % 次级路径: g = (B/A)*u
        g(k) = -a_plant(2:end) * g(k-1:-1:max(1,k-nA))' ...
               + b_plant * u(k-1:-1:max(1,k-Lb))';
        % 残差 = 主噪声 + 抗噪声
        y(k) = g(k) + d_sig(k);
        z(k) = y(k) - g(k);

        % H(z) 滤波误差
        xH = [z(k); xH(1:end-1)];
        zf(k) = H_num(2:end) * xH - H_den(2:end) * yH;
        yH = [zf(k); yH(1:end-1)];

        % NLMS 更新 与 FIR 输出（需要足够长的 z 历史）
        if k > N_adapt && ~clipped(max(1, k-1))
            phi = zf(k-(1:N_adapt))';
            % Use the measured closed-loop residual for the filtered-x
            % gradient. The former one-step identification error ignored
            % the actual controller/plant state and destroyed a good FIR
            % initialisation even on the fixed-frequency case.
            e_nlms = y(k);
            theta = theta + mu * e_nlms * phi / (phi' * phi + 0.01);
            thetaNorm = norm(theta);
            if isfinite(thetaNormLimit) && thetaNorm > thetaNormLimit
                theta = theta * (thetaNormLimit / thetaNorm);
            end
            th(:, k) = theta;

            v = theta' * z(k-(1:N_adapt))';
            xF = [v; xF(1:nF_num-1)];
            uk = F_num * xF - F_den(2:end) * yF(1:nF_den);
            yF = [uk; yF(1:end-1)];
        else
            uk = 0;
        end

        % 启动渐入
        if regulatorOnset == 0
            regulatorOnset = k;
        end
        ramp = min(1, max(0, (k - regulatorOnset) / max(1, rampSamples)));
        u_request = -controlScale * ramp * uk;

        if ~isfinite(u_request)
            u_request = 0;
        end
        [u(k), clipped(k)] = apply_actuator_limit(u_request, actuatorLimit);
        uDemand(k) = u_request;

        % 帧记录
        if mod(k, frameLength) == 0 && frameIndex <= frameCount
            coefficientTime(frameIndex) = k / fs;
            coefficientNorm(frameIndex) = norm(theta);
            frameIndex = frameIndex + 1;
        end
    end

    internalError = zf;
    thetaFinal = theta;
end

function v = polyval_freqz(num, den, w)
    v = sum(num .* exp(-1j*w*(0:length(num)-1))) / ...
        sum(den .* exp(-1j*w*(0:length(den)-1)));
end

function [u_actual, clipped] = apply_actuator_limit(u_demand, limit)
% 归一化执行器硬限幅
    u_actual = min(max(u_demand, -limit), limit);
    clipped = abs(u_actual - u_demand) > 1e-10;
end
