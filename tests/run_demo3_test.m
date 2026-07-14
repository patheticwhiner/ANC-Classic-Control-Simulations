function result = run_demo3_test(signals, test_name, variant, params)
% RUN_DEMO3_TEST  Demo3: Jafari F(z) + FIR(固定) + NLMS(自适应)
%
%   用法:
%     result = run_demo3_test(signals, test_name, variant, params)
%
%   params 专用字段:
%     .N_fir    — 固定 FIR 阶数，默认 120
%     .N_nlms   — NLMS 自适应阶数，默认 64
%     .mu_nlms  — NLMS 步长，默认 0.05

    %% ====================================================================
    % §0. 参数解析与模型加载
    %% ====================================================================

    if nargin < 4, params = struct(); end
    if ~isfield(params, 'N_fir'),   params.N_fir   = 120; end
    if ~isfield(params, 'N_nlms'),  params.N_nlms  = 64; end
    if ~isfield(params, 'mu_nlms'), params.mu_nlms = 0.05; end

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(projectRoot, 'project_init.m'));

    info = DataManager(signals.model);
    A = info.model.A;
    B_poly = info.model.B;
    d = info.orders(4);
    fs = info.fs;
    Ts = 1/fs;

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

    design_key = sprintf('%s_N%d', signals.model, params.N_fir);

    if ~strcmp(ctrl_cache.model, design_key)
        fprintf('--- Demo3: F(z) 频谱展平 + FIR 设计 ---\n');

        % NMP 零点 → 单位圆镜像反射
        z_orig = roots(B);
        z_nmp  = z_orig(abs(z_orig) > 1);
        z_mp   = z_orig(abs(z_orig) <= 1);
        B_tilde = real(poly([z_mp; 1./conj(z_nmp)]));
        B_tilde = B_tilde * (B(1) / B_tilde(1));

        % 一阶低通 L(z)
        rho = exp(-2*pi*2000*Ts);
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
        w_scan = linspace(0.001, pi, 8000);
        G_scan = abs(freqz(B, A, w_scan));
        [~, locs] = findpeaks(G_scan, 'MinPeakHeight', 1.0, 'MinPeakDistance', 50);
        f_peaks = w_scan(locs) * fs / (2*pi);
        f0 = f_peaks(1);
        omega0 = 2*pi*f0 / fs;

        N_fir = params.N_fir;
        zi0 = exp(1j*omega0);
        G_eff_0 = polyval_freqz(H_num, H_den, omega0);
        Aeq = [real(G_eff_0 * zi0.^(-(0:N_fir-1)));
               imag(G_eff_0 * zi0.^(-(0:N_fir-1)))];
        theta_fixed = Aeq \ [1; 0];

        % γ 值
        Nfreq = 2000;
        w_grid = linspace(0, pi, Nfreq);
        G_eff_grid = freqz(H_num, H_den, w_grid).';
        gamma = max(abs(G_eff_grid(:) .* ...
            (exp(-1j*(0:N_fir-1)'*w_grid)' * theta_fixed)));

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
    b_plant = [zeros(1,d), B];
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
            [y_closed, u] = sim_fixed_jafari(d_sig, Tsim, fs, Ts, ...
                N_fir, theta_fixed, a_plant, b_plant, nA, Lb, ...
                nF_num, nF_den, F_num, F_den, k_start);

        case 'adaptive'
            [y_closed, u, theta_final] = sim_nlms_jafari(d_sig, Tsim, fs, Ts, ...
                params.N_nlms, params.mu_nlms, a_plant, b_plant, nA, Lb, ...
                nF_num, nF_den, F_num, F_den, H_num, H_den, k_start);

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

    if strcmp(variant, 'adaptive')
        result.extra.N_nlms  = params.N_nlms;
        result.extra.mu_nlms = params.mu_nlms;
        result.extra.norm_theta = norm(theta_final);
    end

    fprintf('  Demo3 %-9s | %s: %.1f dB | t_conv=%.2fs | max|u|=%.3f\n', ...
        variant, test_name, result.supp_db, result.t_conv_s, result.u_max);
end

%% ====================================================================
% 仿真引擎 (提取自 JafariANC_RealAcoustic.m)
%% ====================================================================

function [y, u] = sim_fixed_jafari(d_sig, Tsim, fs, Ts, N_fir, theta_fixed, ...
        a_plant, b_plant, nA, Lb, nF_num, nF_den, F_num, F_den, k_start)
    Nsim = round(Tsim * fs);
    y = zeros(1, Nsim);
    u = zeros(1, Nsim);
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
            u(k) = -uk;
            if isnan(u(k)) || abs(u(k)) > 1e4
                u(k) = 0;
            end
        end
    end
end

function [y, u, theta_final] = sim_nlms_jafari(d_sig, Tsim, fs, Ts, ...
        N_adapt, mu, a_plant, b_plant, nA, Lb, ...
        nF_num, nF_den, F_num, F_den, H_num, H_den, k_start)
    Nsim = round(Tsim * fs);
    y = zeros(1, Nsim);
    u = zeros(1, Nsim);
    g = zeros(1, Nsim);
    z = zeros(1, Nsim);
    zf = zeros(1, Nsim);
    theta = zeros(N_adapt, 1);
    xF = zeros(nF_num, 1);
    yF = zeros(max(1, nF_den), 1);
    xH = zeros(length(H_num)-1, 1);
    yH = zeros(max(1, length(H_den)-1), 1);

    for k = k_start:Nsim
        % 次级路径
        g(k) = -a_plant(2:end) * g(k-1:-1:max(1,k-nA))' ...
               + b_plant * u(k-1:-1:max(1,k-Lb))';
        y(k) = g(k) + d_sig(k);
        z(k) = y(k) - g(k);

        % 展平滤波: zf = H(z) * z
        xH = [z(k); xH(1:end-1)];
        zf(k) = H_num(2:end) * xH(1:length(H_num)-1) ...
                - H_den(2:end) * yH(1:length(H_den)-1);
        yH = [zf(k); yH(1:end-1)];

        % NLMS 自适应
        if k > N_adapt
            phi = zf(k-(1:N_adapt))';
        else
            phi = zeros(N_adapt, 1);
        end
        e = z(k) - theta' * phi;
        theta = theta + mu * e * phi / (phi'*phi + 0.01);

        % F(z) 重建控制信号
        v = theta' * z(k-(1:N_adapt))';
        xF = [v; xF(1:nF_num-1)];
        uk = F_num * xF - F_den(2:end) * yF(1:nF_den);
        yF = [uk; yF(1:end-1)];
        u(k) = -uk;
        if isnan(u(k)) || abs(u(k)) > 1e4
            u(k) = 0;
        end
    end
    theta_final = theta;
end

function v = polyval_freqz(num, den, w)
    v = sum(num .* exp(-1j*w*(0:length(num)-1))) / ...
        sum(den .* exp(-1j*w*(0:length(den)-1)));
end
