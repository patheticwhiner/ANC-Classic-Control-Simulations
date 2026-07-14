function result = run_demo2_test(signals, test_name, variant, params)
% RUN_DEMO2_TEST  Demo2: LQG (增广/裸) + Youla-Q 自适应
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
%     .Q_lqr, .R_lqr — LQR 权重 (默认 1, 1)
%     .Qn_plant  — 对象状态过程噪声 (默认 1e-4)
%     .Qn_dist   — 扰动状态过程噪声 (默认 1)
%     .Rn        — 测量噪声 (默认 0.01)
%     .nQ, .lambda, .F_init, .bp_band — Q-RLS 参数

    if nargin < 4, params = struct(); end
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

    projectRoot = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(projectRoot, 'project_init.m'));

    %% ====================================================================
    % §1. 模型加载
    %% ====================================================================
    info = DataManager(signals.model);
    fs = info.fs;
    Ts = 1/fs;

    sys_idss = idss(info.model);
    Af = sys_idss.A;
    Bf = sys_idss.B;
    Cf = sys_idss.C;
    Kf = sys_idss.K;
    n_plant = size(Af, 1);

    %% ====================================================================
    % §2. 构建增广扰动模型
    %% ====================================================================
    test_sig = signals.(test_name);

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
        % 裸 LQG: 使用辨识的 innovation Kalman 增益 Kf (已验证稳定)
        A_design = Af;  B_design = Bf;  C_design = Cf;
        G_design = eye(n_plant);
        n_total   = n_plant;  n_dist = 0;

        Qn = params.Qn_plant * eye(n_plant);
        Rn = params.Rn * eye(size(Cf, 1));
        % 标记: 跳过 dlqe, 直接使用 Kf
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

    cache_key = sprintf('%s_%s_aug%d_Q%.0e_R%.0e_Qp%.0e_Qd%.0e_Rn%.0e', ...
        signals.model, test_name, params.augmented, ...
        params.Q_lqr, params.R_lqr, params.Qn_plant, params.Qn_dist, params.Rn);

    if ~strcmp(lqg_cache.key, cache_key)
        fprintf('--- Demo2: LQG (n_plant=%d, n_dist=%d) ---\n', n_plant, n_dist);

        % LQR
        Q_lqr = params.Q_lqr * (C_design' * C_design);
        R_lqr = params.R_lqr;
        [K, ~, ~] = dlqr(A_design, B_design, Q_lqr, R_lqr);

        % Kalman: 增广用 dlqe; 裸用辨识 Kf
        if params.augmented
            [L, ~, ~] = dlqe(A_design, G_design, C_design, Qn, Rn);
        else
            L = Kf;  % 辨识 innovation 增益 (已验证稳定, max|eig(A_cl)|=0.995)
        end

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

    % 使用完整对象动力学 (仅植物部分 + 扰动叠加)
    switch variant
        case 'fixed'
            [y_closed, u] = sim_lqg_augmented(...
                Af, Bf, Cf, n_plant, n_dist, ...
                A_design, B_design, C_design, K, L, A_cl, d_sig);

        case 'adaptive'
            [y_closed, u, theta_final] = sim_lqg_augmented_adaptive(...
                Af, Bf, Cf, n_plant, n_dist, ...
                A_design, B_design, C_design, K, L, A_cl, d_sig, fs, params);

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
    result.extra.max_eig_Acl = max(abs(eig(A_cl)));

    if strcmp(variant, 'adaptive')
        result.extra.nQ         = params.nQ;
        result.extra.lambda     = params.lambda;
        result.extra.norm_theta = norm(theta_final);
    end

    fprintf('  Demo2 %-9s | %s: %.1f dB | t_conv=%.2fs | max|u|=%.3f\n', ...
        variant, test_name, result.supp_db, result.t_conv_s, result.u_max);
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
function [y, u] = sim_lqg_augmented(...
        Af_plant, Bf_plant, Cf_plant, n_plant, n_dist, ...
        A_design, B_design, C_design, K, L, A_cl, d)
% 仿真: 真实扰动 d(k) 直接加到输出 (观测器内部用增广模型估计)
    N = length(d);
    x_plant = zeros(n_plant, 1);
    x_hat   = zeros(n_plant + n_dist, 1);
    y = zeros(1, N);
    u = zeros(1, N);

    for k = 1:N
        % 测量: y = plant_output + REAL disturbance
        y_k = Cf_plant * x_plant + d(k);

        % 控制律 (作用于增广状态估计)
        u_k = -K * x_hat;

        % Kalman 更新 (增广观测器学习分离扰动与对象状态)
        x_hat = A_cl * x_hat + L * y_k;

        % 对象状态更新
        x_plant = Af_plant * x_plant + Bf_plant * u_k;

        y(k) = y_k;
        u(k) = u_k;
    end
end

function [y, u, theta_final] = sim_lqg_augmented_adaptive(...
        Af_plant, Bf_plant, Cf_plant, n_plant, n_dist, ...
        A_design, B_design, C_design, K, L, A_cl, d, fs, params)
% 增广 LQG + Youla-Q RLS (真实扰动 d(k) 直接注入输出)

    n_total = n_plant + n_dist;
    N = length(d);
    x_plant = zeros(n_plant, 1);
    x_hat   = zeros(n_total, 1);
    y = zeros(1, N);
    u = zeros(1, N);

    bp = params.bp_band;
    Wn = bp / (fs/2);
    Wn = max(1e-3, min(Wn, 0.999));
    [b_h, a_h] = butter(4, Wn, 'bandpass');
    buf_in  = zeros(1, length(b_h));
    buf_out = zeros(1, length(a_h)-1);

    nQ = params.nQ;
    lambda = params.lambda;
    theta = zeros(nQ, 1);
    P_mat = params.F_init * eye(nQ);
    r_filt_hist = zeros(1, nQ + 1);
    warmup_samples = round(0.5 * fs);

    for k = 1:N
        y_k = Cf_plant * x_plant + d(k);

        r_k = y_k - C_design * x_hat;

        buf_in = [r_k, buf_in(1:end-1)];
        r_filt = b_h * buf_in' - a_h(2:end) * buf_out';
        buf_out = [r_filt, buf_out(1:end-1)];

        r_filt_hist = [r_filt, r_filt_hist(1:end-1)];
        phi = r_filt_hist(1:nQ)';
        y_Q_k = theta' * phi;
        if k < warmup_samples, y_Q_k = 0; end
        y_Q_k = max(-5, min(5, y_Q_k));

        u_k = -K * x_hat + y_Q_k;

        if k >= warmup_samples
            denom = 1 + phi' * P_mat * phi;
            e_post = r_filt / max(denom, 1e-10);
            theta = theta + P_mat * phi * e_post;
            P_mat = (P_mat - (P_mat * (phi * phi') * P_mat) / max(denom, 1e-10)) / lambda;
            if norm(theta) > 10, theta = theta / norm(theta) * 10; end
        end
        P_mat = P_mat + 1e-8 * eye(nQ);
        if any(isnan(theta)) || any(isinf(theta))
            theta = zeros(nQ, 1);  P_mat = params.F_init * eye(nQ);
        end

        x_hat = A_cl * x_hat + L * y_k + B_design * y_Q_k;
        x_plant = Af_plant * x_plant + Bf_plant * u_k;

        y(k) = y_k;  u(k) = u_k;
    end
    theta_final = theta;
end
