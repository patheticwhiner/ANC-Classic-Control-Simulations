function result = run_demo1_test(signals, test_name, variant, params)
% RUN_DEMO1_TEST  Demo1: Youla 参数化 + Q-RLS 自适应（Phase 3 统一接口）
%
%   用法:
%     result = run_demo1_test(signals, test_name, variant, params)
%
%   输入: 见 tests/controller_interface.m 接口规范
%
%   params 专用字段:
%     .nQ        — Q 滤波器阶数，默认 2
%     .lambda1   — RLS λ₁ (遗忘因子)，默认 0.98
%     .lambda2   — RLS λ₂ (数据权重)，默认 1
%     .F_diag    — 初始协方差矩阵对角元，默认 1e-3
%
%   ⚠ 已知壁垒 (2026-06-13):
%     中央控制器 R₀=0, S₀=1 为开环占位。bezoutd 对 ARMAX(30,30,30,22)
%     @48kHz 存在结构性病态（Sylvester 条件数 ~1e16, S₀ 极点 >1.8）。
%     详见 tests/../demo1_RST/About_RST_Bezout_Barrier.md

    %% ====================================================================
    % §0. 参数解析与模型加载
    %% ====================================================================

    % 默认参数
    if nargin < 4, params = struct(); end
    if ~isfield(params, 'nQ'),      params.nQ      = 2; end
    if ~isfield(params, 'lambda1'), params.lambda1 = 0.98; end
    if ~isfield(params, 'lambda2'), params.lambda2 = 1; end
    if ~isfield(params, 'F_diag'),  params.F_diag  = 1e-3; end

    % 持久化：控制器设计只做一次
    persistent ctrl_cache;
    if isempty(ctrl_cache)
        ctrl_cache = struct('model', '', 'R0', [], 'S0', [], ...
            'HR', [], 'HS', [], 'P0', [], 'B_star', []);
    end

    % 路径设置
    addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'dataset'));
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));

    % 加载模型
    info = DataManager(signals.model);
    A = info.model.A(:).';          % 行向量
    B_full = info.model.B(:).';     % 含 d 步延迟前导零
    d = info.orders(4);
    fs = info.fs;
    Ts = 1/fs;

    B_coeffs = B_full(d+1:end);     % 去除延迟的分子系数

    %% ====================================================================
    % §1. 中央控制器设计
    %% ====================================================================
    %
    %  当前使用开环占位 R₀=0, S₀=1:
    %    - P₀ = A (被控对象开环特征多项式，极点 max|p|=0.9958)
    %    - 固定 variant: 抑制 ≈ 0 dB（无反馈，框架验证用）
    %    - 自适应 variant: Q-RLS 提供全部性能
    %
    %  bezoutd 极点配置不可行的原因:
    %    ARMAX(30,30,30,22) @48kHz → Sylvester 矩阵 (82×82) 条件数 ~1e16
    %    → S₀ 零点模 >1.8 → 控制器内部不稳定
    %    → 即使 P₀ Schur，仿真因 S₀ IIR 不稳定而发散
    %
    %  解决方案 (未来):
    %    1. 模型降阶 (平衡截断 → 10–15 阶)
    %    2. 连续域设计 + 离散化
    %    3. 直接 Q 参数化跳过固定控制器
    %    4. 正则化 Bezout 求解 (当前 bezoutd_reg 同样不稳定)

    B = [zeros(1, d), B_coeffs];     % 含延迟完整分子
    B_star = B(2:end);               % Landau 风格: 去首零，抗噪/观测多项式一致

    design_key = signals.model;
    if ~strcmp(ctrl_cache.model, design_key)
        fprintf('--- Demo1: 中央控制器 (开环占位) ---\n');

        R0 = 0;   HR = 1;
        S0 = 1;   HS = 1;
        P0 = A;   % P₀ = A·S₀ + B·R₀ = A

        ctrl_cache.model  = design_key;
        ctrl_cache.R0     = R0;
        ctrl_cache.S0     = S0;
        ctrl_cache.HR     = HR;
        ctrl_cache.HS     = HS;
        ctrl_cache.P0     = P0;
        ctrl_cache.B_star = B_star;

        fprintf('  R₀=0, S₀=1 (开环占位, P₀=A, max|root|=%.4f)\n', ...
            max(abs(roots(P0))));
    else
        R0     = ctrl_cache.R0;
        S0     = ctrl_cache.S0;
        HR     = ctrl_cache.HR;
        HS     = ctrl_cache.HS;
        P0     = ctrl_cache.P0;
        B_star = ctrl_cache.B_star;
    end

    %% ====================================================================
    % §2. 选择测试信号
    %% ====================================================================
    test_sig = signals.(test_name);
    d_sig    = test_sig.d(:)';      % 行向量
    Nsim     = length(d_sig);
    y_open   = test_sig.y_open(:)';

    %% ====================================================================
    % §3. 仿真
    %% ====================================================================

    switch variant
        case 'fixed'
            % ---- 仅固定 R₀/S₀ 控制器 ----
            [y_closed, u] = sim_fixed_rst(d_sig, A, B_full, R0, S0);
            Q_final = [];

        case 'adaptive'
            % ---- R₀/S₀ + Youla-Q RLS 自适应 (Landau 2017 结构) ----
            fprintf('  Q-RLS: nQ=%d, λ₁=%.3f, λ₂=%.1f, F₀=%.0e\n', ...
                params.nQ, params.lambda1, params.lambda2, params.F_diag);
            [y_closed, u, Q_final] = sim_adaptive_landau(...
                d_sig, A, B, B_star, R0, S0, HR, HS, P0, params);
            fprintf('  Q_final: [%s], ||Q||=%.3f\n', ...
                num2str(Q_final, '%.3f '), norm(Q_final));

        otherwise
            error('Demo1: 未知 variant ''%s''', variant);
    end

    %% ====================================================================
    % §4. 计算标准化指标
    %% ====================================================================
    meta = struct('demo', 'demo1', 'variant', variant, 'test', test_name);
    result = compute_metrics(y_open, y_closed, u, test_sig, meta);

    %% ====================================================================
    % §5. 追加 demo1 专用字段
    %% ====================================================================
    result.extra.deg_R0 = length(R0) - 1;
    result.extra.deg_S0 = length(S0) - 1;
    result.extra.R0_is_placeholder = isequal(R0, 0);

    if strcmp(variant, 'adaptive')
        result.extra.Q_final  = Q_final;
        result.extra.nQ       = params.nQ;
        result.extra.lambda1  = params.lambda1;
        result.extra.lambda2  = params.lambda2;
    end

    fprintf('  Demo1 %-9s | %s: %.1f dB | t_conv=%.2fs | max|u|=%.3f\n', ...
        variant, test_name, result.supp_db, result.t_conv_s, result.u_max);
end

%% ====================================================================
% 辅助函数
%% ====================================================================

function [y, u] = sim_fixed_rst(d, A, B, R, S)
% 固定 RST 控制器仿真
%   被控对象: y(k) = d(k) + (B/A)*u(k)
%   控制律:   S(z⁻¹)·u(k) = -R(z⁻¹)·y(k)
    N = length(d);
    y = zeros(1, N);
    u = zeros(1, N);

    nA = length(A);  nB = length(B);
    nR = length(R);  nS = length(S);
    buf_len = max([nA, nB, nR, nS]) + 1;

    u_buf    = zeros(1, buf_len);
    anti_buf = zeros(1, buf_len);
    y_buf    = zeros(1, buf_len);

    for k = 1:N
        anti_new = FIR(B, u_buf) - FIR(A(2:end), anti_buf);
        y_new = d(k) + anti_new;

        u_new = -(FIR(R, [y_new, y_buf(1:end-1)]) + ...
                  FIR(S(2:end), u_buf)) / S(1);

        anti_buf = [anti_new, anti_buf(1:end-1)];
        u_buf    = [u_new,    u_buf(1:end-1)];
        y_buf    = [y_new,    y_buf(1:end-1)];

        y(k) = y_new;
        u(k) = u_new;
    end
end

function [y, u, Q_final] = sim_adaptive_landau(...
    disturbance, A, B, B_star, Ro, So, HR, HS, Po, params)
% R₀/S₀ + Youla-Q RLS 自适应（移植自 Landau 2017 YKNANC_05 验证结构）
%
%   关键: 抗噪路径用 B(2:end)，扰动观测用 B_star = B(2:end) → 延迟一致
%         P₀ = A·S₀ + B·R₀ 为闭环特征多项式
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
        u_new = -FIR(Ro, Y) - FIR(HR_HS_Q, W) - FIR(So(2:end), U);
        U = [u_new, U(1:end-1)];

        % ---- 7. 回归向量: w2 = (H_R*H_S*B_star/P₀)*w ----
        w2_new = (FIR(HR_HS_Bstar, W) - FIR(Po(2:end), W2)) / Po(1);
        W2 = [w2_new, W2(1:end-1)];
        Phi = [w2_new, Phi(1:end-1)];

        % 记录
        y(k) = y_new;
        u(k) = u_new;
    end

    Q_final = Q_vec;
end

function y = FIR(b, x)
% FIR 点积滤波: y = Σ b(i)·x(i)
    L = min(length(b), length(x));
    y = b(1:L) * x(1:L)';
end
