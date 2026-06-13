function result = run_demo4_test(signals, test_name, variant, params)
% RUN_DEMO4_TEST  Demo4: RST 自动设计 (降阶 + 极点配置)
%
%   策略: balred 降阶 ARMAX 30→reduced_order → Diophantine 极点配置
%   → 全阶模型仿真验证。ε-MOPSO 优化框架延后至调参阶段。
%
%   用法:
%     result = run_demo4_test(signals, test_name, variant, params)
%
%   仅支持 'fixed' variant (RST 为固定控制器)。
%
%   params 专用字段:
%     .reduced_order  — 降阶目标阶数，默认 8
%     .f_design       — 陷波设计频率 (Hz)，默认 334
%     .zeta           — 陷波阻尼，默认 0.01

    if nargin < 4, params = struct(); end
    if ~isfield(params, 'reduced_order'), params.reduced_order = 8; end
    if ~isfield(params, 'f_design'),      params.f_design      = 334; end
    if ~isfield(params, 'zeta'),          params.zeta          = 0.01; end

    if ~strcmp(variant, 'fixed')
        error('Demo4: 仅支持 ''fixed'' variant。ε-MOPSO 优化延后。');
    end

    addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'dataset'));
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));

    %% ====================================================================
    % §1. 模型加载 + 降阶
    %% ====================================================================
    info = DataManager(signals.model);
    A_full = info.model.A;
    B_full = info.model.B;
    d = info.orders(4);
    fs = info.fs;
    Ts = 1/fs;

    %% ====================================================================
    % §2. RST 设计 (在降阶模型上)
    %% ====================================================================
    persistent ctrl_cache;
    if isempty(ctrl_cache)
        ctrl_cache = struct('design_key', '', 'R0', [], 'S0', []);
    end

    design_key = sprintf('%s_red%d_f%.0f_z%.2f', ...
        signals.model, params.reduced_order, params.f_design, params.zeta);

    if ~strcmp(ctrl_cache.design_key, design_key)
        fprintf('--- Demo4: RST 自动设计 (降阶 %d 阶) ---\n', params.reduced_order);

        % balred 降阶
        plant_full = tf(B_full, A_full, Ts);
        n_red = params.reduced_order;
        fprintf('  balred: %d → %d 阶...\n', length(A_full)-1, n_red);
        plant_red = balred(plant_full, n_red);
        [B_red_num, A_red_den] = tfdata(plant_red, 'v');
        A_red = A_red_den(:)';
        B_red = B_red_num(:)';

        % H_S: 陷波 @ f_design
        H_S = second_order_factor(params.f_design, params.zeta, Ts);
        H_R = 1;

        % P_F: 辅助滚降极点
        P_F = 1;
        for i = 1:min(4, n_red)
            P_F = conv(P_F, [1, -0.85]);
        end

        % P_D = A_red: 保留降阶模型的自然模态
        P_D = A_red;
        P_desired = conv(P_D, P_F);

        % 检查度数预算
        nah = length(A_red)-1 + length(H_S)-1;
        nbh = length(B_red)-1 + length(H_R)-1;
        if length(P_desired)-1 > nah + nbh - 1
            % 截断 P_desired
            P_desired = P_desired(1:(nah + nbh));
            fprintf('  P_desired 截断至 %d 阶\n', length(P_desired)-1);
        end

        % 求解 Diophantine (在降阶模型上)
        fprintf('  求解 Diophantine (nA_red=%d, nB_red=%d)...\n', ...
            length(A_red)-1, length(B_red)-1);
        try
            [R_prime, S_prime] = bezoutd(A_red, B_red, H_S, H_R, P_desired);
        catch
            warning('Demo4: Diophantine 失败 — 降阶至更简模型或使用占位 RST');
            R_prime = 0;  S_prime = 1;
        end

        R0 = conv(R_prime, H_R);
        S0 = conv(S_prime, H_S);

        % 降阶闭环检查
        P0_red = addPolynomials(conv(A_red, S0), conv(B_red, R0));
        if ~isschur(P0_red) && ~isequal(R0, 0)
            warning('Demo4: 降阶闭环非 Schur — 使用占位 R₀=0');
            R0 = 0;  S0 = 1;
        end

        ctrl_cache.design_key = design_key;
        ctrl_cache.R0 = R0;
        ctrl_cache.S0 = S0;

        fprintf('  完成: deg(R₀)=%d, deg(S₀)=%d\n', length(R0)-1, length(S0)-1);
    else
        R0 = ctrl_cache.R0;
        S0 = ctrl_cache.S0;
    end

    %% ====================================================================
    % §3. 仿真 (全阶模型)
    %% ====================================================================
    test_sig = signals.(test_name);
    d_sig    = test_sig.d(:)';
    y_open   = test_sig.y_open(:)';

    [y_closed, u] = sim_fixed_rst(d_sig, A_full, B_full, R0, S0);

    %% ====================================================================
    % §4. 指标
    %% ====================================================================
    meta = struct('demo', 'demo4', 'variant', 'fixed', 'test', test_name);
    result = compute_metrics(y_open, y_closed, u, test_sig, meta);

    result.extra.reduced_order = params.reduced_order;
    result.extra.deg_R0 = length(R0) - 1;
    result.extra.deg_S0 = length(S0) - 1;
    result.extra.R0_is_placeholder = isequal(R0, 0);

    fprintf('  Demo4 fixed     | %s: %.1f dB | max|u|=%.3f\n', ...
        test_name, result.supp_db, result.u_max);
end

%% ====================================================================
% 辅助函数 (复制自 demo1)
%% ====================================================================

function fac = second_order_factor(f_hz, zeta, Ts)
    if f_hz == 0, fac = 1; return; end
    wn = 2*pi*f_hz;
    r = exp(-zeta*wn*Ts);
    wd = wn*sqrt(max(0, 1 - zeta^2));
    c1 = -2*r*cos(wd*Ts);
    c2 = r^2;
    fac = [1, c1, c2];
end

function [y, u] = sim_fixed_rst(d, A, B, R, S)
    N = length(d);
    y = zeros(1, N);  u = zeros(1, N);
    nA = length(A);  nB = length(B);
    nR = length(R);  nS = length(S);
    buf_len = max([nA, nB, nR, nS]) + 1;

    u_buf   = zeros(1, buf_len);
    anti_buf = zeros(1, buf_len);
    y_buf    = zeros(1, buf_len);

    for k = 1:N
        anti_new = FIR(B, u_buf) - FIR(A(2:end), anti_buf);
        y_new = d(k) + anti_new;
        u_new = -(FIR(R, [y_new, y_buf(1:end-1)]) + ...
                  FIR(S(2:end), u_buf)) / max(abs(S(1)), 1e-10);
        if isnan(u_new) || abs(u_new) > 1e4, u_new = 0; end

        anti_buf = [anti_new, anti_buf(1:end-1)];
        u_buf    = [u_new,    u_buf(1:end-1)];
        y_buf    = [y_new,    y_buf(1:end-1)];

        y(k) = y_new;  u(k) = u_new;
    end
end

function y = FIR(b, x)
    L = min(length(b), length(x));
    y = b(1:L) * x(1:L)';
end
