function [R, S, T, Syp] = postprocess_RST(theta_opt, sys_info, varargin)
% POSTPROCESS_RST  从优化结果反推RST控制器并进行仿真验证
%
%   输入:
%     theta_opt  - 最优决策变量 (x_1,...,x_{nX}, y_1,...,y_{nY})
%     sys_info   - 系统信息结构体 (同RST_objective)
%     varargin   - 可选参数:
%       'Plot'      - true/false, 是否绘图 (默认true)
%       'SaveFig'   - true/false, 是否保存图片 (默认false)
%       'FigPrefix' - 图片文件名前缀 (默认'RST_')
%
%   输出:
%     R          - R 多项式系数 (降幂)
%     S          - S 多项式系数 (降幂)
%     T          - T 多项式系数 (标量)
%     Syp        - 输出灵敏度函数结构体 (含 .omega, .magnitude 字段)
%
%   工作流程:
%     1. 从 θ 构建 X 和 Y 多项式
%     2. 构建增广对象模型 A*H_S, B*H_R, P_D*X
%     3. 求解 Bezout 方程得到 S 和 R
%     4. 计算 T 多项式（跟踪增益）
%     5. 稳定性检查
%     6. 频域分析（灵敏度函数）
%     7. 频域绘图（含约束边界）
%     8. 时域仿真（阶跃响应）
%     9. 输出 RST 多项式摘要
%
%   参考文献:
%     Landau, I.D. and Zito, G., "Digital Control Systems: Design,
%     Identification and Implementation", Springer, 2006.

    % =====================================================================
    % 输入参数解析
    % =====================================================================
    p = inputParser;
    addRequired(p, 'theta_opt', @isnumeric);
    addRequired(p, 'sys_info', @isstruct);
    addParameter(p, 'Plot', true, @islogical);
    addParameter(p, 'SaveFig', false, @islogical);
    addParameter(p, 'FigPrefix', 'RST_', @ischar);
    parse(p, theta_opt, sys_info, varargin{:});

    plot_flag   = p.Results.Plot;
    save_flag   = p.Results.SaveFig;
    fig_prefix  = p.Results.FigPrefix;

    % 获取当前脚本所在目录，确保路径正确
    scriptDir = fileparts(mfilename('fullpath'));
    if isempty(scriptDir)
        scriptDir = pwd;
    end
    addpath(fullfile(scriptDir, '..', 'functions'));
    addpath(scriptDir);  % 确保 ./utils 中 bezoutd_reg 等文件可解析

    % 确保 theta_opt 为列向量
    theta_opt = theta_opt(:);

    % =====================================================================
    % 步骤1：从 θ 构建 X 和 Y 多项式
    % =====================================================================
    nX = sys_info.nX;
    nY = sys_info.nY;

    x_roots = theta_opt(1:nX);
    y_roots = theta_opt(nX+1:nX+nY);

    % 构建 X 和 Y 多项式 (降幂, 首项为1)
    X = real(poly(x_roots));   % X(z^{-1}) = 1 + x1*z^{-1} + ...
    Y = real(poly(y_roots));   % Y(z^{-1}) = 1 + y1*z^{-1} + ...

    % 确保多项式为非空
    if isempty(X), X = 1; end
    if isempty(Y), Y = 1; end

    fprintf('\n========== 多项式构建 ==========\n');
    fprintf('X(z^{-1}) 阶数: %d, 根均在单位圆内: %d\n', ...
        length(X)-1, all(abs(x_roots) < 1));
    fprintf('Y(z^{-1}) 阶数: %d, 根均在单位圆内: %d\n', ...
        length(Y)-1, all(abs(y_roots) < 1));

    % =====================================================================
    % 步驟2：构建增广对象模型
    % =====================================================================
    % RST 闭环特征多项式:
    %   P_c = A*H_S*S + z^{-d}*B*H_R*R = P_D*X
    %
    % 等价于:
    %   A*H_S*S + B_aug*H_R*R = P_D*X
    % 其中 B_aug = z^{-d}*B (已将延迟并入B)

    A_sys = sys_info.A(:).';
    B_sys = sys_info.B(:).';
    d_sys = sys_info.d;
    P_D   = sys_info.P_D(:).';
    H_R   = sys_info.H_R(:).';
    H_S   = sys_info.H_S(:).';

    % 将延迟并入分子: B_tilde = [zeros(1,d_sys), B_sys]
    B_tilde = [zeros(1, d_sys), B_sys];

    % A*H_S 与 B_tilde*H_R
    A_HS = conv(A_sys, H_S);
    B_HR = conv(B_tilde, H_R);

    % 期望闭环特征多项式: P_D * X
    P_DX = conv(P_D, X);

    fprintf('增广对象 A*H_S 阶数: %d\n', length(A_HS)-1);
    fprintf('增广对象 B*H_R 阶数: %d\n', length(B_HR)-1);
    fprintf('期望闭环多项式 P_D*X 阶数: %d\n', length(P_DX)-1);

    % =====================================================================
    % 步骤3：求解 Bezout 方程
    % =====================================================================
    % 使用 bezoutd_reg（改进版）: [Rp, Sp] = bezoutd_reg(A, B, Hs, Hr, P)
    %   - 原点填充替代 rmin=1e-16，消除大规模系统 Sylvester 矩阵病态
    %   - 内置 Tikhonov 正则化自动回退
    % 注意: 这里的 A, B 是原始对象多项式(不含延迟)

    try
        % 注意: bezoutd_reg 要求 B 包含延迟前导零，使用 B_tilde（已含 z^{-d}）
        [R_raw, S_raw] = bezoutd_reg(A_sys, B_tilde, H_S, H_R, P_DX);

        % 去除前导零 (trim from left)
        S = trimPolynomial(S_raw, "left");
        R = trimPolynomial(R_raw, "left");

        fprintf('\nBezout 方程求解成功 (bezoutd_reg).\n');
    catch ME
        fprintf('\nbezoutd_reg 求解失败: %s\n', ME.message);
        fprintf('尝试 Sylvester 矩阵备选方法...\n');

        [S, R] = solve_bezout_sylvester(A_sys, B_tilde, H_S, H_R, P_DX);

        if isempty(S) || isempty(R)
            error('postprocess_RST:BezoutFailed', ...
                  'Bezout 方程求解失败，Sylvester 备选方法也失败。');
        end
    end

    % 确保多项式为行向量
    S = S(:).';
    R = R(:).';

    fprintf('R 多项式阶数: %d\n', length(R)-1);
    fprintf('S 多项式阶数: %d\n', length(S)-1);

    % =====================================================================
    % 步骤4：计算 T 多项式（跟踪增益）
    % =====================================================================
    % T = P_D(1) * Y(1) / B(1)，用于零稳态误差跟踪阶跃输入
    % 注意: B(1) = sum(B_sys)，因为 z^{-d} 在 z=1 时 = 1

    PD1 = polyval(P_D, 1);       % P_D(z=1)
    Y1  = polyval(Y, 1);         % Y(z=1)
    B1  = polyval(B_sys, 1);     % B(z=1), 不含延迟

    if abs(B1) > 1e-10
        T = PD1 * Y1 / B1;
    else
        % B(1) ≈ 0, 对象含积分作用 → 使用 T = 0 (纯反馈)
        T = 0;
        warning('postprocess_RST:B1NearZero', ...
                'B(1) ≈ 0 (%.2e), 设置 T=0 (纯反馈控制)。', B1);
    end

    fprintf('\n跟踪增益 T = %.6f\n', T);

    % =====================================================================
    % 步骤5：稳定性检查
    % =====================================================================
    % 实际闭环特征多项式: P_cl = A*H_S*S + B_tilde*H_R*R
    % Pad to same length before adding (conv result lengths may differ)
    P_cl1 = conv(A_HS, S);
    P_cl2 = conv(B_HR, R);
    max_len = max(length(P_cl1), length(P_cl2));
    P_cl1(end+1:max_len) = 0;
    P_cl2(end+1:max_len) = 0;
    P_cl = P_cl1 + P_cl2;
    P_cl = trimPolynomial(P_cl, "left");

    cl_roots = roots(P_cl);
    max_root_mag = max(abs(cl_roots));
    is_stable = max_root_mag < 1;

    % 检查与 P_DX 的匹配度
    P_DX_padded = padPolynomial(P_DX, length(P_cl));
    P_cl_padded = padPolynomial(P_cl, length(P_DX));
    bezout_error = max(abs(P_DX_padded - P_cl_padded));

    fprintf('\n========== 稳定性验证 ==========\n');
    fprintf('闭环极点最大模值: %.6f\n', max_root_mag);
    if is_stable
        fprintf('✓ 闭环系统稳定\n');
    else
        fprintf('✗ 警告: 闭环系统不稳定! (max|root| = %.4f)\n', max_root_mag);
    end
    fprintf('Bezout 残差 ||P_DX - P_cl||_∞ = %.2e\n', bezout_error);

    % 检查 R 和 S 的根
    R_roots = roots(R);
    S_roots = roots(S);
    fprintf('R 多项式: 阶数=%d, 根均在单位圆内=%d\n', ...
        length(R)-1, all(abs(R_roots) < 1));
    fprintf('S 多项式: 阶数=%d, 根均在单位圆内=%d\n', ...
        length(S)-1, all(abs(S_roots) < 1));

    % =====================================================================
    % 步骤6：频域分析 — 灵敏度函数
    % =====================================================================
    nFreq = 2048;
    omega = linspace(0, pi, nFreq);

    Syp_mag = zeros(1, nFreq);
    T_mag   = zeros(1, nFreq);

    for k = 1:nFreq
        w = omega(k);
        z_inv = exp(-1j * w);

        % 多项式在 e^{-jw} 处的值
        A_val   = polyval(A_sys, z_inv);
        HS_val  = polyval(H_S, z_inv);
        S_val   = polyval(S, z_inv);
        B_val   = polyval(B_sys, z_inv);
        HR_val  = polyval(H_R, z_inv);
        R_val   = polyval(R, z_inv);
        Pcl_val = polyval(P_cl, z_inv);

        % 考虑延迟: B_delayed = z^{-d} * B
        B_delayed_val = z_inv^d_sys * B_val;

        % 输出灵敏度: S_yp = A*H_S*S / P_cl
        Syp_mag(k) = abs(A_val * HS_val * S_val / (Pcl_val + 1e-12));

        % 互补灵敏度: T_yp = z^{-d}*B*H_R*R / P_cl
        T_mag(k) = abs(B_delayed_val * HR_val * R_val / (Pcl_val + 1e-12));
    end

    % 灵敏度函数输出结构体
    Syp = struct('omega', omega, 'magnitude', Syp_mag);

    % =====================================================================
    % 步骤7：频域绘图
    % =====================================================================
    if plot_flag
        % 计算约束边界
        one_minus_exp = abs(1 - exp(-1j * omega));
        one_minus_exp = max(one_minus_exp, 1e-8);
        inv_diff = 1 ./ one_minus_exp;

        upper_bound = 1 + inv_diff;
        lower_bound = 1 - inv_diff;
        lower_bound(lower_bound < 0) = 0;

        % ---- 图1: 灵敏度函数 Bode 图 (含约束边界) ----
        figure('Name', 'Sensitivity Function');

        semilogx(omega(2:end), 20*log10(Syp_mag(2:end)), 'b-', 'LineWidth', 1.5);
        hold on;
        semilogx(omega(2:end), 20*log10(upper_bound(2:end)), 'r--', 'LineWidth', 1);
        semilogx(omega(2:end), 20*log10(lower_bound(2:end)), 'r--', 'LineWidth', 1);
        hold off;
        xlabel('Frequency (rad/sample)');
        ylabel('|S_{yp}| (dB)');
        title('Output Sensitivity Function with Constraints');
        legend('|S_{yp}|', 'Upper bound (Δτ ≥ T_s)', 'Lower bound (ΔM ≥ 0.5)', ...
            'Location', 'best');
        grid on;
        xlim([omega(2), pi]);

        if save_flag
            saveas(gcf, [fig_prefix 'sensitivity.png']);
        end

        % ---- 图2: 互补灵敏度函数 ----
        figure('Name', 'Complementary Sensitivity');

        semilogx(omega(2:end), 20*log10(T_mag(2:end)), 'b-', 'LineWidth', 1.5);
        xlabel('Frequency (rad/sample)');
        ylabel('|T_{yp}| (dB)');
        title('Complementary Sensitivity Function');
        grid on;
        xlim([omega(2), pi]);

        if save_flag
            saveas(gcf, [fig_prefix 'complementary.png']);
        end
    end

    % =====================================================================
    % 步骤8：时域仿真（阶跃响应）
    % =====================================================================
    if plot_flag
        N_sim = 200;
        Ts_sim = sys_info.Ts;

        % 参考信号 (延迟20步的阶跃)
        r = [zeros(20, 1); ones(N_sim - 20, 1)];

        % --- 使用 filter 进行差分方程仿真 ---
        % 闭环传递函数:
        %   y(k) = T * (z^{-d}*B*H_R*R / P_cl) * r(k)   (跟踪)
        %   y(k) = (A*H_S*S / P_cl) * d(k)               (对输出扰动的响应，此处不激励)
        %
        % 跟踪响应: y_r = filter(num_yr, den_yr, r)

        num_yr = T * conv(B_tilde, conv(H_R, R));
        den_yr = P_cl;

        % 归一化（使分母首项为1）
        num_yr = num_yr / den_yr(1);
        den_yr = den_yr / den_yr(1);

        y = filter(num_yr, den_yr, r);

        % 控制信号传递函数:
        %   u(k) = T * (A*H_S*S / P_cl) * r(k)
        num_ur = T * conv(A_HS, S);
        den_ur = P_cl;
        num_ur = num_ur / den_ur(1);
        den_ur = den_ur / den_ur(1);

        u = filter(num_ur, den_ur, r);

        % 计算性能指标
        steady_state_idx = round(N_sim * 0.7):N_sim;
        settling_error = abs(y(steady_state_idx) - 1);
        max_settling_err = max(settling_error);
        overshoot = (max(y) - 1) * 100;

        % 上升到90%的时间
        rise_idx = find(y >= 0.9, 1, 'first');
        if isempty(rise_idx), rise_idx = N_sim; end
        rise_time = (rise_idx - 20) * Ts_sim;  % 减去参考延迟

        fprintf('\n========== 时域性能指标 ==========\n');
        fprintf('最大超调: %.2f %%\n', overshoot);
        if overshoot < 0
            fprintf('  (无超调，单调收敛)\n');
        end
        fprintf('上升时间 (10%%→90%%): %.3f 采样周期\n', rise_idx - 20);
        fprintf('稳态误差 (70%%→100%%): max=%.4e\n', max_settling_err);

        % ---- 图3: 阶跃响应 ----
        t = (0:N_sim-1) * Ts_sim;

        figure('Name', 'Step Response');

        subplot(2,1,1);
        stairs(t, y, 'b-', 'LineWidth', 1.5);
        hold on;
        stairs(t, r, 'r--', 'LineWidth', 1);
        hold off;
        xlabel('Time (s)');
        ylabel('Output y');
        title(sprintf('Closed-loop Step Response (OS=%.1f%%, t_r=%.2f s)', ...
            overshoot, rise_time));
        legend('y(k)', 'r(k)', 'Location', 'best');
        grid on;
        ylim_padding = 0.1 * (max(y) - min(y)) + 0.05;
        ylim([min(y)-ylim_padding, max(y)+ylim_padding]);

        subplot(2,1,2);
        stairs(t, u, 'b-', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('Control u');
        title(sprintf('Control Signal (max|u| = %.2f)', max(abs(u))));
        grid on;

        if save_flag
            saveas(gcf, [fig_prefix 'step_response.png']);
        end

        % ---- 图4: 零极点图 ----
        figure('Name', 'Pole-Zero Map');

        % 开环零极点
        ol_zeros = roots(B_sys);
        ol_poles = roots(A_sys);
        % 闭环极点
        cl_poles = roots(P_cl);

        h_zplane = zplane([], []);  % 创建空白 zplane
        hold on;

        % 单位圆
        theta_circle = linspace(0, 2*pi, 200);
        plot(cos(theta_circle), sin(theta_circle), 'k-', 'LineWidth', 0.5);

        % 绘制开环极点 (x)
        plot(real(ol_poles), imag(ol_poles), 'rx', 'MarkerSize', 10, ...
            'LineWidth', 1.5, 'DisplayName', 'Open-loop poles');
        % 绘制开环零点 (o)
        plot(real(ol_zeros), imag(ol_zeros), 'ro', 'MarkerSize', 10, ...
            'LineWidth', 1.5, 'DisplayName', 'Open-loop zeros');
        % 绘制闭环极点 (filled circle)
        plot(real(cl_poles), imag(cl_poles), 'b.', 'MarkerSize', 20, ...
            'DisplayName', 'Closed-loop poles');

        hold off;
        title('Pole-Zero Map');
        legend('Location', 'best');
        grid on;
        axis equal;
        axis([-1.5, 1.5, -1.5, 1.5]);

        if save_flag
            saveas(gcf, [fig_prefix 'pzmap.png']);
        end
    end

    % =====================================================================
    % 步骤9：输出 RST 多项式摘要
    % =====================================================================
    fprintf('\n========== RST 控制器多项式 ==========\n');

    % R 多项式
    fprintf('R(z^{-1}) = ');
    print_polynomial(R);

    % S 多项式
    fprintf('S(z^{-1}) = ');
    print_polynomial(S);

    % T
    fprintf('T(z^{-1}) = %.6f\n', T);

    fprintf('\n========== RST 控制律 ==========\n');
    fprintf('S(z^{-1}) * u(k) = T * r(k) - R(z^{-1}) * y(k)\n');
    if length(S) >= 2
        fprintf('即: u(k) = [T*r(k) - R(z^{-1})*y(k)');
        fprintf(' - (%.4f*z^{-1}+...+%.4f*z^{-%d})*u(k)] / %.4f\n', ...
            S(2), S(end), length(S)-1, S(1));
    else
        fprintf('即: u(k) = [T*r(k) - R(z^{-1})*y(k)] / %.4f\n', S(1));
    end

end

% =========================================================================
% 辅助函数: 打印多项式（格式化输出）
% =========================================================================
function print_polynomial(poly)
% PRINT_POLYNOMIAL  格式化打印多项式
%   输入 poly 为降幂排列的系数向量

    poly = poly(:).';
    first = true;
    for i = 1:length(poly)
        coeff = poly(i);
        if i == 1
            fprintf('%.4f', coeff);
            first = false;
        elseif coeff >= 0
            fprintf(' + %.4f z^{-%d}', coeff, i-1);
        else
            fprintf(' - %.4f z^{-%d}', abs(coeff), i-1);
        end
    end
    fprintf('\n');
end

% =========================================================================
% 辅助函数: Sylvester 矩阵法求解 Bezout 方程（备选方案）
% =========================================================================
function [S, R] = solve_bezout_sylvester(A, B, H_S, H_R, P_target)
% SOLVE_BEZOUT_SYLVESTER  使用 Sylvester 矩阵求解 Bezout 方程（备选方案）
%   求解: A*H_S*S + B*H_R*R = P_target
%
%   输入:
%     A, B     - 原始对象多项式系数 (降幂)
%     H_S, H_R - 控制器固定因子 (降幂)
%     P_target - 目标闭环多项式系数 (降幂)
%
%   输出:
%     S    - S 多项式系数 (降幂)
%     R    - R 多项式系数 (降幂)

    % 确保为行向量
    A = A(:).';
    B = B(:).';
    H_S = H_S(:).';
    H_R = H_R(:).';
    P_target = P_target(:).';

    % 计算增广多项式
    A_HS = conv(A, H_S);
    B_HR = conv(B, H_R);

    n_AHS = length(A_HS);
    n_BHR = length(B_HR);
    n_P = length(P_target);

    % 标准解: deg(S) = deg(B_HR) - 1, deg(R) = deg(A_HS) - 1
    nS = n_BHR - 1;
    nR = n_AHS - 1;

    % 如果 P_target 阶数更大，均匀分配额外阶数给 R 和 S
    expected_deg = n_AHS + n_BHR - 2;  % deg(A_HS*S) = n_AHS-1 + nS
    if n_P - 1 > expected_deg
        extra = (n_P - 1) - expected_deg;
        nR = nR + ceil(extra / 2);
        nS = nS + floor(extra / 2);
    end

    % 确保非负
    nS = max(0, nS);
    nR = max(0, nR);

    % 如果 P_target 阶数不足，补零
    if n_P < max(n_AHS + nS, n_BHR + nR)
        P_target = [P_target, zeros(1, max(n_AHS + nS, n_BHR + nR) - n_P)];
        n_P = length(P_target);
    end

    % 构建 Sylvester 矩阵: M * [S_coeff; R_coeff] = P_target'
    % M 尺寸: (n_P) × (nS + nR)
    M = zeros(n_P, nS + nR);

    % A_HS 的卷积列
    for col = 1:nS
        M(col:col + n_AHS - 1, col) = A_HS';
    end

    % B_HR 的卷积列
    for col = 1:nR
        row_start = col;
        M(row_start:row_start + n_BHR - 1, nS + col) = B_HR';
    end

    % 求解
    coeff = M \ P_target(:);

    S = coeff(1:nS)';
    R = coeff(nS+1:end)';

    % 去除前导零
    S = trimPolynomial(S, "left");
    R = trimPolynomial(R, "left");

    fprintf('  Sylvester 解: S 阶数=%d, R 阶数=%d\n', length(S)-1, length(R)-1);
end
