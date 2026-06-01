function f = RST_objective(theta, sys_info)
% RST_OBJECTIVE 基于Zames-Francis积分的RST控制器优化目标函数
%   计算RST多项式控制器参数的优化目标值，基于频域加权灵敏度积分。
%
%   输入:
%     theta    - 决策变量向量 [x_1, ..., x_{nX}, y_1, ..., y_{nY}]
%                其中 x_i 为 X 多项式的根, y_i 为 Y 多项式的根
%     sys_info - 系统信息结构体，包含以下字段:
%       .A        - 分母多项式系数 (降幂: A = [a0, a1, ..., a_nA])
%       .B        - 分子多项式系数 (降幂: B = [b0, b1, ..., b_nB])
%       .d        - 纯延迟步数
%       .P_D      - 主导极点多项式系数 (降幂)
%       .H_R      - 控制器分子固定因子 (如 [1, 1], 降幂)
%       .H_S      - 控制器分母固定因子 (如 [1, -1], 降幂)
%       .nX       - X多项式阶数
%       .nY       - Y多项式阶数
%       .beta     - 不稳定零点 (复数列向量, |beta| > 1)
%       .alpha    - 不稳定极点 (复数列向量, |alpha| > 1)
%       .Ts       - 采样时间 (可选, 默认: 1)
%       .nFreq    - 频率离散点数 (可选, 默认: 1000)
%
%   输出:
%     f         - 目标函数值向量 [f_β1; f_β2; ...; f_delay]
%                 其中 f_βi 对应于每个不稳定零点的Poisson积分匹配误差
%                 f_delay 为延迟积分匹配误差
%
%   参考文献:
%     Zames, G. and Francis, B.A., "Feedback, Minimax Sensitivity, and
%     Optimal Robustness", IEEE Trans. AC, 28(5), pp. 585-601, 1983.

    % =====================================================================
    % 输入参数验证
    % =====================================================================
    if nargin < 2
        error('RST_objective:NotEnoughInputs', ...
              '需要2个输入参数: theta 和 sys_info。');
    end

    % 确保 theta 为列向量
    theta = theta(:);

    % 提取系统信息
    required_fields = {'A', 'B', 'd', 'P_D', 'H_R', 'H_S', 'nX', 'nY', 'beta', 'alpha'};
    for i = 1:length(required_fields)
        if ~isfield(sys_info, required_fields{i})
            error('RST_objective:MissingField', ...
                  'sys_info 缺少必要字段: ''%s''。', required_fields{i});
        end
    end

    % 设置默认值
    if ~isfield(sys_info, 'Ts'),    sys_info.Ts    = 1;    end
    if ~isfield(sys_info, 'nFreq'), sys_info.nFreq = 1000; end

    % 提取参数（便于引用）
    A     = sys_info.A(:).';
    B     = sys_info.B(:).';
    d     = sys_info.d;
    P_D   = sys_info.P_D(:).';
    H_R   = sys_info.H_R(:).';
    H_S   = sys_info.H_S(:).';
    nX    = sys_info.nX;
    nY    = sys_info.nY;
    beta  = sys_info.beta(:);
    alpha = sys_info.alpha(:);
    nFreq = sys_info.nFreq;

    % 验证 theta 维度
    if length(theta) ~= nX + nY
        error('RST_objective:DimensionMismatch', ...
              'theta 长度 (%d) 与 nX+nY (%d) 不匹配。', ...
              length(theta), nX + nY);
    end

    % =====================================================================
    % 步骤1: 拆分决策变量并构建多项式
    % =====================================================================
    % 提取 X 和 Y 多项式的根
    x_roots = theta(1:nX);
    y_roots = theta(nX+1:nX+nY);

    % 使用 poly() 从根构建多项式系数 (降幂排列)
    X_coeff = poly(x_roots);     % 1× (nX+1)
    Y_coeff = poly(y_roots);     % 1× (nY+1)

    % =====================================================================
    % 步骤2: 生成频率向量
    % =====================================================================
    omega = linspace(-pi, pi, nFreq);

    % =====================================================================
    % 步骤3: 计算频域响应 |S_d(e^{jω})|
    % =====================================================================
    [S_d_mag, log_S_d_mag] = compute_Sd_magnitude(A, B, d, P_D, H_R, H_S, ...
                                                   X_coeff, Y_coeff, omega);

    % =====================================================================
    % 步骤4: 对每个不稳定零点计算目标值
    % =====================================================================
    n_beta = length(beta);
    f_beta = zeros(n_beta, 1);

    for i = 1:n_beta
        beta_i = beta(i);

        % --- 4.1 计算 Blaschke 乘积 B_a(beta_i) ---
        B_a_val = blaschke_product(beta_i, alpha);

        % 理论上 |B_a(beta_i)| <= 1 (对于不稳定极点alpha满足|alpha|>1)
        % 避免对数中出现零或负值
        B_a_val = max(abs(B_a_val), 1e-12);

        % --- 4.2 计算 Poisson 核 ---
        K = poisson_kernel(omega, beta_i);

        % --- 4.3 数值积分 ---
        integral_val = trapz(omega, log_S_d_mag .* K);

        % --- 4.4 目标值: 积分与理论值的绝对偏差 ---
        % Zames-Francis 等式: ∫ log|S_d| * K dω = π * log(|B_a^{-1}(beta)|)
        % 即: ∫ log|S_d| * K dω = -π * log(|B_a(beta)|)
        target_val = -pi * log(B_a_val);

        f_beta(i) = abs(integral_val - target_val);
    end

    % =====================================================================
    % 步骤5: 计算延迟目标值
    % =====================================================================
    % Zames-Francis 延迟积分:
    % ∫ log|S_d(e^{jω})| dω = π * Σ log(|alpha_k|)

    integral_delay = trapz(omega, log_S_d_mag);

    % 理论值
    alpha_mag = abs(alpha);
    alpha_mag = max(alpha_mag, 1e-12);  % 防止 log(0)
    target_delay = pi * sum(log(alpha_mag));

    f_delay = abs(integral_delay - target_delay);

    % =====================================================================
    % 步骤6: 组装输出向量
    % =====================================================================
    f = [f_beta; f_delay];

    % =====================================================================
    % 步骤7: Bezout 残差目标 — 确保 X,Y 可反解为合法 R,S
    % =====================================================================
    % 构建增广对象
    A_HS = conv(sys_info.A, sys_info.H_S);
    B_HR = conv(sys_info.B, sys_info.H_R);
    P_DX = conv(sys_info.P_D, X_coeff);

    % 尝试求解 Bezout 方程
    try
        % 包含延迟的 B 多项式
        B_delayed = [zeros(1, sys_info.d), sys_info.B];

        [Rp, Sp] = bezoutd(sys_info.A, B_delayed, sys_info.H_S, sys_info.H_R, P_DX);
        % 实际闭环多项式
        P_cl = conv(A_HS, Sp) + conv(B_HR, Rp);
        % 对齐长度
        max_len = max(length(P_DX), length(P_cl));
        P_DX_pad = [P_DX, zeros(1, max_len - length(P_DX))];
        P_cl_pad = [P_cl, zeros(1, max_len - length(P_cl))];
        % Bezout 残差（归一化）
        f_bezout = norm(P_DX_pad - P_cl_pad) / max(abs(P_DX) + 1e-10);
    catch
        f_bezout = 100;  % 求解失败则赋大值
    end

    % 追加到目标向量
    f = [f; f_bezout];

end

% =========================================================================
% 辅助函数: 计算 |S_d(e^{jω})| 和 log|S_d|
% =========================================================================
function [mag, log_mag] = compute_Sd_magnitude(A, B, d, P_D, H_R, H_S, ...
                                               X_coeff, Y_coeff, omega)
% COMPUTE_SD_MAGNITUDE 计算灵敏度函数 S_d 在频域的幅值
%
%   S_d(e^{jω}) = (A(e^{jω}) / P_D(e^{jω})) * (X(e^{jω}) / Y(e^{jω}))
%
%   输出:
%     mag     - |S_d(e^{jω})| 行向量
%     log_mag - log|S_d(e^{jω})| 行向量

    nFreq = length(omega);
    mag = zeros(1, nFreq);

    % 计算 A 和 P_D 的阶数
    nA = length(A) - 1;
    nPD = length(P_D) - 1;
    nX = length(X_coeff) - 1;
    nY = length(Y_coeff) - 1;

    for k = 1:nFreq
        w = omega(k);

        % 计算 A(e^{jω}) = Σ A(m+1) * exp(-jω*m), m=0,...,nA
        A_val = sum(A .* exp(-1j * w * (0:nA)));

        % 计算 B(e^{jω})
        % B_val = sum(B .* exp(-1j * w * (0:length(B)-1)));

        % 计算 P_D(e^{jω})
        P_D_val = sum(P_D .* exp(-1j * w * (0:nPD)));

        % 计算 X(e^{jω})
        X_val = sum(X_coeff .* exp(-1j * w * (0:nX)));

        % 计算 Y(e^{jω})
        Y_val = sum(Y_coeff .* exp(-1j * w * (0:nY)));

        % 计算 S_d(e^{jω})
        % S_d = (A/P_D) * (X/Y)
        % 注意: 这里假设 H_R 和 H_S 已包含在 X 和 Y 中或另行处理
        % 完整的灵敏度函数需要包含 H_R, H_S 因子，
        % 此处按用户规范计算核心部分

        if abs(Y_val) < 1e-12
            Y_val = 1e-12;  % 避免除零
        end

        S_d_val = (A_val / max(abs(P_D_val), 1e-12)) * (X_val / Y_val);

        mag(k) = abs(S_d_val);
    end

    % 计算对数幅值
    log_mag = log(max(mag, 1e-12));
end

% =========================================================================
% 辅助函数: 计算 Blaschke 乘积
% =========================================================================
function B = blaschke_product(z, poles)
% BLASCHKE_PRODUCT 计算不稳定极点的 Blaschke 乘积在 z 处的值
%   B_a(z) = Π (z - alpha_k) / (1 - conj(alpha_k) * z)
%
%   输入:
%     z     - 评估点 (复数标量)
%     poles - 不稳定极点向量 (|alpha_k| > 1)
%   输出:
%     B     - Blaschke 乘积值

    B = 1;
    for k = 1:length(poles)
        alpha_k = poles(k);
        num = z - alpha_k;
        den = 1 - conj(alpha_k) * z;
        B = B * (num / den);
    end
end

% =========================================================================
% 辅助函数: 计算 Poisson 核
% =========================================================================
function K = poisson_kernel(omega, beta)
% POISSON_KERNEL 计算不稳定零点的 Poisson 核
%   K(ω) = (|β|² - 1) / (1 - 2|β|cos(ω - ∠β) + |β|²)
%
%   输入:
%     omega - 频率向量 (弧度)
%     beta  - 不稳定零点 (复数, |beta| > 1)
%   输出:
%     K     - Poisson 核值 (与 omega 同长度行向量)

    r = abs(beta);
    phi = angle(beta);

    % Poisson 核公式
    K = (r^2 - 1) ./ (1 - 2 * r * cos(omega - phi) + r^2);

    % 确保为行向量
    K = K(:)';
end
