function [g1, g2, max_g1, max_g2] = RST_constraints(theta, sys_info)
% RST_CONSTRAINTS 频域模值裕度和延迟裕度约束计算
%   计算RST控制器的频域鲁棒性约束值，用于约束优化。
%
%   输入:
%     theta    - 决策变量向量 [x_1, ..., x_{nX}, y_1, ..., y_{nY}]
%     sys_info - 系统信息结构体 (同 RST_objective)
%       必要字段: .A, .B, .d, .P_D, .H_R, .H_S, .nX, .nY
%       可选字段: .nFreq (默认: 200), .verbose (默认: false)
%
%   输出:
%     g1      - 下界约束值向量 (nFreq × 1):
%               g1(ω) = 1 - 1/|1 - e^{-jω}| - |S_d(e^{jω})| ≤ 0
%     g2      - 上界约束值向量 (nFreq × 1):
%               g2(ω) = |S_d(e^{jω})| - 1 - 1/|1 - e^{-jω}| ≤ 0
%     max_g1  - max(g1)，用于快速约束违例检查 (>0表示违例)
%     max_g2  - max(g2)，用于快速约束违例检查 (>0表示违例)
%
%   约束物理意义:
%     - g1 ≤ 0: 保证足够的模值裕度（下界）
%     - g2 ≤ 0: 保证足够的延迟裕度（上界）
%     - 两者同时满足即保证系统的鲁棒稳定性
%
%   参考文献:
%     Landau, I.D. and Zito, G., "Digital Control Systems: Design,
%     Identification and Implementation", Springer, 2006.

    % =====================================================================
    % 输入参数验证
    % =====================================================================
    if nargin < 2
        error('RST_constraints:NotEnoughInputs', ...
              '需要2个输入参数: theta 和 sys_info。');
    end

    % 确保 theta 为列向量
    theta = theta(:);

    % 检查必要字段
    required_fields = {'A', 'B', 'd', 'P_D', 'H_R', 'H_S', 'nX', 'nY'};
    for i = 1:length(required_fields)
        if ~isfield(sys_info, required_fields{i})
            error('RST_constraints:MissingField', ...
                  'sys_info 缺少必要字段: ''%s''。', required_fields{i});
        end
    end

    % 设置默认频率点数和详细输出标志
    if ~isfield(sys_info, 'nFreq')
        sys_info.nFreq = 200;
    end
    if ~isfield(sys_info, 'verbose')
        sys_info.verbose = false;
    end

    % 提取参数
    A     = sys_info.A(:).';
    B     = sys_info.B(:).';
    d     = sys_info.d;
    P_D   = sys_info.P_D(:).';
    H_R   = sys_info.H_R(:).';
    H_S   = sys_info.H_S(:).';
    nX    = sys_info.nX;
    nY    = sys_info.nY;
    nFreq = sys_info.nFreq;

    % 验证 theta 维度
    if length(theta) ~= nX + nY
        error('RST_constraints:DimensionMismatch', ...
              'theta 长度 (%d) 与 nX+nY (%d) 不匹配。', ...
              length(theta), nX + nY);
    end

    % =====================================================================
    % 步骤1: 拆分决策变量并构建多项式
    % =====================================================================
    x_roots = theta(1:nX);
    y_roots = theta(nX+1:nX+nY);

    X_coeff = poly(x_roots);
    Y_coeff = poly(y_roots);

    % =====================================================================
    % 步骤2: 生成频率向量（对数间隔在低频加密，提高鲁棒性约束精度）
    % =====================================================================
    half = max(2, floor(nFreq / 2));
    omega_low = logspace(-4, log10(pi), half);
    omega_high = linspace(omega_low(end) + (pi - omega_low(end)) / (nFreq - half), pi, nFreq - half);
    omega = unique([omega_low, omega_high]);
    nFreq_actual = length(omega);

    % =====================================================================
    % 步骤3: 计算 |S_d(e^{jω})|
    % =====================================================================
    S_d_mag = compute_Sd_mag_constraint(A, B, d, P_D, H_R, H_S, ...
                                         X_coeff, Y_coeff, omega, nFreq_actual);

    % =====================================================================
    % 步骤4: 计算约束边界 |1 - e^{-jω}|
    % =====================================================================
    % |1 - e^{-jω}| = |1 - (cos ω - j sin ω)|
    %                = |(1 - cos ω) + j sin ω|
    %                = sqrt((1 - cos ω)² + sin² ω)
    %                = sqrt(1 - 2cos ω + cos² ω + sin² ω)
    %                = sqrt(2 - 2cos ω)
    %                = 2 * |sin(ω/2)|

    one_minus_exp = abs(1 - exp(-1j * omega));
    % 避免除零: 当 ω → 0 时, |1 - e^{-jω}| → 0
    one_minus_exp = max(one_minus_exp, 1e-8);

    inv_one_minus_exp = 1 ./ one_minus_exp;

    % =====================================================================
    % 步骤5: 计算约束值
    % =====================================================================
    % g1(ω) = 1 - 1/|1 - e^{-jω}| - |S_d(e^{jω})| ≤ 0
    g1 = 1 - inv_one_minus_exp - S_d_mag;

    % g2(ω) = |S_d(e^{jω})| - 1 - 1/|1 - e^{-jω}| ≤ 0
    g2 = S_d_mag - 1 - inv_one_minus_exp;

    % 确保输出为列向量
    g1 = g1(:);
    g2 = g2(:);

    % =====================================================================
    % 步骤5.5: 低频段约束松弛 — 考虑水床效应 (Bode积分约束)
    % =====================================================================
    omega_col = omega(:);  % 确保为列向量

    % 低频松弛因子：ω < 0.1π 时逐渐放松约束
    relax = ones(size(omega_col));
    low_freq_idx = omega_col < 0.1 * pi;
    if any(low_freq_idx)
        % 在 ω=0 处完全放松(×0.5)，在 ω=0.1π 处回到正常约束(×1.0)
        relax(low_freq_idx) = 1.0 - 0.5 * (1 - omega_col(low_freq_idx) / (0.1*pi));
    end

    % 对 g2 应用松弛
    g2 = g2 .* relax;

    % g1 低频也适当放松（但程度较轻）
    relax_g1 = ones(size(omega_col));
    if any(low_freq_idx)
        relax_g1(low_freq_idx) = 1.0 - 0.2 * (1 - omega_col(low_freq_idx) / (0.1*pi));
    end
    g1 = g1 .* relax_g1;

    % =====================================================================
    % 步骤6: 计算最大值（用于快速违例检查）
    % =====================================================================
    max_g1 = max(g1);
    max_g2 = max(g2);

    % =====================================================================
    % 步骤7: 约束诊断输出（仅在 verbose 模式或违例严重时）
    % =====================================================================
    if sys_info.verbose || max_g1 > 0 || max_g2 > 0
        n_viol_g1 = sum(g1 > 0);
        n_viol_g2 = sum(g2 > 0);
        if sys_info.verbose
            fprintf('  [约束] max_g1=%.4f (违例点数:%d/%d), max_g2=%.4f (违例点数:%d/%d)\n', ...
                max_g1, n_viol_g1, nFreq_actual, max_g2, n_viol_g2, nFreq_actual);
        end
    end

end

% =========================================================================
% 辅助函数: 计算 |S_d(e^{jω})| （约束版本，使用 [0, π] 频率范围）
% =========================================================================
function mag = compute_Sd_mag_constraint(A, B, d, P_D, H_R, H_S, ...
                                          X_coeff, Y_coeff, omega, nFreq)
% COMPUTE_SD_MAG_CONSTRAINT 计算灵敏度函数幅值
%   S_d(e^{jω}) = (A(e^{jω}) / P_D(e^{jω})) * (X(e^{jω}) / Y(e^{jω}))
%
%   输出:
%     mag - |S_d(e^{jω})| 行向量
%
%   注意: 此处不直接使用 d，S_d 表达式中不含延迟因子 z^{-d}

    if nargin < 9, nFreq = length(omega); end
    mag = zeros(1, nFreq);

    % |S_d| 饱和上限，防止 Y(z) 根靠近单位圆时幅度爆炸
    S_D_MAG_MAX = 100;

    nA = length(A) - 1;
    nPD = length(P_D) - 1;
    nX = length(X_coeff) - 1;
    nY = length(Y_coeff) - 1;

    for k = 1:nFreq
        w = omega(k);

        % A(e^{jω})
        A_val = sum(A .* exp(-1j * w * (0:nA)));

        % P_D(e^{jω})
        P_D_val = sum(P_D .* exp(-1j * w * (0:nPD)));

        % X(e^{jω})
        X_val = sum(X_coeff .* exp(-1j * w * (0:nX)));

        % Y(e^{jω})
        Y_val = sum(Y_coeff .* exp(-1j * w * (0:nY)));

        % 避免除零并钳制 |Y| 最小值，防止 |S_d| 爆炸
        Y_min_mag = 1e-3;
        if abs(Y_val) < Y_min_mag
            Y_val = Y_min_mag * exp(1j * angle(Y_val + 1e-12));
        end
        if abs(P_D_val) < 1e-12
            P_D_val = 1e-12;
        end

        % S_d(e^{jω})，钳制幅值上限
        S_d_val = (A_val / P_D_val) * (X_val / Y_val);
        mag(k) = min(abs(S_d_val), S_D_MAG_MAX);
    end
end
