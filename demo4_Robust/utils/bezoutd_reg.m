function [Rp, Sp, nrp, nsp] = bezoutd_reg(A, B, Hs, Hr, P, varargin)
%BEZOUTD_REG  改进的离散时间 Bezout/Diophantine 方程求解器
%
%   通过系数比较法求解 Bezout 方程：
%      A(z⁻¹)·Hs(z⁻¹)·Sp(z⁻¹) + B(z⁻¹)·Hr(z⁻¹)·Rp(z⁻¹) = P(z⁻¹)
%
%   相比原始 bezoutd.m 的改进：
%     1. 额外极点精确放置在原点 (z=0) — 等价于 P_full = [P, zeros(...)]，
%        彻底消除 rmin=1e-16 的浮点污染（原始 bezoutd 中 poly([...; rmin*e^{jθ}])
%        会产生微小但非零的系数，导致 Sylvester 矩阵不必要的条件数恶化）。
%     2. 条件数检测与详细诊断输出 — 输出系统维度、条件数、残差、最大系数幅值。
%     3. 实验性截断模式 (method='truncated') — 仅保留前 K 个方程 + pinv 最小范
%        数解。⚠ 注意：对于大规模系统（nah+nbh >> np），高阶方程是结构必需的，
%        截断会导致控制器退化和闭环不稳定。
%
% 语法:
%   [Rp, Sp, nrp, nsp] = bezoutd_reg(A, B, Hs, Hr, P)
%   [Rp, Sp, nrp, nsp] = bezoutd_reg(A, B, Hs, Hr, P, 'method', 'truncated')
%   [Rp, Sp, nrp, nsp] = bezoutd_reg(A, B, Hs, Hr, P, 'method', 'truncated', 'nextra_cap', 20)
%
% 输入参数:
%   A   - 系统分母多项式系数向量 [a0 a1 ... aNa]（降幂）
%   B   - 系统分子多项式系数向量 [b0 b1 ... bNb]（降幂，含延迟前导零）
%   Hs  - 控制器分母固定部分（降幂）
%   Hr  - 控制器分子固定部分（降幂）
%   P   - 期望闭环特征多项式（降幂）
%
% 可选参数（名值对）:
%   'method'     - 'standard'  (默认) — 完整 Sylvester 方阵 + 原点填充
%                  'truncated' (实验性) — 截断 + pinv 最小范数解
%   'nextra_cap' - 截断时允许的额外方程数上限 (默认: 2*np_orig)
%                  仅 method='truncated' 时使用
%
% 输出参数:
%   Rp  - 控制器分子多项式系数（降幂）
%   Sp  - 控制器分母多项式系数（降幂）
%   nrp - Rp 多项式的阶数
%   nsp - Sp 多项式的阶数
%
% 设计原理与限制:
%   - 'standard' 模式与原始 bezoutd 产生相同的 S、R 多项式（仅数值实现不同）。
%     原点填充等价于 poly([rootsPdes; zeros(nextra,1)])，物理含义为"辅助动态
%     尽可能快"，在数值上完全无浮点舍入误差。
%   - 对于大规模系统（如 ARMAX(30,30,30,22)），P_DX 仅 7 阶，需要 75 个额外
%     方程全部精确为零。这要求 A_HS*S + B_HR*R 在高阶项上严格对消，天然导致
%     R/S 系数达到 ±10⁴ 量级。这是结构性的，非数值问题。
%   - 'truncated' 模式通过 pinv 产生最小范数解，但实验表明对于此类系统会导
%     致控制器退化（R≡0）和闭环不稳定。仅适用于 nah+nbh ≈ np 的良定系统。
%
% 参考文献:
%   Landau, I.D. & Zito, G. "Digital Control Systems: Design,
%   Identification and Implementation", Springer, 2006.

    %% ──────────── 输入解析 ────────────
    p = inputParser;
    p.addRequired('A', @isnumeric);
    p.addRequired('B', @isnumeric);
    p.addRequired('Hs', @isnumeric);
    p.addRequired('Hr', @isnumeric);
    p.addRequired('P', @isnumeric);
    p.addParameter('method', 'standard', @(x) ismember(x, {'standard', 'truncated'}));
    p.addParameter('nextra_cap', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x) && x >= 0));
    p.parse(A, B, Hs, Hr, P, varargin{:});

    method      = p.Results.method;
    nextra_cap  = p.Results.nextra_cap;

    %% ──────────── 初始化与输入处理 ────────────
    PRECISION = 1e-14;

    % 确保所有输入多项式向量为行向量
    if size(A,1) > 1, A = A'; end
    if size(B,1) > 1, B = B'; end
    if size(Hs,1) > 1, Hs = Hs'; end
    if size(Hs,1) == 0, Hs = 1; end
    if size(Hr,1) > 1, Hr = Hr'; end
    if size(Hr,1) == 0, Hr = 1; end
    if size(P,1) > 1, P = P'; end

    %% ──────────── 计算多项式阶数 ────────────
    np_orig = length(P) - 1;
    nhs = length(Hs) - 1;
    nhr = length(Hr) - 1;

    %% ──────────── 计算复合多项式 ────────────
    if nhs > 0
        Ah = conv(A, Hs);
    else
        Ah = A * Hs;
    end
    nah = length(Ah) - 1;

    if nhr > 0
        Bh = conv(B, Hr);
    else
        Bh = B * Hr;
    end
    nbh = length(Bh) - 1;

    if np_orig > nah + nbh - 1
        warning('bezoutd_reg:TooManyPoles', ...
            '期望极点过多 (np=%d > nah+nbh-1=%d)。', np_orig, nah+nbh-1);
    end

    %% ──────────── 计算额外约束数量 ────────────
    nextra_full = nah + nbh - 1 - np_orig;
    total_dim   = nah + nbh;

    %% ──────────── 控制器多项式阶数 ────────────
    nsp = nbh - 1;
    nrp = nah - 1;

    %% ──────────── 处理期望闭环多项式 P ────────────
    switch method
        case 'standard'
            % ★ 标准模式：原点填充，等价于将所有额外极点放在 z=0
            %   物理含义：辅助极点尽可能快（dead-beat 风格）
            %   数值实现：直接用零填充 P，避免原始 bezoutd 中
            %   poly([rootsPdes; rmin*exp(jθ)]) 的浮点污染。
            if nextra_full > 0
                P_full = [P, zeros(1, nextra_full)];
            else
                P_full = P;
            end
            np = nah + nbh - 1;

            % 构建完整 Sylvester 方阵 (nah+nbh) × (nah+nbh)
            M = zeros(total_dim, total_dim);

            % S 系数对应列（Ah 的 Toeplitz 结构）
            for j = 1:(nsp + 1)
                col = zeros(total_dim, 1);
                col(j : j + nah) = Ah';
                M(:, j) = col;
            end

            % R 系数对应列（Bh 的 Toeplitz 结构）
            for j = 1:(nrp + 1)
                col = zeros(total_dim, 1);
                col(j : j + nbh) = Bh';
                M(:, nsp + 1 + j) = col;
            end

            P_col = P_full(:);
            cond_M = cond(M);

            if cond_M > 1e12
                warning('bezoutd_reg:IllConditioned', ...
                    'Sylvester 矩阵条件数极高 (%.2e)。对于大规模系统(nah+nbh>>np)，这是结构性的，不一定是数值错误。', cond_M);
            end

            X = M \ P_col;
            method_str = 'standard';

        case 'truncated'
            % ★ 实验性截断模式：仅保留前 K 个方程 + 最小范数解
            %   ⚠ 警告：对于 nah+nbh >> np 的系统，截断会破坏必要的
            %   高阶对消，导致控制器退化。

            if isempty(nextra_cap)
                nextra_cap_val = 2 * np_orig;
            else
                nextra_cap_val = nextra_cap;
            end
            nextra_used = max(0, min(nextra_cap_val, nextra_full));

            if nextra_used > 0
                P_full = [P, zeros(1, nextra_used)];
            else
                P_full = P;
            end
            np = length(P_full) - 1;

            n_rows = np + 1;
            M = zeros(n_rows, nsp + 1 + nrp + 1);

            % S 系数列
            for j = 1:(nsp + 1)
                col = zeros(n_rows, 1);
                row_end = min(n_rows, j + nah);
                if row_end >= j
                    n_copy = row_end - j + 1;
                    col(j:row_end) = Ah(1:n_copy)';
                end
                M(:, j) = col;
            end

            % R 系数列
            for j = 1:(nrp + 1)
                col = zeros(n_rows, 1);
                row_end = min(n_rows, j + nbh);
                if row_end >= j
                    n_copy = row_end - j + 1;
                    col(j:row_end) = Bh(1:n_copy)';
                end
                M(:, nsp + 1 + j) = col;
            end

            P_col = P_full(:);
            cond_M = cond(M);

            % pinv 求最小范数解
            X = pinv(M) * P_col;
            method_str = sprintf('truncated(n_rows=%d)', n_rows);

            fprintf('[bezoutd_reg truncated] n_rows=%d (完整: %d), 截断比=%.1f%%\n', ...
                n_rows, total_dim, 100*n_rows/total_dim);
    end

    % 确保系数为实数
    X = real(X);

    %% ──────────── 提取 R 和 S 多项式 ────────────
    Sp = X(1 : nsp+1)';
    Rp = X(nsp+2 : nsp+nrp+2)';

    %% ──────────── 数值验证 ────────────
    % 验证前 np_orig+1 个系数（来自 P_DX 的有意义约束）
    P_cl1 = conv(Ah, Sp);
    P_cl2 = conv(Bh, Rp);
    max_len = max(length(P_cl1), length(P_cl2));
    P_cl1(end+1:max_len) = 0;
    P_cl2(end+1:max_len) = 0;
    P_cl = P_cl1 + P_cl2;

    cmp_len_active = min(np_orig + 1, length(P_cl));
    P_target_active = P(1:cmp_len_active);
    P_cl_active = P_cl(1:cmp_len_active);
    residual_active = max(abs(P_target_active - P_cl_active));

    % 完整残差（如果标准模式）
    if nextra_full > 0 && strcmp(method, 'standard')
        cmp_len_full = min(total_dim, length(P_cl));
        P_target_full = [P, zeros(1, cmp_len_full - np_orig - 1)];
        P_cl_full = P_cl(1:cmp_len_full);
        residual_full = max(abs(P_target_full - P_cl_full));
    else
        residual_full = NaN;
    end

    max_R = max(abs(Rp));
    max_S = max(abs(Sp));
    norm_X = norm(X);

    fprintf('bezoutd_reg: 维度=%d, np_orig=%d, nextra_full=%d, 方法=%s, cond=%.2e\n', ...
        total_dim, np_orig, nextra_full, method_str, cond_M);
    fprintf('  残差(P_DX部分)=%.2e, 残差(完整)=%.2e, max|R|=%.2e, max|S|=%.2e, ||X||=%.2e\n', ...
        residual_active, residual_full, max_R, max_S, norm_X);

    if residual_active > 1e-6
        warning('bezoutd_reg:ActiveResidual', ...
            'P_DX 部分 Bezout 残差较大: %.2e', residual_active);
    end

end
