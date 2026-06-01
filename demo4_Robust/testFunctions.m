function f = testFunctions(x, obj_func)
% TESTFUNCTIONS 多目标标准测试函数（附录A）
%   用于验证 ε-MOPSO 算法的正确性。
%
%   输入:
%     x        - 决策变量向量 (列向量, n_var × 1), 所有分量 ∈ [0, 1]
%     func_num - 函数编号 (1-4)
%
%   输出:
%     f        - 目标函数值向量 (列向量, n_obj × 1)
%
%   F1: 2目标, ≥2维, 含Rastrigin型g函数，Pareto前沿为凸曲线
%   F2: 2目标, 30维, 高维ZDT型函数，验证高维优化能力
%   F3: 3目标, 3维, 球坐标参数化, Pareto前沿为球面第一象限
%   F4: 3目标, 3维, f3含正弦振荡项, Pareto前沿为非连续曲面
%
%   参考文献:
%     Coello Coello, C.A., et al. "Evolutionary Algorithms for Solving
%     Multi-Objective Problems", 2nd Ed., Springer, 2007.

    % =====================================================================
    % 输入参数验证
    % =====================================================================
    if nargin < 2
        error('testFunctions:NotEnoughInputs', ...
              '需要2个输入参数: x 和 func_num。');
    end

    % 确保 x 是列向量
    x = x(:);
    n = length(x);

    % 检查 x 是否在 [0,1] 范围内（允许微小容差）
    tol = 1e-8;
    if any(x < -tol) || any(x > 1 + tol)
        warning('testFunctions:OutOfBounds', ...
                '部分决策变量超出 [0,1] 范围，将被裁剪。');
        x = max(0, min(1, x));
    end

    % =====================================================================
    % 根据函数编号计算目标值
    % =====================================================================
    switch obj_func
        case 1
            % -------------------------------------------------------------
            % F1: 2目标, n维 (n ≥ 2)
            % f1 = x(1)
            % f2 = g * h
            % g  = 1 + 10*(n-1) + Σ[x_i² - 10*cos(4π*x_i)], i=2,...,n
            % h  = 1 - sqrt(f1 / g)
            %
            % Pareto前沿: f2 = 1 - sqrt(f1), f1 ∈ [0, 1]
            % 特点: g函数含Rastrigin多峰项，难以收敛到Pareto前沿
            % -------------------------------------------------------------
            if n < 2
                error('testFunctions:InvalidDimension', ...
                      'F1 要求决策变量维度至少为2，当前维度为 %d。', n);
            end

            f1 = x(1);

            % Rastrigin型 g 函数
            xi = x(2:end);
            g = 1 + 10 * (n - 1) + sum(xi.^2 - 10 * cos(4 * pi * xi));

            h = 1 - sqrt(f1 / g);
            f2 = g * h;

            f = [f1; f2];

        case 2
            % -------------------------------------------------------------
            % F2: 2目标, n维 (n ≥ 2)
            % f1 = x(1)
            % f2 = g * (1 - (f1/g)^2)
            % g  = 1 + 9/(n-1) * Σ x(i), i=2,...,n
            %
            % Pareto前沿: f2 = 1 - f1², f1 ∈ [0, 1]
            % 特点: 高维ZDT型，验证算法在高维空间中的收敛性
            % -------------------------------------------------------------
            if n < 2
                error('testFunctions:InvalidDimension', ...
                      'F2 要求决策变量维度至少为2，当前维度为 %d。', n);
            end

            f1 = x(1);

            % ZDT型 g 函数
            g = 1 + 9 / (n - 1) * sum(x(2:end));

            f2 = g * (1 - (f1 / g)^2);

            f = [f1; f2];

        case 3
            % -------------------------------------------------------------
            % F3: 3目标, 3维
            % 将 x 映射到球坐标:
            %   r     = 2 + 4*(x1 - 0.5)     ∈ [0, 4]
            %   theta = π * x2 / 2           ∈ [0, π/2]
            %   phi   = π * x3 / 2           ∈ [0, π/2]
            % f1 = r * cos(theta) * cos(phi)
            % f2 = r * cos(theta) * sin(phi)
            % f3 = r * sin(theta)
            %
            % Pareto前沿: 球面第一象限部分，即 f1² + f2² + f3² = r²
            % 特点: 球坐标参数化，前沿为连续光滑曲面
            % -------------------------------------------------------------
            if n < 3
                error('testFunctions:InvalidDimension', ...
                      'F3 要求决策变量维度至少为3，当前维度为 %d。', n);
            end

            r     = 2 + 4 * (x(1) - 0.5);       % r ∈ [0, 4]
            theta = pi * x(2) / 2;               % θ ∈ [0, π/2]
            phi   = pi * x(3) / 2;               % φ ∈ [0, π/2]

            f1 = r * cos(theta) * cos(phi);
            f2 = r * cos(theta) * sin(phi);
            f3 = r * sin(theta);

            f = [f1; f2; f3];

        case 4
            % -------------------------------------------------------------
            % F4: 3目标, 3维
            % f1 = x(1)
            % f2 = x(2)
            % f3 = x(3) + sin(5π*x1) * sin(5π*x2)
            %
            % Pareto前沿: f3 = f3_min(f1,f2) (非连续曲面)
            % 特点: f3含正弦振荡项，Pareto前沿不连续，
            %       用于测试算法在复杂前沿上的分布能力
            % -------------------------------------------------------------
            if n < 3
                error('testFunctions:InvalidDimension', ...
                      'F4 要求决策变量维度至少为3，当前维度为 %d。', n);
            end

            f1 = x(1);
            f2 = x(2);
            f3 = x(3) + sin(5 * pi * x(1)) * sin(5 * pi * x(2));

            f = [f1; f2; f3];

        otherwise
            error('testFunctions:InvalidFuncNum', ...
                  '无效的函数编号 %d。有效范围: 1-4。', obj_func);
    end

end
