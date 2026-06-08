function phi = penalty_function(f, g1, g2, Lambda, kappa)
% PENALTY_FUNCTION 静态惩罚函数
%   将约束违例转换为惩罚项，生成增广目标函数值。
%   用于将约束多目标优化问题转换为无约束问题。
%
%   输入:
%     f       - 目标函数值向量 (列向量)
%     g1      - 约束1违例值向量: g1 > 0 表示违例（下界约束）
%     g2      - 约束2违例值向量: g2 > 0 表示违例（上界约束）
%     Lambda  - 惩罚系数 (可选, 默认: 10)
%     kappa   - 惩罚指数 (可选, 默认: 1)
%
%   输出:
%     phi     - 惩罚后的增广目标函数值向量 (与 f 同维度)
%
%   惩罚公式:
%     δ = Lambda * ( mean(max(0, g1)^κ) + mean(max(0, g2)^κ) )
%     φ = f + δ
%
%   其中:
%     - max(0, g) 只对违例部分进行惩罚（g > 0）
%     - 使用 mean() 而非 sum()，使惩罚项与频率采样点数无关
%     - κ = 1: 线性惩罚（温和，推荐用于频域约束）
%     - κ = 2: 二次惩罚（对严重违例更敏感）
%
%   使用建议:
%     - Lambda 应适中，避免惩罚项淹没目标函数
%     - 典型值: Lambda ∈ [1, 100], kappa = 1
%     - 若优化结果违例严重，可适度增大 Lambda
%     - 若 Archive 退化（解数量过少），减小 Lambda
%
%   示例:
%     f = [0.5; 0.3];           % 2个目标
%     g1 = [-0.1; 0.2; 0.05];   % 有2处违例
%     g2 = [-0.3; -0.1; -0.2];  % 无违例
%     phi = penalty_function(f, g1, g2, 10, 1);
%     % phi ≈ [0.5+10*mean(0.2+0.05); 0.3+10*mean(0.2+0.05)]
%     % phi ≈ [0.5+0.83; 0.3+0.83] = [1.33; 1.13]

    % =====================================================================
    % 输入参数验证
    % =====================================================================
    if nargin < 3
        error('penalty_function:NotEnoughInputs', ...
              '至少需要3个输入参数: f, g1, g2。');
    end

    % 设置默认值
    if nargin < 4 || isempty(Lambda)
        Lambda = 10;
    end
    if nargin < 5 || isempty(kappa)
        kappa = 1;
    end

    % 验证参数有效性
    if Lambda <= 0
        error('penalty_function:InvalidLambda', ...
              '惩罚系数 Lambda 必须为正数，当前值为 %g。', Lambda);
    end

    if kappa <= 0
        error('penalty_function:InvalidKappa', ...
              '惩罚指数 kappa 必须为正数，当前值为 %g。', kappa);
    end

    % 确保输入为列向量
    f  = f(:);
    g1 = g1(:);
    g2 = g2(:);

    % =====================================================================
    % 计算惩罚项
    % =====================================================================
    % 只对违例部分 (g > 0) 进行惩罚
    g1_violation = max(0, g1);   % 下界违例
    g2_violation = max(0, g2);   % 上界违例

    % 使用均值替代求和，使惩罚项与频率采样点数无关
    penalty_g1 = mean(g1_violation .^ kappa);
    penalty_g2 = mean(g2_violation .^ kappa);

    delta = Lambda * (penalty_g1 + penalty_g2);

    % =====================================================================
    % 计算增广目标值
    % =====================================================================
    % 惩罚项作为标量加到每个目标分量上
    phi = f + delta;

    % =====================================================================
    % 诊断信息（仅在违例严重时给出警告）
    % =====================================================================
    total_violation = mean(g1_violation) + mean(g2_violation);
    if total_violation > 1
        warning('penalty_function:HighViolation', ...
                '约束违例严重: avg_g1_max=%.3f, avg_g2_max=%.3f, delta=%.3f', ...
                mean(g1_violation), mean(g2_violation), delta);
    end

end
