function design_results = design_observer_controller(Ao, Bo, Co, Do, fs, varargin)
% design_observer_controller - 基于能观标准型设计观测器与状态反馈控制器
%
% 用法:
%   results = design_observer_controller(Ao, Bo, Co, Do, fs)
%   results = design_observer_controller(Ao, Bo, Co, Do, fs, 'Q_obs', 10, 'R_obs', 1, ...
%                                                     'Q_ctrl', 1e-6, 'R_ctrl', 1000)
%
% 输入:
%   Ao, Bo, Co, Do - 能观标准型状态空间矩阵 (离散时间)
%   fs             - 采样频率 (Hz)
%   可选参数:
%     'Q_obs'  - 观测器状态噪声权重 (默认: eye(n)*10)
%     'R_obs'  - 观测器测量噪声权重 (默认: 1)
%     'Q_ctrl' - 控制器状态权重 (默认: eye(n)*1e-6)
%     'R_ctrl' - 控制器输入权重 (默认: 1000)
%
% 输出:
%   design_results - struct，包含:
%     .Ao, .Bo, .Co, .Do       - 状态空间矩阵
%     .L_lqr                   - LQR观测器增益
%     .K_lqr                   - LQR控制器增益
%     .Nbar                    - 前馈增益
%     .observer_poles          - 观测器极点
%     .closed_loop_poles       - 闭环极点
%     .is_observable           - 是否完全能观
%     .is_controllable         - 是否完全能控

    %% 解析参数
    p = inputParser;
    addParameter(p, 'Q_obs',  [], @isnumeric);
    addParameter(p, 'R_obs',  [], @isnumeric);
    addParameter(p, 'Q_ctrl', [], @isnumeric);
    addParameter(p, 'R_ctrl', [], @isnumeric);
    parse(p, varargin{:});

    nA = size(Ao, 1);
    ts = 1/fs;

    %% 观测器设计
    fprintf('\n=== 观测器设计 ===\n');

    obsv_matrix = obsv(Ao, Co);
    rank_obsv = rank(obsv_matrix);
    cond_obsv = cond(obsv_matrix);
    fprintf('能观性矩阵秩: %d (系统阶数: %d), 条件数: %.2e\n', rank_obsv, nA, cond_obsv);

    if cond_obsv > 1e10
        fprintf('警告: 系统几乎不能观 (条件数极大)，极点配置法可能失败。\n');
    end

    if rank_obsv == nA
        fprintf('✓ 系统完全能观，开始设计观测器。\n');

        % LQR/Kalman滤波器方法 (推荐)
        if isempty(p.Results.Q_obs)
            Q_obs = eye(nA) * 10;
        else
            Q_obs = p.Results.Q_obs;
        end
        if isempty(p.Results.R_obs)
            R_obs = 1;
        else
            R_obs = p.Results.R_obs;
        end

        try
            L_lqr = lqr(Ao', Co', Q_obs, R_obs)';
            obs_poles = eig(Ao - L_lqr*Co);
            fprintf('LQR观测器设计成功。极点: ');
            print_poles(obs_poles);
        catch ME
            fprintf('✗ LQR观测器设计失败: %s\n', ME.message);
            L_lqr = [];
            obs_poles = [];
        end
    else
        fprintf('✗ 系统不完全能观，无法设计全阶观测器。\n');
        L_lqr = [];
        obs_poles = [];
    end

    %% 状态反馈控制器设计
    fprintf('\n=== 状态反馈控制器设计 ===\n');

    ctrb_matrix = ctrb(Ao, Bo);
    rank_ctrb = rank(ctrb_matrix);
    fprintf('能控性矩阵秩: %d (系统阶数: %d)\n', rank_ctrb, nA);

    if rank_ctrb == nA
        fprintf('✓ 系统完全能控，开始设计状态反馈控制器。\n');

        if isempty(p.Results.Q_ctrl)
            Q_ctrl = eye(nA) * 1e-6;
        else
            Q_ctrl = p.Results.Q_ctrl;
        end
        if isempty(p.Results.R_ctrl)
            R_ctrl = 1000;
        else
            R_ctrl = p.Results.R_ctrl;
        end

        try
            [K_lqr, ~, cl_poles] = lqr(Ao, Bo, Q_ctrl, R_ctrl);
            fprintf('LQR控制器设计成功。闭环极点: ');
            print_poles(cl_poles);

            % 前馈增益
            Nbar = 1 / (Co * inv(eye(nA) - Ao + Bo*K_lqr) * Bo);
            fprintf('前馈增益 Nbar = %.4f\n', Nbar);
        catch ME
            fprintf('✗ LQR控制器设计失败: %s\n', ME.message);
            K_lqr = [];
            cl_poles = [];
            Nbar = [];
        end
    else
        fprintf('✗ 系统不完全能控，无法设计状态反馈控制器。\n');
        K_lqr = [];
        cl_poles = [];
        Nbar = [];
    end

    %% 输出反馈控制器 (分离原理)
    fprintf('\n=== 输出反馈控制器 (分离原理) ===\n');

    if rank_obsv == nA && rank_ctrb == nA
        fprintf('✓ 系统能观能控，可设计输出反馈控制器。\n');

        L = L_lqr;
        K = K_lqr;

        A_aug = [Ao - Bo*K,        Bo*K;
                 zeros(nA, nA),    Ao - L*Co];
        aug_poles = eig(A_aug);

        fprintf('增广系统极点: ');
        print_poles(aug_poles);
        fprintf('  控制律: u = -K*x_hat + Nbar*r\n');
        fprintf('  观测器: x_hat[k+1] = (Ao - L*Co)*x_hat + Bo*u + L*y\n');
    else
        fprintf('✗ 无法设计输出反馈控制器 (能观=%d, 能控=%d)\n', ...
            rank_obsv == nA, rank_ctrb == nA);
        L = []; K = []; aug_poles = [];
    end

    %% 输出结果
    design_results = struct();
    design_results.Ao = Ao;  design_results.Bo = Bo;
    design_results.Co = Co;  design_results.Do = Do;
    design_results.fs = fs;  design_results.ts = ts;

    if ~isempty(L)
        design_results.L = L;
        design_results.L_lqr = L_lqr;
        design_results.observer_poles = obs_poles;
    end
    if ~isempty(K)
        design_results.K = K;
        design_results.K_lqr = K_lqr;
        design_results.closed_loop_poles = cl_poles;
        design_results.Nbar = Nbar;
    end

    design_results.is_observable = (rank_obsv == nA);
    design_results.is_controllable = (rank_ctrb == nA);

    fprintf('\n=== 设计结果已保存到 design_results 结构体 ===\n');
end

function print_poles(poles)
    for i = 1:length(poles)
        if isreal(poles(i))
            fprintf('%.4f ', poles(i));
        else
            fprintf('%.4f±%.4fj ', real(poles(i)), abs(imag(poles(i))));
        end
    end
    fprintf('\n');
end
