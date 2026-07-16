function frozen = export_frozen_params()
%EXPORT_FROZEN_PARAMS 从冻结 stage 结果收割 Simulink/dSPACE 部署参数
%
%   frozen = export_frozen_params()
%
%   数据流（绝不手抄系数）:
%     tests/output/.../demo1234_stage_results.mat (冻结调参)
%       → 对每个 (demo, test, variant) 重跑 tests/internal/controller_demoN
%         （复用其设计代码路径: RST 扫描 / eMOPSO / dlqr+dlqe / F(z) 设计）
%       → 收割 result.extra 中的设计量为扁平 double struct
%       → dspace/params/anc_frozen_params.mat
%
%   键名: demo1__T1__fixed / demo2__T2__adaptive / ... 以及 plant__secondary。
%   范围: 主实验 T1/T2；demo5 仅 T1（T2 Case A 不适用）；demo5_imc 不迁移
%   （评估需求超限）。
%
%   注意: 依赖本地存在 stage MAT（*.mat 不入 git）。若缺失，先运行
%   tests/run_cylinder1dm_stage()。

dspaceRoot = fileparts(fileparts(mfilename('fullpath')));
projectRoot = fileparts(dspaceRoot);
stageFile = fullfile(projectRoot, 'tests', 'output', 'cylinder1dm_2k', ...
    'demo1234_stage', 'demo1234_stage_results.mat');
if ~isfile(stageFile)
    error('export_frozen_params:MissingStage', ...
        ['未找到冻结 stage 结果: %s\n' ...
         '请先运行 tests/run_cylinder1dm_stage() 生成。'], stageFile);
end
payload = load(stageFile, 'stage');
stage = payload.stage;

signals = load_cylinder1dm_signals('evaluation');
fs = signals.fs;

frozen = struct();
frozen.meta = struct( ...
    'source_stage_mat', stageFile, ...
    'stage_generated_at', char(string(stage.generated_at)), ...
    'exported_at', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), ...
    'git_describe', repo_git_describe(projectRoot), ...
    'fs', fs, ...
    'Ts', 1/fs, ...
    'evaluation_seed', stage.suite.meta.evaluation_seed, ...
    'actuator_limit', stage.suite.meta.actuator_limit, ...
    'tuning_peak_limit', stage.suite.meta.tuning_peak_limit);

% 被控对象（离线模型的次级路径；B 含前导零 = 一拍计算延迟约定）
frozen.plant__secondary = struct( ...
    'A', signals.model_data.A(:).', ...
    'B_full', signals.model_data.B(:).', ...
    'fs', fs, ...
    'path_family', signals.model_data.path_family);

for index = 1:numel(stage.selections)
    sel = stage.selections(index);

    % 迁移范围过滤
    if strcmp(sel.demo, 'demo5_imc')
        fprintf('跳过 %s/%s/%s（评估需求超限，不迁移）\n', ...
            sel.demo, sel.variant, sel.test);
        continue;
    end
    if strcmp(sel.demo, 'demo5') && ~strcmp(sel.test, 'T1')
        fprintf('跳过 %s/%s/%s（Case A 不适用，仅 T1 部署）\n', ...
            sel.demo, sel.variant, sel.test);
        continue;
    end

    fprintf('\n===== 收割 %s / %s / %s =====\n', sel.demo, sel.test, sel.variant);
    runner = local_runner(sel.demo);
    result = runner(signals, sel.test, sel.variant, sel.params);

    % 门禁: 重跑结果必须与冻结 stage 评价逐位一致（同时证明
    % controller_demo2/5 的 extra 附加是纯输出、不影响仿真数值）。
    row = strcmp(stage.summary.demo, sel.demo) ...
        & strcmp(stage.summary.variant, sel.variant) ...
        & strcmp(stage.summary.test, sel.test);
    supp_frozen = stage.summary.evaluation_suppression_db(row);
    if abs(result.supp_db - supp_frozen) > 1e-9
        error('export_frozen_params:SuppMismatch', ...
            '%s/%s/%s 重跑 supp_db=%.12f 与冻结值 %.12f 不一致。', ...
            sel.demo, sel.variant, sel.test, result.supp_db, supp_frozen);
    end

    key = sprintf('%s__%s__%s', sel.demo, sel.test, sel.variant);
    entry = harvest(sel, result, fs);
    entry.supp_db_reference = result.supp_db;
    entry.u_demand_max_reference = result.extra.u_demand_max;
    entry.candidate = sel.candidate;
    entry.params_json = jsonencode(sel.params);   % 等效性验证重跑脚本控制器用
    frozen.(key) = entry;
end

outFile = fullfile(dspaceRoot, 'params', 'anc_frozen_params.mat');
save(outFile, 'frozen');
fprintf('\n冻结参数已导出: %s\n', outFile);
fprintf('键: %s\n', strjoin(setdiff(fieldnames(frozen), {'meta'}), ', '));
end

%% ========================================================================
function runner = local_runner(demo)
switch demo
    case 'demo1', runner = @controller_demo1;
    case 'demo2', runner = @controller_demo2;
    case 'demo3', runner = @controller_demo3;
    case 'demo4', runner = @controller_demo4;
    case 'demo5', runner = @controller_demo5;
    otherwise
        error('export_frozen_params:UnknownDemo', 'Unknown demo: %s', demo);
end
end

function entry = harvest(sel, result, fs)
%HARVEST 收割一个控制器实例的部署参数（扁平 double 字段）
ex = result.extra;
p = sel.params;
entry = struct();
entry.demo = sel.demo;
entry.test = sel.test;
entry.variant = sel.variant;
entry.actuator_limit = ex.actuator_limit;

switch sel.demo
    case {'demo1', 'demo4'}
        entry.R0 = ex.R0(:).';
        entry.S0 = ex.S0(:).';
        entry.P0 = ex.P0(:).';
        entry.A  = ex.A(:).';
        entry.B_full = ex.B(:).';
        entry.B_star = entry.B_full(2:end);
        % demo1 运行时无 control_scale（等价于 1）；demo4 有。
        entry.control_scale = get_or(p, 'control_scale', 1);
        % demo1 adaptive 控制律不做 /S0(1)（S0(1)=1 时两者位级一致），
        % demo4 adaptive 显式 /S0(1)。步函数用 s0_norm 参数区分：
        % 传 1 即精确复现 demo1（IEEE754 除以 1.0 为恒等）。
        if strcmp(sel.demo, 'demo4')
            entry.s0_norm = ex.S0(1);
        else
            entry.s0_norm = 1;
        end
        if strcmp(sel.variant, 'adaptive')
            entry.nQ      = ex.nQ;
            entry.lambda1 = ex.lambda1;
            entry.lambda2 = ex.lambda2;
            entry.F_diag  = get_or(p, 'F_diag', 1e-3);
            entry.Q_final_reference = ex.Q_final(:).';
        end

    case 'demo2'
        entry.Af = ex.Af;  entry.Bf = ex.Bf;  entry.Cf = ex.Cf;
        entry.A_design = ex.A_design;
        entry.B_design = ex.B_design;
        entry.C_design = ex.C_design;
        entry.K = ex.K;   entry.L = ex.L;
        entry.n_plant = ex.n_plant;
        entry.n_dist  = ex.n_dist;
        entry.ramp_samples = ex.ramp_samples;
        if strcmp(sel.variant, 'adaptive')
            entry.A_t12 = ex.A_t12;
            entry.B_t12 = ex.B_t12;
            entry.C_t12 = ex.C_t12;
            entry.nQ = ex.nQ;
            entry.lambda = ex.lambda;
            entry.F_init = ex.F_init;
            entry.adaptation_gain = ex.adaptation_gain;
            entry.method_is_lms = double(strcmp(ex.adaptation_method, 'lms'));
            entry.bp_b = ex.bp_b(:).';
            entry.bp_a = ex.bp_a(:).';
            entry.warmup_samples = ex.warmup_samples;
            entry.theta_norm_limit = ex.theta_norm_limit;
            entry.lms_epsilon = ex.lms_epsilon;
            entry.lms_leakage = ex.lms_leakage;
            entry.rls_regularization = ex.rls_regularization;
        end

    case 'demo3'
        % 运行时被控对象多项式（与引擎一致: b_plant = B_poly(2:end)）
        modelA = evaluation_model_polynomials();
        entry.a_plant = modelA.A(:).';
        entry.b_plant = modelA.B_full(2:end);
        entry.ramp_samples = round(get_or(p, 'ramp_seconds', 0.1) * fs);
        entry.control_scale = get_or(p, 'control_scale', 1);
        entry.theta_norm_limit = get_or(p, 'theta_norm_limit', Inf);
        nA = numel(entry.a_plant) - 1;
        Lb = numel(entry.b_plant);
        if strcmp(sel.variant, 'fixed')
            entry.N_fir = ex.N_fir;
            entry.theta_fixed = ex.theta_fixed(:).';
            entry.F_num = ex.F_num(:).';
            entry.F_den = ex.F_den(:).';
            entry.k_start = max(ex.N_fir, nA + Lb) + 1;
        else
            if ~strcmp(ex.adaptive_structure, 'imc_fxnlms')
                error('export_frozen_params:UnexpectedStructure', ...
                    'demo3 冻结自适应结构应为 imc_fxnlms，实际: %s', ...
                    ex.adaptive_structure);
            end
            entry.N_adapt = ex.N_adapt;
            entry.mu_nlms = ex.mu_nlms;
            entry.theta_init = ex.theta_init(:).';
            entry.k_start = max([ex.N_adapt, nA, Lb]) + 1;
        end

    case 'demo5'
        entry.q = ex.q;
        entry.k_gain = ex.k;
        entry.epsilon = ex.epsilon;
        entry.theta_init = ex.theta_init(:).';
        entry.theta_min = ex.theta_min;
        entry.theta_max = ex.theta_max;
        entry.sign_S = ex.sign_S(:).';
        entry.sign_S0 = ex.sign_S0;
        entry.dc_cancel = double(ex.dc_cancel);
        entry.ramp_samples = ex.ramp_samples;
        entry.k_start = ex.k_start;
        entry.A_sec = ex.A_sec(:).';
        entry.B_poly = ex.B_poly(:).';
        entry.method_is_exact = double(strcmp(ex.method, 'exact'));
        entry.output_is_updated = double(strcmp(ex.output_timing, 'updated'));
        if ~entry.method_is_exact || ~entry.output_is_updated
            error('export_frozen_params:UnexpectedDemo5Config', ...
                'demo5 冻结配置应为 exact/updated，实际: %s/%s', ...
                ex.method, ex.output_timing);
        end

    otherwise
        error('export_frozen_params:UnknownDemo', 'Unknown demo: %s', sel.demo);
end
end

function modelA = evaluation_model_polynomials()
% demo3 的被控对象多项式与全局评价信号一致，直接取 evaluation 信号约定
signals = load_cylinder1dm_signals('evaluation');
modelA = struct('A', signals.model_data.A(:).', ...
    'B_full', signals.model_data.B(:).');
end

function value = get_or(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    value = s.(name);
else
    value = default;
end
end

function description = repo_git_describe(projectRoot)
[status, out] = system(sprintf( ...
    'git -C "%s" describe --always --dirty 2>/dev/null', projectRoot));
if status == 0
    description = strtrim(out);
else
    description = 'unknown';
end
end
