function mdl = slk_new_model(name, kind, Ts)
%SLK_NEW_MODEL 新建统一求解器配置的模型/库
%
%   mdl = slk_new_model('offline_demo1_T1_fixed', 'model', 5e-4)
%   mdl = slk_new_model('lib_anc_controllers', 'library')
%
%   已存在同名已加载模型时先强制关闭；返回模型名。

if nargin < 2, kind = 'model'; end

if bdIsLoaded(name)
    close_system(name, 0);
end

switch kind
    case 'library'
        new_system(name, 'Library');
    case 'model'
        new_system(name, 'Model');
        set_param(name, ...
            'SolverType', 'Fixed-step', ...
            'Solver', 'FixedStepDiscrete', ...
            'FixedStep', num2str(Ts, '%.17g'), ...
            'EnableMultiTasking', 'off', ...
            'AlgebraicLoopMsg', 'error', ...
            'SaveFormat', 'Array', ...
            'ReturnWorkspaceOutputs', 'on');
    otherwise
        error('slk_new_model:UnknownKind', 'Unknown kind: %s', kind);
end
mdl = name;
end
