function build_all()
%BUILD_ALL 生成控制器库 + 全部离线模型 + 全部 RT 模型
%
%   前置: 已运行 export_frozen_params（存在 anc_frozen_params.mat）。
%   离线/RT 矩阵直接由冻结参数键驱动，避免清单漂移。

if ~license('test', 'Simulink')
    error('build_all:NoSimulink', '本机无 Simulink 许可证，无法生成模型。');
end

build_controller_library();

dspaceRoot = fileparts(fileparts(mfilename('fullpath')));
matFile = fullfile(dspaceRoot, 'params', 'anc_frozen_params.mat');
payload = load(matFile, 'frozen');
keys = setdiff(fieldnames(payload.frozen), {'meta', 'plant__secondary'});

for index = 1:numel(keys)
    parts = strsplit(keys{index}, '__');
    build_offline_model(parts{1}, parts{2}, parts{3});
end

% RT 部署矩阵: 全部冻结案例都生成（部署优先级见实验规程 Stage 2/3）
for index = 1:numel(keys)
    parts = strsplit(keys{index}, '__');
    build_rt_model(parts{1}, parts{2}, parts{3});
end

fprintf('\n全部模型生成完毕（models/）。下一步: verify_offline_equivalence\n');
end
