function mdlFile = build_offline_model(demo, test, variant)
%BUILD_OFFLINE_MODEL 生成一个离线闭环模型 models/offline_<demo>_<test>_<variant>.slx
%
%   结构（与 tests/internal 引擎的样本内顺序一致，Unit Delay 断开代数环）:
%
%     d_in ─────────────────────────┐
%     plant(u(k-1)) → anti(k) ──→ Sum → e(k) ─→ [controller] ─→ u(k) ─→ 1/z ─┐
%     enable(=1) ───────────────────────────────↑        └→ 日志          └→ plant
%
%   - 扰动 d 由 From Workspace 注入（变量 anc_d_in = [t, d]），与脚本仿真
%     使用完全同一向量 → 可逐样本对拍。
%   - 控制器与被控对象都是库链接（lib_anc_controllers），参数经模型
%     工作区解析（build 时写入，随 .slx 保存）。
%   - demo2 的被控对象用状态空间块（与其脚本引擎一致），其余用多项式递推。

dspaceRoot = fileparts(fileparts(mfilename('fullpath')));
modelsDir = fullfile(dspaceRoot, 'models');
if ~isfolder(modelsDir), mkdir(modelsDir); end

entry = load_anc_params(demo, test, variant);
plant = load_anc_params('plant', 'secondary');
spec = controller_registry(demo, variant);
Ts = 1 / plant.fs;

libFile = fullfile(modelsDir, 'lib_anc_controllers.slx');
if ~bdIsLoaded('lib_anc_controllers')
    if ~isfile(libFile)
        error('build_offline_model:MissingLibrary', ...
            '未找到 %s，请先运行 build_controller_library。', libFile);
    end
    load_system(libFile);
end

mdl = sprintf('offline_%s_%s_%s', demo, test, variant);
slk_new_model(mdl, 'model', Ts);

%% ---- 块 ----
add_block('simulink/Sources/From Workspace', [mdl '/d_in'], ...
    'VariableName', 'anc_d_in', ...
    'SampleTime', num2str(Ts, '%.17g'), ...
    'Interpolate', 'off', ...
    'ZeroCross', 'off', ...
    'OutputAfterFinalValue', 'Setting to zero', ...
    'Position', [40 96 110 124]);
add_block('simulink/Sources/Constant', [mdl '/enable'], ...
    'Value', '1', 'Position', [40 200 110 230]);
add_block('simulink/Math Operations/Sum', [mdl '/err_sum'], ...
    'Inputs', '++', 'Position', [170 95 200 125]);
add_block(['lib_anc_controllers/' spec.block], [mdl '/controller'], ...
    'Position', [260 80 440 240]);
add_block('simulink/Discrete/Unit Delay', [mdl '/u_delay'], ...
    'SampleTime', num2str(Ts, '%.17g'), 'Position', [500 320 540 360]);

if strcmp(demo, 'demo2')
    add_block('lib_anc_controllers/plant_secondary_ss', [mdl '/plant'], ...
        'Position', [170 300 300 380]);
else
    add_block('lib_anc_controllers/plant_secondary', [mdl '/plant'], ...
        'Position', [170 300 300 380]);
end

logNames = {'e_log', 'u_log', 'ud_log', 'clip_log', 'diag_log'};
for index = 1:numel(logNames)
    add_block('simulink/Sinks/To Workspace', [mdl '/' logNames{index}], ...
        'VariableName', logNames{index}, ...
        'SaveFormat', 'Array', ...
        'Position', [520 60*index+20 590 60*index+50]);
end

%% ---- 连线 ----
% 控制器端口顺序: e[, ref], enable → u, u_demand, clipped, diag_out
add_line(mdl, 'd_in/1', 'err_sum/1');
add_line(mdl, 'plant/1', 'err_sum/2');
add_line(mdl, 'err_sum/1', 'controller/1');
if spec.has_ref
    add_block('simulink/Sources/From Workspace', [mdl '/ref_in'], ...
        'VariableName', 'anc_ref_in', ...
        'SampleTime', num2str(Ts, '%.17g'), ...
        'Interpolate', 'off', ...
        'ZeroCross', 'off', ...
        'OutputAfterFinalValue', 'Setting to zero', ...
        'Position', [40 146 110 174]);
    add_line(mdl, 'ref_in/1', 'controller/2');
    add_line(mdl, 'enable/1', 'controller/3');
else
    add_line(mdl, 'enable/1', 'controller/2');
end
add_line(mdl, 'controller/1', 'u_delay/1');
add_line(mdl, 'u_delay/1', 'plant/1');
add_line(mdl, 'err_sum/1', 'e_log/1');
add_line(mdl, 'controller/1', 'u_log/1');
add_line(mdl, 'controller/2', 'ud_log/1');
add_line(mdl, 'controller/3', 'clip_log/1');
add_line(mdl, 'controller/4', 'diag_log/1');

%% ---- 模型工作区参数（随 .slx 保存）----
mws = get_param(mdl, 'ModelWorkspace');
paramMap = spec.params;
for index = 1:size(paramMap, 1)
    assignin(mws, paramMap{index, 1}, entry.(paramMap{index, 2}));
end
if strcmp(demo, 'demo2')
    assignin(mws, 'PLANT_AF', entry.Af);
    assignin(mws, 'PLANT_BF', entry.Bf);
    assignin(mws, 'PLANT_CF', entry.Cf);
else
    assignin(mws, 'PLANT_A', plant.A);
    assignin(mws, 'PLANT_B_STAR', plant.B_full(2:end));
end

%% ---- 模型说明与保存 ----
note = sprintf(['离线闭环等效性模型（脚本生成，勿手改）\n' ...
    'demo=%s test=%s variant=%s candidate=%s\n' ...
    '参数来源: anc_frozen_params.mat (%s)\n' ...
    '对拍真值: tests/internal/controller_%s.m'], ...
    demo, test, variant, entry.candidate, ...
    load_anc_params('meta').git_describe, demo);
set_param(mdl, 'Description', note);

mdlFile = fullfile(modelsDir, [mdl '.slx']);
save_system(mdl, mdlFile);
fprintf('离线模型已生成: %s\n', mdlFile);
end
