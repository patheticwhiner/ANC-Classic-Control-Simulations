function mdlFile = build_rt_model(demo, test, variant)
%BUILD_RT_MODEL 生成实时部署模型 models/rt_<demo>_<test>_<variant>.slx
%
%   与离线模型共用同一控制器库链接，仅替换 I/O:
%
%     [IO_IN: ADC→V→模型单位] → e ─→ [controller] ─ u ─→ ×MasterGain
%         ─→ ×safety_gate ─→ Sat(±u_limit) ─→ [IO_OUT: 模型单位→V→DAC]
%
%   - rti_available() 为真（dSPACE 主机）: IO 用 RTI1202 ADC/DAC 块，
%     并设置 rti1202.tlc 代码生成目标。
%   - 否则（本 Linux 开发机）: IO_IN/IO_OUT 用 Inport/Outport 占位，
%     端口接口与缩放/饱和链完全一致；在 dSPACE 主机重跑本脚本即换入
%     RTI 块。分支写入模型注释。
%   - 安全层: EnableSwitch / PanicSwitch / LatchReset 常量（ControlDesk
%     可调）+ e 滑动 RMS 看门狗（anc_safety_gate_step）。
%   - ControlDesk 监视信号线命名: e / u / u_demand / clipped / diag /
%     e_rms / watchdog_latched / frame_clock。
%
%   ⚠ 本函数绝不下载/运行到硬件 —— 构建与上板由人工按
%     docs/MICROLABBOX_EXPERIMENT_PROTOCOL.md 执行。

dspaceRoot = fileparts(fileparts(mfilename('fullpath')));
modelsDir = fullfile(dspaceRoot, 'models');
if ~isfolder(modelsDir), mkdir(modelsDir); end

entry = load_anc_params(demo, test, variant);
plant = load_anc_params('plant', 'secondary');
spec = controller_registry(demo, variant);
cal = rt_calibration();
Ts = 1 / plant.fs;
hasRti = rti_available();

libFile = fullfile(modelsDir, 'lib_anc_controllers.slx');
if ~bdIsLoaded('lib_anc_controllers')
    if ~isfile(libFile)
        error('build_rt_model:MissingLibrary', ...
            '未找到 %s，请先运行 build_controller_library。', libFile);
    end
    load_system(libFile);
end

mdl = sprintf('rt_%s_%s_%s', demo, test, variant);
slk_new_model(mdl, 'model', Ts);

%% ---- 输入链 ----
if hasRti
    adcIn = add_rti_adc(mdl, 'adc_error', cal.adc_ch_error, [40 80 120 120]);
    add_block('simulink/Math Operations/Gain', [mdl '/adc_to_volts'], ...
        'Gain', '10', 'Position', [150 85 180 115]);
    add_line(mdl, [adcIn '/1'], 'adc_to_volts/1');
    vSrc = 'adc_to_volts';
else
    add_block('simulink/Sources/In1', [mdl '/e_volts_in'], ...
        'Position', [40 85 80 115]);
    vSrc = 'e_volts_in';
end
add_block('simulink/Math Operations/Gain', [mdl '/e_cal'], ...
    'Gain', num2str(cal.e_per_volt, '%.17g'), 'Position', [210 85 240 115]);
add_line(mdl, [vSrc '/1'], 'e_cal/1');

if spec.has_ref
    if hasRti
        adcRef = add_rti_adc(mdl, 'adc_ref', cal.adc_ch_ref, [40 150 120 190]);
        add_block('simulink/Math Operations/Gain', [mdl '/adc_ref_to_volts'], ...
            'Gain', '10', 'Position', [150 155 180 185]);
        add_line(mdl, [adcRef '/1'], 'adc_ref_to_volts/1');
        refSrc = 'adc_ref_to_volts';
    else
        add_block('simulink/Sources/In1', [mdl '/ref_volts_in'], ...
            'Position', [40 155 80 185]);
        refSrc = 'ref_volts_in';
    end
    add_block('simulink/Math Operations/Gain', [mdl '/ref_cal'], ...
        'Gain', num2str(cal.ref_per_volt, '%.17g'), 'Position', [210 155 240 185]);
    add_line(mdl, [refSrc '/1'], 'ref_cal/1');
end

%% ---- 控制器 + 开关 ----
add_block('simulink/Sources/Constant', [mdl '/EnableSwitch'], ...
    'Value', '1', 'Position', [40 230 100 260]);
add_block('simulink/Sources/Constant', [mdl '/PanicSwitch'], ...
    'Value', '0', 'Position', [40 290 100 320]);
add_block('simulink/Sources/Constant', [mdl '/LatchReset'], ...
    'Value', '0', 'Position', [40 350 100 380]);

add_block(['lib_anc_controllers/' spec.block], [mdl '/controller'], ...
    'Position', [300 80 480 260]);
add_line(mdl, 'e_cal/1', 'controller/1');
if spec.has_ref
    add_line(mdl, 'ref_cal/1', 'controller/2');
    add_line(mdl, 'EnableSwitch/1', 'controller/3');
else
    add_line(mdl, 'EnableSwitch/1', 'controller/2');
end

%% ---- 安全层 ----
sourceFile = fullfile(dspaceRoot, 'controllers', 'anc_safety_gate_step.m');
slk_add_mlfcn(mdl, 'safety_gate', sourceFile, ...
    {'window_samples', 'rms_limit'}, [300 300 440 380]);
add_line(mdl, 'e_cal/1', 'safety_gate/1');
add_line(mdl, 'PanicSwitch/1', 'safety_gate/2');
add_line(mdl, 'LatchReset/1', 'safety_gate/3');

add_block('simulink/Math Operations/Gain', [mdl '/MasterGain'], ...
    'Gain', '1', 'Position', [520 95 560 125]);
add_line(mdl, 'controller/1', 'MasterGain/1');
add_block('simulink/Math Operations/Product', [mdl '/gate_mult'], ...
    'Inputs', '2', 'Position', [590 100 620 130]);
add_line(mdl, 'MasterGain/1', 'gate_mult/1');
add_line(mdl, 'safety_gate/1', 'gate_mult/2');
add_block('simulink/Discontinuities/Saturation', [mdl '/u_sat'], ...
    'UpperLimit', num2str(entry.actuator_limit, '%.17g'), ...
    'LowerLimit', num2str(-entry.actuator_limit, '%.17g'), ...
    'Position', [650 100 680 130]);
add_line(mdl, 'gate_mult/1', 'u_sat/1');

%% ---- 输出链 ----
add_block('simulink/Math Operations/Gain', [mdl '/u_to_volts'], ...
    'Gain', num2str(cal.volt_per_u, '%.17g'), 'Position', [710 100 740 130]);
add_line(mdl, 'u_sat/1', 'u_to_volts/1');
add_block('simulink/Discontinuities/Saturation', [mdl '/dac_volt_sat'], ...
    'UpperLimit', num2str(cal.dac_volt_limit, '%.17g'), ...
    'LowerLimit', num2str(-cal.dac_volt_limit, '%.17g'), ...
    'Position', [770 100 800 130]);
add_line(mdl, 'u_to_volts/1', 'dac_volt_sat/1');
if hasRti
    add_block('simulink/Math Operations/Gain', [mdl '/volts_to_dac'], ...
        'Gain', '0.1', 'Position', [830 100 860 130]);
    add_line(mdl, 'dac_volt_sat/1', 'volts_to_dac/1');
    dacOut = add_rti_dac(mdl, 'dac_control', cal.dac_ch_control, [890 95 970 135]);
    add_line(mdl, 'volts_to_dac/1', [dacOut '/1']);
else
    add_block('simulink/Sinks/Out1', [mdl '/u_volts_out'], ...
        'Position', [830 100 870 130]);
    add_line(mdl, 'dac_volt_sat/1', 'u_volts_out/1');
end

%% ---- 监视信号（帧钟 + Terminator 收尾 + 命名线）----
add_block('simulink/Sources/Digital Clock', [mdl '/frame_clock'], ...
    'SampleTime', num2str(Ts, '%.17g'), 'Position', [40 420 100 450]);
% 未接汇的监视信号先加 Terminator（线存在后才能命名）
termAt = 500;
for sig = {{'controller', 2}, {'controller', 3}, {'controller', 4}, ...
        {'safety_gate', 2}, {'safety_gate', 3}, {'frame_clock', 1}}
    termAt = termAt + 40;
    tname = sprintf('term_%d', termAt);
    add_block('simulink/Sinks/Terminator', [mdl '/' tname], ...
        'Position', [700 termAt 720 termAt+20]);
    add_line(mdl, sprintf('%s/%d', sig{1}{1}, sig{1}{2}), [tname '/1']);
end
name_line(mdl, 'e_cal', 1, 'e');
name_line(mdl, 'controller', 1, 'u');
name_line(mdl, 'controller', 2, 'u_demand');
name_line(mdl, 'controller', 3, 'clipped');
name_line(mdl, 'controller', 4, 'diag');
name_line(mdl, 'safety_gate', 2, 'e_rms');
name_line(mdl, 'safety_gate', 3, 'watchdog_latched');
name_line(mdl, 'frame_clock', 1, 'frame_clock');

%% ---- 模型工作区参数 ----
mws = get_param(mdl, 'ModelWorkspace');
paramMap = spec.params;
for index = 1:size(paramMap, 1)
    assignin(mws, paramMap{index, 1}, entry.(paramMap{index, 2}));
end
assignin(mws, 'window_samples', round(cal.watchdog_window_s * plant.fs));
assignin(mws, 'rms_limit', cal.watchdog_rms_limit);

%% ---- 代码生成目标 ----
if hasRti
    try
        set_param(mdl, 'SystemTargetFile', 'rti1202.tlc');
    catch rtiErr
        warning('build_rt_model:TargetFile', ...
            '设置 rti1202.tlc 失败（%s）—— 请在 dSPACE 主机上确认 RCP 安装。', ...
            rtiErr.message);
    end
    ioNote = 'IO: RTI1202 ADC/DAC（dSPACE 主机构建）';
else
    ioNote = ['IO: Inport/Outport 占位（本机无 RTI）。' ...
        '在 dSPACE 主机重跑 build_rt_model 换入 DS1202 ADC/DAC 块。'];
end

note = sprintf(['实时部署模型（脚本生成，勿手改）\n' ...
    'demo=%s test=%s variant=%s candidate=%s\n%s\n' ...
    '标定增益来自 params/rt_calibration.m（上板 Stage 0/1 实测后更新）\n' ...
    '构建与上板流程: docs/MICROLABBOX_EXPERIMENT_PROTOCOL.md\n' ...
    '参数溯源: anc_frozen_params.mat (%s)'], ...
    demo, test, variant, entry.candidate, ioNote, ...
    load_anc_params('meta').git_describe);
set_param(mdl, 'Description', note);

mdlFile = fullfile(modelsDir, [mdl '.slx']);
save_system(mdl, mdlFile);
fprintf('RT 模型已生成: %s (%s)\n', mdlFile, ...
    pick(hasRti, 'RTI IO', '占位 IO'));
end

%% ========================================================================
function blk = add_rti_adc(mdl, name, channel, position)
% 在 rtilib1202 中定位 Class-1 ADC 块并实例化（块名随 RCP 版本可能变化，
% 用正则搜索并给出清晰错误）。
blk = add_rti_block(mdl, name, sprintf('DS1202ADC_C%d', channel), position);
end

function blk = add_rti_dac(mdl, name, channel, position)
blk = add_rti_block(mdl, name, sprintf('DS1202DAC_C%d', channel), position);
end

function blk = add_rti_block(mdl, name, blockPattern, position)
hits = find_system('rtilib1202', 'LookUnderMasks', 'all', ...
    'Regexp', 'on', 'Name', ['^' blockPattern '$']);
if isempty(hits)
    error('build_rt_model:RtiBlockNotFound', ...
        ['rtilib1202 中未找到块 %s —— 请在 dSPACE 主机上打开 rtilib1202 ' ...
         '确认实际块名并更新 add_rti_block。'], blockPattern);
end
add_block(hits{1}, [mdl '/' name], 'Position', position);
blk = name;
end

function name_line(mdl, srcBlock, srcPort, signalName)
lines = get_param([mdl '/' srcBlock], 'LineHandles');
handle = lines.Outport(srcPort);
if handle ~= -1
    set_param(handle, 'Name', signalName);
end
end

function out = pick(cond, a, b)
if cond, out = a; else, out = b; end
end
