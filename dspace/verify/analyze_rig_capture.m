function result = analyze_rig_capture(captureFile, openLoopFile)
%ANALYZE_RIG_CAPTURE ControlDesk 录波 → 与仿真报告同口径的指标
%
%   result = analyze_rig_capture('rig_demo1_T1_fixed_20260801.mat', ...
%                                'rig_openloop_T1_20260801.mat')
%
%   录波约定（见实验规程 Stage 4）: .mat 内含变量
%     t, e, u, u_demand, clipped, frame_clock  （2 kHz，列向量或行向量）
%   开环基线文件: 同一噪声工况、控制器 disable 时的 e（即 y_open）。
%
%   输出与 tests/internal/compute_metrics 相同的指标结构
%   （supp_db / u_max / 饱和计数 / 需求峰值），外加丢帧检查。
%   文件名约定 rig_<demo>_<test>_<variant>_<date>.mat 用于解析 meta。

cap = load(captureFile);
ol = load(openLoopFile);

required = {'t', 'e', 'u', 'u_demand', 'frame_clock'};
for index = 1:numel(required)
    if ~isfield(cap, required{index})
        error('analyze_rig_capture:MissingVariable', ...
            '录波 %s 缺少变量 %s（检查 ControlDesk 记录配置）。', ...
            captureFile, required{index});
    end
end

% ---- 丢帧检查（frame_clock 应严格步进 Ts）----
fc = cap.frame_clock(:);
dfc = diff(fc);
Ts = median(dfc);
fs = 1 / Ts;
dropped = sum(dfc > 1.5 * Ts);
if dropped > 0
    warning('analyze_rig_capture:DroppedFrames', ...
        '检测到 %d 处帧间隔异常（>1.5 Ts）—— 指标可信度受损，须在报告披露。', ...
        dropped);
end

% ---- 从文件名解析 meta ----
[~, stem] = fileparts(captureFile);
parts = strsplit(stem, '_');
if numel(parts) >= 4
    meta = struct('demo', parts{2}, 'test', parts{3}, 'variant', parts{4});
else
    meta = struct('demo', 'rig', 'test', 'T1', 'variant', 'fixed');
end

% ---- 对齐长度并计算与仿真同口径的指标 ----
n = min(numel(cap.e), numel(ol.e));
y_closed = cap.e(1:n);  y_closed = y_closed(:).';
y_open = ol.e(1:n);     y_open = y_open(:).';
u = cap.u(1:n);         u = u(:).';

test_struct = struct('fs', fs, 'Tsim', n / fs, 'd', y_open(:));
result = compute_metrics(y_open, y_closed, u, test_struct, meta);

result.extra = struct();
result.extra.u_demand_max = max(abs(cap.u_demand(:)));
if isfield(cap, 'clipped')
    result.extra.saturation_count = sum(cap.clipped(:) > 0.5);
else
    result.extra.saturation_count = NaN;
end
result.extra.dropped_frames = dropped;
result.extra.fs_measured = fs;
result.extra.capture_file = captureFile;
result.extra.openloop_file = openLoopFile;

fprintf('%s/%s/%s (台架): supp=%.2f dB, max|u_demand|=%.3f, sat=%d, 丢帧=%d\n', ...
    meta.demo, meta.test, meta.variant, result.supp_db, ...
    result.extra.u_demand_max, result.extra.saturation_count, dropped);
end
