function signals = load_cylinder1dm_signals(splitName)
%LOAD_CYLINDER1DM_SIGNALS Load a stage calibration/evaluation signal set.
%
%   signals = load_cylinder1dm_signals('calibration'|'evaluation'|'suite')
%
%   默认使用统一的 2 kHz LMS FIR 主/次级路径；优先读取带路径族元数据的
%   dataset/cylinder1dm_2k_stage_suite.mat，缓存不匹配时自动重新生成。

if nargin < 1 || isempty(splitName)
    splitName = 'evaluation';
end

projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
dataFile = fullfile(projectRoot, 'dataset', 'cylinder1dm_2k_stage_suite.mat');
expectedSecondaryModel = 'cylinder1dm_2k_secondary_fir';
expectedPrimaryModel = 'cylinder1dm_2k_primary_fir';
needsGeneration = ~isfile(dataFile);

if ~needsGeneration
    payload = load(dataFile, 'suite');
    needsGeneration = ~isfield(payload, 'suite') ...
        || ~isfield(payload.suite, splitName) ...
        || ~isfield(payload.suite, 'meta') ...
        || ~isfield(payload.suite.meta, 'secondary_model') ...
        || ~isfield(payload.suite.meta, 'primary_model') ...
        || ~strcmp(payload.suite.meta.secondary_model, expectedSecondaryModel) ...
        || ~strcmp(payload.suite.meta.primary_model, expectedPrimaryModel);
end

if needsGeneration
    suite = build_stage_suite(projectRoot, dataFile);
else
    suite = payload.suite;
end

if strcmpi(splitName, 'suite')
    signals = suite;
    return;
end
signals = suite.(splitName);
end

function suite = build_stage_suite(projectRoot, outputFile)
secondaryModel = 'cylinder1dm_2k_secondary_fir';
primaryModel = 'cylinder1dm_2k_primary_fir';
run(fullfile(projectRoot, 'project_init.m'));   % DataManager 需要公共路径

calibration = generate_test_signals(secondaryModel, ...
    'PrimaryModel', primaryModel, 'Duration', 8, 'Seed', 42, ...
    'T1Frequency', 357, 'T1Phase', 0, 'ToneNoiseStd', 0.005, ...
    'T2Range', [300 420], 'T2Direction', 'up', 'T3Band', [100 500], ...
    'SaveFile', fullfile(tempdir, 'cylinder1dm_calibration_signals.mat'));

evaluation = generate_test_signals(secondaryModel, ...
    'PrimaryModel', primaryModel, 'Duration', 10, 'Seed', 142, ...
    'T1Frequency', 357, 'T1Phase', pi/3, 'ToneNoiseStd', 0.005, ...
    'T2Range', [300 420], 'T2Direction', 'down', 'T3Band', [100 500], ...
    'SaveFile', fullfile(tempdir, 'cylinder1dm_evaluation_signals.mat'));

suite = struct();
suite.calibration = calibration;
suite.evaluation = evaluation;
suite.meta = struct( ...
    'secondary_model', secondaryModel, ...
    'primary_model', primaryModel, ...
    'path_family', 'measured_lms_fir', ...
    'target_band_hz', [100 500], ...
    'design_frequency_hz', 357, ...
    'calibration_seed', 42, ...
    'evaluation_seed', 142, ...
    'actuator_limit', 5, ...
    'tuning_peak_limit', 4);
save(outputFile, 'suite');
fprintf('Cylinder stage suite saved: %s\n', outputFile);
end

function signals = generate_test_signals(modelName, varargin)
% GENERATE_TEST_SIGNALS  生成统一测试信号并预计算开环响应
%
%   对指定被控对象模型生成三组标准测试信号，RMS 归一化后通过被控对象
%   预计算开环响应 y_open。所有结果保存为单个 .mat 文件。
%
%   测试信号:
%     T1 — 固定频率正弦，10s
%     T2 — 线性扫频，10s
%     T3 — 带限白噪声（频域法+升余弦窗），10s
%
%   统一参数: Tsim=10s, rng(42), RMS(d)=0.8
%
%   输出结构:
%     signals.model       — 模型名
%     signals.fs          — 采样频率 (Hz)
%     signals.rng_seed    — 随机种子
%     signals.norm_target — 'rms_d'
%     signals.norm_value  — 0.8
%     signals.T1 / .T2 / .T3 — 各含 t, d, y_open (, f_inst for T2)

    if nargin < 1 || isempty(modelName)
        modelName = 'cylinder1dm_2k_secondary_fir';
    end

    parser = inputParser;
    addParameter(parser, 'PrimaryModel', '', @(x) ischar(x) || isstring(x));
    addParameter(parser, 'Duration', 10, @(x) isnumeric(x) && isscalar(x) && x > 0);
    addParameter(parser, 'Seed', 42, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'T1Frequency', [], @(x) isempty(x) || ...
        (isnumeric(x) && isscalar(x) && x > 0));
    addParameter(parser, 'T1Phase', 0, @(x) isnumeric(x) && isscalar(x));
    addParameter(parser, 'ToneNoiseStd', 0.02, @(x) isnumeric(x) && ...
        isscalar(x) && x >= 0);
    addParameter(parser, 'T2Range', [], @(x) isempty(x) || ...
        (isnumeric(x) && numel(x) == 2 && all(x > 0)));
    addParameter(parser, 'T2Direction', 'up', @(x) any(strcmpi(char(x), {'up','down'})));
    addParameter(parser, 'T3Band', [], @(x) isempty(x) || ...
        (isnumeric(x) && numel(x) == 2 && all(x > 0)));
    addParameter(parser, 'SaveFile', '', @(x) ischar(x) || isstring(x));
    parse(parser, varargin{:});
    primaryModelName = char(parser.Results.PrimaryModel);

    %% ====================================================================
    % §0. 模型加载
    % ====================================================================

    info = DataManager(modelName);

    primaryInfo = [];
    if ~isempty(primaryModelName)
        primaryInfo = DataManager(primaryModelName);
        if primaryInfo.fs ~= info.fs
            error('generate_test_signals:SampleRateMismatch', ...
                '主通道 %s (%g Hz) 与次通道 %s (%g Hz) 采样率不一致。', ...
                primaryModelName, primaryInfo.fs, modelName, info.fs);
        end
    end

    fs = info.fs;
    Ts = 1/fs;

    fprintf('========== 测试信号生成 ==========\n');
    fprintf('  次级通道: %s, fs=%d Hz\n', modelName, fs);
    if ~isempty(primaryModelName)
        fprintf('  主通道: %s\n', primaryModelName);
    end
    %% ====================================================================
    % §1. 公共参数
    %% ====================================================================

    Tsim   = parser.Results.Duration;
    Nsim   = round(Tsim * fs);
    rngSeed = parser.Results.Seed;
    fprintf('  Tsim=%.1fs, rng(%d), RMS(d)=0.8\n\n', Tsim, rngSeed);

    signals = struct();
    signals.model       = modelName;
    signals.primary_model = primaryModelName;
    signals.fs          = fs;
    modelB = info.model.B(:).';
    modelOrders = info.orders;
    if strcmp(info.type, 'fir_model')
        % The controller loop applies u(k-1) before the current error sample.
        % Keep the identified FIR taps intact and expose that computation
        % delay explicitly as the leading zero in B_poly.
        modelB = [0, modelB];
        modelOrders = [0, numel(modelB) - 1, 0, 1];
    end
    pathFamily = 'measured_armax';
    if strcmp(info.type, 'fir_model')
        pathFamily = 'measured_lms_fir';
    end
    signals.model_data  = struct( ...
        'A', info.model.A(:).', ...
        'B', modelB, ...
        'orders', modelOrders, ...
        'fs', info.fs, ...
        'path_family', pathFamily);
    signals.rng_seed    = rngSeed;
    signals.norm_target = 'rms_d';
    signals.norm_value  = 0.8;

    %% ====================================================================
    % §2. T1 — 固定频率正弦
    %% ====================================================================
    if isempty(parser.Results.T1Frequency)
        if startsWith(modelName, 'cylinder1dm_2k')
            f0 = 357;
        else
            f0 = 334;
        end
    else
        f0 = parser.Results.T1Frequency;
    end
    f0 = min(f0, 0.45 * fs);
    fprintf('--- T1: 固定频率正弦 (%.0f Hz) ---\n', f0);

    t1 = (0:Nsim-1)' * Ts;
    rng(rngSeed);
    source1 = sin(2*pi*f0*t1 + parser.Results.T1Phase) ...
        + parser.Results.ToneNoiseStd*randn(Nsim, 1);
    d1 = apply_primary_path(source1, primaryInfo);

    % RMS 归一化
    scale1 = 0.8 / rms(d1);
    source1 = source1 * scale1;
    d1 = d1 * scale1;

    % 开环响应 = 主噪声本身 (无 ANC 时误差传声器处的声压)
    y_open1 = d1;

    signals.T1 = struct(...
        'type',   'fixed_sine', ...
        'f_hz',   f0, ...
        'phase_rad', parser.Results.T1Phase, ...
        'fs',     fs, ...
        'Tsim',   Tsim, ...
        't',      t1, ...
        'source', source1, ...
        'd',      d1, ...
        'y_open', y_open1 ...
    );

    fprintf('  f₀=%.1f Hz, RMS(d)=%.4f, max|d|=%.4f\n', ...
        f0, rms(d1), max(abs(d1)));

    %% ====================================================================
    % §3. T2 — 线性扫频
    %% ====================================================================
    if isempty(parser.Results.T2Range)
        if startsWith(modelName, 'cylinder1dm_2k')
            chirpRange = [100, 500];
        else
            chirpRange = [300, 500];
        end
    else
        chirpRange = parser.Results.T2Range(:).';
    end
    chirpRange = min(chirpRange, 0.45 * fs);
    if strcmpi(parser.Results.T2Direction, 'down')
        f_start = max(chirpRange);
        f_end = min(chirpRange);
    else
        f_start = min(chirpRange);
        f_end = max(chirpRange);
    end
    fprintf('--- T2: 线性扫频 (%.0f→%.0f Hz) ---\n', f_start, f_end);

    t2 = (0:Nsim-1)' * Ts;
    rng(rngSeed);

    % 瞬时频率: f(t) = f_start + (f_end - f_start) * t / Tsim
    % 相位 = ∫ 2π·f(t) dt = 2π·(f_start·t + (f_end-f_start)·t²/(2·Tsim))
    f_inst2 = f_start + (f_end - f_start) * t2 / Tsim;
    phase2  = 2*pi * (f_start*t2 + (f_end - f_start)*t2.^2 / (2*Tsim));
    source2 = sin(phase2) + parser.Results.ToneNoiseStd*randn(Nsim, 1);
    d2 = apply_primary_path(source2, primaryInfo);

    % RMS 归一化
    scale2 = 0.8 / rms(d2);
    source2 = source2 * scale2;
    d2 = d2 * scale2;

    y_open2 = d2;

    signals.T2 = struct(...
        'type',    'linear_chirp', ...
        'f_range', [f_start, f_end], ...
        'direction', lower(char(parser.Results.T2Direction)), ...
        'fs',      fs, ...
        'Tsim',    Tsim, ...
        't',       t2, ...
        'source',  source2, ...
        'd',       d2, ...
        'f_inst',  f_inst2, ...
        'y_open',  y_open2 ...
    );

    fprintf('  %.0f→%.0f Hz (%.1f Hz/s), RMS(d)=%.4f\n', ...
        f_start, f_end, (f_end-f_start)/Tsim, rms(d2));

    %% ====================================================================
    % §4. T3 — 带限白噪声
    %% ====================================================================
    if isempty(parser.Results.T3Band)
        if startsWith(modelName, 'cylinder1dm_2k')
            noiseBand = [100, 500];
        else
            noiseBand = [100, 800];
        end
    else
        noiseBand = parser.Results.T3Band(:).';
    end
    f_low  = min(noiseBand);
    f_high = min(max(noiseBand), 0.45 * fs);
    fprintf('--- T3: 带限白噪声 (%.0f–%.0f Hz) ---\n', f_low, f_high);

    t3 = (0:Nsim-1)' * Ts;
    rng(rngSeed);

    % 频域法生成带限噪声: 白噪声 → FFT → 升余弦窗 → IFFT
    noise_raw = randn(Nsim, 1);
    Noise_fft = fft(noise_raw);

    % 频率轴 (只取正半轴用于构造窗函数)
    if mod(Nsim, 2) == 0
        f_axis = (0:Nsim/2)' * fs / Nsim;
        n_nyq = Nsim/2 + 1;
    else
        f_axis = (0:(Nsim-1)/2)' * fs / Nsim;
        n_nyq = (Nsim+1)/2;
    end

    % 升余弦窗: 在 [f_low, f_high] 内为 1，两侧用升余弦过渡
    transition_bw = 15;  % Hz，过渡带宽
    window = zeros(n_nyq, 1);

    for i = 1:n_nyq
        f = f_axis(i);
        if f >= f_low && f <= f_high
            window(i) = 1;
        elseif f > f_high && f <= f_high + transition_bw
            % 高频侧升余弦过渡: 1→0
            window(i) = 0.5 * (1 + cos(pi * (f - f_high) / transition_bw));
        elseif f >= f_low - transition_bw && f < f_low
            % 低频侧升余弦过渡: 0→1
            window(i) = 0.5 * (1 - cos(pi * (f - (f_low - transition_bw)) / transition_bw));
        else
            window(i) = 0;
        end
    end

    % 构造对称的双边频谱窗
    if mod(Nsim, 2) == 0
        % N even: DC + N/2-1 positive freqs + Nyquist + N/2-1 negative freqs (conj)
        full_window = [window; flipud(window(2:end-1))];
    else
        % N odd: DC + (N-1)/2 positive freqs + (N-1)/2 negative freqs (conj)
        full_window = [window; flipud(window(2:end))];
    end

    % 施加频域窗
    Noise_fft_shaped = Noise_fft .* full_window;

    % IFFT 回时域（取实部，因对称窗保实信号）
    source3 = real(ifft(Noise_fft_shaped));
    d3 = apply_primary_path(source3, primaryInfo);

    % RMS 归一化
    scale3 = 0.8 / rms(d3);
    source3 = source3 * scale3;
    d3 = d3 * scale3;

    y_open3 = d3;

    signals.T3 = struct(...
        'type',    'bandlimited_noise', ...
        'f_range', [f_low, f_high], ...
        'method',  'freq_domain_raised_cosine', ...
        'fs',      fs, ...
        'Tsim',    Tsim, ...
        't',       t3, ...
        'source',  source3, ...
        'd',       d3, ...
        'y_open',  y_open3 ...
    );

    fprintf('  %.0f–%.0f Hz, 过渡带=%.0f Hz, RMS(d)=%.4f, max|d|=%.4f\n', ...
        f_low, f_high, transition_bw, rms(d3), max(abs(d3)));

    %% ====================================================================
    % §5. 保存
    %% ====================================================================

    % 默认输出到 dataset/（与模型文件同目录）
    projectRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
    if isempty(parser.Results.SaveFile)
        outFile = fullfile(projectRoot, 'dataset', ...
            sprintf('test_signals_%s.mat', modelName));
    else
        outFile = char(parser.Results.SaveFile);
        if ~isfolder(fileparts(outFile))
            mkdir(fileparts(outFile));
        end
    end

    save(outFile, 'signals');
    fprintf('\n测试信号已保存到: %s\n', outFile);
    fprintf('========== 信号生成完成 ==========\n');
end

function y = apply_primary_path(x, primaryInfo)
    if isempty(primaryInfo)
        y = x;
        return;
    end

    y = filter(primaryInfo.model.B, primaryInfo.model.A, x);
end
