function signals = generate_test_signals(modelName)
% GENERATE_TEST_SIGNALS  生成统一测试信号并预计算开环响应
%
%   对指定被控对象模型生成三组标准测试信号，RMS 归一化后通过被控对象
%   预计算开环响应 y_open。所有结果保存为单个 .mat 文件。
%
%   用法:
%     signals = generate_test_signals()               % 默认 ARMAX 模型
%     signals = generate_test_signals('armax_30303022')
%
%   测试信号:
%     T1 — 固定频率正弦，f₀ = 334 Hz（第一共振峰），10s
%     T2 — 线性扫频，300→500 Hz，10s
%     T3 — 带限白噪声，100–800 Hz（频域法+升余弦窗），10s
%
%   统一参数: Tsim=10s, rng(42), RMS(d)=0.8
%
%   输出文件:
%     dataset/test_signals_<modelName>.mat
%
%   输出结构:
%     signals.model       — 模型名
%     signals.fs          — 采样频率 (Hz)
%     signals.rng_seed    — 随机种子
%     signals.norm_target — 'rms_d'
%     signals.norm_value  — 0.8
%     signals.T1 / .T2 / .T3 — 各含 t, d, y_open (, f_inst for T2)

    if nargin < 1 || isempty(modelName)
        modelName = 'armax_30303022';
    end

    %% ====================================================================
    % §0. 路径与模型加载
    % ====================================================================

    addpath(fullfile(fileparts(mfilename('fullpath')), '..'));
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'dataset'));
    addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'functions'));

    info = DataManager(modelName);

    A = info.model.A;
    B_poly = info.model.B;
    d = info.orders(4);
    fs = info.fs;
    Ts = 1/fs;

    % 去除纯延迟的分子多项式
    B = B_poly(d+1:end);

    % 延迟补偿后的完整分子（用于开环响应仿真）
    b_plant = [zeros(1, d), B];
    a_plant = A;

    fprintf('========== 测试信号生成 ==========\n');
    fprintf('  模型: %s, fs=%d Hz, d=%d\n', modelName, fs, d);
    fprintf('  Tsim=10s, rng(42), RMS(d)=0.8\n\n');

    %% ====================================================================
    % §1. 公共参数
    %% ====================================================================

    Tsim   = 10;
    Nsim   = round(Tsim * fs);
    rngSeed = 42;

    signals = struct();
    signals.model       = modelName;
    signals.fs          = fs;
    signals.rng_seed    = rngSeed;
    signals.norm_target = 'rms_d';
    signals.norm_value  = 0.8;

    %% ====================================================================
    % §2. T1 — 固定频率正弦 @ 334 Hz
    %% ====================================================================
    fprintf('--- T1: 固定频率正弦 (334 Hz) ---\n');

    f0 = 334;
    t1 = (0:Nsim-1)' * Ts;
    rng(rngSeed);
    d1 = sin(2*pi*f0*t1) + 0.02*randn(Nsim, 1);

    % RMS 归一化
    d1 = d1 / rms(d1) * 0.8;

    % 开环响应 = 主噪声本身 (无 ANC 时误差传声器处的声压)
    y_open1 = d1;

    signals.T1 = struct(...
        'type',   'fixed_sine', ...
        'f_hz',   f0, ...
        'fs',     fs, ...
        'Tsim',   Tsim, ...
        't',      t1, ...
        'd',      d1, ...
        'y_open', y_open1 ...
    );

    fprintf('  f₀=%.1f Hz, RMS(d)=%.4f, max|d|=%.4f\n', ...
        f0, rms(d1), max(abs(d1)));

    %% ====================================================================
    % §3. T2 — 线性扫频 300→500 Hz
    %% ====================================================================
    fprintf('--- T2: 线性扫频 (300→500 Hz) ---\n');

    f_start = 300;
    f_end   = 500;
    t2 = (0:Nsim-1)' * Ts;
    rng(rngSeed);

    % 瞬时频率: f(t) = f_start + (f_end - f_start) * t / Tsim
    % 相位 = ∫ 2π·f(t) dt = 2π·(f_start·t + (f_end-f_start)·t²/(2·Tsim))
    f_inst2 = f_start + (f_end - f_start) * t2 / Tsim;
    phase2  = 2*pi * (f_start*t2 + (f_end - f_start)*t2.^2 / (2*Tsim));
    d2 = sin(phase2) + 0.02*randn(Nsim, 1);

    % RMS 归一化
    d2 = d2 / rms(d2) * 0.8;

    % 预计算开环响应
    y_open2 = filter(b_plant, a_plant, d2);

    signals.T2 = struct(...
        'type',    'linear_chirp', ...
        'f_range', [f_start, f_end], ...
        'fs',      fs, ...
        'Tsim',    Tsim, ...
        't',       t2, ...
        'd',       d2, ...
        'f_inst',  f_inst2, ...
        'y_open',  y_open2 ...
    );

    fprintf('  %.0f→%.0f Hz (%.1f Hz/s), RMS(d)=%.4f\n', ...
        f_start, f_end, (f_end-f_start)/Tsim, rms(d2));

    %% ====================================================================
    % §4. T3 — 带限白噪声 100–800 Hz
    %% ====================================================================
    fprintf('--- T3: 带限白噪声 (100–800 Hz) ---\n');

    f_low  = 100;
    f_high = 800;
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
    d3 = real(ifft(Noise_fft_shaped));

    % RMS 归一化
    d3 = d3 / rms(d3) * 0.8;

    % 预计算开环响应
    y_open3 = filter(b_plant, a_plant, d3);

    signals.T3 = struct(...
        'type',    'bandlimited_noise', ...
        'f_range', [f_low, f_high], ...
        'method',  'freq_domain_raised_cosine', ...
        'fs',      fs, ...
        'Tsim',    Tsim, ...
        't',       t3, ...
        'd',       d3, ...
        'y_open',  y_open3 ...
    );

    fprintf('  %.0f–%.0f Hz, 过渡带=%.0f Hz, RMS(d)=%.4f, max|d|=%.4f\n', ...
        f_low, f_high, transition_bw, rms(d3), max(abs(d3)));

    %% ====================================================================
    % §5. 保存
    %% ====================================================================

    % 输出到 dataset/（与模型文件同目录）
    rootDir = fileparts(mfilename('fullpath'));
    dsDir = fullfile(rootDir, '..', 'dataset');
    outFile = fullfile(dsDir, sprintf('test_signals_%s.mat', modelName));

    save(outFile, 'signals');
    fprintf('\n测试信号已保存到: %s\n', outFile);
    fprintf('========== 信号生成完成 ==========\n');
end
