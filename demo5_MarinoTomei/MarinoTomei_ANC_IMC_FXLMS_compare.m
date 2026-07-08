%% Marino-Tomei narrowband ANC vs classic IMC-FxLMS
% A more realistic feedback ANC comparison:
%   e(n) = d_p(n) + S(z)u(n)
%
% Both controllers use the error microphone signal and a secondary-path
% model. There is no external reference microphone in this demo.
%
% - Marino-Tomei: narrowband adaptive internal model. It assumes that the
%   dominant disturbance can be represented by a small number of sinusoidal
%   components and that the sign of Re{S(exp(jw))} is known near those tones.
% - IMC-FxLMS: classic feedback ANC baseline. It reconstructs an internal
%   disturbance reference through an IMC structure and adapts an FIR control
%   filter using the filtered-x LMS rule.

clear; close all; clc;

base_cfg = make_base_config();
scenarios = [ ...
    make_mixed_noise_scenario(base_cfg), ...
    make_narrowband_scenario(base_cfg)];

for scenario_idx = 1:numel(scenarios)
    run_comparison_scenario(scenarios(scenario_idx));
end

%% Local functions

function cfg = make_base_config()
    cfg.fs = 8000;
    cfg.Tsim = 35;
    cfg.control_start = 10;
    cfg.secondary_path = make_realistic_secondary_path(cfg.fs);
    cfg.secondary_est = cfg.secondary_path;

    cfg.modeled_freqs_hz = [180; 260];
    cfg.initial_freqs_hz = [165; 245];
    cfg.tone_amps = [0.75; 0.52];
    cfg.unmodeled_freqs_hz = [420; 735];
    cfg.unmodeled_amps = [0.34; 0.18];
    cfg.broadband_rms = 0.22;
    cfg.tone_mod_depth = 0.13;
    cfg.unmodeled_mod_depth = 0.20;
    cfg.rng_seed = 7;

    cfg.marino.k = 0.008;
    cfg.marino.epsilon = 3e-6;
    cfg.marino.q = numel(cfg.modeled_freqs_hz);

    cfg.imc.mu = 0.012;
    cfg.imc.delta = 1e-4;
    cfg.imc.filter_length = 256;

    cfg.title = 'Realistic mixed-noise feedback ANC comparison';
end

function cfg = make_mixed_noise_scenario(base_cfg)
    cfg = base_cfg;
    cfg.title = 'Realistic mixed-noise feedback ANC comparison';
    cfg.description = 'Modeled tones plus unmodeled tones and broadband noise';
end

function cfg = make_narrowband_scenario(base_cfg)
    cfg = base_cfg;
    cfg.title = 'Narrowband-dominant feedback ANC comparison';
    cfg.description = 'Only the modeled tones plus a very small measurement background';
    cfg.initial_freqs_hz = [178; 258];
    cfg.tone_amps = [0.75; 0.52];
    cfg.unmodeled_freqs_hz = [];
    cfg.unmodeled_amps = [];
    cfg.broadband_rms = 0.015;
    cfg.tone_mod_depth = 0.04;
    cfg.unmodeled_mod_depth = 0;
    cfg.rng_seed = 11;

    % The narrowband case matches Marino-Tomei's information structure:
    % few persistent tones, reliable initial frequency guesses, and little
    % broadband contamination in the adaptation signal.
    cfg.marino.k = 0.026;
    cfg.marino.epsilon = 8e-6;
    cfg.marino.q = numel(cfg.modeled_freqs_hz);

    % Keep IMC-FxLMS as a conservative classic baseline. It still sees only
    % e(n), but it must learn a long FIR controller instead of two oscillators.
    cfg.imc.mu = 0.012;
    cfg.imc.filter_length = 256;
end

function run_comparison_scenario(cfg)
    t = (0:round(cfg.Tsim * cfg.fs) - 1)' / cfg.fs;
    d_primary = make_primary_noise(cfg, t);

    open_loop = struct();
    open_loop.name = 'ANC off';
    open_loop.e = d_primary;
    open_loop.u = zeros(size(d_primary));
    open_loop.freq_hat_hz = [];
    open_loop.adaptive_params = 0;
    open_loop.complexity = make_complexity(0, 0, 0, 0);

    marino = sim_marino_tomei_feedback_anc(cfg, d_primary);
    imc = sim_imc_fxlms_feedback_anc(cfg, d_primary);

    results = [open_loop, marino, imc];
    metrics = summarize_results(cfg, results);

    fprintf('\n=== %s ===\n', cfg.title);
    fprintf('%s\n', cfg.description);
    fprintf('fs = %.0f Hz, Tsim = %.1f s, ANC starts at %.1f s\n', ...
        cfg.fs, cfg.Tsim, cfg.control_start);
    fprintf('Modeled tones: %s Hz\n', mat2str(cfg.modeled_freqs_hz(:)'));
    if isempty(cfg.unmodeled_freqs_hz)
        fprintf('Unmodeled tones: none, broadband RMS %.3f\n', cfg.broadband_rms);
    else
        fprintf('Unmodeled tones: %s Hz, broadband RMS %.3f\n', ...
            mat2str(cfg.unmodeled_freqs_hz(:)'), cfg.broadband_rms);
    end
    fprintf('\n%-28s %-12s %-12s %-12s %-12s %-12s %-12s\n', ...
        'Controller', 'RMS pre', 'RMS tail', 'Supp. dB', ...
        'u RMS tail', 'T20dB(s)', 'Adapt params');
    for i = 1:numel(results)
        fprintf('%-28s %-12.4f %-12.4f %-12.2f %-12.4f %-12s %-12d\n', ...
            results(i).name, metrics(i).rms_pre, metrics(i).rms_tail, ...
            metrics(i).suppression_db, metrics(i).u_rms_tail, ...
            format_time_metric(metrics(i).t20_db_s), results(i).adaptive_params);
    end

    fprintf('\nMarino-Tomei estimated frequencies at the end: %s Hz\n', ...
        mat2str(marino.freq_hat_hz(end, :), 4));
    print_complexity_estimate(results);

    plot_results(cfg, t, results, metrics);
end

function c = make_complexity(adaptive_params, memory_states, mults_per_sample, trans_per_sample)
    c = struct();
    c.adaptive_params = adaptive_params;
    c.memory_states = memory_states;
    c.mults_per_sample = mults_per_sample;
    c.trans_per_sample = trans_per_sample;
end

function print_complexity_estimate(results)
    fprintf('\nApproximate online controller complexity ');
    fprintf('(excluding plant simulation, logging, and plotting):\n');
    fprintf('%-28s %-14s %-14s %-14s %-14s\n', ...
        'Controller', 'Adapt params', 'Memory states', 'Mult/sample', 'Trans/sample');
    for i = 1:numel(results)
        c = results(i).complexity;
        fprintf('%-28s %-14d %-14d %-14d %-14d\n', ...
            results(i).name, c.adaptive_params, c.memory_states, ...
            c.mults_per_sample, c.trans_per_sample);
    end
end

function h = make_realistic_secondary_path(fs)
    % A compact acoustic-like FIR: pure delay, early reflection, and a
    % decaying multi-resonance tail. This keeps the demo reproducible while
    % being less idealized than a short textbook transfer function.
    delay = round(0.0022 * fs);
    tail_len = round(0.045 * fs);
    n = (0:tail_len)';

    early = zeros(delay, 1);
    tail_env = exp(-n / (0.010 * fs));
    tail = tail_env .* ( ...
        0.46 * cos(2 * pi * 520 * n / fs) + ...
        0.30 * cos(2 * pi * 1080 * n / fs + 0.55) + ...
        0.16 * cos(2 * pi * 1740 * n / fs + 1.10));
    h = [early; 0.42; tail(2:end)];
    [H, f] = freqz(h, 1, 4096, fs);
    band = f >= 100 & f <= 900;
    h = h * (0.35 / median(abs(H(band))));
end

function d = make_primary_noise(cfg, t)
    rng(cfg.rng_seed);
    d = zeros(size(t));

    for i = 1:numel(cfg.modeled_freqs_hz)
        amp_env = 1 + cfg.tone_mod_depth * ...
            sin(2 * pi * (0.09 + 0.04 * i) * t + 0.6 * i);
        d = d + cfg.tone_amps(i) * amp_env .* ...
            sin(2 * pi * cfg.modeled_freqs_hz(i) * t + 0.8 * i);
    end

    for i = 1:numel(cfg.unmodeled_freqs_hz)
        amp_env = 1 + cfg.unmodeled_mod_depth * ...
            sin(2 * pi * (0.06 + 0.03 * i) * t + 1.2 * i);
        d = d + cfg.unmodeled_amps(i) * amp_env .* ...
            sin(2 * pi * cfg.unmodeled_freqs_hz(i) * t + 0.4 * i);
    end

    if cfg.broadband_rms > 0
        broadband = make_band_limited_noise(numel(t), cfg.fs, 70, 950);
        broadband = broadband / rms(broadband) * cfg.broadband_rms;
        d = d + broadband;
    end
end

function x = make_band_limited_noise(N, fs, f1, f2)
    xw = randn(N, 1);
    X = fft(xw);
    f = (0:N - 1)' * fs / N;
    mask = (f >= f1 & f <= f2) | (f >= fs - f2 & f <= fs - f1);
    X(~mask) = 0;
    x = real(ifft(X));
end

function out = sim_marino_tomei_feedback_anc(cfg, d_primary)
    fs = cfg.fs;
    N = numel(d_primary);
    q = cfg.marino.q;
    h = cfg.secondary_path(:);
    Ls = numel(h);

    theta = (2 * pi * cfg.initial_freqs_hz(:)' / fs) .^ 2;
    u_hist = zeros(Ls, 1);
    w = zeros(2, q);
    e = zeros(N, 1);
    u = zeros(N, 1);
    freq_hat_hz = zeros(N, q);

    sgn = zeros(q, 1);
    h_idx = (0:Ls - 1)';
    for i = 1:q
        omega = 2 * pi * cfg.modeled_freqs_hz(i) / fs;
        Hw = sum(h .* exp(-1j * omega * h_idx));
        sgn(i) = sign(real(Hw));
        if sgn(i) == 0
            sgn(i) = 1;
        end
    end

    for n = 1:N
        y_secondary = h.' * u_hist;
        e(n) = d_primary(n) + y_secondary;

        if n / fs >= cfg.control_start && n > 1
            w_old = w;
            theta_old = theta;

            for i = 1:q
                alpha = sqrt(max(theta_old(i), 1e-12));
                ca = cos(alpha);
                sa = sin(alpha);

                w(1, i) = ca * w_old(1, i) + (sa / alpha) * w_old(2, i) + ...
                    sgn(i) * cfg.marino.k * (sa / alpha) * e(n);
                w(2, i) = -alpha * sa * w_old(1, i) + ca * w_old(2, i) + ...
                    sgn(i) * cfg.marino.k * (ca - 1) * e(n);

                theta(i) = theta_old(i) + ...
                    cfg.marino.epsilon * sgn(i) * w_old(2, i) * e(n);
                theta(i) = min(max(theta(i), 1e-8), pi ^ 2);
            end

            u(n) = -sum(w(1, :));
        end

        u_hist = [u(n); u_hist(1:end - 1)];
        freq_hat_hz(n, :) = sqrt(theta) * fs / (2 * pi);
    end

    out.name = 'Marino-Tomei narrowband ANC';
    out.e = e;
    out.u = u;
    out.freq_hat_hz = freq_hat_hz;
    out.adaptive_params = 3 * q;
    out.complexity = make_complexity(3 * q, 3 * q, 15 * q, 3 * q);
end

function out = sim_imc_fxlms_feedback_anc(cfg, d_primary)
    N = numel(d_primary);
    fs = cfg.fs;
    h = cfg.secondary_path(:);
    h_hat = cfg.secondary_est(:);
    Ls = numel(h);
    Lw = cfg.imc.filter_length;

    u_hist = zeros(Ls, 1);
    u_hat_hist = zeros(numel(h_hat), 1);
    z_hist = zeros(Lw, 1);
    z_for_filter = zeros(numel(h_hat), 1);
    zf_hist = zeros(Lw, 1);
    w = zeros(Lw, 1);

    e = zeros(N, 1);
    u = zeros(N, 1);

    for n = 1:N
        y_secondary = h.' * u_hist;
        y_secondary_hat = h_hat.' * u_hat_hist;
        e(n) = d_primary(n) + y_secondary;
        z = e(n) - y_secondary_hat;

        z_hist = [z; z_hist(1:end - 1)];
        z_for_filter = [z; z_for_filter(1:end - 1)];
        zf = h_hat.' * z_for_filter;
        zf_hist = [zf; zf_hist(1:end - 1)];

        if n / fs >= cfg.control_start
            u(n) = -w.' * z_hist;
            norm_zf = zf_hist.' * zf_hist + cfg.imc.delta;
            w = w + cfg.imc.mu * e(n) * zf_hist / norm_zf;
        end

        u_hist = [u(n); u_hist(1:end - 1)];
        u_hat_hist = [u(n); u_hat_hist(1:end - 1)];
    end

    out.name = 'Classic IMC-FxLMS';
    out.e = e;
    out.u = u;
    out.freq_hat_hz = [];
    out.adaptive_params = Lw;
    out.complexity = make_complexity(Lw, 2 * numel(h_hat) + 3 * Lw, ...
        2 * numel(h_hat) + 3 * Lw, 0);
end

function metrics = summarize_results(cfg, results)
    fs = cfg.fs;
    idx_pre = 1:round(cfg.control_start * fs);
    idx_tail = round((cfg.Tsim - 8) * fs):numel(results(1).e);
    metrics = repmat(struct('rms_pre', 0, 'rms_tail', 0, ...
        'suppression_db', 0, 'u_rms_tail', 0, 't20_db_s', NaN), size(results));

    ref_pre = rms(results(1).e(idx_pre));
    rms_win = max(1, round(0.20 * fs));
    idx_start = round(cfg.control_start * fs) + 1;
    threshold_20_db = ref_pre / 10;
    for i = 1:numel(results)
        metrics(i).rms_pre = rms(results(i).e(idx_pre));
        metrics(i).rms_tail = rms(results(i).e(idx_tail));
        metrics(i).suppression_db = 20 * log10(ref_pre / metrics(i).rms_tail);
        metrics(i).u_rms_tail = rms(results(i).u(idx_tail));

        env = sqrt(movmean(results(i).e .^ 2, rms_win));
        idx_cross = find(env(idx_start:end) <= threshold_20_db, 1, 'first');
        if ~isempty(idx_cross)
            metrics(i).t20_db_s = (idx_cross - 1) / fs;
        end
    end
end

function s = format_time_metric(t)
    if isnan(t)
        s = 'n/a';
    else
        s = sprintf('%.2f', t);
    end
end

function plot_results(cfg, t, results, metrics)
    fs = cfg.fs;
    colors = lines(numel(results));
    idx_plot = 1:20:numel(t);
    rms_win = max(1, round(0.20 * fs));

    figure('Name', cfg.title, 'Color', 'w');
    tiledlayout(3, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile([1 2]); hold on; grid on;
    for i = 1:numel(results)
        plot(t(idx_plot), results(i).e(idx_plot), ...
            'Color', colors(i, :), 'DisplayName', results(i).name);
    end
    xline(cfg.control_start, ':k', 'ANC on', 'LabelOrientation', 'aligned');
    xlabel('Time (s)'); ylabel('e(n)');
    title(sprintf('%s: error microphone signal', cfg.title));
    legend('Location', 'northeast');

    nexttile; hold on; grid on;
    for i = 1:numel(results)
        env = sqrt(movmean(results(i).e .^ 2, rms_win));
        plot(t(idx_plot), env(idx_plot), ...
            'Color', colors(i, :), 'DisplayName', results(i).name);
    end
    xline(cfg.control_start, ':k', 'ANC on', 'LabelOrientation', 'aligned');
    xlabel('Time (s)'); ylabel('moving RMS');
    title('Error moving RMS');

    nexttile; hold on; grid on;
    for i = 1:numel(results)
        plot(t(idx_plot), results(i).u(idx_plot), ...
            'Color', colors(i, :), 'DisplayName', results(i).name);
    end
    xline(cfg.control_start, ':k', 'ANC on', 'LabelOrientation', 'aligned');
    xlabel('Time (s)'); ylabel('u(n)');
    title('Control loudspeaker signal');

    nexttile; hold on; grid on;
    nfft = 4096;
    tail_idx = round((cfg.Tsim - 8) * fs):numel(t);
    for i = 1:numel(results)
        [Pxx, f] = pwelch(results(i).e(tail_idx), hann(nfft), nfft / 2, nfft, fs);
        plot(f, 10 * log10(Pxx + eps), ...
            'Color', colors(i, :), 'DisplayName', results(i).name);
    end
    xlim([0 1000]);
    xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
    title('Tail error spectrum');

    nexttile; hold on; grid on;
    if ~isempty(results(2).freq_hat_hz)
        for i = 1:size(results(2).freq_hat_hz, 2)
            plot(t(idx_plot), results(2).freq_hat_hz(idx_plot, i), ...
                'LineWidth', 1.0, ...
                'DisplayName', sprintf('hat f_%d', i));
            yline(cfg.modeled_freqs_hz(i), '--', ...
                sprintf('true %.0f Hz', cfg.modeled_freqs_hz(i)), ...
                'HandleVisibility', 'off');
        end
    end
    xline(cfg.control_start, ':k', 'ANC on', 'LabelOrientation', 'aligned');
    xlabel('Time (s)'); ylabel('Frequency (Hz)');
    title('Marino-Tomei frequency estimates');
    legend('Location', 'best');

    figure('Name', [cfg.title ' metrics'], 'Color', 'w');
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

    nexttile; grid on;
    bar(categorical({results.name}), [metrics.rms_tail]);
    ylabel('Tail RMS of e(n)');
    title('Residual error');

    nexttile; grid on;
    bar(categorical({results.name}), [metrics.suppression_db]);
    ylabel('Suppression relative to pre-control (dB)');
    title('Overall suppression');
end
