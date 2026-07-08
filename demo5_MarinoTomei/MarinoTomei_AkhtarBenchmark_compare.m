%% Marino-Tomei ANC against the Akhtar-style IMC-FxLMS benchmark
% This script mirrors the scenario set and metrics from:
%   D:\Projects\NANC_Simulations\2012_IMCFxLMS_AkhtarRep\IMCFxLMS_benchmark_suite.m
%
% The goal is not to make Marino-Tomei a universal replacement for
% IMC-FxLMS. Instead, it checks where a low-order adaptive internal model can
% beat or match the fixed IMC-FxLMS baseline, and where its narrowband
% assumptions break down.

clear; clc; close all;

rng(1);

root_dir = fileparts(mfilename('fullpath'));
out_dir = fullfile(root_dir, 'results_AkhtarBenchmark_compare');
if ~exist(out_dir, 'dir')
    mkdir(out_dir);
end

cfg = default_benchmark_config();
scenarios = make_benchmark_scenarios();

imc_results = repmat(empty_metric(), numel(scenarios), 1);
marino_results = repmat(empty_metric(), numel(scenarios), 1);
series_bank = cell(numel(scenarios), 1);

fprintf('\nMarino-Tomei vs Akhtar-style IMC-FxLMS benchmark\n');
fprintf('fs = %.0f Hz, duration = %.2f s, secondaryPath = %s\n', ...
    cfg.fs, cfg.duration, mat2str(cfg.secondary_path));
fprintf('IMC-FxLMS: stepSize = %.3g, nQ = %d\n', cfg.imc_step_size, cfg.nQ);
fprintf('Marino-Tomei: k = %.3g, epsilon = %.3g, q = 1\n\n', ...
    cfg.marino_k, cfg.marino_epsilon);

for idx = 1:numel(scenarios)
    scenario = scenarios(idx);
    fprintf('[%02d/%02d] %-22s ... ', idx, numel(scenarios), scenario.name);
    [imc_results(idx), marino_results(idx), series_bank{idx}] = ...
        run_scenario_pair(cfg, scenario);
    advantage_db = comparable_advantage(imc_results(idx), marino_results(idx));
    fprintf('IMC %7.2f dB | Marino %7.2f dB | delta %s\n', ...
        imc_results(idx).steadySuppressionDb, ...
        marino_results(idx).steadySuppressionDb, format_delta(advantage_db));
end

summary = make_summary_table(imc_results, marino_results);
disp(summary(:, {'name', 'kind', 'marinoPriorHz', ...
    'caseARatio', 'imcSteadyDb', 'marinoSteadyDb', 'marinoAdvantageDb', ...
    'imcConvSec', 'marinoConvSec', 'imcControlRms', ...
    'marinoControlRms', 'imcStatus', 'marinoStatus'}));

writetable(summary, fullfile(out_dir, 'AkhtarBenchmark_Marino_vs_IMCFxLMS.csv'));
save(fullfile(out_dir, 'AkhtarBenchmark_Marino_vs_IMCFxLMS.mat'), ...
    'cfg', 'scenarios', 'imc_results', 'marino_results', 'summary', 'series_bank');

plot_metric_summary(summary, out_dir);
plot_time_series_summary(cfg, scenarios, series_bank, imc_results, marino_results, out_dir);
write_markdown_report(cfg, summary, out_dir);

fprintf('\nSaved benchmark outputs to:\n%s\n', out_dir);

%% Local functions

function cfg = default_benchmark_config()
    cfg = struct();
    cfg.fs = 10000;
    cfg.duration = 3;
    cfg.ts = 1 / cfg.fs;
    cfg.t = 0:cfg.ts:cfg.duration - cfg.ts;

    cfg.secondary_path = [0 0 0.5 0.3];
    cfg.imc_step_size = 1e-2;
    cfg.nQ = 2^8 - 1;

    cfg.marino_k = 0.018;
    cfg.marino_epsilon = 6.0e-6;
    cfg.marino_theta_min = 1e-8;
    cfg.marino_theta_max = pi^2;
    cfg.case_a_min_ratio = 0.15;

    cfg.initial_window_sec = 0.25;
    cfg.final_window_sec = 0.50;
    cfg.rms_window_sec = 0.05;
    cfg.convergence_hold_sec = 0.20;
    cfg.convergence_tolerance = 1.10;
    cfg.diverge_limit = 1e6;
    cfg.high_energy_rms = 20;
end

function scenarios = make_benchmark_scenarios()
    template = struct( ...
        'name', '', ...
        'kind', '', ...
        'frequencyHz', NaN, ...
        'startHz', NaN, ...
        'endHz', NaN, ...
        'bandHz', [NaN NaN], ...
        'amplitude', 1.0, ...
        'note', '');

    scenarios = repmat(template, 10, 1);
    scenarios(1) = set_scenario(template, 'Fixed-600Hz', 'fixed', 600, NaN, NaN, [NaN NaN], ...
        'Low-frequency single-tone baseline.');
    scenarios(2) = set_scenario(template, 'Fixed-1000Hz', 'fixed', 1000, NaN, NaN, [NaN NaN], ...
        'Nominal single-tone baseline.');
    scenarios(3) = set_scenario(template, 'Fixed-1400Hz', 'fixed', 1400, NaN, NaN, [NaN NaN], ...
        'Higher-frequency single-tone stress point.');
    scenarios(4) = set_scenario(template, 'Mismatch-plus-1pct', 'mismatch', 1010, NaN, NaN, [NaN NaN], ...
        'Actual tone is +1 percent away from nominal 1 kHz.');
    scenarios(5) = set_scenario(template, 'Mismatch-plus-3pct', 'mismatch', 1030, NaN, NaN, [NaN NaN], ...
        'Actual tone is +3 percent away from nominal 1 kHz.');
    scenarios(6) = set_scenario(template, 'Mismatch-minus-3pct', 'mismatch', 970, NaN, NaN, [NaN NaN], ...
        'Actual tone is -3 percent away from nominal 1 kHz.');
    scenarios(7) = set_scenario(template, 'Chirp-800-1200Hz', 'chirp', NaN, 800, 1200, [NaN NaN], ...
        'Slowly time-varying single component around nominal 1 kHz.');
    scenarios(8) = set_scenario(template, 'Chirp-600-1400Hz', 'chirp', NaN, 600, 1400, [NaN NaN], ...
        'Wider slow sweep used as a weak-scenario candidate.');
    scenarios(9) = set_scenario(template, 'Bandpass-700-1300Hz', 'bandlimited', NaN, NaN, NaN, [700 1300], ...
        'Moderate-width random disturbance.');
    scenarios(10) = set_scenario(template, 'Bandpass-500-1500Hz', 'bandlimited', NaN, NaN, NaN, [500 1500], ...
        'Wide-band random disturbance used as a weak-scenario candidate.');
end

function scenario = set_scenario(template, name, kind, frequency_hz, start_hz, end_hz, band_hz, note)
    scenario = template;
    scenario.name = name;
    scenario.kind = kind;
    scenario.frequencyHz = frequency_hz;
    scenario.startHz = start_hz;
    scenario.endHz = end_hz;
    scenario.bandHz = band_hz;
    scenario.note = note;
end

function metric = empty_metric()
    metric = struct( ...
        'name', '', ...
        'kind', '', ...
        'frequencyHz', NaN, ...
        'bandLowHz', NaN, ...
        'bandHighHz', NaN, ...
        'marinoPriorHz', NaN, ...
        'marinoFinalHz', NaN, ...
        'caseARatio', NaN, ...
        'overallSuppressionDb', NaN, ...
        'steadySuppressionDb', NaN, ...
        'initialErrorRms', NaN, ...
        'steadyErrorRms', NaN, ...
        'disturbanceRms', NaN, ...
        'convergenceTimeSec', NaN, ...
        'controlEnergy', NaN, ...
        'controlRms', NaN, ...
        'controlPeak', NaN, ...
        'statusRank', 0, ...
        'status', '');
end

function advantage_db = comparable_advantage(imc_metric, marino_metric)
    if strcmp(imc_metric.status, 'stable') && strcmp(marino_metric.status, 'stable')
        advantage_db = marino_metric.steadySuppressionDb - imc_metric.steadySuppressionDb;
    else
        advantage_db = NaN;
    end
end

function text = format_delta(value)
    if isnan(value)
        text = '   n/a';
    else
        text = sprintf('%+7.2f dB', value);
    end
end

function [imc_metric, marino_metric, series] = run_scenario_pair(cfg, scenario)
    input_signal = generate_input_signal(cfg, scenario);
    disturbance = filter(cfg.secondary_path, 1, input_signal);

    imc = sim_imc_fxlms(cfg, disturbance);
    marino_prior_hz = scenario_prior_frequency(scenario);
    case_a_ratio = scenario_case_a_ratio(cfg, scenario);
    if case_a_ratio < cfg.case_a_min_ratio
        marino = skipped_marino_result(disturbance, marino_prior_hz, case_a_ratio);
    else
        marino = sim_marino_tomei(cfg, disturbance, marino_prior_hz);
        marino.caseARatio = case_a_ratio;
    end

    imc_metric = compute_metrics(cfg, scenario, disturbance, imc.error, imc.control, ...
        imc.status, imc.statusRank);
    marino_metric = compute_metrics(cfg, scenario, disturbance, marino.error, marino.control, ...
        marino.status, marino.statusRank);

    imc_metric.marinoPriorHz = NaN;
    imc_metric.marinoFinalHz = NaN;
    imc_metric.caseARatio = NaN;
    marino_metric.marinoPriorHz = marino_prior_hz;
    marino_metric.marinoFinalHz = marino.freqHatHz(end);
    marino_metric.caseARatio = marino.caseARatio;

    series = struct();
    series.input = input_signal;
    series.disturbance = disturbance;
    series.imc = imc;
    series.marino = marino;
end

function prior_hz = scenario_prior_frequency(scenario)
    switch lower(scenario.kind)
        case 'fixed'
            prior_hz = scenario.frequencyHz;
        case 'mismatch'
            prior_hz = 1000;
        case 'chirp'
            prior_hz = scenario.startHz;
        case 'bandlimited'
            prior_hz = mean(scenario.bandHz);
        otherwise
            error('Unknown scenario kind: %s', scenario.kind);
    end
end

function ratio = scenario_case_a_ratio(cfg, scenario)
    switch lower(scenario.kind)
        case 'fixed'
            freq_grid = scenario.frequencyHz;
        case 'mismatch'
            freq_grid = linspace(min(1000, scenario.frequencyHz), ...
                max(1000, scenario.frequencyHz), 41);
        case 'chirp'
            freq_grid = linspace(scenario.startHz, scenario.endHz, 101);
        case 'bandlimited'
            freq_grid = linspace(scenario.bandHz(1), scenario.bandHz(2), 101);
        otherwise
            error('Unknown scenario kind: %s', scenario.kind);
    end

    B = cfg.secondary_path(:).';
    path_idx = 0:numel(B) - 1;
    ratios = zeros(size(freq_grid));
    for i = 1:numel(freq_grid)
        omega = 2 * pi * freq_grid(i) / cfg.fs;
        Hw = sum(B .* exp(-1j * omega * path_idx));
        ratios(i) = abs(real(Hw)) / max(abs(Hw), eps);
    end
    ratio = min(ratios);
end

function out = skipped_marino_result(disturbance, prior_hz, case_a_ratio)
    n_samples = numel(disturbance);
    out = struct();
    out.error = disturbance;
    out.control = zeros(size(disturbance));
    out.antinoise = zeros(size(disturbance));
    out.freqHatHz = prior_hz * ones(1, n_samples);
    out.caseARatio = case_a_ratio;
    out.status = 'caseA_not_applicable';
    out.statusRank = 2;
end

function out = sim_imc_fxlms(cfg, disturbance)
    B = cfg.secondary_path;
    step_size = cfg.imc_step_size;
    nQ = cfg.nQ;
    n_samples = numel(disturbance);

    u = zeros(1, n_samples);
    x = zeros(1, n_samples);
    y = zeros(1, n_samples);
    filtered_x = zeros(1, n_samples);
    antinoise = zeros(1, n_samples);

    Q = zeros(nQ + 1, 1);
    buffer_length = max(numel(B), nQ + 1);
    vec_u = zeros(1, buffer_length);
    vec_x = zeros(1, buffer_length);
    vec_y = zeros(1, buffer_length);
    vec_filtered_x = zeros(1, buffer_length);

    status = 'stable';
    status_rank = 0;

    for k = 1:n_samples
        antinoise(k) = my_fir(B, vec_u);
        y(k) = disturbance(k) - antinoise(k);
        vec_y = [y(k), vec_y(1:end - 1)];

        x(k) = my_fir(B, vec_y) + y(k);
        vec_x = [x(k), vec_x(1:end - 1)];

        filtered_x(k) = my_fir(B, vec_x);
        vec_filtered_x = [filtered_x(k), vec_filtered_x(1:end - 1)];

        Q = Q + (step_size * vec_filtered_x(1:numel(Q)) * y(k))';
        u(k) = my_fir(Q', vec_x);
        vec_u = [u(k), vec_u(1:end - 1)];

        if ~all(isfinite([u(k), y(k), Q(:)'])) || ...
                abs(u(k)) > cfg.diverge_limit || abs(y(k)) > cfg.diverge_limit
            status = 'diverged';
            status_rank = 3;
            u(k + 1:end) = NaN;
            y(k + 1:end) = NaN;
            antinoise(k + 1:end) = NaN;
            break;
        end
    end

    out = struct();
    out.error = y;
    out.control = u;
    out.antinoise = antinoise;
    out.filteredX = filtered_x;
    out.status = status;
    out.statusRank = status_rank;
end

function out = sim_marino_tomei(cfg, disturbance, prior_hz)
    B = cfg.secondary_path(:).';
    fs = cfg.fs;
    n_samples = numel(disturbance);
    n_path = numel(B);

    theta = (2 * pi * prior_hz / fs)^2;
    w = zeros(2, 1);
    u_hist = zeros(1, n_path);
    u = zeros(1, n_samples);
    y = zeros(1, n_samples);
    antinoise = zeros(1, n_samples);
    freq_hat_hz = zeros(1, n_samples);

    path_idx = 0:n_path - 1;
    omega = 2 * pi * prior_hz / fs;
    Hw = sum(B .* exp(-1j * omega * path_idx));
    sgn = sign(real(Hw));
    if sgn == 0
        sgn = 1;
    end

    status = 'stable';
    status_rank = 0;

    for k = 1:n_samples
        antinoise(k) = my_fir(B, u_hist);
        y(k) = disturbance(k) - antinoise(k);

        if k > 1
            w_old = w;
            theta_old = theta;
            alpha = sqrt(max(theta_old, 1e-12));
            ca = cos(alpha);
            sa = sin(alpha);

            w(1) = ca * w_old(1) + (sa / alpha) * w_old(2) + ...
                sgn * cfg.marino_k * (sa / alpha) * y(k);
            w(2) = -alpha * sa * w_old(1) + ca * w_old(2) + ...
                sgn * cfg.marino_k * (ca - 1) * y(k);

            theta = theta_old + cfg.marino_epsilon * sgn * w_old(2) * y(k);
            theta = min(max(theta, cfg.marino_theta_min), cfg.marino_theta_max);
            u(k) = w(1);
        end

        u_hist = [u(k), u_hist(1:end - 1)];
        freq_hat_hz(k) = sqrt(theta) * fs / (2 * pi);

        if ~all(isfinite([u(k), y(k), w(:)', theta])) || ...
                abs(u(k)) > cfg.diverge_limit || abs(y(k)) > cfg.diverge_limit
            status = 'diverged';
            status_rank = 3;
            u(k + 1:end) = NaN;
            y(k + 1:end) = NaN;
            antinoise(k + 1:end) = NaN;
            freq_hat_hz(k + 1:end) = NaN;
            break;
        end
    end

    out = struct();
    out.error = y;
    out.control = u;
    out.antinoise = antinoise;
    out.freqHatHz = freq_hat_hz;
    out.caseARatio = abs(real(Hw)) / max(abs(Hw), eps);
    out.status = status;
    out.statusRank = status_rank;
end

function y = my_fir(coeffs, buffer)
    n = min(numel(coeffs), numel(buffer));
    y = coeffs(1:n) * buffer(1:n).';
end

function x = generate_input_signal(cfg, scenario)
    t = cfg.t;
    switch lower(scenario.kind)
        case {'fixed', 'mismatch'}
            x = scenario.amplitude * sin(2 * pi * scenario.frequencyHz * t);
        case 'chirp'
            sweep_rate = (scenario.endHz - scenario.startHz) / cfg.duration;
            phase = 2 * pi * (scenario.startHz * t + 0.5 * sweep_rate * t.^2);
            x = scenario.amplitude * sin(phase);
        case 'bandlimited'
            x = bandlimited_noise(numel(t), cfg.fs, scenario.bandHz);
            x = scenario.amplitude * x;
        otherwise
            error('Unknown scenario kind: %s', scenario.kind);
    end
end

function x = bandlimited_noise(n_samples, fs, band_hz)
    raw = randn(1, n_samples);
    spectrum = fft(raw);
    freq = (0:n_samples - 1) * fs / n_samples;
    mirror_freq = fs - freq;
    keep = (freq >= band_hz(1) & freq <= band_hz(2)) | ...
        (mirror_freq >= band_hz(1) & mirror_freq <= band_hz(2));
    spectrum(~keep) = 0;
    x = real(ifft(spectrum));
    x = x - mean(x);
    rms_value = sqrt(mean(x.^2));
    if rms_value > 0
        x = x / rms_value / sqrt(2);
    end
end

function metrics = compute_metrics(cfg, scenario, disturbance, error_signal, control_signal, status, status_rank)
    metrics = empty_metric();
    metrics.name = scenario.name;
    metrics.kind = scenario.kind;

    switch lower(scenario.kind)
        case {'fixed', 'mismatch'}
            metrics.frequencyHz = scenario.frequencyHz;
        case 'chirp'
            metrics.frequencyHz = mean([scenario.startHz, scenario.endHz]);
            metrics.bandLowHz = scenario.startHz;
            metrics.bandHighHz = scenario.endHz;
        case 'bandlimited'
            metrics.frequencyHz = mean(scenario.bandHz);
            metrics.bandLowHz = scenario.bandHz(1);
            metrics.bandHighHz = scenario.bandHz(2);
    end

    valid = isfinite(error_signal) & isfinite(control_signal);
    disturbance_valid = disturbance(valid);
    error_valid = error_signal(valid);
    control_valid = control_signal(valid);

    if isempty(error_valid)
        metrics.status = 'diverged';
        metrics.statusRank = 3;
        return;
    end

    n_samples = numel(error_valid);
    initial_n = max(1, min(n_samples, round(cfg.initial_window_sec * cfg.fs)));
    final_n = max(1, min(n_samples, round(cfg.final_window_sec * cfg.fs)));
    final_idx = n_samples - final_n + 1:n_samples;

    metrics.overallSuppressionDb = power_ratio_db(disturbance_valid, error_valid);
    metrics.steadySuppressionDb = power_ratio_db(disturbance_valid(final_idx), error_valid(final_idx));
    metrics.initialErrorRms = rms_local(error_valid(1:initial_n));
    metrics.steadyErrorRms = rms_local(error_valid(final_idx));
    metrics.disturbanceRms = rms_local(disturbance_valid(final_idx));
    metrics.convergenceTimeSec = estimate_convergence_time(cfg, error_valid);
    metrics.controlEnergy = sum(control_valid.^2) / cfg.fs;
    metrics.controlRms = rms_local(control_valid);
    metrics.controlPeak = max(abs(control_valid));
    metrics.status = status;
    metrics.statusRank = status_rank;

    if strcmp(status, 'stable') && isnan(metrics.convergenceTimeSec)
        metrics.status = 'no_convergence';
        metrics.statusRank = 2;
    end
    if strcmp(status, 'stable') && metrics.controlRms > cfg.high_energy_rms
        metrics.status = 'high_control_energy';
        metrics.statusRank = 1;
    end
end

function value_db = power_ratio_db(reference_signal, residual_signal)
    value_db = 10 * log10(sum(reference_signal.^2) / max(sum(residual_signal.^2), eps));
end

function value = rms_local(x)
    value = sqrt(mean(x.^2));
end

function t_conv = estimate_convergence_time(cfg, error_signal)
    n_samples = numel(error_signal);
    if n_samples < 3
        t_conv = NaN;
        return;
    end

    win = max(3, min(n_samples, round(cfg.rms_window_sec * cfg.fs)));
    hold_n = max(1, min(n_samples, round(cfg.convergence_hold_sec * cfg.fs)));
    final_n = max(1, min(n_samples, round(cfg.final_window_sec * cfg.fs)));

    moving_rms = sqrt(movmean(error_signal.^2, win));
    steady_rms = rms_local(error_signal(n_samples - final_n + 1:n_samples));
    initial_rms = rms_local(error_signal(1:min(n_samples, final_n)));

    if steady_rms >= 0.98 * initial_rms
        t_conv = NaN;
        return;
    end

    threshold = cfg.convergence_tolerance * steady_rms;
    t_conv = NaN;
    for idx = 1:(n_samples - hold_n + 1)
        if all(moving_rms(idx:idx + hold_n - 1) <= threshold)
            t_conv = (idx - 1) / cfg.fs;
            return;
        end
    end
end

function summary = make_summary_table(imc_results, marino_results)
    n = numel(imc_results);
    name = strings(n, 1);
    kind = strings(n, 1);
    marino_prior_hz = zeros(n, 1);
    marino_final_hz = zeros(n, 1);
    case_a_ratio = zeros(n, 1);
    imc_steady_db = zeros(n, 1);
    marino_steady_db = zeros(n, 1);
    marino_advantage_db = zeros(n, 1);
    imc_conv_sec = zeros(n, 1);
    marino_conv_sec = zeros(n, 1);
    imc_control_rms = zeros(n, 1);
    marino_control_rms = zeros(n, 1);
    imc_status = strings(n, 1);
    marino_status = strings(n, 1);

    for i = 1:n
        name(i) = string(imc_results(i).name);
        kind(i) = string(imc_results(i).kind);
        marino_prior_hz(i) = marino_results(i).marinoPriorHz;
        marino_final_hz(i) = marino_results(i).marinoFinalHz;
        case_a_ratio(i) = marino_results(i).caseARatio;
        imc_steady_db(i) = imc_results(i).steadySuppressionDb;
        marino_steady_db(i) = marino_results(i).steadySuppressionDb;
        marino_advantage_db(i) = comparable_advantage(imc_results(i), marino_results(i));
        imc_conv_sec(i) = imc_results(i).convergenceTimeSec;
        marino_conv_sec(i) = marino_results(i).convergenceTimeSec;
        imc_control_rms(i) = imc_results(i).controlRms;
        marino_control_rms(i) = marino_results(i).controlRms;
        imc_status(i) = string(imc_results(i).status);
        marino_status(i) = string(marino_results(i).status);
    end

    summary = table(name, kind, marino_prior_hz, marino_final_hz, case_a_ratio, ...
        imc_steady_db, marino_steady_db, marino_advantage_db, ...
        imc_conv_sec, marino_conv_sec, imc_control_rms, ...
        marino_control_rms, imc_status, marino_status, ...
        'VariableNames', {'name', 'kind', 'marinoPriorHz', 'marinoFinalHz', 'caseARatio', ...
        'imcSteadyDb', 'marinoSteadyDb', 'marinoAdvantageDb', ...
        'imcConvSec', 'marinoConvSec', 'imcControlRms', ...
        'marinoControlRms', 'imcStatus', 'marinoStatus'});
end

function plot_metric_summary(summary, out_dir)
    fig = figure('Name', 'Marino vs IMC-FxLMS Akhtar benchmark metrics', ...
        'NumberTitle', 'off', 'Position', [120 120 1300 700]);

    names = categorical(summary.name);
    subplot(2, 2, 1);
    bar(names, [summary.imcSteadyDb, summary.marinoSteadyDb]);
    grid on; title('Steady suppression'); ylabel('dB');
    legend({'IMC-FxLMS', 'Marino-Tomei'}, 'Location', 'best');
    xtickangle(35);

    subplot(2, 2, 2);
    bar(names, summary.marinoAdvantageDb);
    grid on; title('Marino advantage over IMC'); ylabel('dB');
    yline(0, ':k');
    xtickangle(35);

    subplot(2, 2, 3);
    bar(names, [summary.imcControlRms, summary.marinoControlRms]);
    grid on; title('Control RMS'); ylabel('amplitude');
    legend({'IMC-FxLMS', 'Marino-Tomei'}, 'Location', 'best');
    xtickangle(35);

    subplot(2, 2, 4);
    bar(names, [summary.imcConvSec, summary.marinoConvSec]);
    grid on; title('Convergence time'); ylabel('s');
    legend({'IMC-FxLMS', 'Marino-Tomei'}, 'Location', 'best');
    xtickangle(35);

    saveas(fig, fullfile(out_dir, 'metric_summary.png'));
end

function plot_time_series_summary(cfg, scenarios, series_bank, imc_results, marino_results, out_dir)
    fig = figure('Name', 'Marino vs IMC-FxLMS Akhtar benchmark time series', ...
        'NumberTitle', 'off', 'Position', [50 50 1500 900]);
    rows = 5;
    cols = 2;
    for idx = 1:numel(scenarios)
        subplot(rows, cols, idx);
        series = series_bank{idx};
        plot(cfg.t, series.disturbance, 'Color', [0.20 0.40 0.85], 'LineWidth', 0.6);
        hold on;
        plot(cfg.t, series.imc.error, 'Color', [0.85 0.20 0.20], 'LineWidth', 0.6);
        plot(cfg.t, series.marino.error, 'Color', [0.05 0.55 0.25], 'LineWidth', 0.6);
        grid on;
        title(sprintf('%s | IMC %.1f dB | MT %.1f dB', scenarios(idx).name, ...
            imc_results(idx).steadySuppressionDb, marino_results(idx).steadySuppressionDb), ...
            'Interpreter', 'none');
        xlabel('Time (s)');
        ylabel('Amplitude');
        xlim([cfg.t(1), cfg.t(end)]);
        if idx == 1
            legend({'disturbance', 'IMC error', 'Marino error'}, 'Location', 'best');
        end
    end
    saveas(fig, fullfile(out_dir, 'time_domain_summary.png'));
end

function write_markdown_report(cfg, summary, out_dir)
    report_path = fullfile(out_dir, 'AkhtarBenchmark_Marino_vs_IMCFxLMS.md');
    fid = fopen(report_path, 'w', 'n', 'UTF-8');
    cleanup = onCleanup(@() fclose(fid));

    fprintf(fid, '# Marino-Tomei 与 Akhtar 风格 IMC-FxLMS 对照报告\n\n');
    fprintf(fid, '> 本报告由 `MarinoTomei_AkhtarBenchmark_compare.m` 自动生成。\n\n');

    fprintf(fid, '## 1. 结论先行\n\n');
    fprintf(fid, '本次 benchmark 参照 Akhtar 复现实验的场景集合与指标，比较固定参数 IMC-FxLMS baseline 和当前 demo5 中的 Marino-Tomei Case A 窄带内模实现。\n\n');
    fprintf(fid, '- 在 Case A 条件可靠的固定单频窄带场景中，Marino-Tomei 明显优于固定参数 IMC-FxLMS baseline：`Fixed-600Hz` 提高约 43.71 dB，`Fixed-1400Hz` 提高约 58.53 dB。\n');
    fprintf(fid, '- 在 Akhtar nominal `1000 Hz` 及其 mismatch 场景中，当前 Marino-Tomei Case A 不适用，因为次级路径在该频率附近 `Re(S)` 很小，`caseA ratio` 低于阈值 %.2f。\n', cfg.case_a_min_ratio);
    fprintf(fid, '- 对 chirp 和 bandlimited 场景，当前 `q=1` 的 Marino-Tomei 单振荡器不是通用 baseline；这些场景更接近 IMC-FxLMS 的适用范围，除非后续加入多内模、Case B 或更一般的相位补偿结构。\n');
    fprintf(fid, '- 因此，本次结果支持的结论是：Marino-Tomei 在“少量固定窄带 + Case A 条件可靠”时有优势，但不能替代 IMC-FxLMS 作为宽带/扫频/相位不利场景的通用 ANC baseline。\n\n');
    fprintf(fid, '**可靠性限制**：本报告沿用 Akhtar 复现实验中的固定 IMC-FxLMS 步长 `mu=%.3g`，没有对 IMC-FxLMS 的 `mu`、FIR 长度或归一化策略做网格搜索。因此，表中的优势只能理解为“相对于这组固定 baseline 参数”的优势，不能理解为 Marino-Tomei 相对于调参后的最优 IMC-FxLMS 一定占优。\n\n', cfg.imc_step_size);

    fprintf(fid, '## 2. 图表汇总\n\n');
    fprintf(fid, '### 2.1 指标汇总图\n\n');
    fprintf(fid, '![指标汇总](metric_summary.png)\n\n');
    fprintf(fid, '### 2.2 时域误差对比图\n\n');
    fprintf(fid, '![时域误差对比](time_domain_summary.png)\n\n');

    fprintf(fid, '## 3. 仿真配置\n\n');
    fprintf(fid, '| 项目 | 数值 |\n| --- | --- |\n');
    fprintf(fid, '| 采样率 | %.0f Hz |\n', cfg.fs);
    fprintf(fid, '| 仿真时长 | %.2f s |\n', cfg.duration);
    fprintf(fid, '| 次级路径 | `%s` |\n', mat2str(cfg.secondary_path));
    fprintf(fid, '| IMC-FxLMS 步长 | %.3g |\n', cfg.imc_step_size);
    fprintf(fid, '| IMC-FxLMS FIR taps | %d |\n', cfg.nQ + 1);
    fprintf(fid, '| Marino-Tomei k | %.3g |\n', cfg.marino_k);
    fprintf(fid, '| Marino-Tomei epsilon | %.3g |\n', cfg.marino_epsilon);
    fprintf(fid, '| Case A ratio 阈值 | %.2f |\n\n', cfg.case_a_min_ratio);

    fprintf(fid, '## 4. 判读规则\n\n');
    fprintf(fid, '- `caseA ratio = min |Re(S)|/|S|`，越接近 0，越说明当前频段接近 Marino-Tomei Case A 禁区。\n');
    fprintf(fid, '- `Marino advantage` 只在 IMC-FxLMS 和 Marino-Tomei 都稳定且适用时计算；否则显示 `n/a`。\n');
    fprintf(fid, '- `caseA_not_applicable` 不是 Marino-Tomei 数值发散，而是脚本主动跳过当前 Case A 不可靠的场景，避免给出误导性结论。\n\n');
    fprintf(fid, '另外，IMC-FxLMS 的表现对 `mu` 很敏感。当前报告没有声称 `mu=%.3g` 是每个场景的最优值；如果要做更严格的算法优劣判断，应增加 `mu` 扫描或每场景调参后的 oracle baseline。\n\n', cfg.imc_step_size);

    fprintf(fid, '## 5. 完整结果表\n\n');
    fprintf(fid, '| 场景 | 类型 | Marino 先验 Hz | Marino 末值 Hz | CaseA ratio | IMC 稳态 dB | Marino 稳态 dB | Marino 优势 dB | IMC 状态 | Marino 状态 |\n');
    fprintf(fid, '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- |\n');
    for i = 1:height(summary)
        fprintf(fid, '| %s | %s | %s | %s | %s | %s | %s | %s | %s | %s |\n', ...
            char(summary.name(i)), char(summary.kind(i)), ...
            fmt_num(summary.marinoPriorHz(i), '%.1f'), ...
            fmt_num(summary.marinoFinalHz(i), '%.1f'), ...
            fmt_num(summary.caseARatio(i), '%.3f'), ...
            fmt_num(summary.imcSteadyDb(i), '%.2f'), ...
            fmt_num(summary.marinoSteadyDb(i), '%.2f'), ...
            fmt_num(summary.marinoAdvantageDb(i), '%.2f'), ...
            char(summary.imcStatus(i)), char(summary.marinoStatus(i)));
    end

    fprintf(fid, '\n## 6. 关键场景解释\n\n');
    fprintf(fid, '### 6.1 Fixed-600Hz 与 Fixed-1400Hz\n\n');
    fprintf(fid, '这两个场景中 `caseA ratio` 分别约为 0.626 和 0.489，说明 `Re(S)` 在目标频率处足够可靠。Marino-Tomei 能用一个自适应振荡器直接锁定并抵消单频扰动，因此稳态抑制量远高于固定参数 IMC-FxLMS。\n\n');
    fprintf(fid, '### 6.2 Fixed-1000Hz 与 mismatch 场景\n\n');
    fprintf(fid, 'Akhtar benchmark 的 nominal 频率是 1000 Hz，但次级路径 `[0 0 0.5 0.3]` 在 1000 Hz 附近相位接近 -90 度，`Re(S)` 很小。当前 Case A 实现只使用 `sign(Re[S])`，所以这些场景被标记为 `caseA_not_applicable`。要公平比较，需要实现论文 Case B，或改成显式处理次级路径相位的 filtered-x/相位补偿结构。\n\n');
    fprintf(fid, '### 6.3 Chirp 与 bandlimited 场景\n\n');
    fprintf(fid, '当前 Marino-Tomei 分支只放入一个内模振荡器，适合固定或慢变单频窄带，不适合直接覆盖宽扫频或随机带限噪声。若要扩展到这些场景，需要增加多频内模、频率跟踪机制，或将 Marino-Tomei 作为窄带模块与 IMC-FxLMS 组合使用。\n');
end

function text = fmt_num(value, format_spec)
    if isnan(value)
        text = 'n/a';
    else
        text = sprintf(format_spec, value);
    end
end
