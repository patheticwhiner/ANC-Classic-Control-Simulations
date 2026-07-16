function generate_demo1234_reports(stage, figureRoot, testDir, demoFilter, partFilter)
%GENERATE_DEMO1234_REPORTS Create image-first reports for frozen Demo1-5.

if nargin < 4, demoFilter = ''; end
if nargin < 5, partFilter = ''; end

demoNames = {'demo1', 'demo2', 'demo3', 'demo4', 'demo5'};
for demoIndex = 1:numel(demoNames)
    demo = demoNames{demoIndex};
    if ~isempty(demoFilter) && ~strcmp(demo, demoFilter)
        continue;
    end
    demoFigureDir = fullfile(figureRoot, demo);
    if ~isfolder(demoFigureDir), mkdir(demoFigureDir); end

    results = results_for_demo(stage, demo);
    if isempty(partFilter) || strcmp(partFilter, 'overview')
        fprintf('Plot %s overview\n', demo);
        plot_demo_overview(results, ...
            fullfile(demoFigureDir, sprintf('%s_analysis.png', demo)));
    end
    if isempty(partFilter) || strcmp(partFilter, 't1')
        fprintf('Plot %s T1 fixed\n', demo);
        plot_t1_analysis(results.T1_fixed, fullfile(demoFigureDir, 't1_analysis.png'));
        if isfield(results, 'T1_adaptive')
            fprintf('Plot %s T1 adaptive\n', demo);
            plot_t1_analysis(results.T1_adaptive, ...
                fullfile(demoFigureDir, 't1_adaptive_analysis.png'));
        end
    end
    if isempty(partFilter) || strcmp(partFilter, 't2')
        fprintf('Plot %s T2\n', demo);
        plot_t2_analysis(results.T2, fullfile(demoFigureDir, 't2_analysis.png'));
        if strcmp(demo, 'demo2') && isfield(stage, 'demo2_t2_comparison')
            plot_t2_analysis(stage.demo2_t2_comparison.lms.evaluation, ...
                fullfile(demoFigureDir, 't2_lms_analysis.png'));
            plot_demo2_t2_adaptation_comparison(stage.demo2_t2_comparison, ...
                fullfile(demoFigureDir, 't2_adaptation_comparison.png'));
        end
    end
    if strcmp(partFilter, 't3')
        fprintf('Plot %s T3\n', demo);
        plot_t3_analysis(results.T3, fullfile(demoFigureDir, 't3_analysis.png'));
    end
    if isempty(partFilter) || strcmp(partFilter, 'design')
        fprintf('Plot %s design\n', demo);
        plot_design_analysis(stage, demo, results, ...
            fullfile(demoFigureDir, 'design_analysis.png'));
    end
    if isempty(partFilter) || strcmp(partFilter, 'report')
        fprintf('Write %s report\n', demo);

        reportPath = fullfile(testDir, sprintf('%s_CYLINDER1DM_REPORT.md', upper(demo)));
        write_demo_report(reportPath, stage, demo);
    end
end

if (isempty(partFilter) || strcmp(partFilter, 'report')) && isempty(demoFilter)
    write_stage_report(fullfile(testDir, 'DEMO1234_STAGE_REPORT.md'), stage);
end
end

function results = results_for_demo(stage, demo)
results = struct();
for index = 1:numel(stage.selections)
    if strcmp(stage.selections(index).demo, demo)
        selection = stage.selections(index);
        result = stage.evaluation_results{index};
        if strcmp(selection.test, 'T1')
            results.(sprintf('T1_%s', selection.variant)) = result;
            if strcmp(selection.variant, 'fixed')
                results.T1 = result;
            end
        else
            results.(selection.test) = result;
        end
    end
end
end

function plot_t1_analysis(result, outputPath)
t = result.extra.t(:);
yOpen = result.extra.y_open(:);
yClosed = result.extra.y_closed(:);
u = result.extra.u(:);
uDemand = result.extra.u_demand(:);
fs = result.fs;

fig = create_report_figure([80 80 1280 820]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(layout);
windowStart = max(1, numel(t) - round(0.08 * fs));
plot(ax, t(windowStart:end), yOpen(windowStart:end), 'Color', [0.55 0.55 0.55]);
hold(ax, 'on');
plot(ax, t(windowStart:end), yClosed(windowStart:end), 'Color', [0.00 0.36 0.62], ...
    'LineWidth', 1.2);
grid(ax, 'on');
xlabel(ax, 'Time (s)'); ylabel(ax, 'Residual amplitude');
title(ax, sprintf('Steady waveform: %.1f dB (gray open, blue controlled)', result.supp_db));

ax = nexttile(layout);
[fOpen, pOpen] = psd_curve(yOpen, fs);
[fClosed, pClosed] = psd_curve(yClosed, fs);
plot(ax, fOpen, pOpen, 'Color', [0.55 0.55 0.55], 'LineWidth', 1.1);
hold(ax, 'on');
plot(ax, fClosed, pClosed, 'Color', [0.00 0.36 0.62], 'LineWidth', 1.2);
xlim(ax, [100 500]); grid(ax, 'on');
xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'PSD (dB/Hz)');
title(ax, 'Steady-state spectrum (gray open, blue controlled)');

ax = nexttile(layout);
[tDemand, demandPeak] = block_peak_envelope(t, uDemand, 800);
[tApplied, appliedPeak] = block_peak_envelope(t, u, 800);
plot(ax, tApplied, appliedPeak, 'Color', [0.00 0.45 0.30], ...
    'LineWidth', 2.2, 'DisplayName', 'Applied (solid green)');
hold(ax, 'on');
plot(ax, tDemand, demandPeak, '--', 'Color', [0.80 0.33 0.10], ...
    'LineWidth', 1.2, 'DisplayName', 'Raw demand (dashed orange)');
yline(ax, 4, '--', 'Tuning bound', 'HandleVisibility', 'off');
yline(ax, 5, ':', 'Hard limit', 'HandleVisibility', 'off');
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Block peak |control|');
ylim(ax, [0 max(5.5, 1.05 * max(demandPeak))]);
legend(ax, 'Location', 'best');
title(ax, sprintf('Control: orange demand, green applied; max %.3f, clips %d', ...
    result.extra.u_demand_max, result.extra.saturation_count));

ax = nexttile(layout);
windowLength = max(16, round(0.05 * fs));
movingOpen = sqrt(movmean(yOpen.^2, windowLength));
movingClosed = sqrt(movmean(yClosed.^2, windowLength));
plot_reduced(ax, t, movingOpen, 'Color', [0.55 0.55 0.55], ...
    'LineWidth', 1.1, 'DisplayName', 'Open loop');
hold(ax, 'on');
plot_reduced(ax, t, movingClosed, 'Color', [0.00 0.36 0.62], ...
    'LineWidth', 1.2, 'DisplayName', 'Controlled');
set(ax, 'YScale', 'log');
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, '50 ms moving RMS');
title(ax, 'Full-run residual envelope (gray open, blue controlled)');

sgtitle(layout, sprintf('%s %s T1 at 357 Hz', upper(result.demo), result.variant));
save_figure(fig, outputPath);
end

function plot_t2_analysis(result, outputPath)
t = result.extra.t(:);
yOpen = result.extra.y_open(:);
yClosed = result.extra.y_closed(:);
fs = result.fs;

fig = create_report_figure([80 80 1280 820]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(layout);
windowLength = max(16, round(0.1 * fs));
movingOpen = sqrt(movmean(yOpen.^2, windowLength));
movingClosed = sqrt(movmean(yClosed.^2, windowLength));
plot_reduced(ax, t, movingOpen, 'Color', [0.55 0.55 0.55], ...
    'LineWidth', 1.1, 'DisplayName', 'Open loop');
hold(ax, 'on');
plot_reduced(ax, t, movingClosed, 'Color', [0.00 0.36 0.62], ...
    'LineWidth', 1.2, 'DisplayName', 'Controlled');
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, '100 ms moving RMS');
title(ax, sprintf(['Held-out chirp envelope: %.1f dB overall ' ...
    '(gray open, blue controlled)'], result.supp_db));

ax = nexttile(layout);
breakdown = result.supp_breakdown;
plot(ax, breakdown.bin_centers, breakdown.supp_db, '-o', ...
    'Color', [0.00 0.36 0.62], 'LineWidth', 1.4, 'MarkerSize', 4);
yline(ax, 20, '--', '20 dB'); yline(ax, 0, ':');
grid(ax, 'on'); xlim(ax, [300 420]);
xlabel(ax, 'Instantaneous frequency (Hz)'); ylabel(ax, 'Local suppression (dB)');
title(ax, sprintf('Segmented suppression: worst %.1f dB, median %.1f dB', ...
    min(breakdown.supp_db), median(breakdown.supp_db)));

[openSpec, frequency, specTime] = spectrogram_level(yOpen, fs);
[closedSpec, ~, ~] = spectrogram_level(yClosed, fs);
colorLimits = shared_color_limits(openSpec, closedSpec, 80);
ax = nexttile(layout);
plot_spectrogram(ax, openSpec, frequency, specTime, colorLimits, ...
    'Open-loop spectrogram');

ax = nexttile(layout);
plot_spectrogram(ax, closedSpec, frequency, specTime, colorLimits, ...
    'Controlled spectrogram');

sgtitle(layout, sprintf('%s %s T2 tracking', upper(result.demo), adaptation_label(result)));
save_figure(fig, outputPath);
end

function label = adaptation_label(result)
label = result.variant;
if isfield(result, 'extra') && isstruct(result.extra) ...
        && isfield(result.extra, 'adaptation_method')
    label = upper(result.extra.adaptation_method);
end
end

function plot_demo2_t2_adaptation_comparison(comparison, outputPath)
rls = comparison.rls.evaluation;
lms = comparison.lms.evaluation;
fs = rls.fs;

fig = create_report_figure([80 80 1280 820]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(layout);
plot(ax, rls.supp_breakdown.bin_centers, rls.supp_breakdown.supp_db, '-o', ...
    'Color', [0.00 0.36 0.62], 'LineWidth', 1.3);
hold(ax, 'on');
plot(ax, lms.supp_breakdown.bin_centers, lms.supp_breakdown.supp_db, '-s', ...
    'Color', [0.80 0.33 0.10], 'LineWidth', 1.3);
yline(ax, 20, '--', '20 dB'); yline(ax, 0, ':');
grid(ax, 'on'); xlim(ax, [300 420]);
xlabel(ax, 'Instantaneous frequency (Hz)'); ylabel(ax, 'Local suppression (dB)');
title(ax, sprintf('Local suppression (blue Q-RLS, orange LMS; %.1f/%.1f dB overall)', ...
    rls.supp_db, lms.supp_db));

ax = nexttile(layout);
windowLength = max(16, round(0.1 * fs));
rlsEnvelope = sqrt(movmean(rls.extra.y_closed(:).^2, windowLength));
lmsEnvelope = sqrt(movmean(lms.extra.y_closed(:).^2, windowLength));
plot(ax, rls.extra.t, rlsEnvelope, 'Color', [0.00 0.36 0.62]);
hold(ax, 'on');
plot(ax, lms.extra.t, lmsEnvelope, 'Color', [0.80 0.33 0.10]);
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Moving RMS');
title(ax, 'Residual envelope (blue Q-RLS, orange LMS)');

ax = nexttile(layout);
plot(ax, rls.extra.t, rls.extra.theta_norm_history, ...
    'Color', [0.00 0.36 0.62]);
hold(ax, 'on');
plot(ax, lms.extra.t, lms.extra.theta_norm_history, ...
    'Color', [0.80 0.33 0.10]);
yline(ax, 8, '--', 'Norm limit');
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Parameter norm');
title(ax, 'Parameter norm (blue Q-RLS, orange LMS)');

ax = nexttile(layout);
[tRls, peakRls] = block_peak_envelope(rls.extra.t, rls.extra.u_demand, 800);
[tLms, peakLms] = block_peak_envelope(lms.extra.t, lms.extra.u_demand, 800);
plot(ax, tRls, peakRls, 'Color', [0.00 0.36 0.62], ...
    'LineWidth', 1.1, 'DisplayName', 'Q-RLS');
hold(ax, 'on');
plot(ax, tLms, peakLms, 'Color', [0.80 0.33 0.10], ...
    'LineWidth', 1.1, 'DisplayName', 'LMS');
yline(ax, 4, '--', 'Tuning bound', 'HandleVisibility', 'off');
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Block peak |demand|');
ylim(ax, [0 max(4.4, 1.05 * max([peakRls; peakLms]))]);
title(ax, sprintf('Control demand: blue Q-RLS %.2f, orange LMS %.2f max', ...
    max(peakRls), max(peakLms)));

sgtitle(layout, 'Demo2 T2 adaptation comparison');
save_figure(fig, outputPath);
end

function plot_t3_analysis(result, outputPath)
t = result.extra.t(:);
yOpen = result.extra.y_open(:);
yClosed = result.extra.y_closed(:);
uDemand = result.extra.u_demand(:);
fs = result.fs;

fig = create_report_figure([80 80 1280 820]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(layout);
windowLength = max(16, round(0.1 * fs));
movingOpen = sqrt(movmean(yOpen.^2, windowLength));
movingClosed = sqrt(movmean(yClosed.^2, windowLength));
plot_reduced(ax, t, movingOpen, 'Color', [0.55 0.55 0.55], ...
    'LineWidth', 1.1, 'DisplayName', 'Open loop');
hold(ax, 'on');
plot_reduced(ax, t, movingClosed, 'Color', [0.00 0.36 0.62], ...
    'LineWidth', 1.2, 'DisplayName', 'Controlled');
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, '100 ms moving RMS');
title(ax, sprintf(['100-500 Hz residual envelope: %.1f dB overall ' ...
    '(gray open, blue controlled)'], result.supp_db));

ax = nexttile(layout);
[fOpen, pOpen] = psd_curve(yOpen, fs);
[fClosed, pClosed] = psd_curve(yClosed, fs);
plot(ax, fOpen, pOpen, 'Color', [0.55 0.55 0.55]); hold(ax, 'on');
plot(ax, fClosed, pClosed, 'Color', [0.00 0.36 0.62]);
xlim(ax, [100 500]); grid(ax, 'on');
xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'PSD (dB/Hz)');
title(ax, 'Broadband PSD (gray open, blue controlled)');

ax = nexttile(layout);
[tDemand, demandPeak] = block_peak_envelope(t, uDemand, 800);
plot(ax, tDemand, demandPeak, 'Color', [0.80 0.33 0.10], 'LineWidth', 1.1);
hold(ax, 'on');
yline(ax, 4, '--', 'Tuning bound');
yline(ax, 5, ':', 'Hard limit');
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Block peak |demand|');
ylim(ax, [0 max(5.5, 1.05 * max(demandPeak))]);
title(ax, sprintf('Demand max %.2f, clips %d', ...
    result.extra.u_demand_max, result.extra.saturation_count));

ax = nexttile(layout);
movingSuppression = 20 * log10(movingOpen ./ max(movingClosed, 1e-12));
plot_reduced(ax, t, movingSuppression, 'Color', [0.00 0.36 0.62], ...
    'LineWidth', 1.1);
yline(ax, 20, '--', '20 dB'); yline(ax, 0, ':'); grid(ax, 'on');
xlabel(ax, 'Time (s)'); ylabel(ax, 'Moving suppression (dB)');
title(ax, 'Time-varying broadband performance');

sgtitle(layout, sprintf('%s %s T3 broadband limit', upper(result.demo), result.variant));
save_figure(fig, outputPath);
end

function plot_design_analysis(stage, demo, results, outputPath)
switch demo
    case 'demo1'
        plot_demo1_design(stage, results, outputPath);
    case 'demo2'
        plot_demo2_design(stage, results, outputPath);
    case 'demo3'
        plot_demo3_design(stage, results, outputPath);
    case 'demo4'
        plot_demo4_design(stage, results, outputPath);
    case 'demo5'
        plot_demo5_design(stage, results, outputPath);
end
end

function plot_demo5_design(stage, results, outputPath)
r1 = results.T1;
r1Adaptive = results.T1_adaptive;
r2 = results.T2;
fig = create_report_figure([80 80 1280 820]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(layout);
plot(ax, r1.extra.case_a_frequencies_hz, r1.extra.case_a_ratio, ...
    'o-', 'LineWidth', 1.3, 'Color', [0.00 0.36 0.62]);
hold(ax, 'on');
if isfield(r1.extra, 'case_a_raw_ratio')
    plot(ax, r1.extra.case_a_frequencies_hz, r1.extra.case_a_raw_ratio, ...
        '--', 'LineWidth', 1.0, 'Color', [0.45 0.45 0.45], ...
        'DisplayName', 'raw Re(S)');
end
yline(ax, r1.extra.case_a_threshold, '--', 'Case A threshold');
grid(ax, 'on'); xlabel(ax, 'Frequency (Hz)');
ylabel(ax, '|Re(S_{eff})|/|S_{eff}|');
legend(ax, 'Location', 'best');
title(ax, sprintf('Discrete Case A: min effective ratio %.3f', ...
    r1.extra.case_a_min_ratio));

ax = nexttile(layout);
plot(ax, r1Adaptive.extra.t_log, r1Adaptive.extra.f_est_log(1, :), ...
    'Color', [0.80 0.33 0.10], 'LineWidth', 1.2);
hold(ax, 'on');
yline(ax, r1Adaptive.extra.freq_final_hz(1), ':', 'Final estimate');
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Frequency (Hz)');
title(ax, sprintf('T1 adaptive frequency estimate (final %.2f Hz)', ...
    r1Adaptive.extra.freq_final_hz(1)));

ax = nexttile(layout);
plot(ax, r2.extra.t_log(:), r2.extra.f_est_log(1, :).', ...
    'Color', [0.80 0.33 0.10], 'LineWidth', 1.2);
hold(ax, 'on');
truth = interp1(r2.extra.t(:), r2.extra.f_inst(:), r2.extra.t_log(:), ...
    'linear', 'extrap');
plot(ax, r2.extra.t_log(:), truth, 'Color', [0.55 0.55 0.55], ...
    'LineWidth', 1.1);
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Frequency (Hz)');
title(ax, sprintf('T2 estimate; Case A applicable = %d', ...
    r2.extra.case_a_applicable));

ax = nexttile(layout);
plot_tuning_points(ax, stage.tuning, 'demo5');
title(ax, 'Marino-Tomei candidate search (calibration set)');

sgtitle(layout, 'Demo5 Marino-Tomei Case A design and applicability');
save_figure(fig, outputPath);
end

function plot_demo1_design(stage, results, outputPath)
r1 = results.T1;
r2 = results.T2;
A = r1.extra.A;
B = r1.extra.B;
R = r1.extra.R0;
S = r1.extra.S0;
denominator = addPolynomials(conv(A, S), conv(B, R));
frequency = (100:500).';
fprintf('  demo1 design: frequency responses\n');
sensitivity = freqz(conv(A, S), denominator, frequency, r1.fs);
complementary = freqz(conv(B, R), denominator, frequency, r1.fs);

fprintf('  demo1 design: create figure\n');
fig = create_report_figure([80 80 1280 820]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(layout);
plot(ax, frequency, 20*log10(abs(sensitivity)), 'LineWidth', 1.3); hold(ax, 'on');
plot(ax, frequency, 20*log10(abs(complementary)), 'LineWidth', 1.3);
yline(ax, -20, '--', '20 dB suppression'); grid(ax, 'on');
xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'Magnitude (dB)');
title(ax, 'Sensitivity S (blue) and complementary T (orange)');

ax = nexttile(layout);
plot(ax, frequency, -20*log10(abs(sensitivity)), 'LineWidth', 1.4);
yline(ax, 20, '--', '20 dB'); xline(ax, 357, ':', '357 Hz');
grid(ax, 'on'); xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'Predicted rejection (dB)');
title(ax, 'RST rejection curve');

ax = nexttile(layout);
fprintf('  demo1 design: poles\n');
plot_unit_circle(ax); hold(ax, 'on');
scatter(ax, real(roots(denominator)), imag(roots(denominator)), 45, 'x', ...
    'LineWidth', 1.4);
scatter(ax, real(roots(S)), imag(roots(S)), 38, 'o', 'LineWidth', 1.2);
axis(ax, 'equal'); grid(ax, 'on'); xlabel(ax, 'Real'); ylabel(ax, 'Imaginary');
title(ax, 'Poles: closed loop x, controller o, dashed unit circle');

ax = nexttile(layout);
fprintf('  demo1 design: adaptation\n');
plot(ax, r2.extra.t, r2.extra.q_norm_history, 'LineWidth', 1.2); hold(ax, 'on');
yyaxis(ax, 'right');
plot(ax, r2.extra.t, abs(r2.extra.adaptation_error), ...
    'Color', [0.80 0.33 0.10]);
set(ax, 'YScale', 'log');
grid(ax, 'on'); xlabel(ax, 'Time (s)');
title(ax, 'Q-RLS adaptation on held-out chirp');
yyaxis(ax, 'left'); ylabel(ax, '||Q||');
yyaxis(ax, 'right'); ylabel(ax, '|posterior error|');

sgtitle(layout, 'Demo1 low-order RST and Youla-Q design');
fprintf('  demo1 design: export\n');
save_figure(fig, outputPath);
fprintf('  demo1 design: done\n');
end

function plot_demo2_design(stage, results, outputPath)
r1 = results.T1;
r2 = results.T2;
fig = create_report_figure([80 80 1280 820]);
layout = tiledlayout(fig, 2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

ax = nexttile(layout);
plot_unit_circle(ax); hold(ax, 'on');
poles = eig(r1.extra.A_cl);
scatter(ax, real(poles), imag(poles), 45, 'x', 'LineWidth', 1.4);
axis(ax, 'equal'); grid(ax, 'on'); xlabel(ax, 'Real'); ylabel(ax, 'Imaginary');
title(ax, sprintf('LQG poles, max radius %.4f', max(abs(poles))));

ax = nexttile(layout);
stem(ax, r1.extra.K, 'filled'); hold(ax, 'on');
stem(ax, r1.extra.L, 'Marker', 'none');
grid(ax, 'on'); xlabel(ax, 'State index'); ylabel(ax, 'Gain');
title(ax, 'Gains: controller K (blue), estimator L (orange)');

ax = nexttile(layout);
antinoise = r1.extra.y_closed(:) - r1.extra.y_open(:);
plot(ax, r1.extra.t, r1.extra.y_open, 'Color', [0.65 0.65 0.65]); hold(ax, 'on');
plot(ax, r1.extra.t, antinoise, 'Color', [0.00 0.45 0.30]);
plot(ax, r1.extra.t, r1.extra.y_closed, 'Color', [0.00 0.36 0.62]);
xlim(ax, [max(0, r1.Tsim-0.08), r1.Tsim]); grid(ax, 'on');
xlabel(ax, 'Time (s)'); ylabel(ax, 'Amplitude');
title(ax, 'Gray disturbance, green antinoise, blue residual');

ax = nexttile(layout);
plot(ax, r2.extra.t, r2.extra.theta_norm_history, 'LineWidth', 1.3);
yyaxis(ax, 'right');
errorEnvelope = sqrt(movmean(r2.extra.adaptation_error(:).^2, ...
    max(8, round(.05*r2.fs))));
plot(ax, r2.extra.t, errorEnvelope, 'Color', [0.80 0.33 0.10]);
set(ax, 'YScale', 'log'); grid(ax, 'on'); xlabel(ax, 'Time (s)');
    title(ax, 'T2 normalized-LMS convergence');
yyaxis(ax, 'left'); ylabel(ax, '||theta||');
yyaxis(ax, 'right'); ylabel(ax, 'Residual RMS');

ax = nexttile(layout);
[tQ, peakQ] = block_peak_envelope(r2.extra.t, r2.extra.q_output, 700);
[tDemand, peakDemand] = block_peak_envelope( ...
    r2.extra.t, r2.extra.u_demand, 700);
plot(ax, tQ, peakQ, 'LineWidth', 1.1, 'DisplayName', '|u_Q| peak');
hold(ax, 'on');
plot(ax, tDemand, peakDemand, 'LineWidth', 1.1, ...
    'DisplayName', '|raw demand| peak');
yline(ax, 4, 'k--', 'Tuning bound', 'HandleVisibility', 'off');
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Block peak magnitude');
title(ax, 'Adaptive increment (blue) and raw demand (orange)');

ax = nexttile(layout);
plot_tuning_points(ax, stage.tuning, 'demo2');

sgtitle(layout, 'Demo2 fixed LQG and normalized-LMS adaptation');
save_figure(fig, outputPath);
end

function plot_demo3_design(stage, results, outputPath)
r1 = results.T1;
r2 = results.T2;
frequency = (100:500).';
flattening = freqz(r1.extra.F_num, r1.extra.F_den, frequency, r1.fs);
effective = freqz(r1.extra.H_num, r1.extra.H_den, frequency, r1.fs);

fig = create_report_figure([80 80 1280 820]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
ax = nexttile(layout);
plot(ax, frequency, 20*log10(abs(flattening)), 'LineWidth', 1.3); hold(ax, 'on');
plot(ax, frequency, 20*log10(abs(effective)), 'LineWidth', 1.3);
grid(ax, 'on'); xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'Magnitude (dB)');
title(ax, 'F(z) (blue) and effective S(z)F(z) (orange)');

ax = nexttile(layout);
stem(ax, 0:numel(r1.extra.theta_fixed)-1, r1.extra.theta_fixed, 'filled');
grid(ax, 'on'); xlabel(ax, 'Tap'); ylabel(ax, 'Coefficient');
title(ax, sprintf('Fixed FIR controller, N=%d', r1.extra.N_fir));

ax = nexttile(layout);
plot(ax, r2.extra.coefficient_time, r2.extra.coefficient_norm, 'LineWidth', 1.3);
yyaxis(ax, 'right');
errorEnvelope = sqrt(movmean(r2.extra.internal_error(:).^2, max(8, round(.05*r2.fs))));
plot(ax, r2.extra.t, errorEnvelope, 'Color', [0.80 0.33 0.10]);
set(ax, 'YScale', 'log'); grid(ax, 'on'); xlabel(ax, 'Time (s)');
    title(ax, 'IMC-FxNLMS convergence');
yyaxis(ax, 'left'); ylabel(ax, 'Coefficient norm');
yyaxis(ax, 'right'); ylabel(ax, 'Error RMS');

ax = nexttile(layout);
plot_tuning_points(ax, stage.tuning, 'demo3');

sgtitle(layout, 'Demo3 F(z)+FIR and IMC-FxNLMS');
save_figure(fig, outputPath);
end

function plot_demo4_design(stage, results, outputPath)
r1 = results.T1;
r2 = results.T2;
design1 = r1.extra.design;

fig = create_report_figure([80 80 1280 820]);
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% 面板1: T1 输出灵敏度 |Syp| (含 Ms 参考线与陷波频率)
ax = nexttile(layout);
frequencyHz = r1.extra.syp_omega * r1.fs / (2*pi);
sypDb = 20*log10(r1.extra.syp_mag + eps);
plot(ax, frequencyHz, sypDb, 'LineWidth', 1.3, 'Color', [0.00 0.36 0.62]);
hold(ax, 'on');
yline(ax, 6, '--', 'Ms = 6 dB');
if ~isempty(design1.f_notch)
    xline(ax, design1.f_notch, ':', sprintf('%g Hz', design1.f_notch));
end
xlim(ax, [50 700]); grid(ax, 'on');
xlabel(ax, 'Frequency (Hz)'); ylabel(ax, '|S_{yp}| (dB)');
title(ax, sprintf('T1 output sensitivity, Ms %.2f dB, design rejection %.1f dB', ...
    design1.Ms_db, design1.suppression_design_db));

% 面板2: 闭环与控制器极点 (T1 冻结解)
ax = nexttile(layout);
plot_unit_circle(ax); hold(ax, 'on');
scatter(ax, real(roots(r1.extra.P0)), imag(roots(r1.extra.P0)), 45, 'x', ...
    'LineWidth', 1.4);
scatter(ax, real(roots(r1.extra.S0)), imag(roots(r1.extra.S0)), 38, 'o', ...
    'LineWidth', 1.2);
axis(ax, 'equal'); grid(ax, 'on'); xlabel(ax, 'Real'); ylabel(ax, 'Imaginary');
title(ax, sprintf('Poles: closed loop x (%.4f), controller o (%.4f)', ...
    design1.closed_loop_radius, design1.controller_radius));

% 面板3: Pareto 前沿 (T2 冻结中央控制器，前两个目标)
ax = nexttile(layout);
archiveF = r2.extra.archive_f;
archiveDisplay = signed_log1p(archiveF(:, 1:2));
scatter(ax, archiveDisplay(:, 1), archiveDisplay(:, 2), ...
    30, [0.55 0.55 0.55], 'filled', ...
    'MarkerFaceAlpha', 0.6);
hold(ax, 'on');
selectedIndex = r2.extra.design.selected_index;
scatter(ax, archiveDisplay(selectedIndex, 1), archiveDisplay(selectedIndex, 2), 90, ...
    [0.80 0.33 0.10], 'filled');
grid(ax, 'on'); xlabel(ax, 'signed log_{10}(1+|f_1|)');
ylabel(ax, 'signed log_{10}(1+|f_2|)');
title(ax, sprintf('T2 central Pareto archive: %d solutions, %d stable, pick #%d', ...
    r2.extra.design.n_pareto_solutions, r2.extra.design.stable_solutions, ...
    selectedIndex));

% 面板4: 主实验控制需求与约束
ax = nexttile(layout);
tests = {'T1', 'T2'};
colors = {[0.00 0.36 0.62], [0.80 0.33 0.10]};
hold(ax, 'on');
for index = 1:numel(tests)
    result = results.(tests{index});
    [peakTime, peakDemand] = block_peak_envelope( ...
        result.extra.t, result.extra.u_demand, 700);
    plot(ax, peakTime, peakDemand, 'Color', colors{index}, ...
        'LineWidth', 1.1, 'DisplayName', tests{index});
end
yline(ax, 4, '--', 'tuning bound');
yline(ax, 5, ':', 'hard limit');
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, 'Block peak |demand|');
title(ax, 'Control demand: T1 blue, T2 orange');

sgtitle(layout, 'Demo4 eMOPSO sensitivity-shaping RST design');
save_figure(fig, outputPath);
end

function plot_tuning_points(ax, tuningTable, demo)
mask = strcmp(tuningTable.demo, demo) & cellfun(@isempty, tuningTable.error_message);
tests = {'T1', 'T2'};
colors = lines(2);
hold(ax, 'on');
for testIndex = 1:numel(tests)
    rows = mask & strcmp(tuningTable.test, tests{testIndex});
    scatter(ax, tuningTable.control_demand_max(rows), tuningTable.suppression_db(rows), ...
        38, colors(testIndex, :), 'filled', 'DisplayName', tests{testIndex});
    selectedRows = rows & tuningTable.selected;
    scatter(ax, tuningTable.control_demand_max(selectedRows), ...
        tuningTable.suppression_db(selectedRows), 90, colors(testIndex, :), ...
        'o', 'LineWidth', 1.8, 'HandleVisibility', 'off');
end
xline(ax, 4, '--', 'Tuning bound'); yline(ax, 20, '--', '20 dB');
set(ax, 'XScale', 'log');
grid(ax, 'on'); xlabel(ax, 'Peak control demand'); ylabel(ax, 'Calibration suppression (dB)');
title(ax, 'Candidate search: T1/T2 (log demand; rings selected)');
end

function [level, frequency, time] = spectrogram_level(signal, fs)
window = hann(256, 'periodic');
[spectrum, frequency, time] = spectrogram(signal, window, 192, 512, fs);
level = 20 * log10(abs(spectrum) + 1e-8);
end

function limits = shared_color_limits(firstLevel, secondLevel, dynamicRange)
peak = max([firstLevel(:); secondLevel(:)]);
limits = [peak - dynamicRange, peak];
end

function plot_spectrogram(ax, level, frequency, time, colorLimits, plotTitle)
imagesc(ax, time, frequency, level); axis(ax, 'xy');
ylim(ax, [250 470]); colorbar(ax);
clim(ax, colorLimits);
xlabel(ax, 'Time (s)'); ylabel(ax, 'Frequency (Hz)'); title(ax, plotTitle);
end

function [frequency, psdDb] = psd_curve(signal, fs)
[psdValue, frequency] = safeWelch(signal, fs);
psdDb = 10 * log10(psdValue + 1e-14);
end

function [psdValue, frequency] = safeWelch(signal, fs)
% Local robust Welch wrapper for report-only figures.
signal = signal(:);
if numel(signal) < 8
    error('generate_demo1234_reports:ShortPSDSignal', ...
        'PSD plotting requires at least 8 samples.');
end
targetLength = max(8, floor(numel(signal) / 4));
segmentLength = min([1024, numel(signal), 2^floor(log2(targetLength))]);
overlap = floor(segmentLength / 2);
nfft = max(2048, 2^nextpow2(segmentLength));
[psdValue, frequency] = pwelch(signal, hann(segmentLength, 'periodic'), ...
    overlap, nfft, fs);
end

function plot_reduced(ax, time, value, varargin)
sampleStep = max(1, floor(numel(time) / 2500));
plot(ax, time(1:sampleStep:end), value(1:sampleStep:end), varargin{:});
end

function [blockTime, blockPeak] = block_peak_envelope(time, signal, blockCount)
time = time(:);
signal = signal(:);
edges = round(linspace(1, numel(signal) + 1, ...
    min(blockCount, numel(signal)) + 1));
blockTime = zeros(numel(edges) - 1, 1);
blockPeak = zeros(numel(edges) - 1, 1);
for blockIndex = 1:numel(edges) - 1
    samples = edges(blockIndex):edges(blockIndex + 1) - 1;
    blockTime(blockIndex) = mean(time(samples));
    blockPeak(blockIndex) = max(abs(signal(samples)));
end
end

function transformed = signed_log1p(value)
transformed = sign(value) .* log10(1 + abs(value));
end

function plot_unit_circle(ax)
angle = linspace(0, 2*pi, 400);
plot(ax, cos(angle), sin(angle), 'k--', 'LineWidth', 1, ...
    'DisplayName', 'Unit circle');
end

function fig = create_report_figure(position)
options = test_display_options();
fig = figure('Visible', options.figure_visibility, 'Color', 'w', ...
    'Position', position);
end

function save_figure(fig, outputPath)
options = test_display_options();
set(fig, 'PaperPositionMode', 'auto');
print(fig, outputPath, '-dpng', '-r180');
if options.close_after_save
    close(fig);
end
drawnow;
end

function write_demo_report(reportPath, stage, demo)
fid = fopen(reportPath, 'w', 'n', 'UTF-8');
if fid < 0, error('Cannot write report: %s', reportPath); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

demoLabel = upper(demo);
profile = controller_profile(demo);
results = results_for_demo(stage, demo);
fprintf(fid, '# %s Cylinder 1 dm 独立仿真报告\n\n', demoLabel);
if strcmp(demo, 'demo5') && isfield(results.T1_fixed.extra, 'output_timing')
    fprintf(fid, '> 离散输出时序：`%s`。本报告只解释当前冻结结果，不能与其他时序生成的旧 CSV/MAT 混用。\n\n', ...
        results.T1_fixed.extra.output_timing);
end
fprintf(fid, '## 1. 本报告要回答什么\n\n');
fprintf(fid, ['本报告面向需要判断“Cylinder 1 dm 测量 FIR 路径能否支持实际控制器设计”的控制工程读者。' ...
    '目标不是展示脚本能够运行，而是回答以下问题：\n\n']);
fprintf(fid, '1. %s 能否在 357 Hz 主导频点形成稳定且不过载的控制器？\n', profile.method);
fprintf(fid, '2. 当频率在 300–420 Hz 内缓慢变化时，该方法是否仍能维持抑制？\n');
fprintf(fid, ['宽带 T3 不在本报告中，独立见 ' ...
    '`BROADBAND_CYLINDER1DM_REPORT.md`。\n\n']);
fprintf(fid, ['性能结论只使用未参与调参的评价集。成功要求：内部稳定、未限幅控制需求不超过 4、' ...
    '无硬限幅，并达到至少 20 dB 抑制。\n\n']);

fprintf(fid, '## 2. 被比较的控制方法\n\n');
fprintf(fid, '| 项目 | 本 Demo 设置 |\n');
fprintf(fid, '|---|---|\n');
fprintf(fid, '| 控制方法 | %s |\n', profile.method);
fprintf(fid, '| 信息结构 | %s |\n', profile.information);
fprintf(fid, '| 主要自由度 | %s |\n', profile.freedom);
fprintf(fid, '| 预期优势 | %s |\n', profile.strength);
fprintf(fid, '| 预期限制 | %s |\n\n', profile.limitation);

if strcmp(demo, 'demo5')
    fprintf(fid, '## 2.1 实现核对与修复记录\n\n');
    fprintf(fid, ['本次核对将原始 `demo5_MarinoTomei` 与 `tests/internal/controller_demo5` ' ...
        '放在同一条合成次级路径和双音设置下进行对照。内模状态更新、反馈符号和频率自适应律能够复现原始 Demo 的抑制效果，' ...
        '因此没有发现 Marino-Tomei 核心公式被实现成反号或失效版本。\n\n']);
    fprintf(fid, ['此前 `tests/` 的正式候选使用 `output_timing=''previous''`，而原始 Demo 使用更新后的内模状态生成本拍输出；' ...
        '这会引入额外一拍控制滞后。统一测量 FIR 后，正式候选现使用 `output_timing=''updated''`，' ...
        '`previous` 只保留为 ARMAX 离散相位边界诊断。\n\n']);
    fprintf(fid, ['此前可读脚本在没有兼容冻结 MAT 时会回放一个已知无效的 adaptive fallback，导致控制量接近零并产生“没有降噪”的表象；' ...
        'fallback 现已与当前校准候选一致。最终结论仍只使用冻结评价集，不把 fallback 或校准结果冒充为成功设计。\n\n']);
end

fprintf(fid, '## 3. 为什么设置 T1、T2\n\n');
write_scenario_matrix(fid);
fprintf(fid, ['T1 检查“能否设计”，T2 检查“偏离设计点后是否仍有效”。' ...
    '宽带失败边界由独立入口运行，不参与这里的 fixed/adaptive 主比较。\n\n']);

fprintf(fid, '## 4. 独立调参过程\n\n');
fprintf(fid, ['每个场景只在校准集上搜索参数，再把选中的参数原样运行于评价集。' ...
    '完整候选记录位于 `tests/output/cylinder1dm_2k/demo1234_stage/demo1234_tuning.csv`。\n\n']);
for testName = {'T1', 'T2'}
    write_tuning_comparison(fid, stage, demo, testName{1});
end

fprintf(fid, '## 5. 控制器设计与稳定性分析\n\n');
fprintf(fid, '![Design analysis](figures/cylinder1dm_2k/%s/design_analysis.png)\n\n', demo);
write_design_figure_interpretation(fid, stage, demo, results);
fprintf(fid, '%s\n\n', profile.designFigureMeaning);

fprintf(fid, '## 6. 冻结评价跨场景总览\n\n');
write_overview_result(fid, stage, demo);
if strcmp(demo, 'demo5')
    write_demo5_baseline_comparison(fid, stage);
end

fprintf(fid, '## 7. 各场景结果与解释\n\n');
for testName = {'T1', 'T2'}
    write_scenario_result(fid, stage, demo, testName{1});
    if strcmp(testName{1}, 'T1')
        write_t1_adaptive_result(fid, stage, demo);
    end
end

fprintf(fid, '## 8. 冻结评价摘要\n\n');
fprintf(fid, '| Test | Variant | Candidate | Suppression dB | max demand | Clips | Success |\n');
fprintf(fid, '|---|---|---|---:|---:|---:|---|\n');
rows = strcmp(stage.summary.demo, demo);
indices = find(rows).';
for index = indices
    fprintf(fid, '| %s | %s | %s | %.2f | %.3f | %d | %s |\n', ...
        stage.summary.test{index}, stage.summary.variant{index}, ...
        stage.summary.candidate{index}, ...
        stage.summary.evaluation_suppression_db(index), ...
        stage.summary.control_demand_max(index), ...
        stage.summary.saturation_count(index), ...
        yes_no(stage.summary.evaluation_success(index)));
end

fprintf(fid, '\n## 9. 本 Demo 最终说明了什么\n\n');
successRows = rows & stage.summary.evaluation_success;
successTests = strcat(stage.summary.test(successRows), " ", ...
    stage.summary.variant(successRows));
if isempty(successTests)
    fprintf(fid, '- 没有场景同时满足 20 dB、稳定性和控制约束。\n');
else
    fprintf(fid, '- 已证明有效的场景：%s。\n', strjoin(successTests, '、'));
end
failedRows = rows & ~stage.summary.evaluation_success;
failedTests = strcat(stage.summary.test(failedRows), " ", ...
    stage.summary.variant(failedRows));
if ~isempty(failedTests)
    fprintf(fid, '- 尚未证明有效的场景：%s。失败结果用于界定方法适用范围，而不是被隐藏。\n', ...
        strjoin(failedTests, '、'));
end
fprintf(fid, '- %s\n', profile.finalConclusion);
fprintf(fid, '- 当前默认路径是 2 kHz LMS FIR(16) 测量模型；其控制结论仍只代表该辨识路径，不能直接外推到实物系统。\n');
end

function write_demo5_baseline_comparison(fid, stage)
fprintf(fid, '### IMC-FxLMS 纯反馈基线对照\n\n');
fprintf(fid, ['Marino-Tomei 与 IMC-FxLMS 都只使用误差麦克风；后者额外使用次级路径模型重构内部参考，' ...
    '因此这是信息结构相近但自由度不同的对照，不是完全相同算法的排行榜。\n\n']);
fprintf(fid, '| Test | Marino-Tomei dB | Marino demand | Marino clips | IMC-FxLMS dB | IMC demand | IMC clips |\n');
fprintf(fid, '|---|---:|---:|---:|---:|---:|---:|\n');
for testName = {'T1', 'T2'}
    test = testName{1};
    marino = evaluation_result_for(stage, 'demo5', test);
    imc = evaluation_result_for(stage, 'demo5_imc', test);
    fprintf(fid, '| %s | %.2f | %.3f | %d | %.2f | %.3f | %d |\n', ...
        test, marino.supp_db, marino.extra.u_demand_max, ...
        marino.extra.saturation_count, imc.supp_db, ...
        imc.extra.u_demand_max, imc.extra.saturation_count);
end
fprintf(fid, '\n**对照解释**：IMC-FxLMS 的 FIR 自由度不依赖显式恒频内模；Demo5 的 Case A 只在离散有效次级响应符号可靠且扰动近似恒频时适用。\n\n');
end

function write_overview_result(fid, stage, demo)
fixed = evaluation_result_for(stage, demo, 'T1', 'fixed');
adaptive = evaluation_result_for(stage, demo, 'T1', 'adaptive');
chirp = evaluation_result_for(stage, demo, 'T2');
[worstLocal, medianLocal] = local_suppression_summary(chirp);
maxDemand = max([fixed.extra.u_demand_max, ...
    adaptive.extra.u_demand_max, chirp.extra.u_demand_max]);
totalClips = fixed.extra.saturation_count ...
    + adaptive.extra.saturation_count + chirp.extra.saturation_count;

fprintf(fid, '![%s frozen overview](figures/cylinder1dm_2k/%s/%s_analysis.png)\n\n', ...
    upper(demo), demo, demo);
fprintf(fid, ['该图把三个冻结评价放在同一证据面板中：左上直接比较 T1 fixed 与 adaptive 残差，' ...
    '右上显示在线参数状态，左下同时给出 T2 总体和局部分箱抑制，' ...
    '右下用分块峰值展示未限幅控制需求及 4/5 两条约束线。\n\n']);
fprintf(fid, ['**图中结果**：T1 fixed 为 %.2f dB，adaptive 为 %.2f dB；' ...
    'T2 总体为 %.2f dB，但最差/中位局部分箱为 %.2f/%.2f dB；' ...
    '三项评价的最大未限幅需求为 %.3f，合计限幅 %d 次。\n\n'], ...
    fixed.supp_db, adaptive.supp_db, chirp.supp_db, ...
    worstLocal, medianLocal, maxDemand, totalClips);
if adaptive.supp_db > fixed.supp_db
    adaptationText = 'T1 adaptive 相对 fixed 提高了总体抑制';
else
    adaptationText = 'T1 adaptive 相对 fixed 降低了总体抑制';
end
if worstLocal < 20
    fprintf(fid, ['**总览解读**：%s；T2 最差分箱低于 20 dB，' ...
        '因此不能把总体 RMS 成功解释为整个 300–420 Hz 频带均匀成功。'], ...
        adaptationText);
else
    fprintf(fid, ['**总览解读**：T2 各有效分箱均达到 20 dB，' ...
        '总体指标没有掩盖局部失败。']);
end
if maxDemand > 4 || totalClips > 0
    fprintf(fid, [' 同时，控制需求或限幅违反统一约束，' ...
        '图中的抑制改善不构成工程成功。\n\n']);
else
    fprintf(fid, ' 控制需求保持在调参约束内且未发生限幅。\n\n');
end
end

function write_stage_report(reportPath, stage)
fid = fopen(reportPath, 'w', 'n', 'UTF-8');
if fid < 0, error('Cannot write report: %s', reportPath); end
cleanup = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '# Demo1-5 Cylinder 1 dm 阶段报告\n\n');
demo5Index = find(strcmp({stage.selections.demo}, 'demo5'), 1);
if ~isempty(demo5Index) && isfield(stage.evaluation_results{demo5Index}.extra, 'output_timing')
    fprintf(fid, '> Demo5 离散输出时序：`%s`；所有数值结论只适用于本次冻结结果。\n\n', ...
        stage.evaluation_results{demo5Index}.extra.output_timing);
end
fprintf(fid, '## 1. 阶段目标\n\n');
fprintf(fid, ['本阶段要验证的核心命题是：**实测声学路径的 2 kHz LMS FIR 模型，是否足以支持经典控制器设计，' ...
    '并在未参与调参的信号上得到可解释的控制结果。**\n\n']);
fprintf(fid, ['本阶段只回答两个主问题：能否在主导频点设计控制器，以及能否适应频率漂移。' ...
    '100–500 Hz 宽带设计由独立阶段运行并报告，不参与本阶段的 fixed/adaptive 对比。' ...
    'Demo4 使用 epsilon-MOPSO 多目标灵敏度整形离线设计中央 RST，' ...
    '并在其上复用 Demo1 的 Youla-Q RLS 自适应律；Demo5 单列纯反馈窄带频率估计器及 IMC-FxLMS 基线，' ...
    '同时报告 Case A 的适用性边界。\n\n']);

fprintf(fid, '## 2. 统一实验条件\n\n');
fprintf(fid, '| 条件 | 设置 | 作用 |\n');
fprintf(fid, '|---|---|---|\n');
fprintf(fid, '| 被控模型 | 主路径 LMS FIR(16)，次级路径 LMS FIR(16)，2 kHz | 验证统一测量 FIR 路径是否可直接用于设计 |\n');
fprintf(fid, '| 校准集 | 8 s，seed 42 | 仅用于选择参数，不作为最终结论 |\n');
fprintf(fid, '| 评价集 | 10 s，seed 142，T2 反向扫频 | 检查冻结参数是否能泛化到未见波形 |\n');
fprintf(fid, '| 执行器 | 硬限幅 5，调参峰值约束 4 | 防止用不可实现的大控制量换取表面抑制 |\n');
fprintf(fid, '| 成功标准 | 稳定、未限幅需求不超过 4、无饱和、至少 20 dB 抑制 | 区分“仿真运行”与“控制器设计成功” |\n\n');

fprintf(fid, '## 2.1 Demo5 实现核对结论\n\n');
fprintf(fid, ['原始 `demo5_MarinoTomei` 与 `tests/internal/controller_demo5` 在同一条合成次级路径和双音设置下能够得到一致的抑制趋势，' ...
    '说明 Marino-Tomei 内模状态方程、反馈符号和频率自适应律本身没有失效。此前 `tests/` 与参考 Demo 的主要差异是输出时序：' ...
    '正式 FIR 候选曾使用 `previous`，现统一为参考 Demo 的 `updated`；`previous` 仅保留为 ARMAX 离散相位诊断。' ...
    '此外，可读脚本的无 stage fallback 曾回放已知无效的 adaptive 候选，现已与当前校准候选对齐。\n\n']);

fprintf(fid, '## 3. 两种主测试分别验证什么\n\n');
write_scenario_matrix(fid);

fprintf(fid, '## 4. 五个 Demo 的角色与可比性\n\n');
fprintf(fid, '| Demo | 方法 | 信息结构 | 主要验证问题 | 重要限制 |\n');
fprintf(fid, '|---|---|---|---|---|\n');
fprintf(fid, '| Demo1 | 低阶 RST + Youla-Q RLS | 反馈与内部扰动重构 | 低阶多项式设计和慢变频率自适应是否可行 | 中央控制器窄带，Q 自由度有限 |\n');
fprintf(fid, '| Demo2 | 增广 LQG + normalized-LMS（Q-RLS 对照） | 仅输出反馈与扰动状态模型 | 最优控制和反馈自适应能否处理定频与扫频 | 边缘频段局部抑制仍低于 20 dB |\n');
fprintf(fid, '| Demo3 | F(z)+固定 FIR / IMC-FxNLMS | 源参考 + 次级路径模型 | 固定模型逆与源参考自适应是否适合定频和扫频 | 信息结构不同于 Demo1/2/4，必须单独解释 |\n');
fprintf(fid, '| Demo4 | eMOPSO 整形 RST + Youla-Q RLS | 误差输出反馈，自适应律同 Demo1 | 优化设计的中央控制器能否替代手工扫描设计并承载同款自适应 | 每个中央候选需一次完整优化，计算量大 |\n\n');
fprintf(fid, '| Demo5 | Marino-Tomei Case A + 纯反馈 IMC-FxLMS | 误差麦克风；前者用频率内模，后者用次级路径模型重构内部参考 | 少量窄带频率估计结构在当前离散模型上的适用性边界 | Case A 固定符号不能跨越离散有效响应零点；T1 还受控制余量限制 |\n\n');
fprintf(fid, ['因此，横向比较关注的是“不同经典结构在各场景下表现出什么能力和限制”，' ...
    '而不是在完全相同信息条件下给算法排绝对名次。\n\n']);

fprintf(fid, '## 5. 独立报告\n\n');
fprintf(fid, '- [Demo1](DEMO1_CYLINDER1DM_REPORT.md)\n');
fprintf(fid, '- [Demo2](DEMO2_CYLINDER1DM_REPORT.md)\n');
fprintf(fid, '- [Demo3](DEMO3_CYLINDER1DM_REPORT.md)\n');
fprintf(fid, '- [Demo4](DEMO4_CYLINDER1DM_REPORT.md)\n\n');
fprintf(fid, '- [Demo5](DEMO5_CYLINDER1DM_REPORT.md)\n\n');
fprintf(fid, '- [独立宽带边界报告](BROADBAND_CYLINDER1DM_REPORT.md)\n\n');

fprintf(fid, '## 6. 冻结留出结果\n\n');
fprintf(fid, '| Demo | Test | Variant | Suppression dB | max demand | Clips | Success |\n');
fprintf(fid, '|---|---|---|---:|---:|---:|---|\n');
for index = 1:height(stage.summary)
    fprintf(fid, '| %s | %s | %s | %.2f | %.3f | %d | %s |\n', ...
        stage.summary.demo{index}, stage.summary.test{index}, ...
        stage.summary.variant{index}, ...
        stage.summary.evaluation_suppression_db(index), ...
        stage.summary.control_demand_max(index), ...
        stage.summary.saturation_count(index), ...
        yes_no(stage.summary.evaluation_success(index)));
end

fprintf(fid, '\n## 7. 按场景比较得到什么结论\n\n');
demo1T1 = summary_result(stage, 'demo1', 'T1');
demo2T1 = summary_result(stage, 'demo2', 'T1');
demo3T1 = summary_result(stage, 'demo3', 'T1');
demo1T2 = summary_result(stage, 'demo1', 'T2');
demo2T2 = summary_result(stage, 'demo2', 'T2');
demo3T2 = summary_result(stage, 'demo3', 'T2');
demo4T1 = summary_result(stage, 'demo4', 'T1');
demo4T2 = summary_result(stage, 'demo4', 'T2');
demo5T1 = summary_result(stage, 'demo5', 'T1');
demo5T2 = summary_result(stage, 'demo5', 'T2');
fprintf(fid, '### T1：低阶模型能否设计出有效控制器？\n\n');
fprintf(fid, ['四种方法分别达到 %.2f、%.2f、%.2f 和 %.2f dB，控制需求分别为 %.3f、%.3f、%.3f 和 %.3f，' ...
    '限幅次数分别为 %d、%d、%d 和 %d，冻结成功状态依次为 %s、%s、%s、%s。' ...
    '成功判断同时包含内部稳定、20 dB 抑制、未限幅需求不超过 4 和零限幅，不能只按抑制量排名。' ...
    '这里证明的是主导频点上的设计能力，不代表宽带能力。\n\n'], ...
    demo1T1.supp_db, demo2T1.supp_db, demo3T1.supp_db, demo4T1.supp_db, ...
    demo1T1.demand, demo2T1.demand, demo3T1.demand, demo4T1.demand, ...
    demo1T1.clips, demo2T1.clips, demo3T1.clips, demo4T1.clips, ...
    success_word(demo1T1.success), success_word(demo2T1.success), ...
    success_word(demo3T1.success), success_word(demo4T1.success));
if demo5T1.success
    demo5T1Conclusion = ['该冻结结果同时满足 20 dB、未限幅需求和零限幅标准，' ...
        '证明 Case A 固定窄带版本在当前测量 FIR 的 357 Hz 场景有效。'];
else
    demo5T1Conclusion = ['该冻结结果未同时满足 20 dB、未限幅需求和零限幅标准，' ...
        '不能把 Case A 理论适用性误写成当前模型上的成功设计。'];
end
fprintf(fid, ['Demo5 Marino-Tomei T1 为 %.2f dB，需求 %.3f，限幅 %d，判定为 %s；%s\n\n'], ...
    demo5T1.supp_db, demo5T1.demand, demo5T1.clips, ...
    success_word(demo5T1.success), demo5T1Conclusion);
fprintf(fid, '### T2：设计点附近频率变化时谁能继续工作？\n\n');
fprintf(fid, ['Demo1 Q-RLS 为 %.2f dB（%s），Demo3 IMC-FxNLMS 为 %.2f dB（%s），' ...
    'Demo2 normalized-LMS 为 %.2f dB（%s），Demo4 eMOPSO 中央 + Q-RLS 为 %.2f dB（%s）。' ...
    '对比说明，选择与非平稳场景匹配的更新律能够显著改善缓慢频率变化下的性能；' ...
    'Demo2 的 Q-RLS 失败对照与 normalized-LMS 主结果说明更新律是关键解释变量；' ...
    'Demo1 与 Demo4 使用同一条 Q-RLS 自适应律、不同的中央控制器设计，' ...
    '两者的差距可直接归因于中央设计方法（手工扫描 vs 多目标优化）。\n\n'], ...
    demo1T2.supp_db, success_word(demo1T2.success), ...
    demo3T2.supp_db, success_word(demo3T2.success), ...
    demo2T2.supp_db, success_word(demo2T2.success), ...
    demo4T2.supp_db, success_word(demo4T2.success));
    fprintf(fid, ['Demo5 T2 为 %.2f dB，判定为 %s；离散有效次级响应在 300–420 Hz 内跨越零点，' ...
        '因此 Case A 固定符号结构不适用于整个扫频。IMC-FxLMS 基线保留在 Demo5 独立报告中，' ...
    '用于区分窄带内模限制与纯反馈 FIR 自适应能力。\n\n'], ...
    demo5T2.supp_db, success_word(demo5T2.success));
if isfield(stage, 'demo2_t2_comparison')
    comparison = stage.demo2_t2_comparison;
    fprintf(fid, '\n#### Demo2 T2：Q-RLS 与 normalized-LMS 对照\n\n');
    fprintf(fid, ['该对照保持相同反馈信息、相同 LQG 基线和相同评价扫频，' ...
        '只替换在线参数更新律。它用于区分自适应律差异与固定中心频率模型失配，' ...
        '其中 normalized-LMS 是冻结主结果，Q-RLS 作为失败对照保留。\n\n']);
    fprintf(fid, '![Demo2 LMS T2 analysis](figures/cylinder1dm_2k/demo2/t2_lms_analysis.png)\n\n');
    fprintf(fid, '![Demo2 T2 adaptation comparison](figures/cylinder1dm_2k/demo2/t2_adaptation_comparison.png)\n\n');
    [rlsWorst, ~] = local_suppression_summary(comparison.rls.evaluation);
    [lmsWorst, ~] = local_suppression_summary(comparison.lms.evaluation);
    fprintf(fid, ['**图像解读**：收敛图用于判断更新律是否真正改变残差，' ...
        '控制需求图用于检查这种改善是否以超限为代价。Q-RLS 的评价总体抑制 %.2f dB、' ...
        '最差局部 %.2f dB；normalized-LMS 为 %.2f dB、最差局部 %.2f dB。' ...
        '因此 LMS 将整体性能从失效区恢复到通过区，但边缘频段仍未达到 20 dB。\n\n'], ...
        comparison.rls.evaluation.supp_db, rlsWorst, ...
        comparison.lms.evaluation.supp_db, lmsWorst);
    fprintf(fid, '| Method | Candidate | Calibration dB | Evaluation dB | Worst local dB | Demand | Clips |\n');
    fprintf(fid, '|---|---|---:|---:|---:|---:|---:|\n');
    fprintf(fid, '| Q-RLS | %s | %.2f | %.2f | %.2f | %.3f | %d |\n', ...
        comparison.rls.candidate, comparison.rls.calibration.supp_db, ...
        comparison.rls.evaluation.supp_db, rlsWorst, ...
        comparison.rls.evaluation.extra.u_demand_max, ...
        comparison.rls.evaluation.extra.saturation_count);
    fprintf(fid, '| LMS | %s | %.2f | %.2f | %.2f | %.3f | %d |\n\n', ...
        comparison.lms.candidate, comparison.lms.calibration.supp_db, ...
        comparison.lms.evaluation.supp_db, lmsWorst, ...
        comparison.lms.evaluation.extra.u_demand_max, ...
        comparison.lms.evaluation.extra.saturation_count);
    fprintf(fid, ['解释重点：normalized-LMS 恢复了总体扫频抑制，但最差局部分箱仍低于 20 dB；' ...
        '因此主结论是总体有效、边缘频段仍有限，而不是全频带均匀成功。\n\n']);
end
fprintf(fid, '## 8. 总体结论与下一阶段\n\n');
fprintf(fid, '- **已验证**：低阶模型足以支持多种主导频点控制器设计；T2 的成功与失败由冻结评价表逐项给出。\n');
fprintf(fid, '- **独立边界**：100–500 Hz 宽带结论不在本阶段内，见 `BROADBAND_CYLINDER1DM_REPORT.md`。\n');
fprintf(fid, '- **比较边界**：Demo3 使用源参考，Demo1/2/4 使用反馈信息；Demo1 与 Demo4 共享同一自适应律，构成中央设计方法的受控对比。\n');
fprintf(fid, '- **下一步**：针对宽带结构与饱和感知自适应进行修正；Demo4 方向上可扩大 eMOPSO 预算、细化目标频带加权，或在选解评分中显式约束控制需求。\n');
end

function write_scenario_matrix(fid)
fprintf(fid, '| 场景 | 信号设置 | 实验目的 | 主要比较 | 能回答的问题 |\n');
fprintf(fid, '|---|---|---|---|---|\n');
fprintf(fid, '| T1 定频 | 357 Hz 正弦；校准与评价使用不同相位/噪声种子 | 检查控制器能否在模型主导频点实现深抑制 | 开环与闭环波形、PSD、控制需求、稳定性 | 低阶模型是否足以完成基本控制器设计？ |\n');
fprintf(fid, '| T2 扫频 | 校准 300→420 Hz；评价 420→300 Hz | 检查偏离设计点和扫频方向变化后的跟踪能力 | 分段抑制、时频图、固定与自适应结构差异 | 控制器是否只在单一频点有效？ |\n');
fprintf(fid, '\n');
end

function result = summary_result(stage, demo, testName)
index = find(strcmp(stage.summary.demo, demo) & strcmp(stage.summary.test, testName), 1);
if isempty(index)
    error('Missing summary result for %s %s.', demo, testName);
end
result = struct( ...
    'supp_db', stage.summary.evaluation_suppression_db(index), ...
    'demand', stage.summary.control_demand_max(index), ...
    'clips', stage.summary.saturation_count(index), ...
    'success', stage.summary.evaluation_success(index));
end

function [worst, medianValue] = local_suppression_summary(result)
values = result.supp_breakdown.supp_db;
values = values(isfinite(values));
if isempty(values)
    worst = NaN;
    medianValue = NaN;
else
    worst = min(values);
    medianValue = median(values);
end
end

function write_tuning_comparison(fid, stage, demo, testName)
selectionRows = strcmp({stage.selections.demo}, demo) ...
    & strcmp({stage.selections.test}, testName);
indices = find(selectionRows);
for selectionIndex = indices
    selection = stage.selections(selectionIndex);
    rows = strcmp(stage.tuning.demo, demo) ...
        & strcmp(stage.tuning.test, testName) ...
        & strcmp(stage.tuning.variant, selection.variant);
    candidateCount = sum(rows);
    feasibleCount = sum(stage.tuning.feasible(rows));
    selectedRow = rows & stage.tuning.selected;
    fprintf(fid, '### %s %s 参数搜索\n\n', testName, selection.variant);
    fprintf(fid, '- 搜索候选：%d 组；满足稳定性和控制峰值约束：%d 组。\n', ...
        candidateCount, feasibleCount);
    fprintf(fid, '- 选中参数：`%s`，校准集抑制 %.2f dB，未限幅需求 %.3f，限幅 %d 次。\n', ...
        selection.candidate, selection.calibration_result.supp_db, ...
        selection.calibration_result.extra.u_demand_max, ...
        selection.calibration_result.extra.saturation_count);
    if ~selection.feasible
        fprintf(fid, ['- 注意：该场景没有可行候选；当前选择只是惩罚评分最高的失败候选，' ...
            '不能解释为控制器设计成功。\n']);
    end
    if any(selectedRow)
        fprintf(fid, '- 冻结参数：`%s`。\n\n', stage.tuning.params_json{selectedRow});
    else
        fprintf(fid, '\n');
    end
end
end

function write_t1_adaptive_result(fid, stage, demo)
selection = selection_for(stage, demo, 'T1', 'adaptive');
result = evaluation_result_for(stage, demo, 'T1', 'adaptive');
summaryRow = strcmp(stage.summary.demo, demo) ...
    & strcmp(stage.summary.test, 'T1') ...
    & strcmp(stage.summary.variant, 'adaptive');
summaryIndex = find(summaryRow, 1);
fixedResult = evaluation_result_for(stage, demo, 'T1', 'fixed');

fprintf(fid, '### T1 adaptive：定频自适应对照\n\n');
fprintf(fid, ['**作用与目的**：在与固定 T1 相同的 357 Hz 留出信号上启用在线自适应，' ...
    '检查额外自由度是否改善抑制，同时保持控制需求和限幅约束。\n\n']);
fprintf(fid, ['**具体设置**：与 T1 fixed 使用相同评价信号和执行器约束；' ...
    '仅控制器变体切换为冻结的 adaptive 候选。\n\n']);
fprintf(fid, '**比较内容**：与同 Demo 的 T1 fixed 结果直接比较抑制、未限幅需求和限幅次数。\n\n');
fprintf(fid, '![T1 adaptive analysis](figures/cylinder1dm_2k/%s/t1_adaptive_analysis.png)\n\n', demo);
write_t1_figure_interpretation(fid, result);
fprintf(fid, ['**冻结结果**：候选 `%s`；校准集 %.2f dB，评价集 %.2f dB；' ...
    '未限幅需求 %.3f；限幅 %d 次；判定为 **%s**。\n\n'], ...
    selection.candidate, selection.calibration_result.supp_db, result.supp_db, ...
    result.extra.u_demand_max, result.extra.saturation_count, ...
    success_word(stage.summary.evaluation_success(summaryIndex)));
fprintf(fid, ['**结果解释**：相对 fixed 的 %.2f dB，adaptive 为 %.2f dB；' ...
    '该差异必须与控制需求 %.3f 一并判断，不能只按抑制量排名。\n\n'], ...
    fixedResult.supp_db, result.supp_db, result.extra.u_demand_max);
if stage.summary.evaluation_success(summaryIndex)
    conclusion = '该定频自适应候选在留出评价上满足统一成功标准。';
else
    conclusion = '该定频自适应候选未同时满足抑制与控制约束，不能视为有效设计。';
end
fprintf(fid, '**本场景结论**：%s\n\n', conclusion);
end

function write_t1_figure_interpretation(fid, result)
t = result.extra.t(:);
yOpen = result.extra.y_open(:);
yClosed = result.extra.y_closed(:);
uDemand = result.extra.u_demand(:);
uApplied = result.extra.u(:);
windowLength = max(16, round(0.05 * result.fs));
openRms = sqrt(movmean(yOpen.^2, windowLength));
closedRms = sqrt(movmean(yClosed.^2, windowLength));
tailStart = max(1, numel(t) - round(0.1 * numel(t)) + 1);
tailOpen = median(openRms(tailStart:end));
tailClosed = median(closedRms(tailStart:end));
lineDifference = max(abs(uDemand - uApplied));
margin = 4 - result.extra.u_demand_max;

fprintf(fid, '**图像解读**：该 T1 综合图的四个面板分别回答“目标频点是否被压低、频谱峰值是否移动、控制需求是否可实现、全程是否稳定”。');
fprintf(fid, '左上时域面板中，灰色开环扰动与蓝色闭环残差在稳态段的幅值差对应 %.2f dB 总体抑制；', ...
    result.supp_db);
fprintf(fid, '右上 PSD 面板应重点看 357 Hz 主峰，蓝色峰值相对灰色峰值的下降，而不是只看宽频底噪。');
fprintf(fid, '左下控制面板中，绿色实线是实际施加量，橙色虚线是未限幅需求；两者最大差 %.3g，限幅 %d 次，', ...
    lineDifference, result.extra.saturation_count);
if lineDifference < 1e-9 && result.extra.saturation_count == 0
    fprintf(fid, '本例两线重合是“需求全部可施加、没有硬限幅”的证据，不是橙色数据缺失。');
else
    fprintf(fid, '两线的分离部分直接显示限幅压缩发生的时间段。');
end
fprintf(fid, '需求峰值 %.3f，相对调参上限的余量为 %.3f。', ...
    result.extra.u_demand_max, margin);
if isfinite(result.t_conv_s) && result.t_conv_s > 0
    fprintf(fid, '右下移动 RMS 面板显示收敛时间约 %.3f s；末段开环/闭环 RMS 中位数为 %.3g / %.3g。', ...
        result.t_conv_s, tailOpen, tailClosed);
else
    fprintf(fid, '右下移动 RMS 面板显示固定或快速建立后的全程残差包络；末段开环/闭环 RMS 中位数为 %.3g / %.3g。', ...
        tailOpen, tailClosed);
end
fprintf(fid, '\n\n');
end

function write_t2_figure_interpretation(fid, result)
breakdown = result.supp_breakdown;
values = breakdown.supp_db(:);
centers = breakdown.bin_centers(:);
finite = isfinite(values);
values = values(finite);
centers = centers(finite);
if isempty(values)
    worst = NaN; medianValue = NaN; worstFrequency = NaN; belowCount = 0; totalCount = 0;
else
    [worst, worstIndex] = min(values);
    medianValue = median(values);
    worstFrequency = centers(worstIndex);
    belowCount = sum(values < 20);
    totalCount = numel(values);
end

fprintf(fid, '**图像解读**：左上面板是整段留出扫频的开环/闭环 RMS 包络，');
fprintf(fid, '用于判断总体 %.2f dB 是否由整个扫频过程贡献，而不是某个短暂频点。', result.supp_db);
fprintf(fid, ['右上分箱曲线给出局部频率证据：最差分箱约在 %.0f Hz，抑制 %.2f dB，中位数 %.2f dB；' ...
    '%d/%d 个分箱低于 20 dB，因此总体通过不等于全频带均匀通过。'], ...
    worstFrequency, worst, medianValue, belowCount, totalCount);
fprintf(fid, '下方两幅时频图应成对阅读：开环图显示扫频能量脊线，闭环图显示该脊线在哪些频段被压低；');
fprintf(fid, '若闭环脊线在边缘仍清晰，就对应右上曲线中的局部失败。');
fprintf(fid, '控制需求峰值 %.3f、硬限幅 %d 次，是对“频率性能改善是否可实现”的独立约束检查。\n\n', ...
    result.extra.u_demand_max, result.extra.saturation_count);
end

function write_design_figure_interpretation(fid, stage, demo, results)
r1 = results.T1;
r2 = results.T2;
fprintf(fid, '**图像解读**：设计图不是结果截图，而是用来检查控制器结构、稳定性和调参自由度是否解释了 T1/T2 曲线。');
switch demo
    case 'demo1'
        design = r1.extra.design;
        q = r2.extra.q_norm_history(:);
        q = q(isfinite(q));
        feasible = sum(strcmp(stage.tuning.demo, demo) & stage.tuning.feasible);
        total = sum(strcmp(stage.tuning.demo, demo));
        fprintf(fid, '左上/右上显示 RST 灵敏度陷波及其预测抑制 %.2f dB；闭环半径 %.4f、控制器半径 %.4f 均在单位圆内。', ...
            design.suppression_design_db, design.closed_loop_radius, design.controller_radius);
        if isempty(q)
            fprintf(fid, '左下 Q-RLS 历史没有可用有限值；');
        else
            fprintf(fid, '左下 Q-RLS 范数从 %.3g 变化到 %.3g，说明自适应增量的规模；', q(1), q(end));
        end
        fprintf(fid, '右下候选搜索中 %d/%d 个候选满足约束，说明当前选解不是单纯追求抑制峰值。', feasible, total);
    case 'demo2'
        poles = eig(r1.extra.A_cl);
        theta = r2.extra.theta_norm_history(:);
        theta = theta(isfinite(theta));
        feasible = sum(strcmp(stage.tuning.demo, demo) & stage.tuning.feasible);
        total = sum(strcmp(stage.tuning.demo, demo));
        fprintf(fid, '左上 LQG 极点最大半径 %.4f，说明固定基线内部稳定；上中 K/L 增益面板显示控制器与估计器的增益尺度。', max(abs(poles)));
        if isempty(theta)
            fprintf(fid, 'T2 参数范数没有可用有限值；');
        else
            fprintf(fid, 'T2 参数范数从 %.3g 变化到 %.3g，和残差 RMS 的变化共同判断 normalized-LMS 是否在学习而非发散；', theta(1), theta(end));
        end
        fprintf(fid, '候选搜索中 %d/%d 个候选满足约束，右侧控制增量面板用于确认自适应改善没有突破需求上限。', feasible, total);
    case 'demo3'
        frequency = (100:500).';
        flattening = freqz(r1.extra.F_num, r1.extra.F_den, frequency, r1.fs);
        effective = freqz(r1.extra.H_num, r1.extra.H_den, frequency, r1.fs);
        theta = r1.extra.theta_fixed(:);
        coefficientNorm = r2.extra.coefficient_norm(:);
        coefficientNorm = coefficientNorm(isfinite(coefficientNorm));
        feasible = sum(strcmp(stage.tuning.demo, demo) & stage.tuning.feasible);
        total = sum(strcmp(stage.tuning.demo, demo));
        [~, targetIndex] = min(abs(frequency - r1.extra.f_design));
        fprintf(fid, '左上 F(z) 展平后的幅值在设计频带内由 %.2f 到 %.2f dB，有效次级通道在 %.0f Hz 为 %.2f dB；', ...
            min(20*log10(abs(flattening)+eps)), max(20*log10(abs(flattening)+eps)), ...
            r1.extra.f_design, 20*log10(abs(effective(targetIndex))+eps));
        fprintf(fid, '右上固定 FIR 的 tap 范围为 [%.3g, %.3g]，范数 %.3g；', ...
            min(theta), max(theta), norm(theta));
        fprintf(fid, '左下 IMC-FxNLMS 系数范数末值 %.3g，右下 %d/%d 个候选满足约束。', ...
            iff(isempty(coefficientNorm), NaN, coefficientNorm(end)), feasible, total);
    case 'demo4'
        design = r1.extra.design;
        design2 = r2.extra.design;
        feasible = sum(strcmp(stage.tuning.demo, demo) & stage.tuning.feasible);
        total = sum(strcmp(stage.tuning.demo, demo));
        fprintf(fid, '左上灵敏度图显示 Ms=%.2f dB、设计抑制 %.2f dB；右上极点半径为闭环 %.4f、控制器 %.4f，稳定性本身不是失败原因。', ...
            design.Ms_db, design.suppression_design_db, design.closed_loop_radius, design.controller_radius);
        fprintf(fid, '左下 Pareto 图中 T2 选中第 %d 个解，稳定解 %d/%d；右下直接显示 T1/T2 需求均在约束线之上。', ...
            design2.selected_index, design2.stable_solutions, design2.n_pareto_solutions);
        fprintf(fid, '因此图像支持“抑制可以做出来，但控制余量不足”的失败解释，而不是发散解释；%d/%d 个候选满足统一约束。', feasible, total);
    case 'demo5'
        ratio = r1.extra.case_a_ratio(:);
        rawRatio = ratio;
        if isfield(r1.extra, 'case_a_raw_ratio')
            rawRatio = r1.extra.case_a_raw_ratio(:);
        end
        adaptiveT1 = results.T1_adaptive;
        feasible = sum(strcmp(stage.tuning.demo, demo) & stage.tuning.feasible);
        total = sum(strcmp(stage.tuning.demo, demo));
        timing = 'updated';
        if isfield(r1.extra, 'output_timing')
            timing = r1.extra.output_timing;
        end
        fprintf(fid, ['离散有效 Case A 诊断的最小 |Re(S_eff)|/|S_eff| 为 %.3f，' ...
            'raw Re(S) 比值为 %.3f；输出时序为 %s，反馈符号采用 %s 约定。'], ...
            min(ratio), min(rawRatio), timing, r1.extra.feedback_sign);
        fprintf(fid, 'T1 自适应频率估计从 %.2f Hz 变化到 %.2f Hz，T2 估计末值 %.2f Hz，Case A 适用性为 %d。', ...
            adaptiveT1.extra.freq_init_hz(1), adaptiveT1.extra.freq_final_hz(1), ...
            r2.extra.freq_final_hz(1), ...
            r2.extra.case_a_applicable);
        fprintf(fid, ['校准候选中 %d/%d 个满足稳定性、需求和 Case A 约束；' ...
            '结合冻结评价结果，这张图用于区分 T1 fixed 的局部成功、T1 adaptive 的性能不足和 T2 的结构性不适用，' ...
            '而不是把有限抑制笼统地称为成功。'], feasible, total);
end
fprintf(fid, '\n\n');
end

function value = iff(condition, firstValue, secondValue)
if condition
    value = firstValue;
else
    value = secondValue;
end
end

function write_scenario_result(fid, stage, demo, testName)
selection = selection_for(stage, demo, testName);
result = evaluation_result_for(stage, demo, testName);
summaryRow = strcmp(stage.summary.demo, demo) ...
    & strcmp(stage.summary.test, testName) ...
    & strcmp(stage.summary.variant, selection.variant);
summaryIndex = find(summaryRow, 1);
narrative = scenario_narrative(demo, testName);

fprintf(fid, '### %s：%s\n\n', testName, narrative.title);
fprintf(fid, '**作用与目的**：%s\n\n', narrative.purpose);
fprintf(fid, '**具体设置**：%s\n\n', narrative.setup);
fprintf(fid, '**比较内容**：%s\n\n', narrative.comparison);
fprintf(fid, '![%s analysis](figures/cylinder1dm_2k/%s/%s_analysis.png)\n\n', ...
    testName, demo, lower(testName));
if strcmp(testName, 'T1')
    write_t1_figure_interpretation(fid, result);
else
    write_t2_figure_interpretation(fid, result);
end
fprintf(fid, '**冻结结果**：校准集 %.2f dB，评价集 %.2f dB；未限幅需求 %.3f；限幅 %d 次；判定为 **%s**。\n\n', ...
    selection.calibration_result.supp_db, result.supp_db, ...
    result.extra.u_demand_max, result.extra.saturation_count, ...
    success_word(stage.summary.evaluation_success(summaryIndex)));
fprintf(fid, '**结果解释**：%s\n\n', narrative.interpretation);
if strcmp(demo, 'demo5') && strcmp(testName, 'T1')
    if stage.summary.evaluation_success(summaryIndex)
        conclusion = ['当前 Cylinder 1 dm 测量 FIR 与离散内模组合已证明定频 fixed 版本有效；' ...
            '该结论不外推到频率自适应或扫频场景。'];
    else
        conclusion = ['当前 Cylinder 1 dm 测量 FIR 与离散内模组合尚未证明定频 fixed 版本有效；' ...
            '这不是“脚本运行即成功”。'];
    end
else
    conclusion = narrative.conclusion;
end
fprintf(fid, '**本场景结论**：%s\n\n', conclusion);
if strcmp(demo, 'demo2') && strcmp(testName, 'T2') ...
        && isfield(stage, 'demo2_t2_comparison')
    comparison = stage.demo2_t2_comparison;
    fprintf(fid, '#### 更新律对照\n\n');
    fprintf(fid, '![Demo2 LMS T2 analysis](figures/cylinder1dm_2k/demo2/t2_lms_analysis.png)\n\n');
    fprintf(fid, '![Demo2 T2 adaptation comparison](figures/cylinder1dm_2k/demo2/t2_adaptation_comparison.png)\n\n');
    [rlsWorst, ~] = local_suppression_summary(comparison.rls.evaluation);
    [lmsWorst, ~] = local_suppression_summary(comparison.lms.evaluation);
    fprintf(fid, ['**图像解读**：左侧/上方收敛曲线用于判断更新律是否真正改变残差，' ...
        '右侧控制需求曲线用于检查这种改善是否以超限为代价。Q-RLS 的评价总体抑制 %.2f dB、' ...
        '最差局部 %.2f dB；normalized-LMS 为 %.2f dB、最差局部 %.2f dB。' ...
        '因此这里的关键不是“哪条线更高”，而是 LMS 将整体性能从失效区恢复到通过区，' ...
        '同时仍保留边缘频段失效。\n\n'], ...
        comparison.rls.evaluation.supp_db, rlsWorst, ...
        comparison.lms.evaluation.supp_db, lmsWorst);
    fprintf(fid, '| Method | Candidate | Calibration dB | Evaluation dB | Worst local dB | Demand | Clips |\n');
    fprintf(fid, '|---|---|---:|---:|---:|---:|---:|\n');
    fprintf(fid, '| Q-RLS | %s | %.2f | %.2f | %.2f | %.3f | %d |\n', ...
        comparison.rls.candidate, comparison.rls.calibration.supp_db, ...
        comparison.rls.evaluation.supp_db, rlsWorst, ...
        comparison.rls.evaluation.extra.u_demand_max, ...
        comparison.rls.evaluation.extra.saturation_count);
    fprintf(fid, '| LMS | %s | %.2f | %.2f | %.2f | %.3f | %d |\n\n', ...
        comparison.lms.candidate, comparison.lms.calibration.supp_db, ...
        comparison.lms.evaluation.supp_db, lmsWorst, ...
        comparison.lms.evaluation.extra.u_demand_max, ...
        comparison.lms.evaluation.extra.saturation_count);
    fprintf(fid, ['LMS 对照只替换参数更新律；如果局部频率抑制仍集中在中心频段，' ...
        '则拍频式包络的主因仍是中心频率扰动模型失配。\n\n']);
end
end

function selection = selection_for(stage, demo, testName, variant)
if nargin < 4, variant = ''; end
index = find(strcmp({stage.selections.demo}, demo) ...
    & strcmp({stage.selections.test}, testName) ...
    & (isempty(variant) | strcmp({stage.selections.variant}, variant)), 1);
selection = stage.selections(index);
end

function result = evaluation_result_for(stage, demo, testName, variant)
if nargin < 4, variant = ''; end
index = find(strcmp({stage.selections.demo}, demo) ...
    & strcmp({stage.selections.test}, testName) ...
    & (isempty(variant) | strcmp({stage.selections.variant}, variant)), 1);
result = stage.evaluation_results{index};
end

function profile = controller_profile(demo)
switch demo
    case 'demo1'
        profile = struct( ...
            'method', '低阶 RST 中央控制器与 Youla-Q RLS 自适应', ...
            'information', '误差输出反馈，并由对象模型和控制输入重构扰动', ...
            'freedom', 'RST 多项式以及低阶自适应 Q 参数', ...
            'strength', '在指定频点形成深陷波，并通过 Q 参数跟踪有限频率变化', ...
            'limitation', '中央 RST 为窄带设计；宽带能力在独立 T3 阶段评估', ...
            'designFigureMeaning', '设计图同时展示灵敏度陷波、闭环/控制器极点，以及 T2 中 Q-RLS 参数和后验误差的变化，用于确认深抑制不是由不稳定极点造成。', ...
            'finalConclusion', 'Demo1 主实验检验低阶 RST 设计及 Q-RLS 扫频跟踪；宽带结论见独立 T3 报告。');
    case 'demo2'
        profile = struct( ...
            'method', '增广扰动模型 LQG 与 normalized-LMS 自适应（Q-RLS 对照）', ...
            'information', '仅使用误差输出反馈，由 Kalman 估计器重构对象和扰动状态', ...
            'freedom', 'LQG 权重，以及 T2 中在线更新的 normalized-LMS 系数', ...
            'strength', '固定 LQG 处理已知窄带扰动，normalized-LMS 在扫频上保持持续跟踪', ...
            'limitation', '反馈残差缺少源信号预览，边缘频段局部抑制仍低于 20 dB', ...
            'designFigureMeaning', '设计图展示联合闭环极点、LQG 增益、T1 反噪声，以及 T2 的 Q 参数收敛、控制增量和候选搜索，用于区分稳定性、适应增量与控制约束。', ...
            'finalConclusion', 'Demo2 证明低阶状态空间模型足以支持窄带 LQG；T2 主结果采用 normalized-LMS，Q-RLS 作为失败对照保留。');
    case 'demo4'
        profile = struct( ...
            'method', 'epsilon-MOPSO 灵敏度整形中央 RST 与 Youla-Q RLS 自适应', ...
            'information', '误差输出反馈；自适应律与 Demo1 完全相同，仅中央控制器换为多目标优化设计', ...
            'freedom', 'X/Y 整形多项式的根（nX+nY 维决策变量）、可选内模陷波，以及低阶自适应 Q 参数', ...
            'strength', '多目标优化自动在目标频抑制与灵敏度峰值 Ms 之间权衡；与 Demo1 构成"同一自适应律、不同中央设计"的受控对比', ...
            'limitation', '每个中央候选需要一次完整优化，计算量大；宽带能力在独立 T3 阶段评估', ...
            'designFigureMeaning', '设计图展示输出灵敏度 |Syp| 频响（含 Ms 与设计频率标注）、闭环/控制器极点、T2 Pareto 前沿与选中解，以及主实验控制需求，用于确认所选中央解在抑制与鲁棒性之间的位置。', ...
            'finalConclusion', 'Demo4 主实验检验多目标灵敏度整形中央 RST 及其承载 Youla-Q 自适应的能力。');
    case 'demo5'
        profile = struct( ...
            'method', 'Marino-Tomei 自适应频率估计器与纯反馈 IMC-FxLMS 对照', ...
            'information', '两者都只使用误差麦克风；Marino-Tomei 使用 Case A 次级路径符号，IMC-FxLMS 使用次级路径模型重构内部参考', ...
            'freedom', 'Marino-Tomei 的内模频率/振荡器状态；IMC-FxLMS 的 FIR tap', ...
            'strength', '在少量恒频窄带扰动且 Case A 条件可靠时提供低阶频率参数化', ...
            'limitation', '固定符号不能跨越离散有效响应零点的扫频；T1 fixed 的成功不能外推为频率自适应或扫频成功', ...
            'designFigureMeaning', '设计图同时展示 Case A 适用性、频率估计状态、T2 扫频中的符号边界和校准候选搜索；这些证据用于区分理论适用性、参数收敛和实际抑制。', ...
            'finalConclusion', 'Demo5 在当前测量 FIR 上证明了 T1 fixed 的窄带控制能力；T1 adaptive 仍未达到 20 dB，T2 还因 Case A 符号跨零而不适用。');
    otherwise
        profile = struct( ...
            'method', 'F(z)+固定 FIR 与 IMC-FxNLMS 自适应基线', ...
            'information', '使用主噪声源参考信号，并使用次级路径模型进行 filtered-x 更新', ...
            'freedom', '固定 FIR tap 或在线自适应 FIR 系数', ...
            'strength', '源参考和次级路径 filtered-x 更新适合定频及有限扫频跟踪', ...
            'limitation', 'IMC-FxNLMS 与 Jafari F(z) 自适应结构不同，局部频段仍可能失效', ...
            'designFigureMeaning', '设计图展示 F(z) 展平、固定 FIR 系数、IMC-FxNLMS 系数收敛和 T1/T2 参数搜索，用于解释固定与自适应结果的差异。', ...
            'finalConclusion', 'Demo3 主实验采用固定 F(z)+FIR 与 IMC-FxNLMS；原 Jafari 回归量的失败边界见独立 T3 结果。');
end
end

function narrative = scenario_narrative(demo, testName)
switch testName
    case 'T1'
        commonSetup = '357 Hz 定频正弦；校准和评价保持频率一致，但更换相位和随机噪声种子。';
        commonComparison = '比较开环/闭环稳态波形、357 Hz PSD 峰值、控制需求与设计稳定性。';
    case 'T2'
        commonSetup = '校准集采用 300→420 Hz 上扫，评价集采用 420→300 Hz 下扫，避免只记住单一时间轨迹。';
        commonComparison = '比较全时域残差、20 Hz 分箱抑制和开闭环时频图，检查整个扫频过程而非总体 RMS。';
    otherwise
        commonSetup = '100–500 Hz 带限白噪声；校准和评价使用不同随机种子。';
        commonComparison = '比较宽带 PSD、移动窗口抑制、未限幅控制需求和硬限幅次数。';
end

switch demo
    case 'demo1'
        if strcmp(testName, 'T1')
            narrative = make_narrative('固定 RST 的基本设计能力', ...
                '验证低阶多项式对象是否能直接求得稳定 RST，并在目标频点达到深抑制。', ...
                commonSetup, commonComparison, ...
                '候选频点比较显示，控制器需要与该低阶模型的主导频点匹配；冻结参数在评价集保持抑制且无饱和。', ...
                '低阶模型足以完成窄带 RST 控制器设计。');
        elseif strcmp(testName, 'T2')
            narrative = make_narrative('Youla-Q RLS 的频率跟踪能力', ...
                '验证固定 RST 偏离设计点后，自适应 Q 是否能补偿缓慢频率变化。', ...
                commonSetup, commonComparison, ...
                'Q-RLS 在反向留出扫频上仍保持有效抑制，说明性能不是对校准上扫轨迹的简单复现；冻结候选在控制约束内表现最好。', ...
                'Demo1 的自适应增量对有限频率漂移有效。');
        else
            narrative = make_narrative('RST/Youla 结构的宽带边界', ...
                '检验低阶 Q 是否能把窄带中央控制器扩展为 100–500 Hz 宽带控制器。', ...
                commonSetup, commonComparison, ...
                '所有 Q 阶数候选都不可行；冻结结果显示抑制不足且控制需求超过约束并触发限幅。提高 Q 阶数还造成参数膨胀。', ...
                '失败主要是结构自由度和饱和感知不足，不是继续微调 lambda 即可解决。');
        end
    case 'demo2'
        if strcmp(testName, 'T1')
            narrative = make_narrative('增广 LQG 的窄带控制能力', ...
                '验证低阶对象与二阶正弦扰动模型能否构成有效的 LQG 控制器。', ...
                commonSetup, commonComparison, ...
                '激进权重相较默认权重显著改善抑制，同时控制需求仍处于约束内，说明原默认效果差主要是权重未适配。', ...
                '低阶状态空间模型可以支持高性能窄带 LQG。');
        elseif strcmp(testName, 'T2')
            narrative = make_narrative('normalized-LMS 的扫频补偿能力', ...
                '检查在线 normalized-LMS 能否在不破坏固定 LQG 稳定性的前提下持续补偿中心频率扰动模型与连续扫频之间的失配。', ...
                commonSetup, commonComparison, ...
                ['主结果使用 normalized-LMS 生成 filtered-x 回归量；' ...
                '扫频频率逐渐离开中心模型时，反噪声与原扰动的相位差连续变化，' ...
                '两者叠加后形成拍频式残差包络。normalized-LMS 可持续跟踪该失配，' ...
                '但边缘频率仍会出现局部残差。'], ...
                'normalized-LMS 能恢复总体扫频抑制，但边缘局部分箱仍低于 20 dB；Q-RLS 失败对照说明更新律选择仍是重要变量。');
        else
            narrative = make_narrative('低阶随机扰动模型的宽带能力', ...
                '检验四阶带通噪声模型和固定 Kalman/LQR 权重能否抑制随机宽带实现。', ...
                commonSetup, commonComparison, ...
                '控制需求未超限且无饱和，但抑制接近零，说明控制器过于保守且扰动状态模型无法预测具体随机波形。', ...
                '当前 LQG 宽带失败是扰动建模和权重结构问题，不是稳定性问题。');
        end
    case 'demo4'
        if strcmp(testName, 'T1')
            narrative = make_narrative('eMOPSO 灵敏度整形的定频设计能力', ...
                '验证多目标优化能否自动找到带阻尼内模陷波、同时满足 Ms 约束的稳定 RST 参数。', ...
                commonSetup, commonComparison, ...
                'Pareto 存档中按“陷波频率抑制 − Ms 超限惩罚”确定性选解；冻结解的闭环与控制器极点均严格位于单位圆内。', ...
                '冻结解达到抑制阈值但控制需求超过 4，因此当前 eMOPSO T1 设计无效。');
        elseif strcmp(testName, 'T2')
            narrative = make_narrative('eMOPSO 中央控制器上的 Youla-Q 扫频补偿', ...
                '在优化设计的中央 RST 上运行与 Demo1 相同的 Youla-Q 自适应律，隔离比较中央控制器设计方法对扫频跟踪的影响。', ...
                commonSetup, commonComparison, ...
                '自适应律、Q 结构与 RLS 参数网格都与 Demo1 保持一致，唯一差异是中央控制器由 eMOPSO 整形而非手工扫描设计；与 Demo1 T2 的差距即可归因于中央设计。', ...
                '虽然扫频抑制较高，但控制需求超过 4，因此当前 eMOPSO/Youla T2 设计无效。');
        else
            narrative = make_narrative('整形中央控制器 + 低阶 Q 的宽带边界', ...
                '检验 Ms 约束下的整形中央控制器叠加低阶 Q 后，能在 100–500 Hz 上取得多少宽带抑制。', ...
                commonSetup, commonComparison, ...
                'Bode 水床定理与低阶 Q 自由度共同限制宽带反馈抑制；该结果与 Demo1 T3 相互印证：瓶颈在结构自由度，而不是中央控制器的设计方法。', ...
                'T3 的结果界定的是"反馈 + 低阶 Q"这一结构的物理边界，而不是优化算法的缺陷。');
        end
    case 'demo5'
        if strcmp(testName, 'T1')
            narrative = make_narrative('Marino-Tomei 的定频 Case A 设计能力', ...
                '验证误差麦克风反馈、窄带内模和次级路径符号条件能否在 357 Hz 主导频点形成可实现控制。', ...
                commonSetup, commonComparison, ...
                '校准搜索同时记录原始 Re{S}、离散有效反馈符号、频率估计和未限幅需求；最终是否成功由冻结评价的抑制、需求和限幅共同判定。', ...
                '定频 fixed 结论由冻结评价成功状态动态生成，不能从 Case A 适用性单独推断。');
        elseif strcmp(testName, 'T2')
            narrative = make_narrative('Case A 有效响应跨零的扫频适用性边界', ...
                '检查 300–420 Hz 扫频是否满足固定离散有效符号的理论前提，同时保留频率估计轨迹作为诊断证据。', ...
                commonSetup, commonComparison, ...
                '离散有效次级响应在扫频区间内跨过零点，Case A 适用性标记为 false；因此任何有限抑制都不能解释为 Marino-Tomei 对全扫频的成功跟踪。IMC-FxLMS 作为独立基线继续运行。', ...
                '当前 Case A 实现不适用于整个 T2 频带；需要 Case B 或相位补偿结构后才能重新进行公平扫频评价。');
        else
            narrative = make_narrative('Marino-Tomei 的宽带边界', ...
                '检验少量恒频内模是否能覆盖 100–500 Hz 宽带扰动。', ...
                commonSetup, commonComparison, ...
                '该方法没有宽带随机扰动的显式建模自由度，因此宽带结果只用于失败边界记录。', ...
                '宽带能力不由当前 Demo5 主阶段宣称。');
        end
    otherwise
        if strcmp(testName, 'T1')
            narrative = make_narrative('F(z)+固定 FIR 的模型逆设计能力', ...
                '验证低阶次级路径模型是否足以设计一个因果 FIR，在主导频点产生幅相匹配的反噪声。', ...
                commonSetup, commonComparison, ...
                '不同 FIR 阶数候选的性能接近，最终选择最简单的候选；这说明更高阶 FIR 并非该定频场景的必要条件。', ...
                '低阶次级路径模型足以支持紧凑的固定 FIR 设计。');
        elseif strcmp(testName, 'T2')
            narrative = make_narrative('IMC-FxNLMS 的扫频跟踪能力', ...
                '验证具有源参考和次级路径估计的前馈自适应滤波器能否跟踪频率变化。', ...
                commonSetup, commonComparison, ...
                '冻结 IMC-FxNLMS 候选在评价下扫中达到总体成功阈值且无饱和；分段曲线显示边缘频段仍低于 20 dB，因此结论限定为有限频带有效。', ...
                'IMC-FxNLMS 在总体指标上达到 20 dB，但最差局部分箱仍低于 20 dB，结论限定为有限频带有效。');
        else
            narrative = make_narrative('NLMS-on-FIR 的宽带性能与控制代价', ...
                '检验增加 FIR 自由度和源参考后，是否能在控制约束内实现宽带抑制。', ...
                commonSetup, commonComparison, ...
                '图像显示有部分改善，但冻结结果违反控制需求或限幅约束；所有候选均不可行。', ...
                'Demo3 只是比 Demo1/2 获得更多部分抑制，仍未形成可接受的宽带控制器。');
        end
end
end

function value = make_narrative(titleText, purpose, setup, comparison, interpretation, conclusion)
value = struct('title', titleText, 'purpose', purpose, 'setup', setup, ...
    'comparison', comparison, 'interpretation', interpretation, ...
    'conclusion', conclusion);
end

function value = success_word(flag)
if flag
    value = '成功';
else
    value = '失败';
end
end

function value = yes_no(flag)
if flag
    value = 'yes';
else
    value = 'no';
end
end
