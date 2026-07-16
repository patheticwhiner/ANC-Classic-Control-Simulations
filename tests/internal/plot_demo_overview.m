function details = plot_demo_overview(results, outputPath)
%PLOT_DEMO_OVERVIEW Plot the frozen T1/T2 comparison for one Demo.

required = {'T1_fixed', 'T1_adaptive', 'T2'};
for index = 1:numel(required)
    if ~isfield(results, required{index})
        error('plot_demo_overview:MissingResult', ...
            'Missing result field %s.', required{index});
    end
end

fixed = results.T1_fixed;
adaptive = results.T1_adaptive;
chirp = results.T2;
if ~strcmp(fixed.demo, adaptive.demo) || ~strcmp(fixed.demo, chirp.demo)
    error('plot_demo_overview:MixedDemo', ...
        'Overview inputs must belong to the same Demo.');
end

[adaptationTime, adaptationValue, adaptationLabel] = ...
    adaptation_trace(adaptive);
if isempty(adaptationValue) || ~any(isfinite(adaptationValue))
    error('plot_demo_overview:EmptyAdaptationTrace', ...
        '%s has no finite adaptation history.', upper(fixed.demo));
end

options = test_display_options();
fig = figure('Visible', options.figure_visibility, 'Color', 'w', ...
    'Position', [40 40 1440 860]);
layout = tiledlayout(fig, 2, 2, ...
    'TileSpacing', 'compact', 'Padding', 'compact');

% T1 fixed/adaptive comparison.
ax = nexttile(layout);
t = fixed.extra.t(:);
windowLength = max(8, round(0.2 * fixed.fs));
openRms = moving_rms(fixed.extra.y_open(:), windowLength);
fixedRms = moving_rms(fixed.extra.y_closed(:), windowLength);
adaptiveRms = moving_rms(adaptive.extra.y_closed(:), windowLength);
plot_reduced(ax, t, openRms, ...
    'Color', [0.65 0.65 0.65], 'DisplayName', 'Open loop');
hold(ax, 'on');
plot_reduced(ax, t, fixedRms, ...
    'Color', [0.00 0.36 0.62], 'LineWidth', 1.2, ...
    'DisplayName', 'T1 fixed');
plot_reduced(ax, t, adaptiveRms, ...
    'Color', [0.80 0.33 0.10], 'LineWidth', 1.2, ...
    'DisplayName', 'T1 adaptive');
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, '200 ms moving RMS');
title(ax, sprintf('T1 357 Hz: fixed %.1f dB | adaptive %.1f dB', ...
    fixed.supp_db, adaptive.supp_db));
legend(ax, 'Location', 'best');

% Controller-specific adaptation state. Missing fields fail loudly above.
ax = nexttile(layout);
plot(ax, adaptationTime, adaptationValue, ...
    'Color', [0.80 0.33 0.10], 'LineWidth', 1.3);
hold(ax, 'on');
if isfinite(adaptive.t_conv_s) ...
        && adaptive.t_conv_s >= adaptationTime(1) ...
        && adaptive.t_conv_s <= adaptationTime(end)
    xline(ax, adaptive.t_conv_s, '--', 't_{conv}', ...
        'LabelOrientation', 'horizontal');
end
grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, adaptationLabel);
title(ax, sprintf('T1 adaptation state (t_{conv} = %.2f s)', ...
    adaptive.t_conv_s));

% T2 local evidence exposes edge-band failures hidden by overall RMS.
ax = nexttile(layout);
breakdown = chirp.supp_breakdown;
values = breakdown.supp_db(:);
finiteValues = values(isfinite(values));
if isempty(finiteValues)
    error('plot_demo_overview:EmptyT2Breakdown', ...
        '%s T2 has no finite local-suppression values.', upper(fixed.demo));
end
plot(ax, breakdown.bin_centers, values, '-o', ...
    'Color', [0.00 0.36 0.62], 'LineWidth', 1.3, 'MarkerSize', 5);
hold(ax, 'on'); yline(ax, 20, '--', '20 dB'); yline(ax, 0, ':');
grid(ax, 'on'); xlabel(ax, 'Instantaneous frequency (Hz)');
ylabel(ax, 'Local suppression (dB)'); xlim(ax, [300 420]);
title(ax, sprintf('T2: %.1f dB overall, %.1f dB worst local', ...
    chirp.supp_db, min(finiteValues)));

% Raw demand is reduced to block peaks so thresholds remain readable.
ax = nexttile(layout);
hold(ax, 'on');
plot_demand_peak(ax, fixed, [0.00 0.36 0.62], 'T1 fixed');
plot_demand_peak(ax, adaptive, [0.80 0.33 0.10], 'T1 adaptive');
plot_demand_peak(ax, chirp, [0.00 0.45 0.30], 'T2 adaptive');
yline(ax, 4, '--', 'Tuning bound', 'HandleVisibility', 'off');
yline(ax, 5, ':', 'Hard limit', 'HandleVisibility', 'off');
maxDemand = max([fixed.extra.u_demand_max, ...
    adaptive.extra.u_demand_max, chirp.extra.u_demand_max]);
ylim(ax, [0 max(5.5, 1.05 * maxDemand)]);
grid(ax, 'on'); xlabel(ax, 'Time (s)');
ylabel(ax, 'Block peak |raw control demand|');
title(ax, sprintf('Control demand: max %.3f, total clips %d', ...
    maxDemand, fixed.extra.saturation_count ...
    + adaptive.extra.saturation_count + chirp.extra.saturation_count));
legend(ax, 'Location', 'best');

sgtitle(layout, sprintf('%s: T1 fixed %.1f -> adaptive %.1f dB | T2 %.1f dB', ...
    upper(fixed.demo), fixed.supp_db, adaptive.supp_db, chirp.supp_db));

outputDir = fileparts(outputPath);
if ~isempty(outputDir) && ~isfolder(outputDir), mkdir(outputDir); end
drawnow;
exportgraphics(fig, outputPath, 'Resolution', 180);
if options.close_after_save, close(fig); end

details = struct( ...
    'demo', fixed.demo, ...
    'adaptation_sample_count', numel(adaptationValue), ...
    'adaptation_min', min(adaptationValue(isfinite(adaptationValue))), ...
    'adaptation_max', max(adaptationValue(isfinite(adaptationValue))), ...
    't2_worst_local_db', min(finiteValues), ...
    'max_control_demand', maxDemand);
end

function [time, value, label] = adaptation_trace(result)
switch result.demo
    case {'demo1', 'demo4'}
        value = required_vector(result.extra, 'q_norm_history');
        label = '||Q||';
        time = history_time(result, value);
    case 'demo2'
        value = required_vector(result.extra, 'theta_norm_history');
        label = '||theta||';
        time = history_time(result, value);
    case 'demo3'
        value = required_vector(result.extra, 'coefficient_norm');
        time = required_vector(result.extra, 'coefficient_time');
        label = '||theta||';
        if numel(time) ~= numel(value)
            error('plot_demo_overview:HistorySizeMismatch', ...
                'Demo3 coefficient time and norm lengths differ.');
        end
    case 'demo5'
        history = result.extra.f_est_log;
        if isempty(history) || size(history, 1) < 1
            error('plot_demo_overview:InvalidHistoryField', ...
                'Demo5 frequency estimate history is empty.');
        end
        value = history(1, :);
        time = result.extra.t_log;
        label = 'Estimated frequency (Hz)';
    otherwise
        error('plot_demo_overview:UnknownDemo', ...
            'Unsupported Demo %s.', result.demo);
end
time = time(:);
value = value(:);
end

function value = required_vector(extra, fieldName)
if ~isfield(extra, fieldName)
    error('plot_demo_overview:MissingHistoryField', ...
        'Missing adaptation history field extra.%s.', fieldName);
end
value = extra.(fieldName);
if isempty(value) || ~isvector(value)
    error('plot_demo_overview:InvalidHistoryField', ...
        'Adaptation history extra.%s must be a nonempty vector.', fieldName);
end
end

function time = history_time(result, value)
if isfield(result.extra, 't') && numel(result.extra.t) == numel(value)
    time = result.extra.t;
else
    time = linspace(0, result.Tsim, numel(value));
end
end

function plot_reduced(ax, x, y, varargin)
step = max(1, floor(numel(x) / 6000));
plot(ax, x(1:step:end), y(1:step:end), varargin{:});
end

function value = moving_rms(signal, windowLength)
signal = signal(:);
value = sqrt(movmean(signal .^ 2, windowLength));
end

function plot_demand_peak(ax, result, color, label)
[time, peak] = block_peak_envelope( ...
    result.extra.t(:), result.extra.u_demand(:), 500);
plot(ax, time, peak, 'Color', color, 'LineWidth', 1.1, ...
    'DisplayName', label);
end

function [time, peak] = block_peak_envelope(t, signal, blockCount)
edges = unique(round(linspace(1, numel(signal) + 1, blockCount + 1)));
time = zeros(numel(edges) - 1, 1);
peak = zeros(numel(edges) - 1, 1);
for index = 1:numel(edges) - 1
    samples = edges(index):edges(index + 1) - 1;
    time(index) = mean(t(samples));
    peak(index) = max(abs(signal(samples)));
end
end
