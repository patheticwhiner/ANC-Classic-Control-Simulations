function fig = plot_anc_migration_result(migration, axes_handles)
%PLOT_ANC_MIGRATION_RESULT Plot held-out residual and tuning evidence.

if nargin < 2 || isempty(axes_handles)
    fig = figure('Name', 'ANC model migration result', 'Color', 'w');
    layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact');
    axes_handles = gobjects(1,4);
    for index = 1:4, axes_handles(index) = nexttile(layout); end
else
    fig = ancestor(axes_handles(1), 'figure');
end

r = migration.evaluation.result;
signal = migration.evaluation.signals.T1;
t = signal.t(:);
open = signal.y_open(:);
closed = r.extra.y_closed(:);
step = max(1, ceil(numel(t)/20000));
idx = 1:step:numel(t);

ax = axes_handles(1); cla(ax);
plot(ax, t(idx), open(idx), 'Color', [0.55 0.55 0.55]); hold(ax, 'on');
plot(ax, t(idx), closed(idx), 'b'); hold(ax, 'off'); grid(ax, 'on');
xlabel(ax, 'Time (s)'); ylabel(ax, 'Residual');
title(ax, sprintf('Held-out: %.2f dB', r.supp_db));
legend(ax, 'Open', 'Closed', 'Location', 'best');

ax = axes_handles(2); cla(ax);
[p_open, f] = pwelch(open, [], [], [], migration.models.fs);
[p_closed, ~] = pwelch(closed, [], [], [], migration.models.fs);
plot(ax, f, 10*log10(p_open+eps), 'Color', [0.55 0.55 0.55]); hold(ax, 'on');
plot(ax, f, 10*log10(p_closed+eps), 'b'); hold(ax, 'off'); grid(ax, 'on');
xlim(ax, [0 min(migration.models.fs/2, max(500, 1.2*migration.noise.high_hz))]);
xlabel(ax, 'Frequency (Hz)'); ylabel(ax, 'PSD (dB/Hz)'); title(ax, 'Held-out PSD');

records = migration.calibration.candidates;
ax = axes_handles(3); cla(ax);
scores = [records.score];
bar(ax, scores); grid(ax, 'on');
xticks(ax, 1:numel(records)); xticklabels(ax, {records.name});
xtickangle(ax, 25); ylabel(ax, 'Selection score'); title(ax, 'Calibration candidates');

ax = axes_handles(4); cla(ax);
if isfield(r.extra, 'u_demand')
    demand = r.extra.u_demand(:);
    plot(ax, t(idx), demand(idx)); hold(ax, 'on');
    yline(ax, migration.config.design.tuning_demand_limit, '--r');
    yline(ax, -migration.config.design.tuning_demand_limit, '--r');
    hold(ax, 'off'); grid(ax, 'on'); xlabel(ax, 'Time (s)'); ylabel(ax, 'u demand');
else
    text(ax, 0.5, 0.5, 'No demand trace', 'HorizontalAlignment', 'center');
end
title(ax, sprintf('Demand %.3g | clips %g | pass %d', ...
    migration.evaluation.diagnostics.demand_max, ...
    migration.evaluation.diagnostics.saturation_count, ...
    migration.evaluation.passed));
end
