function fig = anc_migration_designer_gui(varargin)
%ANC_MIGRATION_DESIGNER_GUI Model, scenario, controller and tuning editor.
%
%   anc_migration_designer_gui
%   fig = anc_migration_designer_gui('Visible','off')
%   fig = anc_migration_designer_gui('ControllerIds', ...
%       {'demo1_rst_fixed','demo4_emopso_rst'}, ...
%       'DeploymentExporter', @my_exporter)

parser = inputParser;
addParameter(parser, 'Visible', 'on', @(x) any(strcmpi(x, {'on','off'})));
addParameter(parser, 'ModelRoot', '', @(x) ischar(x) || isstring(x));
addParameter(parser, 'ControllerIds', {}, ...
    @(x) iscell(x) || isstring(x) || ischar(x));
addParameter(parser, 'WindowName', 'ANC Model Migration Designer', ...
    @(x) ischar(x) || isstring(x));
addParameter(parser, 'DeploymentExporter', [], ...
    @(x) isempty(x) || isa(x, 'function_handle'));
addParameter(parser, 'DeploymentDirectory', '', ...
    @(x) ischar(x) || isstring(x));
parse(parser, varargin{:});

tool_dir = fileparts(mfilename('fullpath'));
addpath(tool_dir);
defaults = anc_migration_config();
if ~isempty(parser.Results.ModelRoot)
    defaults.model_root = char(parser.Results.ModelRoot);
end
catalog = scan_anc_migration_models(defaults.model_root);
pairs = pair_anc_migration_models(catalog, defaults.models.sample_rate_tolerance);
if isempty(pairs)
    error('anc_migration_designer_gui:noPairs', ...
        'No compatible primary/secondary pairs were found under %s.', ...
        defaults.model_root);
end
controllers = select_controllers(anc_controller_catalog(), ...
    parser.Results.ControllerIds);
if ~ismember(defaults.design.controller_id, {controllers.id})
    defaults.design.controller_id = controllers(1).id;
end
default_pair = select_default_pair(pairs);
deployment_exporter = parser.Results.DeploymentExporter;
deployment_directory = char(parser.Results.DeploymentDirectory);
has_deployment_export = ~isempty(deployment_exporter);

screen = get(groot, 'ScreenSize');
window_width = min(1480, max(1100, screen(3)-100));
window_height = min(920, max(720, screen(4)-140));
position = [max(20, (screen(3)-window_width)/2), ...
    max(30, (screen(4)-window_height)/2), window_width, window_height];
fig = uifigure('Name', char(parser.Results.WindowName), ...
    'Position', position, 'Visible', parser.Results.Visible, ...
    'Tag', 'anc_migration_designer_gui');
fig.UserData = struct('result', [], 'pairs', pairs, 'catalog', catalog);

outer = uigridlayout(fig, [3 1]);
outer.RowHeight = {326, 72, '1x'};
outer.Padding = [8 8 8 8];
outer.RowSpacing = 6;

top_sections = uigridlayout(outer, [1 2]);
top_sections.Layout.Row = 1;
top_sections.ColumnWidth = {'1x', '1x'};
top_sections.Padding = [0 0 0 0];
top_sections.ColumnSpacing = 8;
experiment_panel = uipanel(top_sections, 'Title', 'Experiment setup', ...
    'Tag', 'migration_experiment_panel');
parameter_panel = uipanel(top_sections, 'Title', 'Controller parameters', ...
    'Tag', 'migration_parameter_panel');

%% Compact experiment section.
form = uigridlayout(experiment_panel, [10 4]);
form.ColumnWidth = {105, '1x', 105, '1x'};
form.RowHeight = repmat({'1x'}, 1, 10);
form.Padding = [6 5 6 5];
form.RowSpacing = 4;
form.ColumnSpacing = 6;

add_label(form, 'Model pair', 1, 1);
pair_drop = uidropdown(form, 'Items', {pairs.label}, ...
    'ItemsData', {pairs.id}, 'Value', pairs(default_pair).id, ...
    'Tag', 'migration_model_pair', 'ValueChangedFcn', @model_changed);
place(pair_drop, 1, [2 4]);

add_label(form, 'Primary model', 2, 1);
primary_text = uieditfield(form, 'text', 'Editable', 'off', ...
    'Tag', 'migration_primary'); place(primary_text, 2, [2 4]);
add_label(form, 'Secondary model', 3, 1);
secondary_text = uieditfield(form, 'text', 'Editable', 'off', ...
    'Tag', 'migration_secondary'); place(secondary_text, 3, [2 4]);

add_label(form, 'Controller scheme', 4, 1);
controller_drop = uidropdown(form, 'Items', {controllers.label}, ...
    'ItemsData', {controllers.id}, 'Value', defaults.design.controller_id, ...
    'Tag', 'migration_controller', 'ValueChangedFcn', @controller_changed);
place(controller_drop, 4, 2);
add_label(form, 'Tuning profile', 4, 3);
profile_drop = uidropdown(form, 'Items', {'Quick','Standard'}, ...
    'ItemsData', {'quick','standard'}, 'Value', defaults.design.profile, ...
    'Tag', 'migration_profile', 'ValueChangedFcn', @design_context_changed);
place(profile_drop, 4, 4);

add_label(form, 'Noise / scenario', 5, 1);
noise_drop = uidropdown(form, ...
    'Items', {'Fixed sine','Linear chirp','Band-limited noise', ...
        'White noise','Multisine','Audio file'}, ...
    'ItemsData', {'fixed_sine','linear_chirp','bandlimited_noise', ...
        'white_noise','multisine','file'}, ...
    'Value', defaults.noise.type, 'Tag', 'migration_noise_type', ...
    'ValueChangedFcn', @noise_changed);
place(noise_drop, 5, 2);
add_label(form, 'Sample rate', 5, 3);
fs_text = uieditfield(form, 'text', 'Editable', 'off', ...
    'Tag', 'migration_fs'); place(fs_text, 5, 4);

add_label(form, 'Duration (s)', 6, 1);
duration_field = uieditfield(form, 'numeric', 'Value', defaults.noise.duration, ...
    'Limits', [0.2 Inf], 'Tag', 'migration_duration', ...
    'ValueChangedFcn', @design_context_changed); place(duration_field, 6, 2);
add_label(form, 'Disturbance RMS', 6, 3);
amplitude_field = uieditfield(form, 'numeric', ...
    'Value', defaults.noise.amplitude_rms, 'Limits', [eps Inf], ...
    'Tag', 'migration_amplitude'); place(amplitude_field, 6, 4);

frequency_label = add_label(form, 'Design/tone Hz', 7, 1);
frequency_field = uieditfield(form, 'numeric', ...
    'Value', defaults.noise.frequency_hz, 'Limits', [eps Inf], ...
    'Tag', 'migration_frequency', 'ValueChangedFcn', @design_context_changed);
place(frequency_field, 7, 2);
low_label = add_label(form, 'Band low Hz', 7, 3);
low_field = uieditfield(form, 'numeric', 'Value', defaults.noise.low_hz, ...
    'Limits', [eps Inf], 'Tag', 'migration_low', ...
    'ValueChangedFcn', @design_context_changed); place(low_field, 7, 4);

frequencies_label = add_label(form, 'Frequencies Hz', 8, 1);
frequencies_field = uieditfield(form, 'text', ...
    'Value', strtrim(sprintf('%g ', defaults.noise.frequencies_hz)), ...
    'Tag', 'migration_frequencies', 'ValueChangedFcn', @design_context_changed);
place(frequencies_field, 8, 2);
high_label = add_label(form, 'Band high Hz', 8, 3);
high_field = uieditfield(form, 'numeric', 'Value', defaults.noise.high_hz, ...
    'Limits', [eps Inf], 'Tag', 'migration_high', ...
    'ValueChangedFcn', @design_context_changed); place(high_field, 8, 4);

file_label = add_label(form, 'Audio file', 9, 1);
file_field = uieditfield(form, 'text', 'Tag', 'migration_audio_file');
place(file_field, 9, [2 3]);
file_button = uibutton(form, 'Text', 'Browse...', ...
    'ButtonPushedFcn', @browse_audio); place(file_button, 9, 4);

add_label(form, 'Output MAT', 10, 1);
output_field = uieditfield(form, 'text', 'Value', '', ...
    'Tag', 'migration_output'); place(output_field, 10, [2 3]);
output_button = uibutton(form, 'Text', 'Choose...', ...
    'ButtonPushedFcn', @choose_output); place(output_button, 10, 4);

%% Compact controller-parameter section; the table scrolls for long schemas.
parameter_layout = uigridlayout(parameter_panel, [3 1]);
parameter_layout.RowHeight = {34, '1x', 48};
parameter_layout.Padding = [6 5 6 5];
parameter_layout.RowSpacing = 4;
parameter_toolbar = uigridlayout(parameter_layout, [1 4]);
parameter_toolbar.ColumnWidth = {190, 125, '1x', 120};
parameter_toolbar.Padding = [0 0 0 0];
auto_tune_box = uicheckbox(parameter_toolbar, ...
    'Text', 'Auto tune candidate sweep', 'Value', defaults.design.auto_tune, ...
    'Tag', 'migration_auto_tune', 'ValueChangedFcn', @design_context_changed);
parameter_count_label = uilabel(parameter_toolbar, 'Text', '', ...
    'Tag', 'migration_parameter_count');
uilabel(parameter_toolbar, ...
    'Text', 'Tick Override to apply an edited value.', ...
    'HorizontalAlignment', 'center');
uibutton(parameter_toolbar, 'Text', 'Reset overrides', ...
    'Tag', 'migration_reset_parameters', 'ButtonPushedFcn', @reset_parameters);

parameter_table = uitable(parameter_layout, ...
    'Tag', 'migration_parameter_table', ...
    'ColumnName', {'Override','Parameter','Value','Unit','Implementation meaning'}, ...
    'ColumnEditable', [true false true false false], ...
    'ColumnWidth', {75, 165, 140, 95, 'auto'}, ...
    'CellEditCallback', @parameter_edited);
parameter_footer = uitextarea(parameter_layout, 'Editable', 'off', ...
    'Tag', 'migration_parameter_help', ...
    'Value', {'Manual overrides are validated before simulation.'});

%% Results remain visible on the same page.
plot_grid = uigridlayout(outer, [2 2]);
plot_grid.Layout.Row = 3;
plot_grid.Padding = [0 0 0 0];
axes_handles = gobjects(1,4);
for axis_index = 1:4
    axes_handles(axis_index) = uiaxes(plot_grid, ...
        'Tag', sprintf('migration_ax%d', axis_index));
    grid(axes_handles(axis_index), 'on');
end
title(axes_handles(1), 'Scenario preview');
title(axes_handles(2), 'Spectrum');
title(axes_handles(3), 'Calibration candidates');
title(axes_handles(4), 'Control demand');

%% Persistent action/status strip.
if has_deployment_export
    action_columns = 6;
    action = uigridlayout(outer, [2 action_columns]);
    action.ColumnWidth = {'1x','1x','1.4x','1x','1x','1x'};
else
    action_columns = 5;
    action = uigridlayout(outer, [2 action_columns]);
    action.ColumnWidth = {'1x','1x','1.4x','1x','1x'};
end
action.Layout.Row = 2;
action.RowHeight = {30, 30};
preview_button = uibutton(action, 'Text', 'Preview scenario', ...
    'Tag', 'migration_preview', 'ButtonPushedFcn', @preview_scenario);
preview_button.Layout.Column = 1;
run_button = uibutton(action, 'Text', 'Auto tune + held-out evaluation', ...
    'Tag', 'migration_run', 'ButtonPushedFcn', @run_migration);
run_button.Layout.Column = [2 3];
save_button = uibutton(action, 'Text', 'Save result copy', ...
    'Tag', 'migration_save', 'Enable', 'off', 'ButtonPushedFcn', @save_copy);
save_button.Layout.Column = 4;
if has_deployment_export
    deployment_button = uibutton(action, 'Text', 'Export fixed R/S', ...
        'Tag', 'migration_export_deployment', 'Enable', 'off', ...
        'ButtonPushedFcn', @export_deployment);
    deployment_button.Layout.Column = 5;
else
    deployment_button = gobjects(0);
end
overwrite_box = uicheckbox(action, 'Text', 'Allow overwrite', ...
    'Value', false, 'Tag', 'migration_overwrite');
overwrite_box.Layout.Column = action_columns;
status = uitextarea(action, 'Editable', 'off', 'Value', {'Ready.'}, ...
    'Tag', 'migration_status');
status.Layout.Row = 2; status.Layout.Column = [1 action_columns];

manual_store = struct();
current_controller_id = char(controller_drop.Value);
parameter_schema = struct([]);
apply_pair(pairs(default_pair));
update_noise_controls();
refresh_parameter_table(false);
update_compatibility();

    function model_changed(~,~)
        save_current_parameter_table();
        apply_pair(selected_pair());
        refresh_parameter_table(true);
        invalidate_result();
    end

    function pair = selected_pair()
        pair_index = find(strcmp({pairs.id}, char(pair_drop.Value)), 1);
        pair = pairs(pair_index);
    end

    function apply_pair(pair)
        primary_text.Value = pair.primary.file;
        secondary_text.Value = pair.secondary.file;
        fs_text.Value = sprintf('%.6g Hz (%s)', pair.fs, upper(pair.kind));
        nyquist_guard = 0.45*pair.fs;
        high_field.Value = min(high_field.Value, nyquist_guard);
        frequency_field.Value = min(frequency_field.Value, nyquist_guard);
        frequencies = parse_numeric_list(frequencies_field.Value);
        frequencies = frequencies(frequencies < pair.fs/2);
        if isempty(frequencies), frequencies = min(357, nyquist_guard); end
        frequencies_field.Value = strtrim(sprintf('%g ', frequencies));
        status.Value = {sprintf('Selected %s model pair at %.6g Hz.', ...
            upper(pair.kind), pair.fs); ...
            sprintf('Primary fit %.2f%% | secondary fit %.2f%%', ...
            pair.primary.fit_sim, pair.secondary.fit_sim)};
    end

    function noise_changed(~,~)
        save_current_parameter_table();
        update_noise_controls();
        refresh_parameter_table(true);
        update_compatibility();
        invalidate_result();
    end

    function update_noise_controls()
        kind = char(noise_drop.Value);
        enable(frequency_label, frequency_field, 'on');
        enable(low_label, low_field, 'off');
        enable(high_label, high_field, 'off');
        enable(frequencies_label, frequencies_field, 'off');
        enable(file_label, file_field, 'off');
        file_button.Enable = 'off';
        switch kind
            case {'linear_chirp','bandlimited_noise'}
                enable(low_label, low_field, 'on');
                enable(high_label, high_field, 'on');
            case 'white_noise'
                enable(frequency_label, frequency_field, 'off');
            case 'multisine'
                enable(frequencies_label, frequencies_field, 'on');
            case 'file'
                enable(file_label, file_field, 'on');
                file_button.Enable = 'on';
        end
    end

    function controller_changed(~,~)
        save_current_parameter_table();
        current_controller_id = char(controller_drop.Value);
        refresh_parameter_table(true);
        update_compatibility();
        invalidate_result();
    end

    function design_context_changed(~,~)
        save_current_parameter_table();
        refresh_parameter_table(true);
        update_compatibility();
        invalidate_result();
    end

    function update_compatibility()
        controller = controllers(strcmp({controllers.id}, char(controller_drop.Value)));
        kind = char(noise_drop.Value);
        compatible = ismember(kind, controller.noise_types);
        run_button.Enable = on_off(compatible);
        if auto_tune_box.Value
            run_button.Text = 'Auto tune + held-out evaluation';
        else
            run_button.Text = 'Run selected parameters + held-out evaluation';
        end
        override_count = count_overrides();
        if compatible
            expense = '';
            if controller.expensive
                expense = ' Offline eMOPSO may take minutes.';
            end
            status.Value = {sprintf('%s. Manual overrides: %d.%s', ...
                controller.information, override_count, expense)};
        else
            status.Value = {sprintf('%s is incompatible with %s. Supported: %s.', ...
                controller.label, kind, strjoin(controller.noise_types, ', '))};
        end
    end

    function refresh_parameter_table(preserve)
        controller_id = char(controller_drop.Value);
        if preserve && isfield(manual_store, matlab.lang.makeValidName(controller_id))
            stored = manual_store.(matlab.lang.makeValidName(controller_id));
        else
            stored = table();
        end
        pair = selected_pair();
        design = defaults.design;
        design.profile = char(profile_drop.Value);
        design.auto_tune = false;
        design.manual_params = struct();
        parameter_schema = anc_controller_parameter_schema( ...
            controller_id, pair.fs, current_noise(), design);
        row_count = numel(parameter_schema);
        Override = false(row_count,1);
        Parameter = strings(row_count,1);
        Value = strings(row_count,1);
        Unit = strings(row_count,1);
        Description = strings(row_count,1);
        for row = 1:row_count
            Parameter(row) = string(parameter_schema(row).name);
            Value(row) = format_parameter_value(parameter_schema(row).value);
            Unit(row) = string(parameter_schema(row).unit);
            Description(row) = string(parameter_schema(row).description);
            if ~isempty(stored) && any(strcmp(stored.Parameter, Parameter(row)))
                source_row = find(strcmp(stored.Parameter, Parameter(row)), 1);
                if stored.Override(source_row)
                    Override(row) = true;
                    Value(row) = stored.Value(source_row);
                end
            end
        end
        parameter_table.Data = table(Override, Parameter, Value, Unit, Description);
        parameter_count_label.Text = sprintf('%d editable parameters', row_count);
        parameter_footer.Value = { ...
            sprintf('%s @ %.6g Hz; profile %s.', controller_id, pair.fs, char(profile_drop.Value)); ...
            'Override is applied to every auto-tune candidate. With auto tune off, one sampling-rate-aware default candidate is run.'};
        save_current_parameter_table();
    end

    function save_current_parameter_table()
        if isempty(current_controller_id) || isempty(parameter_table.Data), return; end
        manual_store.(matlab.lang.makeValidName(current_controller_id)) = ...
            parameter_table.Data;
    end

    function parameter_edited(~, event)
        if event.Indices(2) == 3
            data = parameter_table.Data;
            data.Override(event.Indices(1)) = true;
            parameter_table.Data = data;
        end
        save_current_parameter_table();
        update_compatibility();
        invalidate_result();
    end

    function reset_parameters(~,~)
        key = matlab.lang.makeValidName(char(controller_drop.Value));
        if isfield(manual_store, key), manual_store = rmfield(manual_store, key); end
        refresh_parameter_table(false);
        update_compatibility();
        invalidate_result();
    end

    function count = count_overrides()
        data = parameter_table.Data;
        if isempty(data), count = 0; else, count = nnz(data.Override); end
    end

    function manual = collect_manual_parameters()
        manual = struct();
        data = parameter_table.Data;
        for row = 1:height(data)
            if ~data.Override(row), continue; end
            name = char(data.Parameter(row));
            definition = parameter_schema(strcmp({parameter_schema.name}, name));
            manual.(name) = parse_anc_parameter_value(data.Value(row), definition);
        end
    end

    function noise = current_noise()
        noise = defaults.noise;
        noise.type = char(noise_drop.Value);
        noise.duration = duration_field.Value;
        noise.amplitude_rms = amplitude_field.Value;
        noise.frequency_hz = frequency_field.Value;
        noise.low_hz = low_field.Value;
        noise.high_hz = high_field.Value;
        noise.frequencies_hz = parse_numeric_list(frequencies_field.Value);
        noise.file = strtrim(char(file_field.Value));
    end

    function cfg = read_config()
        pair = selected_pair();
        cfg = anc_migration_config();
        cfg.model_root = defaults.model_root;
        cfg.models.primary_file = pair.primary.file;
        cfg.models.secondary_file = pair.secondary.file;
        cfg.noise = current_noise();
        cfg.design.controller_id = char(controller_drop.Value);
        cfg.design.profile = char(profile_drop.Value);
        cfg.design.auto_tune = auto_tune_box.Value;
        cfg.design.manual_params = collect_manual_parameters();
        cfg.output.file = strtrim(char(output_field.Value));
        cfg.output.overwrite = overwrite_box.Value;
    end

    function preview_scenario(~,~)
        try
            cfg = read_config();
            models = import_anc_migration_models(cfg.models, ...
                cfg.models.sample_rate_tolerance);
            signals = prepare_anc_migration_signals(models, cfg.noise, ...
                'evaluation', cfg.models.control_delay_samples);
            plot_preview(signals);
            status.Value = {sprintf('%s: %.3f s, %.4g RMS, %d samples.', ...
                signals.T1.type, signals.T1.Tsim, rms(signals.T1.d), ...
                numel(signals.T1.d))};
        catch exception
            show_error(exception, 'Preview failed');
        end
    end

    function run_migration(src,~)
        src.Enable = 'off';
        cleanup = onCleanup(@() restore_run_button(src)); %#ok<NASGU>
        try
            cfg = read_config();
            cfg.output.plot = false;
            cfg.design.progress_callback = @progress;
            status.Value = {'Starting calibration sweep...'};
            drawnow;
            migration = run_anc_migration(cfg);
            state = fig.UserData;
            state.result = migration;
            fig.UserData = state;
            save_button.Enable = 'on';
            if has_deployment_export
                deployment_button.Enable = on_off(migration.summary.passed);
            end
            plot_anc_migration_result(migration, axes_handles);
            status.Value = {sprintf(['Complete: %s, held-out %.2f dB, demand %.3f, ' ...
                'clips %g, pass=%d.'], migration.calibration.best_name, ...
                migration.summary.suppression_db, migration.summary.demand_max, ...
                migration.summary.saturation_count, migration.summary.passed); ...
                ['Saved: ' migration.output_file]};
        catch exception
            show_error(exception, 'Migration failed');
        end
    end

    function progress(fraction, message)
        status.Value = {sprintf('%3.0f%%: %s', 100*fraction, message)};
        drawnow limitrate;
    end

    function restore_run_button(button)
        controller = controllers(strcmp({controllers.id}, char(controller_drop.Value)));
        button.Enable = on_off(ismember(char(noise_drop.Value), controller.noise_types));
    end

    function save_copy(~,~)
        state = fig.UserData;
        if isempty(state.result), return; end
        [name, folder] = uiputfile({'*.mat','MATLAB result (*.mat)'}, ...
            'Save migration result copy', fullfile(pwd, 'anc_migration_result.mat'));
        if isequal(name, 0), return; end
        migration = state.result; %#ok<NASGU>
        save(fullfile(folder, name), 'migration', '-v7.3');
        status.Value = {['Saved result copy: ' fullfile(folder, name)]};
    end

    function export_deployment(~,~)
        state = fig.UserData;
        if isempty(state.result), return; end
        migration = state.result;
        if ~migration.summary.passed
            show_error(MException('anc_migration_designer_gui:notDeployable', ...
                'Held-out evaluation did not pass; deployment export is disabled.'), ...
                'Controller is not deployable');
            return;
        end
        if isempty(deployment_directory)
            initial_directory = pwd;
        else
            initial_directory = deployment_directory;
        end
        if abs(migration.models.fs-48000)/48000 <= 0.01
            suggested_name = 'Opfilter4.mat';
        else
            suggested_name = sprintf('Opfilter4_fs%.0f.mat', migration.models.fs);
        end
        [name, folder] = uiputfile({'*.mat','dSPACE fixed controller (*.mat)'}, ...
            'Export validated fixed R/S controller', ...
            fullfile(initial_directory, suggested_name));
        if isequal(name, 0), return; end
        output_file = fullfile(folder, name);
        overwrite = false;
        if isfile(output_file)
            if strcmp(fig.Visible, 'on')
                answer = uiconfirm(fig, sprintf('Overwrite %s?', output_file), ...
                    'Confirm deployment overwrite', ...
                    'Options', {'Overwrite','Cancel'}, 'DefaultOption', 2, ...
                    'CancelOption', 2);
                if ~strcmp(answer, 'Overwrite'), return; end
                overwrite = true;
            else
                show_error(MException('anc_migration_designer_gui:outputExists', ...
                    'Deployment file already exists: %s', output_file), ...
                    'Deployment export failed');
                return;
            end
        end
        try
            deployment_exporter(migration, output_file, 'Overwrite', overwrite);
            status.Value = {['Exported validated fixed R/S: ' output_file]};
        catch exception
            show_error(exception, 'Deployment export failed');
        end
    end

    function plot_preview(signals)
        test = signals.T1;
        N = numel(test.d);
        step = max(1, ceil(N/20000));
        idx = 1:step:N;
        cla(axes_handles(1));
        plot(axes_handles(1), test.t(idx), test.source(idx), 'Color', [0.6 0.6 0.6]);
        hold(axes_handles(1), 'on');
        plot(axes_handles(1), test.t(idx), test.d(idx), 'b'); hold(axes_handles(1), 'off');
        grid(axes_handles(1), 'on'); xlabel(axes_handles(1), 'Time (s)');
        legend(axes_handles(1), 'Source', 'Primary-path disturbance');
        title(axes_handles(1), sprintf('%s evaluation preview', test.type));
        [pxx, f] = pwelch(test.d, [], [], [], signals.fs);
        cla(axes_handles(2)); plot(axes_handles(2), f, 10*log10(pxx+eps));
        grid(axes_handles(2), 'on'); xlabel(axes_handles(2), 'Frequency (Hz)');
        ylabel(axes_handles(2), 'PSD (dB/Hz)'); title(axes_handles(2), 'Disturbance PSD');
        cla(axes_handles(3)); cla(axes_handles(4));
    end

    function browse_audio(~,~)
        [name, folder] = uigetfile({'*.wav;*.flac;*.mp3','Audio files'}, ...
            'Select source-noise audio');
        if ~isequal(name, 0), file_field.Value = fullfile(folder, name); end
    end

    function choose_output(~,~)
        [name, folder] = uiputfile({'*.mat','MATLAB result (*.mat)'}, ...
            'Choose automatic migration output', 'anc_migration_result.mat');
        if ~isequal(name, 0), output_field.Value = fullfile(folder, name); end
    end

    function invalidate_result()
        state = fig.UserData;
        state.result = [];
        fig.UserData = state;
        save_button.Enable = 'off';
        if has_deployment_export
            deployment_button.Enable = 'off';
        end
    end

    function show_error(exception, title_text)
        status.Value = {exception.message};
        if strcmp(fig.Visible, 'on')
            uialert(fig, exception.message, title_text, 'Icon', 'error');
        end
    end
end

function controllers = select_controllers(catalog, requested)
if isempty(requested)
    controllers = catalog;
    return;
end
if ischar(requested) || (isstring(requested) && isscalar(requested))
    requested = cellstr(string(requested));
else
    requested = cellstr(string(requested(:)));
end
known = {catalog.id};
unknown = requested(~ismember(requested, known));
if ~isempty(unknown)
    error('anc_migration_designer_gui:unknownControllerFilter', ...
        'Unknown requested controller ids: %s.', strjoin(unknown, ', '));
end
controllers = catalog(ismember(known, requested));
if isempty(controllers)
    error('anc_migration_designer_gui:emptyControllerFilter', ...
        'ControllerIds must select at least one controller.');
end
end

function index = select_default_pair(pairs)
score = zeros(numel(pairs),1);
for k = 1:numel(pairs)
    score(k) = 100*strcmp(pairs(k).kind, 'fir') ...
        - abs(pairs(k).fs-2000)/100 ...
        + 10*contains(pairs(k).primary.file, '20260718');
end
[~, index] = max(score);
end

function values = parse_numeric_list(text_value)
cleaned = regexprep(char(text_value), '[\[\],;]', ' ');
values = sscanf(cleaned, '%f').';
end

function text_value = format_parameter_value(value)
if islogical(value)
    if value, text_value = "true"; else, text_value = "false"; end
elseif isnumeric(value)
    if isempty(value)
        text_value = "[]";
    else
        text_value = string(strtrim(sprintf('%.9g ', value)));
    end
else
    text_value = string(value);
end
end

function enable(label, field, state)
label.Enable = state;
field.Enable = state;
end

function value = on_off(condition)
if condition, value = 'on'; else, value = 'off'; end
end

function label = add_label(parent, text_value, row, column)
label = uilabel(parent, 'Text', text_value, 'HorizontalAlignment', 'right');
place(label, row, column);
end

function place(component, row, column)
component.Layout.Row = row;
component.Layout.Column = column;
end
