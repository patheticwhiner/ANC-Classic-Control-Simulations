function result = test_anc_model_migration()
%TEST_ANC_MODEL_MIGRATION Backend and invisible-GUI regression test.

tool_dir = fileparts(mfilename('fullpath'));
addpath(tool_dir);
cfg = anc_migration_config();
catalog = scan_anc_migration_models(cfg.model_root);
assert(numel(catalog) == 8, 'Expected all eight 20260718 model artifacts.');
pairs = pair_anc_migration_models(catalog, cfg.models.sample_rate_tolerance);
assert(numel(pairs) == 4, 'Expected FIR/ARMAX pairs at two sample rates.');
for pair_index = 1:numel(pairs)
    pair_cfg = cfg.models;
    pair_cfg.primary_file = pairs(pair_index).primary.file;
    pair_cfg.secondary_file = pairs(pair_index).secondary.file;
    pair_models = import_anc_migration_models( ...
        pair_cfg, cfg.models.sample_rate_tolerance);
    pair_noise = cfg.noise;
    pair_noise.duration = 0.5;
    pair_signals = prepare_anc_migration_signals( ...
        pair_models, pair_noise, 'evaluation', 1);
    assert(abs(pair_signals.fs-pairs(pair_index).fs) < 1e-9);
    assert(all(isfinite(pair_signals.T1.d)));
end

index = find(strcmp({pairs.kind}, 'fir') & abs([pairs.fs]-2000) < 1, 1);
assert(~isempty(index), 'The 2 kHz FIR pair is missing.');
cfg.models.primary_file = pairs(index).primary.file;
cfg.models.secondary_file = pairs(index).secondary.file;
cfg.noise.duration = 1.2;
cfg.noise.type = 'fixed_sine';
cfg.design.controller_id = 'demo3_fir_fixed';
cfg.design.profile = 'quick';
cfg.design.candidate_limit = 1;
cfg.output.save = false;
cfg.output.plot = false;

models = import_anc_migration_models(cfg.models, cfg.models.sample_rate_tolerance);
assert(models.fs == 2000 && strcmp(models.kind, 'fir'));
signals = prepare_anc_migration_signals(models, cfg.noise, 'evaluation', 1);
assert(abs(rms(signals.T1.d)-cfg.noise.amplitude_rms) < 1e-10);
assert(signals.model_data.orders(4) == 1);

controllers = anc_controller_catalog();
controller_results = nan(numel(controllers),1);
for k = 1:numel(controllers)
    local_noise = cfg.noise;
    if strcmp(controllers(k).id, 'demo5_marino_adaptive')
        local_noise.type = 'linear_chirp';
    elseif strcmp(controllers(k).id, 'demo5_marino_fixed')
        local_noise.type = 'fixed_sine';
    end
    design = cfg.design;
    design.controller_id = controllers(k).id;
    candidates = build_anc_migration_candidates( ...
        controllers(k).id, models.fs, local_noise, design);
    assert(~isempty(candidates), 'No candidates for %s.', controllers(k).id);
    smoke_cfg = cfg;
    smoke_cfg.noise = local_noise;
    smoke_cfg.noise.duration = 0.8;
    smoke_cfg.design = design;
    smoke_cfg.design.profile = 'quick';
    smoke_cfg.design.candidate_limit = 1;
    smoke_cfg.output.save = false;
    smoke_cfg.output.plot = false;
    smoke = run_anc_migration(smoke_cfg);
    controller_results(k) = smoke.evaluation.result.supp_db;
    assert(isfinite(controller_results(k)), ...
        'Nonfinite result for %s.', controllers(k).id);
    clear_controller_caches();
end

% Manual parameters must override every auto-tune candidate, while turning
% auto tuning off must retain exactly one sampling-rate-aware candidate.
manual_design = cfg.design;
manual_design.controller_id = 'imc_fxlms';
manual_design.profile = 'standard';
manual_design.auto_tune = true;
manual_design.manual_params = struct('N_fir', 20, 'mu', 0.004, ...
    'delta', 2e-4);
manual_candidates = build_anc_migration_candidates( ...
    'imc_fxlms', models.fs, cfg.noise, manual_design);
assert(numel(manual_candidates) > 1);
assert(all(arrayfun(@(x) x.params.N_fir == 20, manual_candidates)));
assert(all(arrayfun(@(x) x.params.mu == 0.004, manual_candidates)));
manual_design.auto_tune = false;
manual_candidate = build_anc_migration_candidates( ...
    'imc_fxlms', models.fs, cfg.noise, manual_design);
assert(numel(manual_candidate) == 1);
assert(manual_candidate.params.N_fir == 20 && manual_candidate.params.mu == 0.004);
manual_cfg = cfg;
manual_cfg.noise.duration = 0.8;
manual_cfg.design = manual_design;
manual_cfg.design.candidate_limit = 1;
manual_cfg.output.save = false;
manual_cfg.output.plot = false;
manual_run = run_anc_migration(manual_cfg);
assert(manual_run.calibration.best_params.N_fir == 20);
assert(manual_run.calibration.best_params.mu == 0.004);
clear_controller_caches();

for k = 1:numel(controllers)
    schema = anc_controller_parameter_schema(controllers(k).id, ...
        models.fs, cfg.noise, cfg.design);
    assert(~isempty(schema), 'Missing parameter schema for %s.', controllers(k).id);
end

migration = run_anc_migration(cfg);
assert(strcmp(migration.controller_id, cfg.design.controller_id));
assert(isfield(migration, 'calibration') && isfield(migration, 'evaluation'));
assert(isfinite(migration.evaluation.result.supp_db));

fig = anc_migration_designer_gui('Visible', 'off', 'ModelRoot', cfg.model_root);
cleanup = onCleanup(@() close_if_valid(fig)); %#ok<NASGU>
assert(strcmp(fig.Tag, 'anc_migration_designer_gui'));
pair_control = findall(fig, 'Tag', 'migration_model_pair');
controller_control = findall(fig, 'Tag', 'migration_controller');
noise_control = findall(fig, 'Tag', 'migration_noise_type');
assert(numel(pair_control.Items) == 4);
assert(numel(controller_control.Items) == numel(controllers));
assert(numel(noise_control.Items) == 6);
assert(~isempty(findall(fig, 'Tag', 'migration_run')));
assert(~isempty(findall(fig, 'Tag', 'migration_experiment_panel')));
assert(~isempty(findall(fig, 'Tag', 'migration_parameter_panel')));
parameter_control = findall(fig, 'Tag', 'migration_parameter_table');
assert(height(parameter_control.Data) > 0);
assert(~isempty(findall(fig, 'Tag', 'migration_auto_tune')));
assert(isfile(fullfile(tool_dir, 'docs', ...
    'anc_model_migration_algorithm_guide.tex')));

% Exercise the real Preview callback. Creating the figure alone does not
% catch nested-workspace name collisions during table-to-config parsing.
preview_control = findall(fig, 'Tag', 'migration_preview');
feval(preview_control.ButtonPushedFcn, preview_control, []);
status_control = findall(fig, 'Tag', 'migration_status');
status_text = strjoin(cellstr(string(status_control.Value)), ' ');
assert(~contains(status_text, 'table'), ...
    'Preview failed while reading the controller parameter table.');
assert(~isempty(findall(fig, 'Tag', 'migration_ax1').Children), ...
    'Preview callback did not draw the scenario signal.');

result = struct('catalog_count', numel(catalog), 'pair_count', numel(pairs), ...
    'controller_count', numel(controllers), ...
    'controller_suppression_db', controller_results, ...
    'smoke_suppression_db', migration.evaluation.result.supp_db, ...
    'gui_created', true);
fprintf('test_anc_model_migration passed: %d models, %d pairs, %.2f dB smoke result.\n', ...
    result.catalog_count, result.pair_count, result.smoke_suppression_db);
end

function clear_controller_caches()
clear controller_demo1 controller_demo2 controller_demo3 controller_demo4 ...
    controller_demo5 controller_imc_fxlms;
end

function close_if_valid(fig)
if ~isempty(fig) && isvalid(fig), close(fig); end
end
