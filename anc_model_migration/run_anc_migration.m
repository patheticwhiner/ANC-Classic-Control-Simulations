function migration = run_anc_migration(cfg)
%RUN_ANC_MIGRATION Tune on calibration data and evaluate on held-out data.

if nargin < 1 || isempty(cfg), cfg = anc_migration_config(); end
validate_config(cfg);
configure_paths(cfg.project_root);

models = import_anc_migration_models(cfg.models, ...
    cfg.models.sample_rate_tolerance);
calibration = prepare_anc_migration_signals(models, cfg.noise, ...
    'calibration', cfg.models.control_delay_samples);
evaluation = prepare_anc_migration_signals(models, cfg.noise, ...
    'evaluation', cfg.models.control_delay_samples);
validate_compatibility(cfg.design.controller_id, calibration.T1.type);

candidates = build_anc_migration_candidates( ...
    cfg.design.controller_id, models.fs, cfg.noise, cfg.design);
limit = min(numel(candidates), floor(cfg.design.candidate_limit));
if isinf(cfg.design.candidate_limit), limit = numel(candidates); end
candidates = candidates(1:limit);
if isempty(candidates)
    error('run_anc_migration:noCandidates', 'No controller candidates were generated.');
end

records = repmat(empty_record(), numel(candidates), 1);
for index = 1:numel(candidates)
    notify_progress(cfg, index-1, numel(candidates), candidates(index).name);
    records(index).name = candidates(index).name;
    records(index).params = candidates(index).params;
    try
        value = run_anc_migration_controller(cfg.design.controller_id, ...
            calibration, candidates(index).params);
        [score, feasible, diagnostics] = score_result(value, cfg.design);
        records(index).completed = true;
        records(index).result = value;
        records(index).score = score;
        records(index).feasible = feasible;
        records(index).diagnostics = diagnostics;
    catch exception
        records(index).error_identifier = exception.identifier;
        records(index).error_message = exception.message;
        records(index).score = -Inf;
    end
end

completed = find([records.completed]);
if isempty(completed)
    messages = strjoin({records.error_message}, newline);
    error('run_anc_migration:allCandidatesFailed', ...
        'All controller candidates failed:%s%s', newline, messages);
end
feasible = completed([records(completed).feasible]);
if isempty(feasible), pool = completed; else, pool = feasible; end
[~, local_best] = max([records(pool).score]);
best_index = pool(local_best);
best = records(best_index);

notify_progress(cfg, numel(candidates), numel(candidates), 'held-out evaluation');
evaluation_result = run_anc_migration_controller( ...
    cfg.design.controller_id, evaluation, best.params);
[evaluation_score, evaluation_feasible, evaluation_diagnostics] = ...
    score_result(evaluation_result, cfg.design);
passed = evaluation_feasible ...
    && evaluation_result.supp_db >= cfg.design.success_suppression_db;

migration = struct();
migration.format_version = 1;
migration.created = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
migration.models = models;
migration.noise = cfg.noise;
migration.controller_id = cfg.design.controller_id;
migration.controller = find_controller(cfg.design.controller_id);
migration.config = cfg;
migration.config.design.progress_callback = [];
migration.calibration = struct('signals', calibration, 'candidates', records, ...
    'best_index', best_index, 'best_name', best.name, ...
    'best_params', best.params, 'had_feasible_candidate', ~isempty(feasible));
migration.evaluation = struct('signals', evaluation, 'result', evaluation_result, ...
    'score', evaluation_score, 'feasible', evaluation_feasible, ...
    'diagnostics', evaluation_diagnostics, 'passed', passed);
migration.summary = make_summary(migration);

if cfg.output.save
    output_file = resolve_output_file(cfg, migration);
    if isfile(output_file) && ~cfg.output.overwrite
        error('run_anc_migration:outputExists', ...
            'Output exists and overwrite is disabled: %s', output_file);
    end
    output_dir = fileparts(output_file);
    if ~isfolder(output_dir), mkdir(output_dir); end
    migration.output_file = output_file;
    save(output_file, 'migration', '-v7.3');
    write_tuning_csv(output_file, records);
else
    migration.output_file = '';
end

if cfg.output.plot
    plot_anc_migration_result(migration);
end
notify_progress(cfg, numel(candidates), numel(candidates), 'complete');
end

function validate_config(cfg)
required = {'project_root','models','noise','design','output'};
for index = 1:numel(required)
    if ~isfield(cfg, required{index})
        error('run_anc_migration:missingSetting', ...
            'cfg.%s is required.', required{index});
    end
end
model_fields = {'primary_file','secondary_file','sample_rate_tolerance','control_delay_samples'};
if ~all(isfield(cfg.models, model_fields))
    error('run_anc_migration:missingModelSetting', ...
        'Model selection and delay settings are incomplete.');
end
validateattributes(cfg.design.candidate_limit, {'numeric'}, ...
    {'scalar','positive'});
end

function configure_paths(project_root)
run(fullfile(project_root, 'project_init.m'));
addpath(fullfile(project_root, 'tests'));
addpath(fullfile(project_root, 'tests', 'internal'));
addpath(fullfile(project_root, 'demo1_RST'));
addpath(fullfile(project_root, 'demo2_LQG'));
addpath(fullfile(project_root, 'demo3_Robust'));
addpath(fullfile(project_root, 'demo4_Robust'));
addpath(fullfile(project_root, 'demo4_Robust', 'utils'));
end

function validate_compatibility(controller_id, noise_type)
catalog = anc_controller_catalog();
match = find(strcmp({catalog.id}, controller_id), 1);
if isempty(match)
    error('run_anc_migration:unknownController', ...
        'Unknown controller id: %s', controller_id);
end
if ~ismember(noise_type, catalog(match).noise_types)
    error('run_anc_migration:incompatibleNoise', ...
        '%s does not support %s. Supported types: %s.', ...
        catalog(match).label, noise_type, ...
        strjoin(catalog(match).noise_types, ', '));
end
end

function record = empty_record()
record = struct('name', '', 'params', struct(), 'completed', false, ...
    'result', struct(), 'score', -Inf, 'feasible', false, ...
    'diagnostics', struct(), 'error_identifier', '', 'error_message', '');
end

function [score, feasible, diagnostics] = score_result(result, design)
demand = Inf;
saturation = Inf;
numeric_failures = 0;
if isfield(result, 'extra')
    if isfield(result.extra, 'u_demand_max'), demand = result.extra.u_demand_max; end
    if isfield(result.extra, 'saturation_count'), saturation = result.extra.saturation_count; end
    if isfield(result.extra, 'numeric_failure_count')
        numeric_failures = result.extra.numeric_failure_count;
    end
end
stable = infer_stability(result);
finite = isfinite(result.supp_db) && isfinite(demand) ...
    && numeric_failures == 0;
feasible = finite && stable && demand <= design.tuning_demand_limit ...
    && saturation == 0;
score = result.supp_db ...
    - 20*max(0, demand-design.tuning_demand_limit) ...
    - 50*double(saturation > 0) ...
    - 1000*double(~finite) ...
    - 100*double(~stable);
diagnostics = struct('suppression_db', result.supp_db, ...
    'demand_max', demand, 'saturation_count', saturation, ...
    'numeric_failure_count', numeric_failures, 'stable', stable, ...
    'finite', finite);
end

function stable = infer_stability(result)
stable = true;
if ~isfield(result, 'extra'), return; end
extra = result.extra;
if isfield(extra, 'design') && isstruct(extra.design)
    if isfield(extra.design, 'closed_loop_radius')
        stable = stable && extra.design.closed_loop_radius < 1;
    end
    if isfield(extra.design, 'controller_radius')
        stable = stable && extra.design.controller_radius < 1;
    end
end
if isfield(extra, 'max_eig_Acl')
    stable = stable && isfinite(extra.max_eig_Acl) && extra.max_eig_Acl < 1.05;
end
end

function controller = find_controller(id)
catalog = anc_controller_catalog();
index = find(strcmp({catalog.id}, id), 1);
controller = catalog(index);
end

function summary = make_summary(migration)
r = migration.evaluation.result;
d = migration.evaluation.diagnostics;
summary = struct('model_id', migration.models.id, ...
    'model_kind', migration.models.kind, 'fs', migration.models.fs, ...
    'controller_id', migration.controller_id, ...
    'candidate', migration.calibration.best_name, ...
    'suppression_db', r.supp_db, 'demand_max', d.demand_max, ...
    'saturation_count', d.saturation_count, ...
    'stable', d.stable, 'passed', migration.evaluation.passed);
end

function output_file = resolve_output_file(cfg, migration)
if ~isempty(cfg.output.file)
    output_file = char(cfg.output.file);
    return;
end
safe_model = matlab.lang.makeValidName(migration.models.id);
safe_controller = matlab.lang.makeValidName(migration.controller_id);
stamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
output_file = fullfile(cfg.project_root, 'anc_model_migration', 'output', ...
    safe_model, safe_controller, ['migration_' stamp '.mat']);
end

function write_tuning_csv(output_file, records)
[folder, stem] = fileparts(output_file);
name = strings(numel(records),1);
completed = false(numel(records),1);
feasible = false(numel(records),1);
score = -inf(numel(records),1);
suppression_db = nan(numel(records),1);
demand_max = nan(numel(records),1);
saturation_count = nan(numel(records),1);
error_message = strings(numel(records),1);
for index = 1:numel(records)
    name(index) = records(index).name;
    completed(index) = records(index).completed;
    feasible(index) = records(index).feasible;
    score(index) = records(index).score;
    error_message(index) = records(index).error_message;
    if records(index).completed
        suppression_db(index) = records(index).diagnostics.suppression_db;
        demand_max(index) = records(index).diagnostics.demand_max;
        saturation_count(index) = records(index).diagnostics.saturation_count;
    end
end
table_value = table(name, completed, feasible, score, suppression_db, ...
    demand_max, saturation_count, error_message);
writetable(table_value, fullfile(folder, [stem '_tuning.csv']));
end

function notify_progress(cfg, index, count, message)
callback = cfg.design.progress_callback;
if isempty(callback), return; end
try
    callback(index/max(1,count), message);
catch
    % A UI callback must never invalidate a deterministic simulation run.
end
end
