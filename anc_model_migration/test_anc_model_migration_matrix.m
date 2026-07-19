function matrix = test_anc_model_migration_matrix()
%TEST_ANC_MODEL_MIGRATION_MATRIX Run every discovered pair/controller entry.
%
% This is the comprehensive integration test.  It verifies execution and
% evidence capture, not that every method reaches the 20 dB success gate.

tool_dir = fileparts(mfilename('fullpath'));
addpath(tool_dir);
base = anc_migration_config();
pairs = pair_anc_migration_models(scan_anc_migration_models(base.model_root));
controllers = anc_controller_catalog();

pair_id = strings(numel(pairs)*numel(controllers),1);
kind = strings(size(pair_id));
fs = nan(size(pair_id));
controller_id = strings(size(pair_id));
ran = false(size(pair_id));
suppression_db = nan(size(pair_id));
passed = false(size(pair_id));
error_message = strings(size(pair_id));
row = 0;

for pair_index = 1:numel(pairs)
    for controller_index = 1:numel(controllers)
        row = row+1;
        pair_id(row) = pairs(pair_index).id;
        kind(row) = pairs(pair_index).kind;
        fs(row) = pairs(pair_index).fs;
        controller_id(row) = controllers(controller_index).id;
        cfg = base;
        cfg.models.primary_file = pairs(pair_index).primary.file;
        cfg.models.secondary_file = pairs(pair_index).secondary.file;
        cfg.noise.duration = 0.8;
        cfg.noise.type = 'fixed_sine';
        if strcmp(controllers(controller_index).id, 'demo5_marino_adaptive')
            cfg.noise.type = 'linear_chirp';
        end
        cfg.design.controller_id = controllers(controller_index).id;
        cfg.design.profile = 'quick';
        cfg.design.candidate_limit = 1;
        cfg.output.save = false;
        cfg.output.plot = false;
        try
            evalc('migration = run_anc_migration(cfg);');
            ran(row) = true;
            suppression_db(row) = migration.summary.suppression_db;
            passed(row) = migration.summary.passed;
        catch exception
            error_message(row) = string(exception.identifier) + ": " + ...
                string(exception.message);
        end
        clear_controller_caches();
    end
end

matrix = table(pair_id, kind, fs, controller_id, ran, ...
    suppression_db, passed, error_message);
assert(all(matrix.ran), 'At least one pair/controller integration failed.');
fprintf('test_anc_model_migration_matrix passed: %d/%d combinations ran.\n', ...
    nnz(matrix.ran), height(matrix));
end

function clear_controller_caches()
clear controller_demo1 controller_demo2 controller_demo3 controller_demo4 ...
    controller_demo5 controller_imc_fxlms;
end
