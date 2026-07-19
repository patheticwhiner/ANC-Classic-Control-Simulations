function schema = anc_controller_parameter_schema(controller_id, fs, noise, design)
%ANC_CONTROLLER_PARAMETER_SCHEMA Resolve sample-rate-aware GUI defaults.

definitions = anc_controller_parameter_definitions(controller_id);
local_design = design;
local_design.auto_tune = false;
local_design.manual_params = struct();
candidates = build_anc_migration_candidates(controller_id, fs, noise, local_design);
defaults = candidates(1).params;

schema = definitions;
for index = 1:numel(schema)
    name = schema(index).name;
    if isfield(defaults, name)
        schema(index).value = defaults.(name);
    else
        schema(index).value = schema(index).default_value;
    end
end
end
