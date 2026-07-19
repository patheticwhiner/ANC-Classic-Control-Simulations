function models = import_anc_migration_models(model_files, relative_tolerance)
%IMPORT_ANC_MIGRATION_MODELS Load a matched primary/secondary model pair.

if nargin < 2 || isempty(relative_tolerance), relative_tolerance = 0.01; end
required = {'primary_file','secondary_file'};
if ~isstruct(model_files) || ~all(isfield(model_files, required))
    error('import_anc_migration_models:invalidInput', ...
        'model_files must contain primary_file and secondary_file.');
end

models = struct();
models.primary = import_one(model_files.primary_file, 'primary');
models.secondary = import_one(model_files.secondary_file, 'secondary');
rate_gap = abs(models.primary.fs-models.secondary.fs) / models.secondary.fs;
if rate_gap > relative_tolerance
    error('import_anc_migration_models:sampleRateMismatch', ...
        'Primary %.6g Hz and secondary %.6g Hz exceed the %.3f%% tolerance.', ...
        models.primary.fs, models.secondary.fs, 100*relative_tolerance);
end
models.fs = models.secondary.fs;
models.kind = models.secondary.kind;
models.id = make_model_id(models.primary.file, models.secondary.file);
end

function model = import_one(file_path, role)
file_path = char(file_path);
if ~isfile(file_path)
    error('import_anc_migration_models:notFound', ...
        '%s model not found: %s', role, file_path);
end
payload = load(file_path);
extra = struct();

if isfield(payload, 's') && isstruct(payload.s) && isfield(payload.s, 'w')
    b = double(payload.s.w(:).');
    a = 1;
    fs = double(payload.s.fs);
    kind = 'fir';
    orders = [0, numel(b)-1, 0, 0];
    if isfield(payload.s, 'mu'), extra.identification_mu = double(payload.s.mu); end
elseif isfield(payload, 'ARMAXmodel') && isstruct(payload.ARMAXmodel) ...
        && isfield(payload.ARMAXmodel, 'model')
    packed = payload.ARMAXmodel;
    a = double(packed.model.A(:).');
    b = double(packed.model.B(:).');
    fs = double(packed.fs);
    kind = 'armax';
    if isfield(packed, 'orders')
        orders = double(packed.orders(:).');
    else
        orders = [numel(a)-1, numel(b)-1, 0, leading_zeros(b)];
    end
    if isfield(packed, 'fit_sim'), extra.fit_sim = double(packed.fit_sim); end
else
    error('import_anc_migration_models:unsupportedFormat', ...
        'Unsupported model format: %s', file_path);
end

if a(1) == 0 || any(~isfinite(a)) || any(~isfinite(b))
    error('import_anc_migration_models:invalidCoefficients', ...
        'Invalid coefficients in %s', file_path);
end
b = b/a(1); a = a/a(1);
poles = roots(a);
stable = isempty(poles) || all(abs(poles) < 1);
if ~stable
    error('import_anc_migration_models:unstableModel', ...
        '%s model is not Schur stable: %s', role, file_path);
end

extra.identification_nk = read_delay_hint(file_path);
model = struct('role', role, 'file', file_path, 'kind', kind, ...
    'A', a, 'B', b, 'orders', orders, 'fs', fs, ...
    'stable', stable, 'poles', poles, 'extra', extra);
end

function count = leading_zeros(values)
tolerance = max(1, max(abs(values))) * 1e-12;
first = find(abs(values) > tolerance, 1, 'first');
if isempty(first), count = numel(values); else, count = first - 1; end
end

function nk = read_delay_hint(file_path)
nk = NaN;
[folder, stem] = fileparts(file_path);
doc = fullfile(folder, [stem '.md']);
if ~isfile(doc), return; end
text = fileread(doc);
token = regexp(text, '\*\*nk\*\*:\s*(\d+)', 'tokens', 'once');
if ~isempty(token), nk = str2double(token{1}); end
end

function id = make_model_id(primary_file, secondary_file)
[~, primary] = fileparts(primary_file);
[~, secondary] = fileparts(secondary_file);
id = matlab.lang.makeValidName([primary '__' secondary], ...
    'ReplacementStyle', 'underscore');
end
