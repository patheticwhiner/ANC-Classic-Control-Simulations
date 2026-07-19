function catalog = scan_anc_migration_models(model_root)
%SCAN_ANC_MIGRATION_MODELS Inspect supported FIR and ARMAX MAT artifacts.

if nargin < 1 || isempty(model_root)
    cfg = anc_migration_config();
    model_root = cfg.model_root;
end
model_root = char(model_root);
if ~isfolder(model_root)
    error('scan_anc_migration_models:notFound', ...
        'Model directory does not exist: %s', model_root);
end

files = dir(fullfile(model_root, '**', '*.mat'));
template = empty_entry();
catalog = repmat(template, 0, 1);
for index = 1:numel(files)
    file_path = fullfile(files(index).folder, files(index).name);
    try
        entry = inspect_file(file_path);
    catch exception
        entry = template;
        entry.file = file_path;
        entry.name = files(index).name;
        entry.valid = false;
        entry.message = exception.message;
    end
    if ~strcmp(entry.role, 'unknown')
        catalog(end+1, 1) = entry; %#ok<AGROW>
    end
end
end

function entry = inspect_file(file_path)
payload = load(file_path);
entry = empty_entry();
entry.file = file_path;
[~, entry.name, extension] = fileparts(file_path);
entry.name = [entry.name extension];
entry.role = infer_role(entry.name);
entry.batch = string(fileparts(file_path));

if isfield(payload, 's') && isstruct(payload.s) && isfield(payload.s, 'w')
    b = double(payload.s.w(:).');
    a = 1;
    entry.kind = 'fir';
    entry.fs = read_fs(payload.s);
    entry.orders = numel(b);
    entry.order_text = sprintf('FIR(%d)', numel(b));
    if isfield(payload.s, 'mu'), entry.identification_mu = double(payload.s.mu); end
elseif isfield(payload, 'ARMAXmodel') && isstruct(payload.ARMAXmodel) ...
        && isfield(payload.ARMAXmodel, 'model')
    [b, a] = read_ab(payload.ARMAXmodel.model);
    entry.kind = 'armax';
    entry.fs = read_fs(payload.ARMAXmodel);
    if isfield(payload.ARMAXmodel, 'orders')
        entry.orders = double(payload.ARMAXmodel.orders(:).');
        entry.order_text = sprintf('ARMAX%s', mat2str(entry.orders));
    else
        entry.orders = [numel(a)-1, numel(b)-1, 0, leading_zeros(b)];
        entry.order_text = sprintf('ARMAX%s', mat2str(entry.orders));
    end
    if isfield(payload.ARMAXmodel, 'fit_sim')
        entry.fit_sim = double(payload.ARMAXmodel.fit_sim);
    end
else
    error('scan_anc_migration_models:unsupportedFormat', ...
        'Unsupported MAT model format: %s', file_path);
end

if isempty(entry.fs) || ~isscalar(entry.fs) || ~isfinite(entry.fs) || entry.fs <= 0
    error('scan_anc_migration_models:invalidSampleRate', ...
        'Model has no valid sample rate: %s', file_path);
end
if isempty(b) || isempty(a) || any(~isfinite(b)) || any(~isfinite(a))
    error('scan_anc_migration_models:invalidCoefficients', ...
        'Model contains invalid coefficients: %s', file_path);
end

entry.stable = isempty(roots(a)) || all(abs(roots(a)) < 1);
entry.valid = entry.stable;
entry.message = 'ok';
entry.delay_hint = read_delay_hint(file_path);
if ~isfinite(entry.fit_sim)
    entry.fit_sim = read_fit_hint(file_path);
end
entry.label = sprintf('%s | %s | %.6g Hz | %s', ...
    upper(entry.role), upper(entry.kind), entry.fs, entry.order_text);
end

function entry = empty_entry()
entry = struct( ...
    'file', '', 'name', '', 'role', 'unknown', 'kind', 'unknown', ...
    'batch', "", 'fs', NaN, 'orders', [], 'order_text', '', ...
    'fit_sim', NaN, 'stable', false, 'valid', false, ...
    'delay_hint', NaN, 'identification_mu', NaN, ...
    'label', '', 'message', 'not inspected');
end

function role = infer_role(name)
value = lower(name);
if any(contains(value, {'pripath','primpath','primary','pri-'}))
    role = 'primary';
elseif any(contains(value, {'secpath','secondary','sec-'}))
    role = 'secondary';
else
    role = 'unknown';
end
end

function fs = read_fs(value)
fs = [];
if isstruct(value) && isfield(value, 'fs')
    fs = double(value.fs);
elseif isstruct(value) && isfield(value, 'Fs')
    fs = double(value.Fs);
end
end

function [b, a] = read_ab(model)
try
    a = double(model.A(:).');
    b = double(model.B(:).');
catch
    error('scan_anc_migration_models:missingAB', ...
        'ARMAX model does not expose numeric A/B coefficients.');
end
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

function fit = read_fit_hint(file_path)
fit = NaN;
[folder, stem] = fileparts(file_path);
doc = fullfile(folder, [stem '.md']);
if ~isfile(doc), return; end
text = fileread(doc);
token = regexp(text, '\*\*sim Fit\*\*:\s*([\d.]+)%', 'tokens', 'once');
if ~isempty(token), fit = str2double(token{1}); end
end
