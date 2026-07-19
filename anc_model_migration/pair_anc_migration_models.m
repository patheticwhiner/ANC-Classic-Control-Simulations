function pairs = pair_anc_migration_models(catalog, relative_tolerance)
%PAIR_ANC_MIGRATION_MODELS Form traceable primary/secondary model pairs.

if nargin < 2 || isempty(relative_tolerance), relative_tolerance = 0.01; end
if isempty(catalog)
    pairs = repmat(empty_pair(), 0, 1);
    return;
end

primary = catalog(strcmp({catalog.role}, 'primary') & [catalog.valid]);
secondary = catalog(strcmp({catalog.role}, 'secondary') & [catalog.valid]);
pairs = repmat(empty_pair(), 0, 1);
for p = 1:numel(primary)
    for s = 1:numel(secondary)
        same_kind = strcmp(primary(p).kind, secondary(s).kind);
        same_batch = strcmp(char(primary(p).batch), char(secondary(s).batch));
        rate_gap = abs(primary(p).fs-secondary(s).fs) / primary(p).fs;
        if ~same_kind || ~same_batch || rate_gap > relative_tolerance
            continue;
        end
        pair = empty_pair();
        pair.primary = primary(p);
        pair.secondary = secondary(s);
        pair.fs = secondary(s).fs;
        pair.kind = secondary(s).kind;
        pair.id = make_pair_id(primary(p), secondary(s));
        pair.label = sprintf('%s @ %.6g Hz | %s + %s', ...
            upper(pair.kind), pair.fs, strip_extension(primary(p).name), ...
            strip_extension(secondary(s).name));
        pair.valid = true;
        pairs(end+1, 1) = pair; %#ok<AGROW>
    end
end

if ~isempty(pairs)
    [~, order] = sortrows(table( ...
        string({pairs.kind}).', [pairs.fs].', string({pairs.id}).'), [1 2 3]);
    pairs = pairs(order);
end
end

function pair = empty_pair()
pair = struct('id', '', 'label', '', 'kind', '', 'fs', NaN, ...
    'primary', struct(), 'secondary', struct(), 'valid', false);
end

function id = make_pair_id(primary, secondary)
[~, p] = fileparts(primary.file);
[~, s] = fileparts(secondary.file);
raw = sprintf('%s__%s', p, s);
id = matlab.lang.makeValidName(raw, 'ReplacementStyle', 'underscore');
end

function value = strip_extension(name)
[~, value] = fileparts(name);
end
