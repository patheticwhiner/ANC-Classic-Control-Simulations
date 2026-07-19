function value = parse_anc_parameter_value(text_value, definition)
%PARSE_ANC_PARAMETER_VALUE Parse and validate one GUI/manual parameter.

if isnumeric(text_value) || islogical(text_value)
    value = text_value;
else
    text_value = strtrim(char(string(text_value)));
    switch definition.type
        case {'scalar','integer'}
            value = str2double(text_value);
        case 'optional_scalar'
            if isempty(text_value) || strcmp(text_value, '[]')
                value = [];
            else
                value = str2double(text_value);
            end
        case {'vector','vector2','nonempty_vector'}
            if isempty(text_value) || strcmp(text_value, '[]')
                value = [];
            else
                cleaned = regexprep(text_value, '[\[\],;]', ' ');
                value = sscanf(cleaned, '%f').';
            end
        case 'logical'
            if any(strcmpi(text_value, {'true','on','yes','1'}))
                value = true;
            elseif any(strcmpi(text_value, {'false','off','no','0'}))
                value = false;
            else
                value = [];
            end
        case 'choice'
            value = text_value;
        otherwise
            error('parse_anc_parameter_value:unknownType', ...
                'Unknown parameter type %s.', definition.type);
    end
end

name = definition.name;
switch definition.type
    case 'scalar'
        valid = isnumeric(value) && isscalar(value) && ~isnan(value) ...
            && value >= definition.minimum && value <= definition.maximum;
    case 'optional_scalar'
        valid = isnumeric(value) && (isempty(value) || ...
            (isscalar(value) && ~isnan(value) && value >= definition.minimum ...
            && value <= definition.maximum));
    case 'integer'
        valid = isnumeric(value) && isscalar(value) && isfinite(value) ...
            && value == round(value) && value >= definition.minimum ...
            && value <= definition.maximum;
    case 'vector'
        valid = isnumeric(value) && (isempty(value) || isvector(value)) ...
            && all(isfinite(value)) && all(value >= definition.minimum) ...
            && all(value <= definition.maximum);
    case 'vector2'
        valid = isnumeric(value) && isvector(value) && numel(value) == 2 ...
            && all(isfinite(value)) && all(value >= definition.minimum) ...
            && all(value <= definition.maximum) && value(2) > value(1);
    case 'nonempty_vector'
        valid = isnumeric(value) && isvector(value) && ~isempty(value) ...
            && all(isfinite(value)) && all(value >= definition.minimum) ...
            && all(value <= definition.maximum);
    case 'logical'
        valid = islogical(value) && isscalar(value);
    case 'choice'
        valid = (ischar(value) || (isstring(value) && isscalar(value))) ...
            && ismember(char(value), definition.choices);
    otherwise
        valid = false;
end
if ~valid
    bounds = '';
    if any(strcmp(definition.type, {'scalar','optional_scalar','integer', ...
            'vector','vector2','nonempty_vector'}))
        bounds = sprintf(' in [%g, %g]', definition.minimum, definition.maximum);
    elseif strcmp(definition.type, 'choice')
        bounds = ['; choose ' strjoin(definition.choices, ', ')];
    end
    error('parse_anc_parameter_value:invalidValue', ...
        'Invalid value for %s (%s%s).', name, definition.type, bounds);
end
if strcmp(definition.type, 'choice'), value = char(value); end
end
