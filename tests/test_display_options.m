function options = test_display_options(mode)
%TEST_DISPLAY_OPTIONS Global display policy for tests/ figures.
%
%   test_display_options()             return current policy
%   test_display_options('silent')     save figures without showing them
%   test_display_options('show')       show and save figures, keep windows open
%   test_display_options('reset')      restore default (show) mode

appKey = 'ANC_CLASSIC_TEST_DISPLAY_OPTIONS';
if isappdata(0, appKey)
    options = getappdata(0, appKey);
else
    options = make_options('show');
    setappdata(0, appKey, options);
end

if nargin == 0 || isempty(mode)
    return;
end

mode = lower(char(mode));
switch mode
    case {'silent', 'save', 'save_only'}
        options = make_options('silent');
    case {'show', 'visible', 'interactive'}
        options = make_options('show');
    case 'reset'
        options = make_options('show');
    otherwise
        error('test_display_options:UnknownMode', ...
            'Use ''silent'', ''show'', or ''reset''; got ''%s''.', mode);
end
setappdata(0, appKey, options);
fprintf('Test figure display mode: %s\n', options.mode);
end

function options = make_options(mode)
options = struct( ...
    'mode', mode, ...
    'figure_visibility', ternary(strcmp(mode, 'show'), 'on', 'off'), ...
    'close_after_save', ~strcmp(mode, 'show'));
end

function value = ternary(condition, trueValue, falseValue)
if condition
    value = trueValue;
else
    value = falseValue;
end
end
