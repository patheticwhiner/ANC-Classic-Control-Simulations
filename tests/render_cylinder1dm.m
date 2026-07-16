function render_cylinder1dm(scopeOrDemo, demo, part)
%RENDER_CYLINDER1DM Render frozen main or broadband figures without retuning.
%
%   render_cylinder1dm()                         Render all main T1/T2 output.
%   render_cylinder1dm('broadband')              Render all independent T3 output.
%   render_cylinder1dm('main', 'demo4', 't1')    Render one main subset.
%   render_cylinder1dm('demo4', 't1')             Backward-compatible main form.

if nargin < 1, scopeOrDemo = ''; end
if nargin < 2, demo = ''; end
if nargin < 3, part = ''; end
[scope, demoFilter, partFilter] = parse_render_arguments(nargin, scopeOrDemo, demo, part);

testDir = fileparts(mfilename('fullpath'));
run(fullfile(testDir, 'startup.m'));
outputRoot = fullfile(testDir, 'output', 'cylinder1dm_2k');
if strcmp(scope, 'broadband')
    stageFile = fullfile(outputRoot, 'broadband_stage', 'broadband_stage_results.mat');
else
    stageFile = fullfile(outputRoot, 'demo1234_stage', 'demo1234_stage_results.mat');
end

if isfile(stageFile)
    payload = load(stageFile, 'stage');
    stage = payload.stage;
elseif strcmp(scope, 'broadband')
    legacyFile = fullfile(outputRoot, 'demo1234_stage', 'demo1234_stage_results.mat');
    if ~isfile(legacyFile)
        error('render_cylinder1dm:MissingStage', ...
            'Run run_cylinder1dm_stage(''broadband'') before rendering T3 figures.');
    end
    payload = load(legacyFile, 'stage');
    stage = payload.stage;
else
    error('render_cylinder1dm:MissingStage', ...
        'Run run_cylinder1dm_stage() before rendering main figures.');
end

stage = filter_stage(stage, strcmp(scope, 'broadband'));
figureDir = fullfile(testDir, 'figures', 'cylinder1dm_2k');
if strcmp(scope, 'broadband')
    generate_demo1234_reports(stage, figureDir, testDir, demoFilter, 't3');
    generate_broadband_report(stage, figureDir, testDir, 'report');
else
    generate_demo1234_reports(stage, figureDir, testDir, demoFilter, partFilter);
end
end

function [scope, demoFilter, partFilter] = parse_render_arguments(argumentCount, first, second, third)
scope = 'main';
demoFilter = '';
partFilter = '';
if argumentCount == 0
    return;
end

first = char(first);
if any(strcmpi(first, {'main', 'broadband'}))
    scope = lower(first);
    if argumentCount >= 2 && ~isempty(second), demoFilter = char(second); end
    if argumentCount >= 3 && ~isempty(third), partFilter = char(third); end
else
    demoFilter = first;
    if argumentCount >= 2 && ~isempty(second), partFilter = char(second); end
end
end

function stage = filter_stage(stage, keepBroadband)
isBroadband = strcmp({stage.selections.test}, 'T3');
keep = isBroadband == keepBroadband;
stage.selections = stage.selections(keep);
stage.evaluation_results = stage.evaluation_results(keep);
stage.tuning = stage.tuning(strcmp(stage.tuning.test, 'T3') == keepBroadband, :);
stage.summary = stage.summary(strcmp(stage.summary.test, 'T3') == keepBroadband, :);
if keepBroadband
    stage.scope = 'broadband';
    if isfield(stage, 'demo2_t2_comparison')
        stage = rmfield(stage, 'demo2_t2_comparison');
    end
else
    stage.scope = 'main';
end
end
