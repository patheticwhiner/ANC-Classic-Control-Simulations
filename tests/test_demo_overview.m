function test_demo_overview()
%TEST_DEMO_OVERVIEW Regression test for frozen overview adaptation panels.

testDir = fileparts(mfilename('fullpath'));
run(fullfile(testDir, 'startup.m'));
test_display_options('silent');

stageFile = fullfile(testDir, 'output', 'cylinder1dm_2k', ...
    'demo1234_stage', 'demo1234_stage_results.mat');
payload = load(stageFile, 'stage');
stage = payload.stage;

outputDir = tempname;
mkdir(outputDir);
cleanup = onCleanup(@() rmdir(outputDir, 's'));

for demo = {'demo1', 'demo2', 'demo3', 'demo4', 'demo5'}
    name = demo{1};
    results = struct( ...
        'T1_fixed', evaluation_result(stage, name, 'T1', 'fixed'), ...
        'T1_adaptive', evaluation_result(stage, name, 'T1', 'adaptive'), ...
        'T2', evaluation_result(stage, name, 'T2', 'adaptive'));
    outputPath = fullfile(outputDir, sprintf('%s_analysis.png', name));
    details = plot_demo_overview(results, outputPath);

    assert(details.adaptation_sample_count > 1, ...
        '%s adaptation panel has no samples.', name);
    assert(details.adaptation_max > details.adaptation_min, ...
        '%s adaptation panel has no visible range.', name);
    assert(isfile(outputPath) && dir(outputPath).bytes > 0, ...
        '%s overview image was not written.', name);
end

fprintf('test_demo_overview: all five overview figures passed.\n');
end

function result = evaluation_result(stage, demo, testName, variant)
index = find(strcmp({stage.selections.demo}, demo) ...
    & strcmp({stage.selections.test}, testName) ...
    & strcmp({stage.selections.variant}, variant), 1);
if isempty(index)
    error('Missing frozen result for %s/%s/%s.', demo, testName, variant);
end
result = stage.evaluation_results{index};
end
