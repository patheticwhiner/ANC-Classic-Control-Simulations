% STARTUP  添加测试所需全部路径
%
%   用法: run('tests/startup.m')  或  cd tests; startup;
%
%   添加: tests/, project_init.m 提供的公共路径，以及各 demo 目录

rootDir = fileparts(mfilename('fullpath'));
projDir = fullfile(rootDir, '..');

run(fullfile(projDir, 'project_init.m'));
addpath(rootDir);
outputDir = fullfile(rootDir, 'output');
if isfolder(outputDir)
    addpath(outputDir);
end
addpath(fullfile(projDir, 'demo1_RST'));
addpath(fullfile(projDir, 'demo2_LQG'));
addpath(fullfile(projDir, 'demo3_Robust'));
addpath(fullfile(projDir, 'demo4_Robust'));
addpath(fullfile(projDir, 'demo4_Robust', 'utils'));

fprintf('测试路径已配置。\n');
fprintf('  信号: generate_test_signals(''armax_30303022'')\n');
fprintf('  测试: run_demoN_test(signals, ''T1'', ''fixed'')\n');
