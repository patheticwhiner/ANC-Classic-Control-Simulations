% STARTUP  添加测试所需全部路径
%
%   用法: run('tests/startup.m')  或  cd tests; startup;
%
%   添加: tests/, project_init.m 提供的公共路径，以及各 demo 目录

rootDir = fileparts(mfilename('fullpath'));
projDir = fullfile(rootDir, '..');

run(fullfile(projDir, 'project_init.m'));
addpath(rootDir);
addpath(fullfile(rootDir, 'internal'));
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
fprintf('  可读脚本: demo1_cylinder1dm / demo2_cylinder1dm / demo3_cylinder1dm / demo4_cylinder1dm / demo5_cylinder1dm\n');
fprintf('  完整阶段: run_cylinder1dm_stage() / run_cylinder1dm_stage(''broadband'')\n');
fprintf('  重绘报告: render_cylinder1dm() / render_cylinder1dm(''broadband'')\n');
