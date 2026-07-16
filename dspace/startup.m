%% DSPACE/STARTUP  初始化 dspace/ 子项目路径
%
% 用法（任意工作目录）:
%   run(fullfile('<repo>', 'dspace', 'startup.m'))
%
% 依赖 tests/ 的控制器与信号生成代码（等效性对拍的真值来源），
% 以及 dataset/ 的 DataManager（FIR 路径系数）。

dspaceRoot = fileparts(mfilename('fullpath'));
projectRoot = fileparts(dspaceRoot);

run(fullfile(projectRoot, 'project_init.m'));
run(fullfile(projectRoot, 'tests', 'startup.m'));

addpath(dspaceRoot);
addpath(fullfile(dspaceRoot, 'controllers'));
addpath(fullfile(dspaceRoot, 'params'));
addpath(fullfile(dspaceRoot, 'build'));
addpath(fullfile(dspaceRoot, 'build', 'internal'));
addpath(fullfile(dspaceRoot, 'verify'));

clear dspaceRoot projectRoot;
