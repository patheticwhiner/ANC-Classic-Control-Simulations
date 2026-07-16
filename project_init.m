%% PROJECT_INIT  初始化 ANC-Classic-Control-Simulations 的公共路径
%
% 该脚本可从任意当前工作目录运行。各 demo 或工具入口应使用：
%   run(fullfile(fileparts(mfilename('fullpath')), '..', 'project_init.m'))
%
% 公共运行期代码集中在 shared/，模型加载和模型分析集中在 dataset/。
% demo 自己的 utils/ 目录由对应入口显式加入，不在这里全局展开。

projectRoot = fileparts(mfilename('fullpath'));

addpath(projectRoot);
addpath(fullfile(projectRoot, 'dataset'));
addpath(fullfile(projectRoot, 'dataset', 'analysis'));
addpath(genpath(fullfile(projectRoot, 'shared')));

% 系统辨识工具箱（通过 mpm 安装，MATLAB 不会自动加入路径）
identRoot = '/usr/local/MATLAB/R2025b/toolbox/ident';
if exist(identRoot, 'dir') == 7
    addpath(genpath(identRoot));
end
