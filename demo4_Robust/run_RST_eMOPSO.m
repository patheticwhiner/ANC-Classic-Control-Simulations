% RUN_RST_EMOPSO  Thin wrapper — 调用统一入口 run_eMOPSO('RST_toy')
%
%   本文件保留为向后兼容层，内部直接委托给:
%     run_eMOPSO('RST_toy')
%
%   如需使用实测 ARMAX 模型或标准测试函数，请直接调用:
%     run_eMOPSO('RST_armax')
%     run_eMOPSO('F1')  ...  ('F4')
%
%   参考文献:
%     About_RST_eMOPSO_spec.md — RST 控制器整定数学规范
%     About_MOEA_algorithms.md — ε-MOPSO 算法说明

clear; close all; clc;
% run_eMOPSO('RST_armax');
run_eMOPSO('F3');
