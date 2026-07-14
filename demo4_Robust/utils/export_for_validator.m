function export_for_validator(results_file, output_dir)
% EXPORT_FOR_VALIDATOR  将 ε-MOPSO 优化结果导出为 controller_validator 兼容格式
%
%   用法:
%     export_for_validator('output/RST_eMOPSO_results.mat')
%     export_for_validator('output/RST_eMOPSO_results.mat', 'output')
%
%   输入:
%     results_file - RST_eMOPSO_results.mat 的路径
%     output_dir   - 输出目录 (默认: 'output')
%
%   输出文件 (在 output_dir 下):
%     RST_model.mat      - 被控对象模型 (A, B, Ts)，B 含延迟前导零
%     RST_controller.mat - RST 控制器 (R, S, T, Am, Bm)

if nargin < 2, output_dir = 'output'; end

% 确保路径可达（与 postprocess_RST 一致的路径解析）
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
projectRoot = fileparts(fileparts(scriptDir));
run(fullfile(projectRoot, 'project_init.m'));
addpath(scriptDir);

% 加载优化结果
data = load(results_file);
sys_info = data.sys_info;
theta_opt = data.theta_opt;

fprintf('========== 导出 RST 控制器验证文件 ==========\n');

% 调用 postprocess_RST 获取 R, S, T (不绘图)
[R, S, T, ~] = postprocess_RST(theta_opt, sys_info, 'Plot', false);

% 跟踪模型（默认：纯延迟1步）
Am = 1;
Bm = 1;

% 将延迟合并到 B 中
B_delayed = [zeros(1, sys_info.d), sys_info.B];

% 验证闭环稳定性
A_HS = conv(sys_info.A, sys_info.H_S);
B_HR = conv(B_delayed, sys_info.H_R);
P_cl = conv(A_HS, S) + conv(B_HR, R);
% 对齐长度后相加
max_len = max(length(conv(A_HS, S)), length(conv(B_HR, R)));
P_cl1 = [conv(A_HS, S), zeros(1, max_len - length(conv(A_HS, S)))];
P_cl2 = [conv(B_HR, R), zeros(1, max_len - length(conv(B_HR, R)))];
P_cl = P_cl1 + P_cl2;
cl_roots = roots(P_cl);
max_mag = max(abs(cl_roots));
fprintf('闭环极点最大模值: %.4f', max_mag);
if max_mag < 1
    fprintf(' (稳定)\n');
else
    fprintf(' (不稳定!)\n');
end

% 步骤 2：保存模型文件 (B 已含延迟)
if ~exist(output_dir, 'dir'), mkdir(output_dir); end

A = sys_info.A;
B = B_delayed;
Ts = sys_info.Ts;
save(fullfile(output_dir, 'RST_model.mat'), 'A', 'B', 'Ts');
fprintf('\n模型已保存: %s/RST_model.mat\n', output_dir);
fprintf('  A = [%s]\n', num2str(A, '%.4f '));
fprintf('  B = [%s]  (含 d=%d 步延迟)\n', num2str(B, '%.4f '), sys_info.d);
fprintf('  Ts = %.4f\n', Ts);

% 步骤 3：保存控制器文件
save(fullfile(output_dir, 'RST_controller.mat'), 'R', 'S', 'T', 'Am', 'Bm');
fprintf('控制器已保存: %s/RST_controller.mat\n', output_dir);
fprintf('  R = [%s]\n', num2str(R, '%.4f '));
fprintf('  S = [%s]\n', num2str(S, '%.4f '));
fprintf('  T = %.4f\n', T);

fprintf('\n========== 导出完成 ==========\n');
fprintf('在 controller_validator 中:\n');
fprintf('  1. 选择模型文件: %s/RST_model.mat\n', output_dir);
fprintf('  2. 选择控制器文件: %s/RST_controller.mat\n', output_dir);
end
