function libName = build_controller_library()
%BUILD_CONTROLLER_LIBRARY 生成控制器库 models/lib_anc_controllers.slx
%
%   每个控制器 = 一个 MATLAB Function 块（源码注入自 dspace/controllers/），
%   离线与 RT 模型都通过库链接实例化，保证换 I/O 时控制器结构逐位相同。
%
%   源文件缺失的条目会被跳过并提示（支持按阶段渐进移植）。

dspaceRoot = fileparts(fileparts(fileparts(mfilename('fullpath'))));
controllersDir = fullfile(dspaceRoot, 'controllers');
modelsDir = fullfile(dspaceRoot, 'models');
if ~isfolder(modelsDir), mkdir(modelsDir); end

libName = 'lib_anc_controllers';
slk_new_model(libName, 'library');

% (demo, variant) 全矩阵 → 去重后的块清单
cases = { ...
    'demo1', 'fixed'; ...
    'demo1', 'adaptive'; ...
    'demo2', 'fixed'; ...
    'demo2', 'adaptive'; ...
    'demo3', 'fixed'; ...
    'demo3', 'adaptive'; ...
    'demo5', 'fixed'};
seen = {};
row = 0;
for index = 1:size(cases, 1)
    spec = controller_registry(cases{index, 1}, cases{index, 2});
    if any(strcmp(seen, spec.block))
        continue;
    end
    seen{end+1} = spec.block; %#ok<AGROW>

    sourceFile = fullfile(controllersDir, spec.source);
    if ~isfile(sourceFile)
        fprintf('  [跳过] %-20s 源文件尚未移植: %s\n', spec.block, spec.source);
        continue;
    end

    row = row + 1;
    position = [40, 40 + 120*(row-1), 220, 130 + 120*(row-1)];
    slk_add_mlfcn(libName, spec.block, sourceFile, spec.params(:, 1)', position);
    fprintf('  [完成] %-20s ← %s\n', spec.block, spec.source);
end

% 离线被控对象步函数也入库，保证所有离线模型共用同一实现
% （demo1/3/4/5 用多项式递推；demo2 用状态空间，与各自脚本引擎逐位一致）
plantEntries = { ...
    'plant_secondary',    'anc_plant_secondary_step.m', {'PLANT_A', 'PLANT_B_STAR'}; ...
    'plant_secondary_ss', 'anc_plant_ss_step.m', {'PLANT_AF', 'PLANT_BF', 'PLANT_CF'}};
for index = 1:size(plantEntries, 1)
    row = row + 1;
    position = [40, 40 + 120*(row-1), 220, 130 + 120*(row-1)];
    slk_add_mlfcn(libName, plantEntries{index, 1}, ...
        fullfile(controllersDir, plantEntries{index, 2}), ...
        plantEntries{index, 3}, position);
    fprintf('  [完成] %-20s ← %s\n', plantEntries{index, 1}, plantEntries{index, 2});
end

libFile = fullfile(modelsDir, [libName '.slx']);
save_system(libName, libFile);
fprintf('库已生成: %s\n', libFile);
end
