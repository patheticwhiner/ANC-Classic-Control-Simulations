function clearMlxOutputs()
%CLEARMLXOUTPUTS 清除项目中所有 .mlx 实时脚本的输出
%
% 原理：.mlx 是 ZIP 包，代码输出存储在 matlab/output.xml 中。
%       本函数解包 → 删除 output.xml → 更新 [Content_Types].xml → 重新打包。
%       不需要任何内部 API，兼容所有 R2016a+ 版本。
%
% 用法：
%   >> cd tools; clearMlxOutputs        % 在 tools/ 目录下运行
%   >> addpath('tools'); clearMlxOutputs % 任意目录下运行

    % 定位项目根目录（本文件的上级目录）
    scriptDir = fileparts(mfilename('fullpath'));
    projectRoot = fileparts(scriptDir);

    % 递归查找所有 .mlx 文件
    mlxFiles = dir(fullfile(projectRoot, '**', '*.mlx'));
    if isempty(mlxFiles)
        fprintf('未找到 .mlx 文件。\n');
        return;
    end

    fprintf('找到 %d 个 .mlx 文件，开始清除输出...\n\n', length(mlxFiles));

    nCleared = 0;
    for i = 1:length(mlxFiles)
        filePath = fullfile(mlxFiles(i).folder, mlxFiles(i).name);
        relPath = strrep(filePath, [projectRoot filesep], '');
        fprintf('[%d/%d] %s', i, length(mlxFiles), relPath);

        if clearSingleMlxOutput(filePath)
            fprintf('  ✓ 已清除\n');
            nCleared = nCleared + 1;
        else
            fprintf('  (无输出，跳过)\n');
        end
    end

    fprintf('\n完成！共扫描 %d 个文件，清除 %d 个文件的输出。\n', ...
        length(mlxFiles), nCleared);
end

function wasCleared = clearSingleMlxOutput(mlxPath)
%CLEARSINGLEMLXOUTPUT 清除单个 .mlx 文件的输出
%   wasCleared: true 表示清除了输出, false 表示本来就没有输出

    wasCleared = false;

    % 创建临时目录
    tmpDir = tempname;
    mkdir(tmpDir);
    cleaner = onCleanup(@() rmdir(tmpDir, 's'));

    % 解包 .mlx（即 ZIP）
    try
        unzip(mlxPath, tmpDir);
    catch
        warning('无法打开 .mlx 文件: %s', mlxPath);
        return;
    end

    % 核心：删除 matlab/output.xml
    outputFile = fullfile(tmpDir, 'matlab', 'output.xml');
    if ~isfile(outputFile)
        return;  % 本来就没有输出
    end
    delete(outputFile);

    % 同时删除 output/ 目录（某些旧版本 MATLAB 使用）
    outputDir = fullfile(tmpDir, 'output');
    if exist(outputDir, 'dir')
        rmdir(outputDir, 's');
    end

    % 更新 [Content_Types].xml：移除 matlab/output.xml 的声明
    ctFile = fullfile(tmpDir, '[Content_Types].xml');
    if isfile(ctFile)
        xmlText = fileread(ctFile);
        xmlText = regexprep(xmlText, ...
            '<Override[^<>]*PartName="/matlab/output\.xml"[^<>]*/>\s*', '');
        fid = fopen(ctFile, 'w');
        if fid > 0
            fwrite(fid, xmlText);
            fclose(fid);
        end
    end

    % 收集所有剩余文件
    allFiles = dir(fullfile(tmpDir, '**', '*'));
    allFiles = allFiles(~[allFiles.isdir]);
    relPaths = cell(length(allFiles), 1);
    for j = 1:length(allFiles)
        absPath = fullfile(allFiles(j).folder, allFiles(j).name);
        relPaths{j} = strrep(absPath, [tmpDir filesep], '');
    end

    % 重新打包为 .mlx
    tmpZip = [tempname '.zip'];
    zip(tmpZip, relPaths, tmpDir);
    movefile(tmpZip, mlxPath, 'f');

    wasCleared = true;
end
