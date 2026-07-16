function blk = slk_add_mlfcn(parent, name, sourceFile, paramNames, position)
%SLK_ADD_MLFCN 添加 MATLAB Function 块并注入源码、设定参数作用域
%
%   blk = slk_add_mlfcn(parent, name, sourceFile, {'R0','S0'}, [x y w h])
%
%   - 源码整体来自 dspace/controllers/*.m（生成器不内联代码文本，
%     保证提交的 .m 文件是唯一事实来源）。
%   - paramNames 列出的函数参数改为 Parameter 作用域（非可调 →
%     编译期常量，从模型工作区按同名变量解析，persistent 数组
%     尺寸随之在编译期确定，满足 dSPACE codegen 约束）。

blk = add_block('simulink/User-Defined Functions/MATLAB Function', ...
    [parent '/' name], 'Position', position);

script = fileread(sourceFile);
cfg = get_param(blk, 'MATLABFunctionConfiguration');
cfg.FunctionScript = script;

chart = find(sfroot, '-isa', 'Stateflow.EMChart', 'Path', getfullname(blk));
if isempty(chart)
    error('slk_add_mlfcn:ChartNotFound', ...
        '未找到 %s 的 EMChart —— MATLABFunctionConfiguration 注入后图未解析。', ...
        getfullname(blk));
end

for index = 1:numel(paramNames)
    data = find(chart, '-isa', 'Stateflow.Data', 'Name', paramNames{index});
    if isempty(data)
        error('slk_add_mlfcn:DataNotFound', ...
            '块 %s 中未找到参数符号 %s（源: %s）。', ...
            getfullname(blk), paramNames{index}, sourceFile);
    end
    data.Scope = 'Parameter';
    try
        data.Tunable = false;   % 非可调 → 编译期常量
    catch
        % 某些版本 Parameter 数据默认非可调且无 Tunable 属性——忽略
    end
end
end
