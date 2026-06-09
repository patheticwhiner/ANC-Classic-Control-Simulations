function plotComparisonTable(results, columns, title_str)
% plotComparisonTable  打印并列对比性能表 + 柱状图
%
% 输入:
%   results   - struct array, 每个元素代表一种解法
%   columns   - cell array, 每个元素 {显示名, 字段名, 格式, 可选: group_name}
%               e.g. {'抑制(dB)', 'supp_db', '%.1f'}
%   title_str - 表格标题

nMethods = length(results);

fprintf('\n========== %s ==========\n', title_str);

% 打印表头
fprintf('%-15s', '解法');
for c = 1:length(columns)
    fprintf(' | %-12s', columns{c}{1});
end
fprintf('\n');
fprintf('%s', repmat('-', 1, 15 + 15*length(columns)));
fprintf('\n');

% 打印数据行
for i = 1:nMethods
    fprintf('%-15s', results(i).name);
    for c = 1:length(columns)
        field = columns{c}{2};
        fmt = columns{c}{3};
        if isfield(results(i), field)
            val = results(i).(field);
            if isempty(val) || (isnumeric(val) && any(isnan(val)))
                fprintf(' | %-12s', 'N/A');
            else
                fprintf([' | ' fmt], val);
            end
        else
            fprintf(' | %-12s', 'N/A');
        end
    end
    fprintf('\n');
end
fprintf('\n');

% 为数值列生成柱状图 (分组)
num_cols = {};
for c = 1:length(columns)
    field = columns{c}{2};
    % 检查所有方法是否都有该数值字段
    all_numeric = true;
    vals = zeros(1, nMethods);
    for i = 1:nMethods
        if isfield(results(i), field) && ~isempty(results(i).(field)) ...
                && isnumeric(results(i).(field)) && isscalar(results(i).(field)) ...
                && ~isnan(results(i).(field))
            vals(i) = results(i).(field);
        else
            all_numeric = false;
            break;
        end
    end
    if all_numeric && ~all(vals == vals(1))  % 跳过全部相同的列
        num_cols{end+1} = struct('name', columns{c}{1}, 'vals', vals);
    end
end

if ~isempty(num_cols)
    figure('Name', title_str, 'Position', [100 100 min(400+150*length(num_cols), 1200) 500]);
    nCols = length(num_cols);
    for c = 1:nCols
        subplot(1, nCols, c);
        bar(num_cols{c}.vals);
        set(gca, 'XTickLabel', {results.name});
        ylabel(num_cols{c}.name);
        title(num_cols{c}.name);
        grid on;
    end
    sgtitle(title_str);
end

fprintf('对比表打印完成。\n');
end
