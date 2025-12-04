function print_analysis_table(poles, zeros, pole_data, zero_data, Fs)
    % 打印分析结果表格
    
    fprintf('系统阶数: %d (分母)\n', length(poles));
    fprintf('零点个数: %d\n', length(zeros));
    fprintf('极点个数: %d\n\n', length(poles));
    
    % 打印极点表格
    fprintf('==================== 极点分析 ====================\n');
    fprintf('%-20s %-15s %-15s %-20s\n', '零/极点', '频率(Hz)', '阻尼系数ξ', '类型');
    fprintf('%s\n', repmat('-', 1, 70));
    
    for i = 1:length(pole_data)
        root_val = pole_data(i).root;
        freq_hz = pole_data(i).freq_hz;
        xi = pole_data(i).damping;
        type_str = pole_data(i).type;
        
        if isreal(root_val)
            root_str = sprintf('%.4f', root_val);
        else
            root_str = sprintf('%.3f%+.3fi', real(root_val), imag(root_val));
        end
        
        if isinf(xi)
            if xi > 0
                xi_str = '∞ (过阻尼)';
            else
                xi_str = '-∞ (不稳定)';
            end
        else
            xi_str = sprintf('%.4f', xi);
        end
        
        fprintf('%-20s %-15.4f %-15s %-20s\n', root_str, freq_hz, xi_str, type_str);
    end
    
    fprintf('\n');
    
    % 打印零点表格
    fprintf('==================== 零点分析 ====================\n');
    fprintf('%-20s %-15s %-15s %-20s\n', '零/极点', '频率(Hz)', '阻尼系数ξ', '类型');
    fprintf('%s\n', repmat('-', 1, 70));
    
    for i = 1:length(zero_data)
        root_val = zero_data(i).root;
        freq_hz = zero_data(i).freq_hz;
        xi = zero_data(i).damping;
        type_str = zero_data(i).type;
        
        if isreal(root_val)
            root_str = sprintf('%.4f', root_val);
        else
            root_str = sprintf('%.3f%+.3fi', real(root_val), imag(root_val));
        end
        
        if isinf(xi)
            if xi > 0
                xi_str = '∞ (过阻尼)';
            else
                xi_str = '-∞ (不稳定)';
            end
        else
            xi_str = sprintf('%.4f', xi);
        end
        
        fprintf('%-20s %-15.4f %-15s %-20s\n', root_str, freq_hz, xi_str, type_str);
    end
    
    fprintf('\n注释:\n');
    fprintf('- ξ = 0: 无阻尼振荡\n');
    fprintf('- 0 < ξ < 1: 欠阻尼\n');
    fprintf('- ξ = 1: 临界阻尼\n');
    fprintf('- ξ > 1: 过阻尼\n');
    fprintf('- ξ < 0: 不稳定\n');
    fprintf('- 频率为0表示非振荡模式\n');
end