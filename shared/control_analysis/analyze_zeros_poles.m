%% 函数定义部分
function [poles, zeros, pole_data, zero_data] = analyze_zeros_poles(A_coeffs, B_coeffs, Fs)
    % 分析零点和极点的函数
    % 输入: A_coeffs - 分母系数, B_coeffs - 分子系数, Fs - 采样频率
    % 输出: poles, zeros - 根的值, pole_data, zero_data - 分析数据
    
    % 计算零点和极点
    poles = roots(A_coeffs);
    zeros = roots(B_coeffs);
    
    % 分析每个极点 - 使用结构体数组避免类型混合
    pole_data = struct('root', {}, 'freq_hz', {}, 'damping', {}, 'type', {});
    for i = 1:length(poles)
        [freq_hz, damping, type_str] = analyze_single_root(poles(i), Fs);
        pole_data(i).root = poles(i);
        pole_data(i).freq_hz = freq_hz;
        pole_data(i).damping = damping;
        pole_data(i).type = type_str;
    end
    
    % 分析每个零点
    zero_data = struct('root', {}, 'freq_hz', {}, 'damping', {}, 'type', {});
    for i = 1:length(zeros)
        [freq_hz, damping, type_str] = analyze_single_root(zeros(i), Fs);
        zero_data(i).root = zeros(i);
        zero_data(i).freq_hz = freq_hz;
        zero_data(i).damping = damping;
        zero_data(i).type = type_str;
    end
end



