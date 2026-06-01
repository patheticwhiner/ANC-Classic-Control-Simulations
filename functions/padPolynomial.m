function poly_padded = padPolynomial(poly, target_length, str)
    % 补零函数，将多项式前面补零，使其长度与目标长度相同
    % 输入：
    %   poly - 原始多项式系数向量
    %   target_length - 目标长度
    %   str   - 可选, 'left' (左侧补零, 默认) 或 'right' (右侧补零)
    % 输出：
    %   poly_padded - 补零后的多项式系数向量
    if nargin < 3
        str = 'left';
    end
    current_length = length(poly);
    if current_length < target_length
        if strcmp(str, 'left') || strcmp(str, "left")
            poly_padded = [zeros(1, target_length - current_length), poly];
        elseif strcmp(str, 'right') || strcmp(str, "right")
            poly_padded = [poly, zeros(1, target_length - current_length)];
        else
            poly_padded = [zeros(1, target_length - current_length), poly];
        end
    else
        poly_padded = poly;
    end
end