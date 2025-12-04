function poly_padded = padPolynomial(poly, target_length, str)
    % 补零函数，将多项式前面补零，使其长度与目标长度相同
    % 输入：
    %   poly - 原始多项式系数向量
    %   target_length - 目标长度
    % 输出：
    %   poly_padded - 补零后的多项式系数向量
    % str:
    %   left - 在左侧补零
    %   right - 在右侧补零
    current_length = length(poly);
    if current_length < target_length
        switch str
            case "left"
                poly_padded = [zeros(1, target_length - current_length), poly];
            case "right"
                poly_padded = [poly, zeros(1, target_length - current_length)];
        end
    else
        poly_padded = poly;
    end
end