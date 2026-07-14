function poly_trimmed = trimPolynomial(poly, str)
    % 去掉多项式前面的零
    % 输入：
    %   poly - 原始多项式系数向量
    % 输出：
    %   poly_trimmed - 去掉前面零后的多项式系数向量
    % str:
    %   left - 在左侧补零
    %   right - 在右侧补零
    switch str
        case "left"
            % Remove trailing zeros if any
            non_zero_index = find(poly, 1, 'first');
        case "right"
            % Remove trailing zeros if any
            non_zero_index = find(poly, 1, 'last');
    end
    if ~isempty(non_zero_index)
        switch str
            case "left"
                poly_trimmed = poly(non_zero_index:end);
            case "right"
                poly_trimmed = poly(1:non_zero_index);
        end
    else
        poly_trimmed = poly; % 如果多项式全为零，则返回原多项式
    end
end