function result = subPolynomials(A, B)
    % SUBPOLYNOMIALS 计算两个多项式的差 A-B
    %
    % 用法:
    %   result = subPolynomials(A, B)
    %
    % 输入:
    %   A - 第一个多项式系数向量
    %   B - 第二个多项式系数向量
    %
    % 输出:
    %   result - 结果多项式系数向量 A-B
    
    % 使两个多项式长度相等
    lenA = length(A);
    lenB = length(B);
    
    if lenA > lenB
        % 如果A更长，用零填充B
        B_padded = [B, zeros(1, lenA - lenB)];
        A_padded = A;
    elseif lenB > lenA
        % 如果B更长，用零填充A
        A_padded = [A, zeros(1, lenB - lenA)];
        B_padded = B;
    else
        % 长度相等，无需填充
        A_padded = A;
        B_padded = B;
    end
    
    % 计算差
    result = A_padded - B_padded;
    
    % 移除尾部的零
    while length(result) > 1 && abs(result(end)) < 1e-10
        result = result(1:end-1);
    end
end