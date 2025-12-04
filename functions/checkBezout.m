function checkBezout(A,B,S1,R1,HS,HR,P)
    % 检查求解结果是否满足Bezout方程
    AHS_Sp = conv(A, conv(HS, S1)); BHR_Rp = conv(B, conv(HR,R1));
    
    % 确保两个向量长度相同后再相加
    max_len = max(length(AHS_Sp), length(BHR_Rp));
    AHS_Sp_padded = [AHS_Sp, zeros(1, max_len - length(AHS_Sp))];
    BHR_Rp_padded = [BHR_Rp, zeros(1, max_len - length(BHR_Rp))];
    left_side = AHS_Sp_padded + BHR_Rp_padded;
    
    % 确保左侧和右侧向量长度相同
    max_len2 = max(length(left_side), length(P));
    left_side_padded = [left_side, zeros(1, max_len2 - length(left_side))];
    P_padded = [P, zeros(1, max_len2 - length(P))];
    
    % 比较左侧和右侧(P)
    error_norm = norm(left_side_padded - P_padded) / norm(P_padded);
    fprintf('  验证结果: 相对误差 = %.10e\n', error_norm);
    if error_norm < 1e-10
        fprintf('  Bezout方程解验证成功!\n');
    else
        fprintf('  警告: Bezout方程解验证误差较大!\n');
    end
end