% filepath: e:\Code\MATLAB_ANC\YK_Landau\demos\interpolatedControllerBlondel\checkBezoutSolvability.m
function [isSolvable, message] = checkBezoutSolvability(A, B, Hs, Hr, P)
%CHECKBEZOUTSOLVABILITY 检查Bezout方程是否可解
%   检查给定的多项式A, B, Hs, Hr, P是否满足Bezout方程的可解条件
%
% 输入:
%   A  - 系统分母多项式
%   B  - 系统分子多项式
%   Hs - 控制器分母固定部分
%   Hr - 控制器分子固定部分
%   P  - 期望闭环特征多项式
%
% 输出:
%   isSolvable - 布尔值，表示方程是否可解
%   message    - 字符串，包含诊断信息或错误原因

    % 默认为可解
    isSolvable = true;
    message = '方程可解';
    
    % 确保输入为行向量
    if size(A, 1) > 1, A = A'; end
    if size(B, 1) > 1, B = B'; end
    if size(Hs, 1) > 1, Hs = Hs'; end
    if size(Hr, 1) > 1, Hr = Hr'; end
    if size(P, 1) > 1, P = P'; end
    
    % 如果Hs或Hr为空，设置默认值
    if isempty(Hs), Hs = 1; end
    if isempty(Hr), Hr = 1; end
    
    % 计算多项式阶数
    na = length(A) - 1;
    nb = length(B) - 1;
    nhs = length(Hs) - 1;
    nhr = length(Hr) - 1;
    np = length(P) - 1;
    
    % 打印各多项式阶数信息
    fprintf('============= Bezout方程阶数分析 =============\n');
    fprintf('A(z)阶数: %d\n', na);
    fprintf('B(z)阶数: %d\n', nb);
    fprintf('Hs(z)阶数: %d\n', nhs);
    fprintf('Hr(z)阶数: %d\n', nhr);
    fprintf('P(z)阶数: %d\n', np);
    
    % 计算复合多项式阶数
    if (nhs > 0), Ah = conv(A, Hs); else, Ah = A * Hs; end
    nah = length(Ah) - 1;
    
    if (nhr > 0), Bh = conv(B, Hr); else, Bh = B * Hr; end
    nbh = length(Bh) - 1;
    
    fprintf('A*Hs阶数: %d\n', nah);
    fprintf('B*Hr阶数: %d\n', nbh);
    
    % 检查B的前导零（系统延迟）
    if abs(B(1)) < eps
        firstNonZero = find(abs(B) > eps, 1, 'first');
        systemDelay = firstNonZero - 1;
        fprintf('系统具有 %d 步延迟\n', systemDelay);
    else
        systemDelay = 0;
    end

    
    % 打印Bezout可解条件
    fprintf('\n============= Bezout方程可解条件 =============\n');
    fprintf('对于Diophantine方程 A*Hs*S + B*Hr*R = P\n');
    fprintf('可解条件: P的阶数 ≤ (A*Hs的阶数 + B*Hr的阶数 - 1)\n');
    fprintf('当前条件: %d ≤ %d\n', np, nah + nbh - 1);
    
    % 检查期望多项式阶数是否符合要求
    if np > nah + nbh - 1
        isSolvable = false;
        message = sprintf(['期望特征多项式阶数过高。\n'...
                          'P的阶数: %d\n'...
                          '最大允许阶数: %d\n'...
                          '建议: 降低P的阶数或增加Hs/Hr的阶数'], ...
                          np, nah + nbh - 1);
        
        % 增加具体建议
        missingDegree = np - (nah + nbh - 1);
        fprintf('解决方案建议: \n');
        fprintf('1. 降低期望特征多项式P的阶数至少%d阶, 或\n', missingDegree);
        fprintf('2. 增加Hs的阶数至少%d阶, 或\n', missingDegree);
        fprintf('3. 增加Hr的阶数至少%d阶, 或\n', missingDegree);
        fprintf('4. Hs和Hr组合增加总计至少%d阶\n', missingDegree);
        
        if nhs == 0 && nhr == 0
            fprintf('例如: 设置Hs = [1, %s] (增加%d阶)\n', repmat('a, ', 1, missingDegree), missingDegree);
        end
        
        % 打印方程可解性检查结果
        fprintf('============= 方程可解性检查 ============= \n');
        disp(message);
    end
    
%     % 检查系统可控性
%     [Ac, Bc] = tf2ss(fliplr(B),fliplr(A));
%     if rank(ctrb(Ac, Bc)) < size(Ac, 1)
%         warning('系统可能不完全可控，无法任意放置所有极点');
%     end
    
    % 如果可解，计算预期的控制器阶数
    if isSolvable
        nsp = nbh - 1;
        nrp = nah - 1;
        message = sprintf(['方程可解。\n'...
                          '预期控制器分母阶数: %d\n'...
                          '预期控制器分子阶数: %d'], ...
                          nsp, nrp);
                      
        % 打印方程可解性检查结果
        fprintf('\n============= 方程可解性检查 =============\n');
        disp(message);
    end
end