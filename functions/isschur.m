function [varargout] = isschur(den)
    % ISSCHUR 检查多项式是否具有Schur稳定性（所有极点在单位圆内）
    %
    % 用法:
    %   msg = isschur(den)      % 返回稳定性结果并绘制零极点图
    %   [msg] = isschur(den)    % 返回稳定性结果并绘制零极点图
    %   [msg, p] = isschur(den) % 返回稳定性结果和极点，不绘制图像
    %
    % 输入:
    %   den - 分母多项式系数向量
    %
    % 输出:
    %   msg - 逻辑值，1表示稳定，0表示不稳定
    %   p - （可选）多项式根（极点）

    p = roots(den);
    % 验证极点是否在单位圆内
    if max(abs(p)) < 1
        msg = 1;
    else
        msg = 0;
    end
    % 根据输出参数数量决定是否显示信息和绘图
    if nargout < 1
        if msg
            disp('所有极点都在单位圆内，系统稳定。');
        else
            disp('至少有一个极点不在单位圆内，系统不稳定。');
        end
        % 只有一个输出参数，显示信息并绘图
        figure;
        zplane(1, den); % 绘制零极点图
        title('零极点分布图');
        grid on;
    end
    
    % 设置输出参数
    varargout{1} = msg;
    if nargout > 1
        varargout{2} = p;
    end
end