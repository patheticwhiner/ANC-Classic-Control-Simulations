function [y_out, u_hist_new, y_hist_new] = IIR_filter(num, den, u_in, u_hist, y_hist)
    % 实现IIR滤波器: H(z) = num(z)/den(z)
    % 差分方程: den(z)Y(z) = num(z)U(z)
    % 即: a0*y[n] + a1*y[n-1] + ... = b0*u[n] + b1*u[n-1] + ...
    
    % 更新输入历史（新样本加入队首，长度与输入一致）
    u_hist_new = [u_in, u_hist(1:end-1)];
    % FIR计算时临时补零或截断
    if length(u_hist_new) < length(num)
        u_hist_pad = [u_hist_new, zeros(1, length(num)-length(u_hist_new))];
    else
        u_hist_pad = u_hist_new(1:length(num));
    end
    if length(y_hist) < length(den)-1
        y_hist_pad = [y_hist, zeros(1, length(den)-1-length(y_hist))];
    else
        y_hist_pad = y_hist(1:length(den)-1);
    end
    % 计算输出：y[n] = (num部分 - den[1:end]部分) / den[0]
    num_part = FIR(num, u_hist_pad);
    den_part = FIR(den(2:end), y_hist_pad);  % 不包括den(1)，因为那是当前输出项
    y_out = (num_part - den_part) / den(1);
    % 更新输出历史（长度与输入一致）
    y_hist_new = [y_out, y_hist(1:end-1)];
end