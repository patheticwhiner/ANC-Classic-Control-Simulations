function [freq_hz, xi, type_str] = analyze_single_root(root_val, Fs)
    % 分析单个根的特性
    % 输入: root_val - 根值, Fs - 采样频率
    % 输出: freq_hz - 频率(Hz), xi - 阻尼系数, type_str - 类型字符串
    
    if isreal(root_val)
        % 实数根
        if abs(root_val) < 1e-10
            freq_hz = 0;
            xi = Inf; % 积分环节
            type_str = '积分环节';
        elseif abs(root_val) >= 1
            % 在单位圆外，不稳定
            freq_hz = 0; % 实数根没有振荡频率
            xi = -Inf; % 不稳定
            type_str = '实数(不稳定)';
        else
            % 在单位圆内，稳定
            freq_hz = 0; % 实数根没有振荡频率
            xi = Inf; % 过阻尼
            type_str = '实数(稳定)';
        end
    else
        % 复数根
        magnitude = abs(root_val);
        angle_val = angle(root_val);
        
        % 计算自然频率 (Hz)
        % 对于离散系统: wn_digital = angle, wn_continuous = angle * Fs / (2*pi)
        freq_hz = abs(angle_val) * Fs / (2*pi);
        
        % 计算阻尼系数 ξ
        % 对于离散系统: s = ln(z)/T, 其中T = 1/Fs
        % ξ = -real(s) / abs(s)
        if magnitude > 1e-10 % 避免除零错误
            s_equiv = log(root_val) * Fs; % 等效连续域极点
            
            if abs(s_equiv) > 1e-10
                xi = -real(s_equiv) / abs(s_equiv);
            else
                xi = 0;
            end
        else
            xi = 0;
        end
        
        % 稳定性判断
        if magnitude < (1 - 1e-10)
            type_str = '复数(稳定)';
        elseif abs(magnitude - 1) < 1e-10
            type_str = '复数(临界稳定)';
        else
            type_str = '复数(不稳定)';
            xi = -abs(xi); % 不稳定时阻尼系数为负
        end
        
        % 修正边界情况
        if xi > 1
            xi = min(xi, 1); % 限制阻尼系数范围
        end
    end
end