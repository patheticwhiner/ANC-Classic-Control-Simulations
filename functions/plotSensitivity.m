function [S_sys, T_sys, S_freq, T_freq, w_out] = plotSensitivity(plant, controller, w, opts)
%PLOTSENSITIVITY 计算并绘制灵敏度S和互补灵敏度T。
%
%   [S_sys, T_sys] = plotSensitivity(P, K)
%       使用默认频率范围绘制奇异值图。
%
%   [S_sys, T_sys, S_freq, T_freq, w_out] = plotSensitivity(P, K, w, opts)
%       使用指定的频率向量w (rad/s)和绘图选项opts。
%
%   输入:
%       plant       - 被控对象 (P), LTI模型 (tf/ss/zpk) 或 FRD。
%       controller  - 控制器 (K), LTI模型或FRD。
%       w           - (可选) 频率向量 (rad/s)。
%       opts        - (可选) 结构体，字段:
%                     .plotType: 'singular' (默认) | 'element'
%                     .channels: [i, j] (当 plotType='element' 时)
%
%   输出:
%       S_sys, T_sys - 灵敏度和互补灵敏度的LTI系统对象 (如果可计算)。
%       S_freq, T_freq - 在各频率点的频域响应矩阵。
%       w_out        - 使用的频率向量 (rad/s)。

    % --- FIX: 检查可选参数是否存在 ---
    if nargin < 4
        opts = struct(); % 如果未提供opts，则创建一个空的结构体
    end
    if nargin < 3
        w = []; % 如果未提供w，则设置为空，由parseInputs处理
    end
    % --- END FIX ---

    % 1. 解析输入参数和设置默认值
    [w, opts, w_out] = parseInputs(w, opts);

    % 2. 计算开环增益 L = P * K 的频率响应
    L_freq = computeLoopGain(plant, controller, w);

    % 3. 根据开环增益计算灵敏度函数
    [S_freq, T_freq] = calculateSensitivities(L_freq, w);

    % 4. (可选) 尝试创建S和T的LTI系统对象
    [S_sys, T_sys] = createLtiSystems(plant, controller);

    % 5. 根据选项绘图
    plotResponses(w, S_freq, T_freq, L_freq, opts);
end

%% 辅助函数

function [w_out, opts_out, w_used] = parseInputs(w_in, opts_in)
    % 解析输入并设置默认值
    if nargin < 1 || isempty(w_in)
        w_out = logspace(-2, 3, 500); % 默认频率向量
    else
        w_out = w_in;
    end
    w_used = w_out;

    if nargin < 2 || isempty(opts_in)
        opts_out = struct();
    else
        opts_out = opts_in;
    end

    if ~isfield(opts_out, 'plotType'), opts_out.plotType = 'singular'; end
    if ~isfield(opts_out, 'channels'), opts_out.channels = []; end
end

function L_freq = computeLoopGain(plant, controller, w)
    % 计算开环增益 L = P * K 的频率响应
    try
        P_freq = freqresp(plant, w);
        K_freq = freqresp(controller, w);
    catch ME
        error('plotSensitivity:freqrespFailed', ...
            '无法从输入对象获取频率响应。确保P和K是LTI或FRD模型。\nMATLAB错误: %s', ME.message);
    end

    % 维度检查
    [ny, nu, ~] = size(P_freq);
    if size(K_freq, 1) ~= nu || size(K_freq, 2) ~= ny
        error('plotSensitivity:dimMismatch', ...
            '维度不匹配: P是%dx%d系统，但K的频响尺寸不匹配 (期望 %dx%d)。', ny, nu, nu, ny);
    end

    % 逐点计算矩阵乘法 L = P * K
    L_freq = zeros(ny, ny, length(w));
    for i = 1:length(w)
        L_freq(:,:,i) = P_freq(:,:,i) * K_freq(:,:,i);
    end
end

function [S_freq, T_freq] = calculateSensitivities(L_freq, w)
    % 计算灵敏度 S 和 T 的频率响应
    [ny, ~, nw] = size(L_freq);
    S_freq = zeros(ny, ny, nw);
    T_freq = zeros(ny, ny, nw);
    I_ny = eye(ny);

    for i = 1:nw
        Lj = L_freq(:,:,i);
        M = I_ny + Lj;
        
        % 使用条件数检查矩阵是否病态
        if rcond(M) > 1e-12
            Sj = M \ I_ny; % 更高效的 inv(M)
        else
            Sj = pinv(M); % 使用伪逆作为后备
            warning('plotSensitivity:illConditioned', ...
                '在频率 %.2f rad/s 处矩阵 (I+L) 病态，使用伪逆。', w(i));
        end
        
        S_freq(:,:,i) = Sj;
        T_freq(:,:,i) = Lj * Sj;
    end
end

function [S_sys, T_sys] = createLtiSystems(plant, controller)
    % 尝试创建S和T的LTI系统对象
    S_sys = [];
    T_sys = [];
    try
        if isa(plant, 'tf') || isa(plant, 'ss') || isa(plant, 'zpk')
            if isa(controller, 'tf') || isa(controller, 'ss') || isa(controller, 'zpk')
                L_sys = plant * controller;
                S_sys = feedback(1, L_sys); % S = (I + L)^-1
                T_sys = feedback(L_sys, 1); % T = L * (I + L)^-1
            end
        end
    catch
        % 构造失败则忽略，保持S_sys/T_sys为空
    end
end

function plotResponses(w, S_freq, T_freq, L_freq, opts)
    % 绘制灵敏度函数的图形，使用用户指定的Bode图风格

    % --- 检查系统是否为SISO ---
    [ny, nu, ~] = size(L_freq);
    if ny > 1 || nu > 1
        warning('plotSensitivity:MIMOUnsupported', ...
            '此绘图风格仅为SISO系统设计。对于MIMO系统，将只绘制第一个通道(1,1)。');
    end

    % --- 数据准备 ---
    % 转换频率单位从 rad/s 到 Hz
    f_hz = w / (2 * pi);

    % 提取SISO或第一个通道的频率响应
    L_response = squeeze(L_freq(1,1,:));
    S_response = squeeze(S_freq(1,1,:));
    T_response = squeeze(T_freq(1,1,:));

    % 计算幅值 (dB)
    mag_L_dB = 20 * log10(abs(L_response));
    mag_S_dB = 20 * log10(abs(S_response));
    mag_T_dB = 20 * log10(abs(T_response));

    % 计算相位 (度) 并进行卷绕处理
    phase_L = mod(angle(L_response) * 180/pi + 180, 360) - 180;
    phase_S = mod(angle(S_response) * 180/pi + 180, 360) - 180;
    phase_T = mod(angle(T_response) * 180/pi + 180, 360) - 180;

    % --- 绘图 ---
    figure('Name', '灵敏度函数分析 (L, S, T)', 'Position', [100, 100, 800, 600]);

    % 幅频响应
    subplot(2, 1, 1);
    semilogx(f_hz, mag_L_dB, 'k-', 'LineWidth', 1.5); hold on;
    semilogx(f_hz, mag_S_dB, 'b--', 'LineWidth', 1.5);
    semilogx(f_hz, mag_T_dB, 'r:', 'LineWidth', 1.5);
    hold off;
    grid on;
    title('幅频响应', 'FontSize', 14);
    xlabel('频率 (Hz)', 'FontSize', 12);
    ylabel('幅度 (dB)', 'FontSize', 12);
    legend('开环增益 |L|', '灵敏度 |S|', '互补灵敏度 |T|','location','best');
    if ~isempty(f_hz), xlim([f_hz(1), f_hz(end)]); end

    % 相频响应
    subplot(2, 1, 2);
    semilogx(f_hz, phase_L, 'k-', 'LineWidth', 1.5); hold on;
    semilogx(f_hz, phase_S, 'b--', 'LineWidth', 1.5);
    semilogx(f_hz, phase_T, 'r:', 'LineWidth', 1.5);
    hold off;
    grid on;
    title('相频响应', 'FontSize', 14);
    xlabel('频率 (Hz)', 'FontSize', 12);
    ylabel('相位 (度)', 'FontSize', 12);
    legend('开环增益 L', '灵敏度 S', '互补灵敏度 T','location','best');
    if ~isempty(f_hz), xlim([f_hz(1), f_hz(end)]); end
end