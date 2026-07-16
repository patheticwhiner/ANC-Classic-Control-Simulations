function result = compute_metrics(y_open, y_closed, u, test_struct, meta)
% COMPUTE_METRICS  统一测试指标计算
%
%   对单次测试运行计算所有标准化评估指标。
%
%   用法:
%     result = compute_metrics(y_open, y_closed, u, test_struct, meta)
%
%   输入:
%     y_open      — [N×1] 开环输出（来自 test_signals.T?.y_open）
%     y_closed    — [N×1] 闭环输出（demo 仿真结果）
%     u           — [N×1] 控制信号
%     test_struct — 测试信号子结构 (signals.T1 / .T2 / .T3)
%     meta        — struct，字段:
%                     .demo    — 'demo1'|'demo2'|'demo3'|'demo4'
%                     .variant — 'fixed'|'adaptive'
%                     .test    — 'T1'|'T2'|'T3'
%
%   输出:
%     result — 标准化 struct:
%       .demo, .variant, .test, .fs, .Tsim
%       .supp_db         — 稳态抑制 (dB)，标量
%       .supp_breakdown  — T2 专用：频率分箱抑制向量
%       .t_conv_s        — 收敛时间 (s)，固定控制器 = 0
%       .u_max, .u_rms   — 控制信号统计
%       .y_open_rms      — 开环 RMS (稳态窗口内)
%       .y_closed_rms    — 闭环 RMS (稳态窗口内)
%       .d_rms           — 扰动 RMS (稳态窗口内)
%       .extra           — []，demo 特定字段由调用方追加

    fs   = test_struct.fs;
    Tsim = test_struct.Tsim;

    %% ---- 稳态窗口定义 ----
    % 固定控制器：跳过滤波器初始瞬态 (0.5s)
    % 自适应控制器：留出收敛时间 (3s)
    % 统一取该窗口后 80% 作为稳态段

    switch meta.variant
        case 'fixed'
            t_win_start = 0.5;
        case 'adaptive'
            t_win_start = 3;
        otherwise
            error('Unknown variant: ''%s''. Use ''fixed'' or ''adaptive''.', meta.variant);
    end

    idx_start = min(round(t_win_start * fs) + 1, ...
        max(1, length(y_open) - round(0.5 * fs)));
    idx_total = idx_start:length(y_open);
    n_total   = length(idx_total);
    % 后 80%
    idx_ss    = idx_total(max(1, round(0.2 * n_total)):end);

    y_open_ss   = y_open(idx_ss);
    y_closed_ss = y_closed(idx_ss);
    d_ss        = test_struct.d(idx_ss);  % 扰动信号在同一稳态窗口

    rms_open   = rms(y_open_ss);
    rms_closed = rms(y_closed_ss);

    %% ---- 稳态抑制 ----
    supp_db = 20 * log10(rms_open / max(rms_closed, 1e-12));

    %% ---- T2 频率分箱抑制 ----
    supp_breakdown = [];
    if strcmp(meta.test, 'T2') && isfield(test_struct, 'f_inst')
        f_inst  = test_struct.f_inst;
        f_start = min(test_struct.f_range);
        f_end   = max(test_struct.f_range);
        bin_width = 20;  % Hz

        bin_edges = f_start:bin_width:f_end;
        n_bins = length(bin_edges) - 1;

        supp_breakdown = struct(...
            'bin_centers', zeros(n_bins, 1), ...
            'supp_db',     zeros(n_bins, 1), ...
            'n_samples',   zeros(n_bins, 1) ...
        );

        for i = 1:n_bins
            f_lo = bin_edges(i);
            f_hi = bin_edges(i+1);
            in_bin = (f_inst >= f_lo) & (f_inst < f_hi);
            % 最后一个 bin 包含上边界
            if i == n_bins
                in_bin = (f_inst >= f_lo) & (f_inst <= f_hi);
            end

            if sum(in_bin) < fs * 0.05  % 少于 50ms 数据，跳过
                supp_breakdown.bin_centers(i) = (f_lo + f_hi) / 2;
                supp_breakdown.supp_db(i)     = NaN;
                supp_breakdown.n_samples(i)   = sum(in_bin);
                continue;
            end

            rms_o_bin = rms(y_open(in_bin));
            rms_c_bin = rms(y_closed(in_bin));

            supp_breakdown.bin_centers(i) = (f_lo + f_hi) / 2;
            supp_breakdown.supp_db(i)     = 20 * log10(rms_o_bin / max(rms_c_bin, 1e-12));
            supp_breakdown.n_samples(i)   = sum(in_bin);
        end
    end

    %% ---- 收敛时间 (仅自适应) ----
    t_conv_s = 0;
    if strcmp(meta.variant, 'adaptive')
        % 用滑动窗计算抑制的移动平均，检测 50% 稳态值穿越点
        win_len   = round(0.2 * fs);  % 200ms 滑动窗
        win_step  = round(0.05 * fs); % 50ms 步进
        n_windows = floor((length(y_open) - win_len) / win_step);

        supp_moving = zeros(n_windows, 1);
        t_moving    = zeros(n_windows, 1);

        for i = 1:n_windows
            i0 = (i-1)*win_step + 1;
            i1 = i0 + win_len - 1;
            rms_o = rms(y_open(i0:i1));
            rms_c = rms(y_closed(i0:i1));
            supp_moving(i) = 20 * log10(rms_o / max(rms_c, 1e-12));
            t_moving(i)    = (i0 + win_len/2) / fs;
        end

        % 稳态抑制 (正 = 有抑制) 的 50% 阈值
        target_supp = 0.5 * supp_db;

        % 找到第一个持续超过阈值的窗口
        sustained_count = 0;
        for i = 1:n_windows
            if supp_moving(i) >= target_supp
                sustained_count = sustained_count + 1;
                if sustained_count >= 3  % 连续 3 窗 (~150ms) 确认收敛
                    t_conv_s = t_moving(i - 2);
                    break;
                end
            else
                sustained_count = 0;
            end
        end

        % 如果始终未达到阈值，报告最接近的时间
        if t_conv_s == 0
            [~, i_best] = max(supp_moving);
            if supp_moving(i_best) >= 0
                t_conv_s = t_moving(i_best);
            else
                t_conv_s = Tsim;  % 未收敛
            end
        end
    end

    %% ---- 控制能量 ----
    u_max = max(abs(u));
    u_rms = rms(u);

    %% ---- 组装输出 ----
    result = struct();
    result.demo           = meta.demo;
    result.variant        = meta.variant;
    result.test           = meta.test;
    result.fs             = fs;
    result.Tsim           = Tsim;
    result.supp_db        = supp_db;
    result.supp_breakdown = supp_breakdown;
    result.t_conv_s       = t_conv_s;
    result.u_max          = u_max;
    result.u_rms          = u_rms;
    result.y_open_rms     = rms_open;
    result.y_closed_rms   = rms_closed;
    result.d_rms          = rms(d_ss);
    result.extra          = [];
end
