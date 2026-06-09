function plotAVCResults(results, cfg)
% plotAVCResults  并列对比 AVC 控制器的仿真结果
%
% 输入:
%   results  - struct array, 每个元素包含:
%              .name         解法名称
%              .t            时间向量
%              .y            系统输出
%              .u            控制输入
%              .d_total      扰动信号
%              .supp_db       抑制 (dB)
%              .rms_y         输出 RMS
%              .rms_u         控制 RMS
%              .gamma         H∞ 范数 (如适用, NaN 表示不适用)
%              .theta_norm    参数范数
%              .sim_time      耗时
%              .theta_hist    参数历史 (可选, N×Nsim)
%              .z_hist        辅佐信号 (可选)
%   cfg      - struct, 包含:
%              .fs        采样频率
%              .t_on      控制器开启时间
%              .t_disturb 扰动变化时间 (可选)
%              .N_fir     FIR 阶数 (可选)

nMethods = length(results);
colors = lines(nMethods);
styles = {'-', '--', ':', '-.'};

fs = cfg.fs;
t_on = cfg.t_on;
t_all_max = max(arrayfun(@(r) r.t(end), results));

% ---- 图1: 时域输出与控制 ----
figure('Name','AVC 时域对比','Position',[50 50 1100 700]);

% 输出对比
subplot(2,2,1);
plot(results(1).t, results(1).d_total, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2); hold on;
legend_entries = {'扰动 d'};
for i = 1:nMethods
    plot(results(i).t, results(i).y, 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth', 1.2);
    legend_entries{end+1} = [results(i).name ' y'];
end
xline(t_on, 'g--', 'LineWidth', 1.5);
if isfield(cfg, 't_disturb') && ~isempty(cfg.t_disturb)
    xline(cfg.t_disturb, 'r--', 'LineWidth', 1.5);
end
hold off; grid on;
legend(legend_entries, 'Location', 'best');
title('系统输出对比'); xlim([0 t_all_max]);

% 控制输入
subplot(2,2,2); hold on;
for i = 1:nMethods
    plot(results(i).t, results(i).u, 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth', 1);
end
xline(t_on, 'g--', 'LineWidth', 1.5);
if isfield(cfg, 't_disturb') && ~isempty(cfg.t_disturb)
    xline(cfg.t_disturb, 'r--', 'LineWidth', 1.5);
end
hold off; grid on;
legend({results.name}, 'Location', 'best');
title('控制输入 u');

% RMS 滑动窗
subplot(2,2,3); hold on;
win = round(fs * 0.5);
rms_d = sqrt(movmean(results(1).d_total.^2, win));
plot(results(1).t, rms_d, 'Color', [0.7 0.7 0.7], 'LineWidth', 1.2);
legend_entries = {'扰动 RMS'};
for i = 1:nMethods
    rms_y = sqrt(movmean(results(i).y.^2, win));
    plot(results(i).t, rms_y, 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth', 1.2);
    legend_entries{end+1} = [results(i).name ' RMS'];
end
xline(t_on, 'g--', 'LineWidth', 1.5);
if isfield(cfg, 't_disturb') && ~isempty(cfg.t_disturb)
    xline(cfg.t_disturb, 'r--', 'LineWidth', 1.5);
end
hold off; grid on;
legend(legend_entries, 'Location', 'best');
title('滑动 RMS');

% 参数演化（如果有）
subplot(2,2,4); hold on;
has_theta = false;
for i = 1:nMethods
    r = results(i);
    if isfield(r, 'theta_hist') && ~isempty(r.theta_hist)
        has_theta = true;
        nth = size(r.theta_hist, 1);
        step = max(1, floor(nth / 5));
        idx = 1:step:nth;
        for j = 1:min(3, length(idx))
            plot(r.t, r.theta_hist(idx(j),:), 'Color', colors(i,:)*(1-0.3*(j-1)/3), 'LineWidth', 0.8);
        end
    end
end
if ~has_theta
    % 用参数范数代替
    for i = 1:nMethods
        if isfield(results(i), 'theta_norm') && ~isempty(results(i).theta_norm)
            theta_norm = results(i).theta_norm;
            if isscalar(theta_norm)
                yline(theta_norm, 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth', 1.5);
            end
        end
    end
end
xline(t_on, 'g--', 'LineWidth', 1.5);
hold off; grid on;
legend({results.name}, 'Location', 'best');
title('参数演化');

% ---- 图2: 性能对比柱状图 + 频谱 ----
figure('Name','AVC 性能对比','Position',[100 100 1100 700]);

% 抑制效果
subplot(2,3,1);
supp_vals = zeros(1, nMethods);
for i = 1:nMethods
    supp_vals(i) = results(i).supp_db;
end
bar(supp_vals);
set(gca, 'XTickLabel', {results.name});
ylabel('抑制 (dB)');
title('噪声抑制效果'); grid on;

% 输出 RMS
subplot(2,3,2);
rms_y_vals = zeros(1, nMethods);
for i = 1:nMethods
    rms_y_vals(i) = results(i).rms_y;
end
bar(rms_y_vals);
set(gca, 'XTickLabel', {results.name});
ylabel('输出 RMS');
title('稳态输出 RMS'); grid on;

% 控制 RMS
subplot(2,3,3);
rms_u_vals = zeros(1, nMethods);
for i = 1:nMethods
    rms_u_vals(i) = results(i).rms_u;
end
bar(rms_u_vals);
set(gca, 'XTickLabel', {results.name});
ylabel('控制 RMS');
title('稳态控制 RMS'); grid on;

% H∞ 范数 (如适用)
subplot(2,3,4);
gamma_vals = zeros(1, nMethods);
has_gamma = false;
for i = 1:nMethods
    if isfield(results(i), 'gamma') && ~isempty(results(i).gamma) && ~isnan(results(i).gamma)
        gamma_vals(i) = results(i).gamma;
        has_gamma = true;
    end
end
if has_gamma
    bar(gamma_vals);
    set(gca, 'XTickLabel', {results.name});
    ylabel('\gamma');
    title('H∞ 范数'); grid on;
else
    text(0.5, 0.5, '不适用', 'HorizontalAlignment', 'center', 'FontSize', 12);
    title('H∞ 范数');
end

% 频谱
subplot(2,3,5);
Nsim = length(results(1).t);
final_idx = max(1, round(0.7*Nsim)):Nsim;
Nfft = min(2048, length(final_idx));
fs_idx = final_idx(end-Nfft+1:end);
[Pdd, f_p] = pwelch(results(1).d_total(fs_idx).*hann(Nfft)', hann(512), 256, 512, fs);
semilogy(f_p(f_p<=50), Pdd(f_p<=50), 'Color', [0.7 0.7 0.7], 'LineWidth', 1.5); hold on;
for i = 1:nMethods
    [Pyy, ~] = pwelch(results(i).y(fs_idx).*hann(Nfft)', hann(512), 256, 512, fs);
    semilogy(f_p(f_p<=50), Pyy(f_p<=50), 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth', 1.2);
end
hold off; grid on; xlim([0 50]);
xlabel('Hz'); title('稳态 PSD 频谱');
legend_str = ['扰动', {results.name}];
legend(legend_str, 'Location', 'best');

% 耗时
subplot(2,3,6);
time_vals = zeros(1, nMethods);
for i = 1:nMethods
    time_vals(i) = results(i).sim_time;
end
bar(time_vals);
set(gca, 'XTickLabel', {results.name});
ylabel('耗时 (s)');
title('仿真耗时'); grid on;

fprintf('AVC 绘图完成。\n');
end
