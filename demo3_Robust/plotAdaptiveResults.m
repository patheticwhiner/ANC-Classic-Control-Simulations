function plotAdaptiveResults(results, cfg)
% plotAdaptiveResults  并列对比多个自适应控制器的仿真结果
%
% 输入:
%   results  - struct array, 每个元素包含:
%              .name         解法名称 (string)
%              .t            时间向量
%              .y            系统输出
%              .u            控制输入
%              .d            扰动信号
%              .plant_out    G 输出
%              .z_hist       辅佐信号
%              .theta_hist   参数历史 (Nparam × Nsim)
%              .P_tr_hist    P 矩阵迹
%              .phi_norm_hist  回归向量范数
%              .pred_err_hist  预测误差
%              .epsilon_hist   归一化误差
%              .m_sq_hist      m_s²
%              .theta_norm_hist ||θ||
%              .Kz_hist        K(s,θ)z
%              .supp_db        抑制效果 (dB)
%              .rms_dist       扰动 RMS
%              .rms_y          输出 RMS
%              .rms_u          控制 RMS
%              .theta_final    最终 θ
%              .P_final        最终 P
%              .sim_time       仿真耗时
%              .final_idx      稳态索引范围
%   cfg      - struct, 包含:
%              .fs        采样频率
%              .omega     扰动频率 (rad/s)
%              .lambda_val Λ 的 λ 参数 (可选)
%              .Nparam    参数个数 (可选)

nMethods = length(results);
colors = lines(nMethods);
styles = {'-', '--', ':', '-.'};

% ---- 图1: 关键时域对比 (多方法叠加) ----
figure('Name','时域对比','Position',[50 50 1100 750]);

% 1a: 扰动信号 (共用)
subplot(3,2,1);
plot(results(1).t, results(1).d, 'Color',[0.7 0.7 0.7], 'LineWidth',1);
hold on;
legend_entries = {'扰动 d(t)'};
for i = 1:nMethods
    plot(results(i).t, results(i).y, 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth',1.2);
    legend_entries{end+1} = [results(i).name ' y(t)'];
end
hold off; grid on;
legend(legend_entries, 'Location','best');
title('系统输出对比');

% 1b: 控制输入
subplot(3,2,2); hold on;
for i = 1:nMethods
    plot(results(i).t, results(i).u, 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth',1);
end
hold off; grid on;
legend({results.name}, 'Location','best');
title('控制信号 u(t)');

% 1c: 辅佐信号
subplot(3,2,3); hold on;
for i = 1:nMethods
    plot(results(i).t, results(i).z_hist, 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth',1);
end
hold off; grid on;
legend({results.name}, 'Location','best');
title('辅佐信号 z(t)');

% 1d: 参数演化
subplot(3,2,4); hold on;
for i = 1:nMethods
    r = results(i);
    if isfield(r, 'theta_hist') && ~isempty(r.theta_hist)
        % 每隔 4 个参数画一条
        nth = size(r.theta_hist, 1);
        idx = 1:max(1, floor(nth/4)):nth;
        for j = 1:length(idx)
            plot(r.t, r.theta_hist(idx(j),:), 'Color', colors(i,:)*(1-0.3*(j-1)/length(idx)), 'LineWidth',0.8);
        end
    end
end
hold off; grid on;
title('参数演化 \theta(t)');
legend_str = {};
for i = 1:nMethods
    legend_str{end+1} = results(i).name;
end
legend(legend_str, 'Location','best');

% 1e: P 矩阵迹
subplot(3,2,5); hold on;
for i = 1:nMethods
    semilogy(results(i).t, max(results(i).P_tr_hist, 1e-10), 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth',1);
end
hold off; grid on;
legend({results.name}, 'Location','best');
title('P 矩阵迹 trace(P)');

% 1f: 回归向量范数
subplot(3,2,6); hold on;
for i = 1:nMethods
    plot(results(i).t, results(i).phi_norm_hist, 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth',1);
end
hold off; grid on;
legend({results.name}, 'Location','best');
title('回归向量范数 ||\phi(t)||');

% ---- 图2: 误差与收敛对比 (2×2) ----
figure('Name','误差与收敛','Position',[100 100 1000 650]);

% 2a: 预测误差 & 归一化误差
subplot(2,2,1); hold on;
for i = 1:nMethods
    plot(results(i).t, abs(results(i).pred_err_hist), 'Color', colors(i,:), 'LineStyle', '-', 'LineWidth',1);
    plot(results(i).t, abs(results(i).epsilon_hist), 'Color', colors(i,:), 'LineStyle', '--', 'LineWidth',1);
end
hold off; grid on;
title('误差对比 (实线:|预误| 虚线:|ε|)');

% 2b: 归一化因子
subplot(2,2,2); hold on;
for i = 1:nMethods
    plot(results(i).t, results(i).m_sq_hist, 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth',1);
end
hold off; grid on;
legend({results.name}, 'Location','best');
title('归一化因子 m_s^2');

% 2c: 参数范数收敛
subplot(2,2,3); hold on;
for i = 1:nMethods
    plot(results(i).t, results(i).theta_norm_hist, 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth',1.2);
end
hold off; grid on;
legend({results.name}, 'Location','best');
title('参数范数 ||\theta(t)||');

% 2d: K(s,θ)z 与 φ 范数 (取第一个方法)
subplot(2,2,4);
yyaxis left;
plot(results(1).t, results(1).Kz_hist, 'b-', 'LineWidth',1); ylabel('K(s,\theta)z');
yyaxis right;
plot(results(1).t, results(1).phi_norm_hist, 'r-', 'LineWidth',1); ylabel('||\phi(t)||');
title([results(1).name ': K(s,\theta)z 与 \phi 范数']);
grid on;

% ---- 图3: 性能指标总览 ----
figure('Name','性能指标','Position',[150 150 1000 700]);

% 3a: RMS 对比柱状图
subplot(2,3,1);
rms_dist = results(1).rms_dist;
rms_vals_d = rms_dist * ones(1, nMethods);
rms_vals_y = zeros(1, nMethods);
rms_vals_u = zeros(1, nMethods);
for i = 1:nMethods
    rms_vals_y(i) = results(i).rms_y;
    rms_vals_u(i) = results(i).rms_u;
end
bar_data = [rms_vals_d; rms_vals_y; rms_vals_u]';
b = bar(bar_data);
set(gca, 'XTickLabel', {results.name});
ylabel('RMS Value');
legend('扰动RMS', '输出RMS', '控制RMS', 'Location','best');
title('RMS 对比'); grid on;

% 3b: 抑制效果柱状图
subplot(2,3,2);
supp_vals = zeros(1, nMethods);
for i = 1:nMethods
    supp_vals(i) = results(i).supp_db;
end
bar(supp_vals);
set(gca, 'XTickLabel', {results.name});
ylabel('抑制 (dB)');
title('噪声抑制效果'); grid on;

% 3c: 关键自适应指标
subplot(2,3,3);
theta_final_vals = zeros(1, nMethods);
for i = 1:nMethods
    theta_final_vals(i) = norm(results(i).theta_final);
end
bar(theta_final_vals);
set(gca, 'XTickLabel', {results.name});
ylabel('||\theta||');
title('最终参数范数'); grid on;

% 3d: 仿真耗时
subplot(2,3,4);
time_vals = zeros(1, nMethods);
for i = 1:nMethods
    time_vals(i) = results(i).sim_time;
end
bar(time_vals);
set(gca, 'XTickLabel', {results.name});
ylabel('耗时 (s)');
title('仿真耗时'); grid on;

% 3e: 输出频谱 (稳态段，多方法叠加)
subplot(2,3,5);
fs = cfg.fs;
final_idx = results(1).final_idx;
Nfft = min(2048, length(final_idx));
fs_idx = final_idx(end-Nfft+1:end);
[Pdd, f_p] = pwelch(results(1).d(fs_idx).*hann(Nfft)', hann(512), 256, 512, fs);
if isfield(cfg, 'omega')
    f_dist = cfg.omega / (2*pi);
else
    f_dist = [];
end
semilogy(f_p(f_p<=100), Pdd(f_p<=100), 'Color',[0.7 0.7 0.7], 'LineWidth',1.5); hold on;
for i = 1:nMethods
    [Pyy, ~] = pwelch(results(i).y(fs_idx).*hann(Nfft)', hann(512), 256, 512, fs);
    semilogy(f_p(f_p<=100), Pyy(f_p<=100), 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth',1.2);
end
for j = 1:length(f_dist)
    xline(f_dist(j), 'r--');
end
hold off; grid on; xlim([0 50]);
xlabel('Hz'); title('稳态 PSD 频谱');
legend_str = ['扰动', {results.name}];
legend(legend_str, 'Location','best');

% 3f: 参数范数收敛轨迹
subplot(2,3,6); hold on;
for i = 1:nMethods
    plot(results(i).t, results(i).theta_norm_hist, 'Color', colors(i,:), 'LineStyle', styles{min(i,end)}, 'LineWidth',1.5);
end
hold off; grid on;
legend({results.name}, 'Location','best');
xlabel('Time (s)'); ylabel('||\theta(t)||');
title('参数收敛对比');

% ---- 图4: Bode 图 (CC 系列, 仅当有 Lambda_tf 和 K_tf) ----
if isfield(results(1), 'Lambda_tf') && isfield(results(1), 'K_tf')
    figure('Name','Bode 分析','Position',[50 50 1000 700]);
    fs = cfg.fs;
    f_plot = logspace(0, log10(fs/2), 800);
    w_plot = 2*pi*f_plot;

    subplot(2,1,1); hold on;
    for i = 1:nMethods
        if ~isempty(results(i).K_tf)
            [magK, ~] = bode(results(i).K_tf, w_plot);
            semilogx(f_plot, 20*log10(squeeze(magK)), 'Color', colors(i,:), ...
                'LineStyle', styles{min(i,end)}, 'LineWidth',2.2);
        end
    end
    % 画 Lambda 基函数 (仅第一个方法)
    if ~isempty(results(1).Lambda_tf)
        Nparam = length(results(1).Lambda_tf);
        for j = 1:max(1, floor(Nparam/4)):Nparam
            [magL, ~] = bode(results(1).Lambda_tf{j}, w_plot);
            semilogx(f_plot, 20*log10(squeeze(magL)), 'k--', 'LineWidth',0.5);
        end
    end
    hold off; grid on;
    legend([{results.name}], 'Location','best');
    xlabel('频率 (Hz)'); ylabel('幅值 (dB)');
    title('\Lambda_i(s) 与 K(s) 幅频响应对比');

    subplot(2,1,2); hold on;
    for i = 1:nMethods
        if ~isempty(results(i).K_tf)
            [~, phaseK] = bode(results(i).K_tf, w_plot);
            semilogx(f_plot, squeeze(phaseK), 'Color', colors(i,:), ...
                'LineStyle', styles{min(i,end)}, 'LineWidth',2.2);
        end
    end
    hold off; grid on;
    legend({results.name}, 'Location','best');
    xlabel('频率 (Hz)'); ylabel('相位 (deg)');
    title('\Lambda_i(s) 与 K(s) 相频响应对比');
end

fprintf('绘图完成。\n');
end
