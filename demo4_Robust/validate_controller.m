function results = validate_controller(results_file, varargin)
% VALIDATE_CONTROLLER  验证 ε-MOPSO 优化的 RST 控制器性能
%
%   用法:
%     results = validate_controller()
%     results = validate_controller('output/RST_eMOPSO_results.mat')
%     results = validate_controller('output/RST_eMOPSO_results.mat', 'Plot', true)
%
%   输出结构体 results 包含:
%     .R, .S, .T        - 控制器多项式
%     .stable           - 闭环稳定性
%     .max_pole_mag     - 最大极点模值
%     .Ms               - 灵敏度峰值 (dB)
%     .deltaM           - 模裕度
%     .GM               - 增益裕度 (dB)
%     .PM               - 相位裕度 (deg)
%     .deltaTau         - 延迟裕度 (采样周期)
%     .overshoot        - 阶跃响应超调 (%)
%     .settling_time    - 调节时间 (采样周期)
%     .Syp_mag          - 灵敏度函数幅频响应
%     .omega            - 频率向量

% 确保路径可达
scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end
addpath(fullfile(scriptDir, '..', 'functions'));
addpath(scriptDir);

% 默认参数
if nargin < 1 || isempty(results_file)
    results_file = fullfile(scriptDir, 'output', 'RST_eMOPSO_results.mat');
end

plot_flag = false;
for i = 1:2:length(varargin)
    if strcmpi(varargin{i}, 'Plot')
        plot_flag = varargin{i+1};
    end
end

% 加载数据
if ~exist(results_file, 'file')
    error('结果文件不存在: %s\n请先运行 run_RST_eMOPSO', results_file);
end

data = load(results_file);
sys_info = data.sys_info;
theta_opt = data.theta_opt;

fprintf('========== RST 控制器性能验证 ==========\n\n');

% 调用 postprocess_RST 获取 R, S, T
[R, S, T, Syp] = postprocess_RST(theta_opt, sys_info, 'Plot', false);

results.R = R;
results.S = S;
results.T = T;

% === 稳定性验证 ===
B_delayed = [zeros(1, sys_info.d), sys_info.B];
A_HS = conv(sys_info.A, sys_info.H_S);
B_HR = conv(B_delayed, sys_info.H_R);
P_cl = conv(A_HS, S) + conv(B_HR, R);
cl_roots = roots(P_cl);
max_pole_mag = max(abs(cl_roots));
stable = max_pole_mag < 1;

results.stable = stable;
results.max_pole_mag = max_pole_mag;

fprintf('【稳定性】\n');
fprintf('  闭环极点最大模值: %.4f', max_pole_mag);
if stable
    fprintf(' ✓ 稳定\n');
else
    fprintf(' ✗ 不稳定!\n');
end

% === 频域分析 ===
nFreq = 2048;
omega = linspace(0, pi, nFreq);

% 构建传递函数用于频域计算
G = tf(B_delayed, sys_info.A, sys_info.Ts, 'Variable', 'z^-1');
K = tf(R, S, sys_info.Ts, 'Variable', 'z^-1');
L = series(G, K);
Syp_tf = feedback(1, L);

[Syp_mag, ~] = bode(Syp_tf, omega);
Syp_mag = squeeze(Syp_mag);

results.Syp_mag = Syp_mag;
results.omega = omega;

% 灵敏度峰值 Ms
Ms = max(20*log10(Syp_mag));
results.Ms = Ms;

% 模裕度 ΔM
deltaM = min(abs(1 + squeeze(bode(L, omega))));
results.deltaM = deltaM;

% 增益/相位裕度
[GM, PM, ~, ~] = margin(L);
if isinf(GM), GM_dB = Inf; else, GM_dB = 20*log10(GM); end
results.GM = GM_dB;
results.PM = PM;

% 延迟裕度
if ~isempty(PM) && PM > 0
    [~, idx] = min(abs(squeeze(bode(L, omega)) - 1));
    w_cp = omega(max(1, idx));
    if w_cp > 0
        deltaTau = (PM * pi/180) / w_cp;
    else
        deltaTau = Inf;
    end
else
    deltaTau = 0;
end
results.deltaTau = deltaTau;

fprintf('\n【频域性能】\n');
fprintf('  灵敏度峰值 Ms: %.2f dB\n', Ms);
fprintf('  模裕度 ΔM: %.4f\n', deltaM);
fprintf('  增益裕度 GM: %.2f dB\n', GM_dB);
fprintf('  相位裕度 PM: %.2f°\n', PM);
fprintf('  延迟裕度 Δτ: %.2f Ts\n', deltaTau);

% === 时域分析 ===
N_sim = 200;
r = [zeros(10, 1); ones(N_sim-10, 1)];  % 阶跃（延迟10步）

% 闭环传递函数（从 r 到 y）
num_yr = T * B_delayed;
den_yr = P_cl;
% 归一化
num_yr = num_yr / den_yr(1);
den_yr = den_yr / den_yr(1);

y = filter(num_yr, den_yr, r);

% 超调量
final_val = y(end);
peak_val = max(y(10:end));
if abs(final_val) > 1e-6
    overshoot = (peak_val - abs(final_val)) / abs(final_val) * 100;
else
    overshoot = 0;
end
results.overshoot = overshoot;

% 调节时间（进入 ±5% 范围内）
settling_idx = find(abs(y - final_val) > 0.05 * abs(final_val) * ones(size(y)), 1, 'last');
if isempty(settling_idx)
    settling_time = 0;
else
    settling_time = settling_idx + 1;
    if settling_time > N_sim, settling_time = N_sim; end
end
results.settling_time = settling_time;

fprintf('\n【时域性能】\n');
fprintf('  超调量: %.2f%%\n', overshoot);
fprintf('  调节时间 (±5%%): %d Ts\n', settling_time);

% === 综合评估 ===
fprintf('\n【综合评估】\n');
pass_count = 0;
total_checks = 5;

checks = {
    stable,           '闭环稳定';
    GM_dB > 6,        sprintf('GM=%.1f dB > 6 dB', GM_dB);
    PM > 30,          sprintf('PM=%.1f° > 30°', PM);
    Ms < 6,           sprintf('Ms=%.1f dB < 6 dB', Ms);
    deltaM > 0.5,     sprintf('ΔM=%.3f > 0.5', deltaM);
};

for i = 1:size(checks, 1)
    if checks{i, 1}
        fprintf('  ✓ %s\n', checks{i, 2});
        pass_count = pass_count + 1;
    else
        fprintf('  ✗ %s\n', checks{i, 2});
    end
end

fprintf('\n通过率: %d/%d\n', pass_count, total_checks);

% === 可选绘图 ===
if plot_flag
    figure('Name', 'Controller Validation', 'Position', [100, 100, 1200, 800]);

    % 子图1: 灵敏度函数 Bode 图
    subplot(2,3,1);
    semilogx(omega(2:end), 20*log10(Syp_mag(2:end)), 'b-', 'LineWidth', 1.5);
    hold on;
    yline(6, 'r--', 'Ms=6dB', 'LineWidth', 1);
    hold off;
    xlabel('\omega (rad/sample)'); ylabel('|S_{yp}| (dB)');
    title(sprintf('Sensitivity (Ms=%.1f dB)', Ms));
    grid on; xlim([omega(2), pi]);

    % 子图2: 闭环极零图
    subplot(2,3,2);
    pzmap(Syp_tf);
    title(sprintf('Pole-Zero Map (max|pole|=%.3f)', max_pole_mag));

    % 子图3: Nyquist 图
    subplot(2,3,3);
    [re, im] = nyquist(L, omega);
    re = squeeze(re); im = squeeze(im);
    plot(re, im, 'b-', 'LineWidth', 1); hold on;
    plot(-1, 0, 'rx', 'MarkerSize', 10, 'LineWidth', 2);
    % 模裕度圆
    theta_circle = linspace(0, 2*pi, 100);
    plot(-1 + deltaM*cos(theta_circle), deltaM*sin(theta_circle), 'g--', 'LineWidth', 1);
    hold off;
    xlabel('Re'); ylabel('Im');
    title(sprintf('Nyquist (\\DeltaM=%.3f)', deltaM));
    axis equal; grid on;

    % 子图4: 阶跃响应
    subplot(2,3,4);
    t = (0:N_sim-1) * sys_info.Ts;
    stairs(t, y, 'b-', 'LineWidth', 1.5); hold on;
    stairs(t, r, 'r--', 'LineWidth', 1);
    hold off;
    xlabel('Time (s)'); ylabel('y');
    title(sprintf('Step Response (OS=%.1f%%, Ts=%.1f)', overshoot, settling_time));
    legend('y(k)', 'r(k)', 'Location', 'southeast');
    grid on;

    % 子图5: S+T 验证
    subplot(2,3,5);
    T_tf = feedback(L, 1);
    [T_mag, ~] = bode(T_tf, omega);
    T_mag = squeeze(T_mag);
    ST_sum = Syp_mag + T_mag;
    semilogx(omega(2:end), 20*log10(abs(ST_sum(2:end))), 'b-', 'LineWidth', 1.5);
    xlabel('\omega (rad/sample)'); ylabel('|S+T| (dB)');
    title('S+T ≡ 1 Verification');
    grid on; xlim([omega(2), pi]);
    ylim([-1, 1]);

    % 子图6: 综合指标雷达图
    subplot(2,3,6);
    metrics = [min(GM_dB/20*100, 100), min(PM/90*100, 100), ...
               max(0, 100-Ms/12*100), min(deltaM/1.0*100, 100), ...
               max(0, 100-overshoot/50*100)];
    labels = {'GM', 'PM', 'Ms^{-1}', '\DeltaM', 'OS^{-1}'};
    spider_plot(metrics, labels, '综合性能');
end

fprintf('\n========== 验证完成 ==========\n');
end

% 简单的雷达图函数
function spider_plot(values, labels, title_str)
    n = length(values);
    angles = linspace(0, 2*pi, n+1);
    angles = angles(1:n);

    values = max(0, min(100, values));
    values_plot = [values, values(1)];
    angles_plot = [angles, angles(1)];

    polarplot(angles_plot, values_plot, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 8);
    hold on;
    polarplot(angles_plot, 60*ones(1,n+1), 'g--', 'LineWidth', 1);
    polarplot(angles_plot, 80*ones(1,n+1), 'r--', 'LineWidth', 1);
    hold off;

    rlim([0, 100]);
    title(title_str);

    % 手动添加标签
    for i = 1:n
        text(1.15*values(i)*cos(angles(i)), 1.15*values(i)*sin(angles(i)), ...
            labels{i}, 'FontSize', 8, 'HorizontalAlignment', 'center');
    end
end
