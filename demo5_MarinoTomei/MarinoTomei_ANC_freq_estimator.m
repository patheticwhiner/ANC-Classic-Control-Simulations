function results = MarinoTomei_ANC_freq_estimator(run_mode)
%% Marino & Tomei (2016) — ANC 窄带噪声自适应控制
% 本脚本从 ANC 问题出发: 主噪声 d_p 在误差麦克风处与次级路径输出叠加,
% 控制器通过扬声器输入 u 生成反噪声, 目标是压低误差麦克风信号 e:
%
%   e(n) = d_p(n) + S(z)u(n)
%
% Marino-Tomei 自适应内模用于在线估计窄带噪声频率, 并生成补偿信号 u。
% 本脚本提供以下 ANC 对照运行模式:
%   run_mode = 'minimal'          : 最小尺度修复版 (显式 Euler 内模)
%   run_mode = 'discrete'         : 精确离散内模版 (含 300 Hz 未建模噪声)
%   run_mode = 'all'              : 运行诊断版 + 真正 ANC 候选版 (默认)
%   run_mode = 'discrete-modeled' : 理想校验, 仅用于调试算法正确性

close all; clc;

%% ====== 用户设置 ======
if nargin < 1
    run_mode = 'all';
end
run_mode = char(run_mode);

%% ====== ANC 系统与噪声参数 ======
fs   = 4000;             % 采样频率 [Hz]
Ts   = 1 / fs;           % 采样时间 [s]
Tsim = 40;               % 仿真时长 [s]
N    = Tsim * fs;        % 仿真步数
Nlog = 800;              % 绘图日志点数
control_start = 10;      % 前10秒关闭ANC, 作为控制前基线

% 次级路径: 从控制扬声器输入 u 到误差麦克风的传递函数
% S(z) = z^-2 + 1.5*z^-3 - z^-4
S_num = [0, 0, 1, 1.5, -1];
S_len = numel(S_num);

% 主噪声中的建模窄带分量与初始频率估计
freqs     = [50; 120];   % 建模窄带噪声频率 [Hz]
freq_init = [45; 110];   % 初始频率估计 [Hz]

% 主噪声中的未建模窄带分量
unmodeled_freq = 300;    % 未建模噪声频率 [Hz]
unmodeled_amp  = 0.3;    % 未建模噪声幅值

% 数字频率参数: theta = (omega*Ts)^2 = (2*pi*f/fs)^2
theta_true = (2*pi*freqs/fs).^2;
theta_min  = (2*pi*20/fs)^2;
theta_max  = (2*pi*500/fs)^2;

%% ====== sign{Re[S(e^{jωTs})]} 预计算 ======
sgn0 = sign(sum(S_num));
sgn  = zeros(numel(freqs), 1);

for i = 1:numel(freqs)
    omega = 2*pi*freqs(i);
    S_w = S_num * exp(-1j * omega*Ts * (0:S_len-1)');
    sgn(i) = sign(real(S_w));
end

fprintf('=== Marino-Tomei ANC narrowband controller ===\n');
fprintf('error mic: e(n) = d_p(n) + S(z)u(n)\n');
fprintf('secondary path S(z) = z^-2 + 1.5 z^-3 - z^-4\n');
fprintf('sign[S(1)]       = %+d\n', sgn0);
fprintf('sign[Re S(50Hz)] = %+d\n', sgn(1));
fprintf('sign[Re S(120Hz)]= %+d\n', sgn(2));
fprintf('theta true       = [%.6f %.6f]  (digital rad/sample)^2\n', theta_true);
fprintf('ANC control starts at t = %.1f s\n', control_start);

%% ====== 构造实验参数 ======
cfg = struct();
cfg.fs             = fs;
cfg.Ts             = Ts;
cfg.Tsim           = Tsim;
cfg.N              = N;
cfg.Nlog           = Nlog;
cfg.control_start  = control_start;
cfg.S_num          = S_num;
cfg.S_len          = S_len;
cfg.freqs          = freqs;
cfg.freq_init      = freq_init;
cfg.unmodeled_freq = unmodeled_freq;
cfg.unmodeled_amp  = unmodeled_amp;
cfg.theta_true     = theta_true;
cfg.theta_min      = theta_min;
cfg.theta_max      = theta_max;
cfg.sgn0           = sgn0;
cfg.sgn            = sgn;

switch run_mode
    case 'minimal'
        cases = make_minimal_case();

    case 'discrete'
        cases = make_discrete_case(true);

    case 'discrete-modeled'
        cases = make_discrete_case(false);

    case 'all'
        cases = [make_minimal_case(), ...
                 make_discrete_case(true)];

    otherwise
        error('Unknown run_mode: %s. Use ''all'', ''minimal'', ''discrete'', or ''discrete-modeled''.', run_mode);
end

%% ====== 仿真执行 ======
for i = 1:numel(cases)
    result_i = run_case(cfg, cases(i));
    if i == 1
        results = repmat(result_i, numel(cases), 1);
    else
        results(i) = result_i;
    end
end

%% ====== 绘图 ======
plot_results(cfg, results);

end

%% ========================================================================
%  实验配置: 最小尺度修复版
%% ========================================================================
function c = make_minimal_case()
% 保留显式 Euler 内模，但改为按采样索引更新。
c = struct();
c.name              = 'Diagnostic: minimal scale fix';
c.method            = 'minimal';
c.include_unmodeled = true;
c.k                 = 0.08;
c.epsilon           = 1.0e-5;
end

%% ========================================================================
%  实验配置: 精确离散内模版
%% ========================================================================
function c = make_discrete_case(include_unmodeled)
% 对归一化内模振荡器做单采样周期精确积分。
c = struct();
c.name              = 'ANC candidate: exact discrete IM';
c.method            = 'exact';
c.include_unmodeled = include_unmodeled;
c.k                 = 0.10;
c.epsilon           = 1.0e-4;

if ~include_unmodeled
    c.name = 'Sanity check: modeled noise only';
end
end

%% ========================================================================
%  单组实验仿真
%% ========================================================================
function result = run_case(cfg, c)
m = numel(cfg.freqs);

theta = (2*pi*cfg.freq_init/cfg.fs).^2;
w0    = 0;
w     = zeros(2, m);       % 第1行: 正弦估计; 第2行: 正交状态
S_buf = zeros(cfg.S_len, 1);

di        = max(1, round(cfg.N / cfg.Nlog));
log_count = floor(cfg.N / di) + 1;
t_log     = zeros(log_count, 1);
f_log     = zeros(m, log_count);
e_log     = zeros(log_count, 1);
u_log     = zeros(log_count, 1);
li        = 1;

e_all       = zeros(cfg.N, 1);
d_primary   = zeros(cfg.N, 1);
u_all       = zeros(cfg.N, 1);
t_all       = zeros(cfg.N, 1);

for n = 1:cfg.N
    t = (n-1) * cfg.Ts;

    primary_noise = sin(2*pi*cfg.freqs(1)*t) + sin(2*pi*cfg.freqs(2)*t);
    if c.include_unmodeled
        primary_noise = primary_noise ...
            + cfg.unmodeled_amp * sin(2*pi*cfg.unmodeled_freq*t);
    end

    anti_noise_at_mic = cfg.S_num * S_buf;
    e                = primary_noise + anti_noise_at_mic;
    y                = e;

    control_on = (t >= cfg.control_start);

    if n > 1 && control_on
        w0 = w0 + cfg.sgn0 * c.k * y;

        w_old     = w;
        theta_old = theta;

        for i = 1:m
            switch c.method
                case 'minimal'
                    w(1,i) = w_old(1,i) + w_old(2,i) + cfg.sgn(i) * c.k * y;
                    w(2,i) = w_old(2,i) - theta_old(i) * w_old(1,i);

                case 'exact'
                    alpha = sqrt(max(theta_old(i), 1e-12));
                    ca    = cos(alpha);
                    sa    = sin(alpha);

                    % 对 x' = A(theta)x + [k*y; 0] 做单步 ZOH 精确积分。
                    w(1,i) = ca*w_old(1,i) + (sa/alpha)*w_old(2,i) ...
                           + cfg.sgn(i) * c.k * (sa/alpha) * y;
                    w(2,i) = -alpha*sa*w_old(1,i) + ca*w_old(2,i) ...
                           + cfg.sgn(i) * c.k * (ca - 1) * y;

                otherwise
                    error('Unknown method: %s.', c.method);
            end

            theta(i) = theta_old(i) + c.epsilon * cfg.sgn(i) * w_old(2,i) * y;
            theta(i) = min(max(theta(i), cfg.theta_min), cfg.theta_max);
        end
    end

    if control_on
        u = -(w0 + sum(w(1,:)));
    else
        u = 0;
    end
    S_buf = [u; S_buf(1:end-1)];

    t_all(n)     = t;
    e_all(n)     = e;
    d_primary(n) = primary_noise;
    u_all(n)     = u;

    if mod(n, di) == 0 || n == 1
        t_log(li)   = t;
        f_log(:,li) = sqrt(theta) * cfg.fs / (2*pi);
        e_log(li)   = e;
        u_log(li)   = u;
        li          = li + 1;
    end
end

t_log = t_log(1:li-1);
f_log = f_log(:,1:li-1);
e_log = e_log(1:li-1);
u_log = u_log(1:li-1);

before_idx = t_all < cfg.control_start;
after_idx  = t_all >= cfg.control_start;
tail_idx   = max(1, floor(cfg.N/2)):cfg.N;

result.name              = c.name;
result.method            = c.method;
result.include_unmodeled = c.include_unmodeled;
result.k                 = c.k;
result.epsilon           = c.epsilon;
result.theta_final       = theta;
result.freq_final        = sqrt(theta) * cfg.fs / (2*pi);
result.d_rms             = sqrt(mean(d_primary.^2));
result.e_rms_before      = sqrt(mean(e_all(before_idx).^2));
result.e_rms_after       = sqrt(mean(e_all(after_idx).^2));
result.e_rms_all         = sqrt(mean(e_all.^2));
result.e_rms_tail        = sqrt(mean(e_all(tail_idx).^2));
result.u_rms_tail        = sqrt(mean(u_all(tail_idx).^2));
result.t_all             = t_all;
result.e_all             = e_all;
result.u_all             = u_all;
result.t_log             = t_log;
result.f_log             = f_log;
result.e_log             = e_log;
result.u_log             = u_log;

fprintf('\n--- %s ---\n', c.name);
fprintf('primary noise includes unmodeled 300Hz: %s\n', bool_text(c.include_unmodeled));
fprintf('k = %.4g, epsilon = %.4g\n', c.k, c.epsilon);
fprintf('theta final: [%.6f %.6f], true: [%.6f %.6f]\n', theta, cfg.theta_true);
fprintf('freq final : [%.3f %.3f] Hz, true: [%.0f %.0f] Hz\n', ...
    result.freq_final, cfg.freqs);
fprintf('RMS primary %.4f | error before %.4f | error after %.4f | error tail %.4f | u tail %.4f\n', ...
    result.d_rms, result.e_rms_before, result.e_rms_after, ...
    result.e_rms_tail, result.u_rms_tail);
end

%% ========================================================================
%  绘图
%% ========================================================================
function plot_results(cfg, results)
colors = lines(numel(results));

figure('Name', 'ANC modeled noise frequency estimates', 'Position', [100, 80, 980, 420]);
for i = 1:2
    subplot(1,2,i);
    hold on;
    for r = 1:numel(results)
        plot(results(r).t_log, results(r).f_log(i,:), ...
            'LineWidth', 1.4, 'Color', colors(r,:));
    end
    xline(cfg.control_start, 'k:', 'ANC on', 'LineWidth', 1.1);
    yline(cfg.freqs(i), 'k--', 'LineWidth', 1.0);
    title(sprintf('Modeled noise frequency f_%d', i));
    xlabel('Time (s)');
    ylabel('Hz');
    grid on;
end
legend({results.name}, 'Location', 'best');

figure('Name', 'ANC error microphone and control signal', 'Position', [140, 120, 980, 520]);

subplot(2,1,1);
hold on;
for r = 1:numel(results)
    plot(results(r).t_all, results(r).e_all, 'LineWidth', 0.7, 'Color', colors(r,:));
end
xline(cfg.control_start, 'k:', 'ANC on', 'LineWidth', 1.1);
title('Error microphone signal e(t)');
xlabel('Time (s)');
grid on;
legend({results.name}, 'Location', 'best');

subplot(2,1,2);
hold on;
for r = 1:numel(results)
    plot(results(r).t_all, results(r).u_all, 'LineWidth', 0.7, 'Color', colors(r,:));
end
xline(cfg.control_start, 'k:', 'ANC on', 'LineWidth', 1.1);
title('Control loudspeaker input u(t)');
xlabel('Time (s)');
grid on;

figure('Name', 'ANC RMS before and after control', 'Position', [180, 160, 980, 360]);
if numel(results) == 1
    rms_values = [results.e_rms_before, results.e_rms_after, results.e_rms_tail];
    bar(rms_values);
    set(gca, 'XTickLabel', {'Before ANC', 'After ANC', 'Tail window'});
else
    rms_values = [[results.e_rms_before]' [results.e_rms_after]' [results.e_rms_tail]'];
    bar(rms_values);
    set(gca, 'XTickLabel', {results.name});
    xtickangle(15);
    legend('Before ANC (0-10s)', 'After ANC', 'Tail window', 'Location', 'best');
end
ylabel('Error microphone RMS');
title('ANC error reduction');
grid on;
end

%% ========================================================================
%  辅助函数
%% ========================================================================
function text = bool_text(value)
if value
    text = 'yes';
else
    text = 'no';
end
end
