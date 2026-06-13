close all; clc;

%% 定义参数 - 修改这些参数即可改变整个仿真
noise_cov = 0.001;    % 噪声协方差（功率）
sample_time = 0.0001; % 采样时间(秒)
random_seed = 23341; % 随机数种子
sim_time = 10^4*sample_time;      % 仿真时间(秒)

%% 打开模型和配置
model_name = 'simBandLtdWN';
open_system(model_name);

% 设置Band-Limited White Noise模块参数
set_param([model_name '/Band-Limited White Noise'], 'Cov', num2str(noise_cov));
set_param([model_name '/Band-Limited White Noise'], 'Ts', num2str(sample_time));
set_param([model_name '/Band-Limited White Noise'], 'Seed', num2str(random_seed));

% 设置其他模块参数
set_param(model_name, 'StopTime', num2str(sim_time));

%% 运行仿真
out = sim(model_name);

%% 提取数据和可视化
% 提取数据
noise_data = out.get('whitenoise');
t = 0:sample_time:sim_time; % 时间轴（与仿真步长一致）

% 去除趋势项
noise_data = detrend(noise_data);

% 绘制时域波形
figure;
subplot(3,1,1);
plot(t, noise_data);
grid on;
title(['Sample Time = ' num2str(sample_time) 's, Cov = ' num2str(noise_cov)]);
xlabel('Time (s)'); ylabel('Amplitude');

% 计算自相关函数
[acf, lags] = xcorr(noise_data, 'coeff');
subplot(3,1,2);
plot(lags*sample_time, acf); % 时滞转换为秒
grid on;
title('自相关函数');
xlabel('Lag (s)'); % 横轴单位为秒
ylabel('Normalized Correlation'); % 纵轴为归一化相关值
xlim([-0.1, 0.1]);

% 计算功率谱密度（PSD）
fs = 1/sample_time; % 采样频率
signal_length = length(noise_data);

% 自适应窗口长度
window_length = min(1024, floor(signal_length/2));
noverlap = floor(window_length/2);
nfft = max(1024, 2^nextpow2(signal_length));

% 使用自适应参数计算PSD
[pxx, f] = pwelch(noise_data, window_length, noverlap, nfft, fs);

% 计算理论PSD值（白噪声的理论PSD是常数）
theoretical_psd = 2 * noise_cov;

% 计算均值
avg_measured_psd = mean(pxx);

% 绘制实际PSD和理论PSD
subplot(3,1,3);
hold on;
semilogy(f, pxx, 'b-', 'LineWidth', 0.7);
semilogy([f(1) f(end)], [avg_measured_psd avg_measured_psd], 'b-.', 'LineWidth', 0.7);
semilogy([f(1) f(end)], [theoretical_psd theoretical_psd], 'r--', 'LineWidth', 1.2);
hold off;
grid on;
title(['PSD (fs = ' num2str(fs) ' Hz, 窗长 = ' num2str(window_length) ')']);
xlabel('Frequency (Hz)'); ylabel('Power/Frequency (V^2/Hz)');
legend('测量PSD', 'PSD均值', '理论PSD', 'Location', 'southwest');

% 添加理论值标注
text(f(end)*0.7, theoretical_psd*1.2, ['理论值: ' num2str(theoretical_psd, '%.2e') ' V^2/Hz'], ...
     'Color', 'r', 'FontWeight', 'bold');
 
%%
fprintf('估计RMS：%.2f\r\n', sqrt(noise_cov*fs));
fprintf('估计幅值：%.2f\r\n', sqrt(2*noise_cov*fs));
fprintf('计算RMS：%.2f\r\n', rms(noise_data));