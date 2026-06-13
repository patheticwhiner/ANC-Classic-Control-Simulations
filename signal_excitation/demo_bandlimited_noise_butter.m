clear; close all; clc;
%% 设计带通滤波器
fs = 8e3;                   % 采样频率 (Hz)
f_low = 200;                  % 低截止频率 (Hz)
f_high = 500;                % 高截止频率 (Hz)
order = 4;                   % 滤波器阶数

% 将实际频率转换为归一化截止频率（相对于奈奎斯特频率fs/2）
% butter函数要求归一化频率范围为[0,1]
fc = [f_low f_high]/(fs/2);  

% 生成Butterworth带通滤波器
[b, a] = butter(order, fc, 'bandpass');

% 转换为状态空间模型
[Aw, Bw, Cw, Dw] = tf2ss(b, a);
state_dim = size(Aw, 1);

% 生成白噪声输入
N = 10000;                   % 数据点数
n = randn(N, 1);             % 高斯白噪声

% 仿真状态空间模型
w = zeros(N, 1);
x = zeros(state_dim, 1);     % 初始化状态
for k = 1:N
    x = Aw * x + Bw * n(k);
    w(k) = Cw * x + Dw * n(k);
end

% 绘制时域波形
figure; subplot(2,1,1);
t = (0:N-1)/fs;  % 时间向量（秒）
plot(t, w); grid on;
xlabel('Time (s)'); ylabel('Amplitude');
title('Bandlimited Noise (Time Domain)');
% 计算信号的频谱
W_fft = fft(w); % 计算FFT
frequencies = (0:N/2-1)*(fs/N); % 频率轴 (Hz)
W_magnitude = abs(W_fft(1:N/2)); % 取前半部分幅值
% 绘制频谱
subplot(2,1,2);
plot(frequencies, 20*log10(W_magnitude));
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Frequency Spectrum of Bandlimited Noise');
grid on;

% 计算并绘制频谱（使用实际频率Hz）
figure; subplot(2,1,1);
[h, f] = freqz(b, a, 1024, fs);  % 注意这里传入fs参数
mag = 20*log10(abs(h));
plot(f, mag); grid on;
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Filter Frequency Response');
grid on;
% 标注截止频率
hold on;
max_mag = max(mag);
yline(max_mag-3, '--r', '-3dB', 'LabelHorizontalAlignment', 'left');
xline(f_low, '--g', [num2str(f_low) ' Hz'], 'LabelVerticalAlignment', 'bottom');
xline(f_high, '--g', [num2str(f_high) ' Hz'], 'LabelVerticalAlignment', 'bottom');
hold off;

% 计算信号的频谱
subplot(2,1,2);
[Pxx, f] = pwelch(w, hanning(512), 256, 1024, fs);
plot(f, 10*log10(Pxx));
xlabel('Frequency (Hz)');
ylabel('Power/Frequency (dB/Hz)');
title('Power Spectral Density of Bandlimited Noise');
grid on; % ylim([-300 0]);
% 标注截止频率
hold on;
xline(f_low, '--g', [num2str(f_low) ' Hz'], 'LabelVerticalAlignment', 'bottom');
xline(f_high, '--g', [num2str(f_high) ' Hz'], 'LabelVerticalAlignment', 'bottom');
hold off;

save('demo_bandlimited_noise_butter.mat', 'Aw', 'Bw', 'Cw', 'Dw');

%% test:根据传函估算截止频率
sys = ss(Aw, Bw, Cw, Dw, 1/fs);  % 创建离散系统模型
[H, w] = freqresp(sys);  % w实际上是角频率(rad/s)
% 当sys包含采样时间时，freqresp返回的是实际角频率，需要除以2π转换为Hz
f = w/(2*pi);  % 将角频率(rad/s)转换为Hz
H_mag = squeeze(20*log10(abs(H)));
% 找到 -3dB 截止频率
max_mag = max(H_mag);
cutoff_idx = find(H_mag >= max_mag - 3);
f_low_actual = f(min(cutoff_idx));   % 实际低截止频率
f_high_actual = f(max(cutoff_idx));  % 实际高截止频率
% 绘制频率响应并标注实际截止频率
figure;
plot(f, H_mag, 'b', 'LineWidth', 1.5); % 绘制频率响应
hold on;
% 标注 -3dB 线
yline(max_mag - 3, '--r', '-3dB', 'LabelHorizontalAlignment', 'left');
% 标注实际低、高截止频率
xline(f_low_actual, '--g', ['f_{low} = ' num2str(f_low_actual, '%.2f') ' Hz'], 'LabelVerticalAlignment', 'bottom');
xline(f_high_actual, '--g', ['f_{high} = ' num2str(f_high_actual, '%.2f') ' Hz'], 'LabelVerticalAlignment', 'bottom');
% 图像设置
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Estimated Filter Frequency Response with Cutoff Frequencies');
grid on;
legend('Frequency Response', '-3dB Line', 'Estimated Cutoff Frequencies');
hold off;