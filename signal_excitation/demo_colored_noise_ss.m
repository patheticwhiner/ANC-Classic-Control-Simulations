%% 带限噪声生成与分析
% 本程序通过状态空间模型对白噪声进行滤波，生成带限噪声
% 用于主动噪声控制(ANC)系统的激励信号

%% 系统参数定义
Fs = 1000;                     % 采样频率（Hz），可根据实际应用调整
N = 2000;                   % 信号长度
noise_var = 1; % 噪声方差

%% 定义状态空间模型
A_w2 = [0.5 0.1; 0.2 0.8];  % 状态矩阵（特征值需在单位圆内以保证系统稳定）
B_w2 = [1; 0.5];            % 输入矩阵
C_w2 = [1 0];               % 输出矩阵（选择第一个状态作为输出）
D_w2 = 0;                   % 直通项（设为0表示无直通影响）

%% 分析系统特性频率
% 计算系统矩阵特征值
eig_A = eig(A_w2);
% 验证系统稳定性
if max(abs(eig(A_w2))) >= 1
    warning('系统不稳定! 特征值: %s', mat2str(eig(A_w2)));
end
% 特征值的虚部对应系统的固有频率
continuous_frequencies = abs(imag(eig_A)) / (2 * pi);
discrete_frequencies = continuous_frequencies * Fs;
% 显示特性频率信息
disp('连续时间系统特性频率 (Hz):');
disp(continuous_frequencies);
disp('离散时间系统特性频率 (Hz):');
disp(discrete_frequencies);

%% 生成白噪声输入
n = randn(N, noise_var);            % 生成均值为0、方差为1的高斯白噪声

%% 通过状态空间模型计算带限噪声
x = [0; 0];                 % 初始化状态向量
w2 = zeros(N, 1);           % 预分配输出向量内存

% 状态空间迭代计算
for k = 1:N
    x = A_w2 * x + B_w2 * n(k);  % 状态更新方程
    w2(k) = C_w2 * x + D_w2 * n(k); % 输出方程
end

%% 频域分析
frequencies = linspace(0, Fs/2, N/2+1); % 频率向量（0到奈奎斯特频率）
W2_fft = fft(w2);           % 计算快速傅里叶变换
W2_magnitude = abs(W2_fft(1:N/2+1)); % 提取幅值谱（仅保留正频率部分）

% 计算功率谱密度
[psd_w2, f_psd] = pwelch(w2, [], [], [], Fs);
figure; plot(f_psd, 10*log10(psd_w2)); grid on;
xlabel('频率 (Hz)'); ylabel('功率/频率 (dB/Hz)');
title('带限噪声功率谱密度');

%% 可视化结果
figure;
% 时域信号
subplot(2,1,1); 
plot(w2); grid on;
xlabel('采样点'); ylabel('幅值');
title('生成的带限噪声信号');

% 频域表示
subplot(2,1,2); 
plot(frequencies, W2_magnitude);
xlabel('频率 (Hz)');
ylabel('幅值');
title('带限噪声的频谱');
grid on;