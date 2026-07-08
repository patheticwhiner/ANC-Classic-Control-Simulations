% 比较级联和并联结构的传递函数
clear; clc; close all; 

% 参数设置
q = 4; % 阶数
zeta = linspace(0.2, 0.5, q);           % 阻尼比
omega = linspace(2, 5, q);              % 固有频率
theta = linspace(100, 1000, q);             % 极点参数，间隔拉大以便谐振峰分开
z0 = 1; % 零点

% 并联结构
sys_parallel = parallel_tf(q, zeta, omega, theta, z0);

% 级联结构
sys_cascade = cascade_tf(q, zeta, omega, theta, z0);

% 脉冲响应功率谱分析（pwelch），并标记谐振峰频率
t = 0:1e-4:10; % 时间向量，步长减小以提高频率分辨率
y_parallel = impulse(sys_parallel, t);
y_cascade = impulse(sys_cascade, t);

% 响应对比
figure;
subplot(2,2,1);
plot(t, y_parallel);
title('并联结构脉冲响应');
subplot(2,2,2);
plot(t, y_cascade);
title('级联结构脉冲响应');

subplot(2,2,3);
fs = 1/(t(2)-t(1)); % 采样频率
nfft = 2^nextpow2(length(t)); % FFT长度
Y_parallel = fft(y_parallel, nfft);
f_fft = (0:nfft-1)*(fs/nfft); % 频率轴
mag_parallel = abs(Y_parallel)/max(abs(Y_parallel)); % 幅值归一化
half_n = floor(nfft/2)+1;
plot(f_fft(1:half_n), mag_parallel(1:half_n));
xlabel('频率 (Hz)');
ylabel('归一化幅值');
title('并联结构脉冲响应幅度频谱 (FFT)');
xlim([0, fs/2]);
grid on;
hold on;
for i = 1:length(theta)
    xline(sqrt(theta(i)), '--r', ['\omega_0=', num2str(sqrt(theta(i)), '%.2f')], 'LabelOrientation','horizontal', 'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
end
xlim([0,10]);
hold off;

subplot(2,2,4);
Y_cascade = fft(y_cascade, nfft);
f_fft_c = (0:nfft-1)*(fs/nfft);
mag_cascade = abs(Y_cascade)/max(abs(Y_cascade));
plot(f_fft_c(1:half_n), mag_cascade(1:half_n));
xlabel('频率 (Hz)');
ylabel('归一化幅值');
title('级联结构脉冲响应幅度频谱 (FFT)');
xlim([0, fs/2]);
grid on;
hold on;
for i = 1:length(theta)
    xline(sqrt(theta(i)), '--r', ['\omega_0=', num2str(sqrt(theta(i)), '%.2f')], 'LabelOrientation','horizontal', 'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
end
xlim([0,10]);
hold off;

figure;
bode(sys_parallel, 'b', sys_cascade, 'r.');
legend('并联结构','级联结构'); grid on;
title('幅频特性对比');
hold on;
for i = 1:length(theta)
    xline(sqrt(theta(i)), '--k', ['\theta_', num2str(i)], 'LabelOrientation','horizontal', 'LabelVerticalAlignment','bottom', 'HandleVisibility','off');
end
hold off;


function sys = parallel_tf(q, zeta, omega, theta, z0)
    % 并联结构传递函数
    sys = 0;
    s = tf('s');
    for i = 1:q
        sys = sys + (s^2 + 2*zeta(i)*omega(i)*s + theta(i)) / (s^2 + theta(i));
    end
    sys = sys + (s + z0)/s;
end

function sys = cascade_tf(q, zeta, omega, theta, z0)
    % 级联结构传递函数
    s = tf('s');
    sys = 1;
    for i = 1:q
        sys = sys * ((s^2 + 2*zeta(i)*omega(i)*s + theta(i)) / (s^2 + theta(i)));
    end
    sys = sys * ((s + z0)/s);
end