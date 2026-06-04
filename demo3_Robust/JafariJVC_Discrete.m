%[text] # 未知周期干扰自适应抑制的鲁棒稳定性与性能（Jafari & Ioannou, 2015）
%%
%[text] ## 仿真系统设置
% Example 6.1 仿真：最优K(z,theta)设计，满足 ||G0(z)F(z)K(z,theta)||_inf=1
clear; clc; close all;

% 离散系统G0(z)
Ts = 1/480;
z = tf('z', Ts);
fs = 1/Ts;
G0 = (-0.00146*(z-0.1438)*(z-1)) / ((z - 0.7096)*(z^2 - 0.04369*z + 0.01392));

% F(z) = 1
F = tf(1,1,Ts);

% 扰动参数
omega0 = 0.0521; % rad/sample, 25rad/sec
Nsim = 480*50;   % 50秒
t = (0:Nsim-1)*Ts;
d = sin(omega0*t*fs) + 0.02*randn(1,Nsim);

figure;
plot(t, d); grid on;
xlabel('Time(s)'); ylabel('Magnitude');
xlim([0,5]);
% --------- G0(z) 幅频响应（Hz横轴，freqz实现） ---------
figure;
[H, f_bode] = freqz(G0.num{1}, G0.den{1}, 1000, fs); % H:复数响应, f_bode:Hz
subplot(2,1,1);
semilogx(f_bode, 20*log10(abs(H)), 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on;
title('Magnitude of G_0(z) (freqz)');
subplot(2,1,2);
semilogx(f_bode, angle(H)*180/pi, 'LineWidth', 1.5);
xlabel('Frequency (Hz)'); ylabel('Phase (deg)'); grid on;
title('Phase of G_0(z) (freqz)');
%%
%[text] ## Example 1 固定滤波器设计

%%
%[text] ### （1）阶次对控制效果的影响分析

%%
%[text] ### （2）控制器的频率范围鲁棒性分析

%%
%[text] ## Example 2 自适应滤波器设计


%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"inline","rightPanelPercent":28.6}
%---
