%[text] # 真实声学次级通路 ANC — 完整仿真脚本
%[text] 被控对象: ARMAX\_SYSID\_30303022 (30阶离散模型, d=22, fs=48kHz)
%[text] 方法: Jafari F(z) + K(z,θ) — 固定FIR(H∞插值) + 自适应(NLMS)
%[text] 测试: 固定频率 / 线性扫频 / 非线性调频
%[text] 结构: Part 0 公共基础 → Parts 1/2/3 独立测试案例
clear; close all; clc;
%%
%[text] ## Part 0 · 问题总述与公共基础
%[text] 被控对象: 实测声学管道 ARMAX(30,30,30,22) @ 48kHz
%[text] 9个非最小相位零点 → F(z) 频谱展平 → 固定 FIR(H∞插值)
%[text] 自适应: NLMS (μ=0.05, N=64) — 替代RLS以避免P矩阵在非平稳环境衰减
%% 0.1 模型加载
modelFile = fullfile('..', 'dataset', 'armax_30303022_2026-01-20.mat');
load(modelFile, 'ARMAXmodel');

A = ARMAXmodel.model.A;             B_poly = ARMAXmodel.model.B;
d = ARMAXmodel.orders(4);           fs = ARMAXmodel.fs;
Ts = 1/fs;
B = B_poly(d+1:end);
nA = length(A)-1;  nB = length(B)-1;

fprintf('========== 真实声学路径 ANC 仿真 ==========\n');
fprintf('  模型: ARMAX(%d,%d,%d,%d), fs=%d Hz\n', ...
    ARMAXmodel.orders(1), ARMAXmodel.orders(2), ARMAXmodel.orders(3), d, fs);
fprintf('  极点: max|p|=%.4f,  零点: %d个NMP (max|z|=%.4f)\n', ...
    max(abs(roots(A))), sum(abs(roots(B))>1), max(abs(roots(B))));

%% 0.2 F(z) 频谱展平
%[text] 步骤: 提取NMP零点 → 单位圆镜像反射(z→1/z̄) → 重构B̃(z)
%[text] F(z) = κ₀·L(z)·A(z⁻¹)/B̃(z⁻¹),  L(z) = (1-ρ)/(1-ρz⁻¹) 一阶低通
z_orig = roots(B);
z_nmp  = z_orig(abs(z_orig)>1);
z_mp   = z_orig(abs(z_orig)<=1);
B_tilde = real(poly([z_mp; 1./conj(z_nmp)]));
B_tilde = B_tilde * (B(1)/B_tilde(1));
rho = exp(-2*pi*2000*Ts);
L_num = 1-rho;  L_den = [1, -rho];

F_num = A;
F_den = conv(L_den, B_tilde);
F_den = F_den / F_den(1);
nF_num = length(F_num);  nF_den = length(F_den);

H_num = conv([zeros(1,d), B], F_num);
H_den = conv(A, F_den);

fprintf('  F(z): %d个NMP零点反射+1阶LPF, 阶次 num=%d/den=%d\n', ...
    length(z_nmp), nF_num-1, nF_den-1);
%% 0.3 固定控制器 (设计频率 = 第一共振峰)
%[text] FIR K(z,θ)=Σθᵢz⁻ⁱ, N=120. 插值约束: G\_eff(e^{jω₀})·K(e^{jω₀})=1
w_scan = linspace(0.001, pi, 8000);
G_scan = abs(freqz(B, A, w_scan));
[pks, locs] = findpeaks(G_scan, 'MinPeakHeight',1.0, 'MinPeakDistance',50);
f_peaks = w_scan(locs)*fs/(2*pi);
f0 = f_peaks(1);
omega0 = 2*pi*f0/fs;

N_fir = 120;
zi0 = exp(1j*omega0);
G_eff_0 = polyval_freqz(H_num, H_den, omega0);
Aeq = [real(G_eff_0*zi0.^(-(0:N_fir-1))); imag(G_eff_0*zi0.^(-(0:N_fir-1)))];
theta_fixed = Aeq \ [1; 0];

Nfreq = 2000;  w_grid = linspace(0,pi,Nfreq);
G_eff_grid = freqz(H_num, H_den, w_grid).';
gamma_fixed = max(abs(G_eff_grid(:).*(exp(-1j*(0:N_fir-1)'.*w_grid)'*theta_fixed)));

fprintf('  固定FIR: f0=%.1fHz, N=%d, ||θ||=%.2f, γ=%.4f\n', ...
    f0, N_fir, norm(theta_fixed), gamma_fixed);
%[text] 被控对象差分方程 + 固定/自适应仿真句柄
%% 0.4 仿真引擎
a_plant = A;  b_plant = [zeros(1,d), B];  Lb = length(b_plant);
k_start = max(N_fir, nA+Lb) + 1;

rms_win = @(y,t1,t2) rms(y(round(t1*fs):round(t2*fs)));
supp_db = @(y_o,y_c,t1,t2) 20*log10(rms_win(y_o,t1,t2)/max(rms_win(y_c,t1,t2),1e-10));

sim_fixed = @(d_sig,Tsim) sim_fixed_impl(d_sig,Tsim,fs,Ts,N_fir,theta_fixed, ...
    a_plant,b_plant,nA,Lb,nF_num,nF_den,F_num,F_den,k_start);
sim_nlms  = @(d_sig,Tsim,N_adapt,mu) sim_nlms_impl(d_sig,Tsim,fs,Ts,N_adapt,mu, ...
    a_plant,b_plant,nA,Lb,nF_num,nF_den,F_num,F_den,H_num,H_den,k_start);

fprintf('\n========== 公共基础完成 ==========\n\n');
%%
%[text] ## Part 1 · 固定频率扰动 — 稳态性能验证
%[text] 问题: 单频正弦扰动 @ 共振峰频率 f₀≈334.6Hz, |G₀|=+2.5dB
%[text] 目标: 验证 F(z)+K(z) 在理想条件下的降噪能力
%[text] 对比: 开环 vs 固定FIR vs 自适应NLMS
%% 1.1 扰动设计
Tsim1 = 10;  Nsim1 = round(Tsim1*fs);  t1 = (0:Nsim1-1)*Ts;
rng(42);
d1 = 0.6*sin(2*pi*f0*t1) + 0.02*randn(1,Nsim1);
d1 = d1/max(abs(d1))*0.8;
y_open1 = filter(b_plant, a_plant, d1);

%% 1.2 固定控制器

y_fix1 = sim_fixed(d1, Tsim1);
supp_f1 = supp_db(y_open1, y_fix1, 2, Tsim1);

%% 1.3 自适应 NLMS
[y_nlms1, th1] = sim_nlms(d1, Tsim1, 64, 0.05);
supp_a1 = supp_db(y_open1, y_nlms1, 2, Tsim1);

%% 1.4 报告与图解
fprintf('--- Part 1: 固定频率 %.0f Hz ---\n', f0);
fprintf('  固定(N=%d):       %.1f dB\n', N_fir, supp_f1);
fprintf('  自适应(N=%d,μ=%.2f): %.1f dB  ||θ||=%.1f\n', 64, 0.05, supp_a1, norm(th1(:,end)));

figure('Name','Part1: 固定频率');  subplot(2,1,1);
plot(t1, y_open1, 'Color',[.7 .7 .7]); hold on;
plot(t1, y_fix1, 'r-', 'LineWidth',1);  plot(t1, y_nlms1, 'b-', 'LineWidth',1);
legend('开环',sprintf('固定 %.1fdB',supp_f1),sprintf('NLMS %.1fdB',supp_a1));
title(sprintf('固定频率 %.0f Hz', f0));  xlim([4 4.5]);  grid on;
subplot(2,1,2);
[Po, fp] = pwelch(y_open1(end-4096:end),hann(1024),512,2048,fs);
[Pc, ~]  = pwelch(y_fix1(end-4096:end),hann(1024),512,2048,fs);
[Pa, ~]  = pwelch(y_nlms1(end-4096:end),hann(1024),512,2048,fs);
semilogy(fp,Po,'Color',[.7 .7 .7]); hold on;
semilogy(fp,Pc,'r-');  semilogy(fp,Pa,'b-');
xlim([0 1000]);  legend('开环','固定','NLMS');
title('稳态 PSD');  xlabel('Hz');  grid on;
sgtitle('Part 1: 固定频率扰动 — 稳态性能验证');
fprintf('\n');
%%
%[text] ## Part 2 · 线性扫频扰动 — 跟踪能力验证
%[text] 问题: f(t) = 300 + 70·t/Tsim Hz, 线性扫过 300→370 Hz
%[text] 关键: 固定控制器设计在 f₀≈334.6Hz — 偏离时性能退化
%[text] NLMS 应跟踪频率变化，维持全频段抑制

%% 2.1 扰动设计
Tsim2 = 8;  Nsim2 = round(Tsim2*fs);  t2 = (0:Nsim2-1)*Ts;
f_inst2 = 300 + 70*t2/Tsim2;
phase2  = 2*pi*(300*t2 + 35*t2.^2/Tsim2);
rng(42);  d2 = 0.6*sin(phase2) + 0.02*randn(1,Nsim2);
d2 = d2/max(abs(d2))*0.8;
y_open2 = filter(b_plant, a_plant, d2);

%% 2.2 固定控制器
y_fix2 = sim_fixed(d2, Tsim2);

%% 2.3 自适应 NLMS
[y_nlms2, ~] = sim_nlms(d2, Tsim2, 64, 0.05);

%% 2.4 分段对比
fprintf('--- Part 2: 线性扫频 300→370 Hz ---\n');
fprintf('  %8s %10s %12s %12s\n', 't(s)','f(Hz)','固定(dB)','NLMS(dB)');
t_chk = 3:0.5:Tsim2-0.5;
supp_f2 = zeros(size(t_chk));  supp_a2 = zeros(size(t_chk));
for i = 1:length(t_chk)
    tc = t_chk(i);
    supp_f2(i) = supp_db(y_open2, y_fix2,  tc-0.25, tc+0.25);
    supp_a2(i) = supp_db(y_open2, y_nlms2, tc-0.25, tc+0.25);
    fprintf('  %8.1f %10.1f %12.1f %12.1f\n', tc, 300+70*tc/Tsim2, supp_f2(i), supp_a2(i));
end
fprintf('  早期(3-4s): 固定 %.1fdB  NLMS %.1fdB\n', mean(supp_f2(1:3)), mean(supp_a2(1:3)));
fprintf('  晚期(6-7s): 固定 %.1fdB  NLMS %.1fdB\n', mean(supp_f2(end-2:end)), mean(supp_a2(end-2:end)));

figure('Name','Part2: 线性扫频');
subplot(2,1,1);
plot(t2, y_open2,'Color',[.7 .7 .7]); hold on;
plot(t2, y_fix2,'r-','LineWidth',.8);  plot(t2, y_nlms2,'b-','LineWidth',1);
legend('开环','固定','NLMS');  xlim([2 Tsim2]);
title(sprintf('扫频 %.0f→%.0f Hz (%.1f Hz/s)',300,370,70/Tsim2));  grid on;
subplot(2,1,2);
plot(t_chk, supp_f2,'r-o','LineWidth',1.3);  hold on;
plot(t_chk, supp_a2,'b-s','LineWidth',1.3);  yline(0,'k-');
xlabel('t(s)');  ylabel('抑制(dB)');  legend('固定','NLMS');  grid on;
title(sprintf('分段抑制 (设计频率 %.0f Hz)', f0));
sgtitle('Part 2: 线性扫频扰动 — 跟踪能力验证');
fprintf('\n');

%%
%[text] ## Part 3 · 非线性调频扰动 — 频率反转测试
%[text] 问题: f(t) = 335 + 30·sin(2π·1.5·t) Hz, 正弦调制
%[text] 关键: 频率方向周期性反转 — 最严苛的NLMS跟踪测试

%% 3.1 扰动设计
Tsim3 = 10;  Nsim3 = round(Tsim3*fs);  t3 = (0:Nsim3-1)*Ts;
f0_fm = 335;  df_fm = 30;  fm = 1.5;
f_inst3 = f0_fm + df_fm*sin(2*pi*fm*t3);
phase3  = 2*pi*f0_fm*t3 - (df_fm/fm)*cos(2*pi*fm*t3);
rng(42);  d3 = 0.6*sin(phase3) + 0.02*randn(1,Nsim3);
d3 = d3/max(abs(d3))*0.8;
y_open3 = filter(b_plant, a_plant, d3);

%% 3.2 固定控制器
y_fix3 = sim_fixed(d3, Tsim3);

%% 3.3 自适应 NLMS
[y_nlms3, ~] = sim_nlms(d3, Tsim3, 64, 0.05);

%% 3.4 分段对比
fprintf('--- Part 3: 正弦FM 335±30 Hz ---\n');
fprintf('  %8s %10s %12s %12s\n', 't(s)','f(Hz)','固定(dB)','NLMS(dB)');
t_chk3 = 3:0.5:Tsim3-0.5;
for i = 1:length(t_chk3)
    tc = t_chk3(i);
    sf = supp_db(y_open3, y_fix3,  tc-0.25, tc+0.25);
    sa = supp_db(y_open3, y_nlms3, tc-0.25, tc+0.25);
    fprintf('  %8.1f %10.1f %12.1f %12.1f\n', tc, f0_fm+df_fm*sin(2*pi*fm*tc), sf, sa);
end

%% 3.5 图解
figure('Name','Part3: 非线性FM');
subplot(2,2,[1 2]);
plot(t3, y_open3,'Color',[.7 .7 .7]); hold on;
plot(t3, y_fix3,'r-','LineWidth',.7);  plot(t3, y_nlms3,'b-','LineWidth',1);
legend('开环','固定','NLMS');  xlim([2 Tsim3]);
title(sprintf('正弦FM: %.0f±%.0f Hz @ %.1f Hz 调制',f0_fm,df_fm,fm));  grid on;

subplot(2,2,3);
plot(t3, f_inst3,'k-','LineWidth',1);  hold on;
yline(f0,'r--','设计频率');
xlabel('t(s)');  ylabel('f(Hz)');  title('瞬时频率');  grid on;  xlim([2 Tsim3]);

subplot(2,2,4);
t_rev = 1/fm*0.75 + 2;  idx = round((t_rev-0.05)*fs):round((t_rev+0.05)*fs);
if idx(end)<=Nsim3
    plot(t3(idx),y_open3(idx),'Color',[.7 .7 .7]); hold on;
    plot(t3(idx),y_fix3(idx),'r-');  plot(t3(idx),y_nlms3(idx),'b-','LineWidth',1.2);
    legend('开环','固定','NLMS');
end
title(sprintf('频率反转处 (t≈%.1fs)', t_rev));  grid on;
sgtitle('Part 3: 非线性调频扰动 — 频率反转测试');

fprintf('\n========== 全部测试完成 ==========\n');


%% 辅助函数
%%
%[text] ## 辅助函数
%[text] sim\_fixed\_impl: 固定FIR控制器仿真循环 (IMC + F(z))
%[text] sim\_nlms\_impl:  NLMS自适应控制器仿真循环 (IMC + F(z) + G₀F滤波回归)

function y = sim_fixed_impl(d_sig, Tsim, fs, Ts, N_fir, theta_fixed, ...
        a_plant, b_plant, nA, Lb, nF_num, nF_den, F_num, F_den, k_start)
    Nsim = round(Tsim*fs);  t = (0:Nsim-1)*Ts;
    y=zeros(1,Nsim); u=zeros(1,Nsim); g=zeros(1,Nsim); z=zeros(1,Nsim);
    xF=zeros(nF_num,1);  yF=zeros(nF_den-1,1);
    for k = k_start:Nsim
        g(k) = -a_plant(2:end)*g(k-1:-1:k-nA)' + b_plant*u(k-1:-1:k-Lb)';
        y(k) = g(k)+d_sig(k);  z(k) = y(k)-g(k);
        if t(k)>=2
            v = theta_fixed'*z(k-(1:N_fir))';
            xF=[v; xF(1:nF_num-1)];
            uk = F_num*xF - F_den(2:end)*yF;
            yF=[uk; yF(1:nF_den-2)];
            u(k)=-uk;
            if isnan(u(k))||abs(u(k))>1e4, u(k)=0; end
        end
    end
end

function [y, th] = sim_nlms_impl(d_sig, Tsim, fs, Ts, N_adapt, mu, ...
        a_plant, b_plant, nA, Lb, nF_num, nF_den, F_num, F_den, H_num, H_den, k_start)
    Nsim = round(Tsim*fs);  t = (0:Nsim-1)*Ts;
    y=zeros(1,Nsim); u=zeros(1,Nsim); g=zeros(1,Nsim); z=zeros(1,Nsim);
    zf=zeros(1,Nsim);  th=zeros(N_adapt,Nsim);  theta=zeros(N_adapt,1);
    xF=zeros(nF_num,1);  yF=zeros(nF_den-1,1);
    xH=zeros(length(H_num)-1,1);  yH=zeros(length(H_den)-1,1);
    for k = k_start:Nsim
        g(k) = -a_plant(2:end)*g(k-1:-1:k-nA)' + b_plant*u(k-1:-1:k-Lb)';
        y(k) = g(k)+d_sig(k);  z(k) = y(k)-g(k);
        xH=[z(k); xH(1:end-1)];
        zf(k) = H_num(2:end)*xH - H_den(2:end)*yH;
        yH=[zf(k); yH(1:end-1)];
        if t(k)<2, continue; end
        if k>N_adapt, phi=zf(k-(1:N_adapt))'; else phi=zeros(N_adapt,1); end
        e = z(k)-theta'*phi;
        theta = theta + mu*e*phi/(phi'*phi+0.01);
        th(:,k) = theta;
        v = theta'*z(k-(1:N_adapt))';
        xF=[v; xF(1:nF_num-1)];
        uk = F_num*xF - F_den(2:end)*yF;
        yF=[uk; yF(1:nF_den-2)];
        u(k) = -uk;
        if isnan(u(k))||abs(u(k))>1e4, u(k)=0; end
    end
end

function v = polyval_freqz(num, den, w)
    v = sum(num.*exp(-1j*w*(0:length(num)-1))) / sum(den.*exp(-1j*w*(0:length(den)-1)));
end

%[appendix]{"version":"1.0"}
%[metadata:view]
