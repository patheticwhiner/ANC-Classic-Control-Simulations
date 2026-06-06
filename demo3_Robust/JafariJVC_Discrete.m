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
%[text] ## Example 1 固定 FIR 滤波器设计
%[text] 设计 FIR 滤波器 $K(z,\theta)=\sum_{k=1}^N \theta_k z^{-k}$，使 $G_0(e^{j\omega_0})K(e^{j\omega_0},\theta)=1$，
%[text] 同时最小化 $\|G_0(z)K(z,\theta)\|_\infty$。

%% ====== 1.1 频率网格与插值约束 ======
N_fir = 30;                        % FIR 阶数
w_interp = 0.0521;                 % 扰动频率 (rad/sample)
z_interp = exp(1j * w_interp);

w_grid = linspace(0, pi, 2000);    % H∞ 约束频率网格
z_grid = exp(1j * w_grid);
G_grid = squeeze(freqresp(G0, w_grid)).';

% 构造 Phi 矩阵: Phi(k,:) = [z_k^{-1}, z_k^{-2}, ..., z_k^{-N}]
Phi_grid = zeros(length(w_grid), N_fir);
for k = 1:N_fir
    Phi_grid(:,k) = z_grid(:) .^ (-k);
end
Phi_interp = z_interp .^ (-(1:N_fir));
G_interp = squeeze(freqresp(G0, w_interp));

% 插值方程: G0(w0) * Phi_interp * theta = 1
Aeq = [real(G_interp * Phi_interp); imag(G_interp * Phi_interp)];
beq = [1; 0];  % 实数约束：实部=1, 虚部=0

%% ====== 1.2 最小二乘解（基准） ======
theta_ls = Aeq \ beq;  % 最小二乘解（欠定时取最小范数）
mag_ls = abs(G_grid(:) .* (Phi_grid * theta_ls));
gamma_ls = max(mag_ls);
fprintf('=== 固定 FIR 滤波器 (LS) ===\n');
fprintf('N=%d, H∞范数 γ=%.4f, 插值残差=%.2e\n', N_fir, gamma_ls, norm(Aeq*theta_ls - beq));

%% ====== 1.3 H∞ 优化 (YALMIP+MOSEK, 若可用) ======
if exist('yalmip','file') == 2
    yalmip('clear');
    theta_opt = sdpvar(N_fir, 1, 'full');
    gamma_opt = sdpvar(1, 1);

    Constraints = [Aeq * theta_opt == beq];
    for k = 1:length(w_grid)
        Fk = Phi_grid(k,:) * theta_opt;
        Constraints = [Constraints, ...
            norm([real(G_grid(k)*Fk), imag(G_grid(k)*Fk)], 2) <= gamma_opt];
    end

    options = sdpsettings('solver','mosek','verbose',0);
    sol = optimize(Constraints, gamma_opt, options);

    if sol.problem == 0
        theta_hinf = value(theta_opt);
        gamma_hinf = value(gamma_opt);
        fprintf('\n=== 固定 FIR 滤波器 (H∞优化) ===\n');
        fprintf('γ=%.4f, ||θ||=%.4f\n', gamma_hinf, norm(theta_hinf));
    else
        fprintf('H∞优化未收敛，使用 LS 解。\n');
        theta_hinf = theta_ls;
        gamma_hinf = gamma_ls;
    end
else
    fprintf('YALMIP 未安装，使用 LS 解。\n');
    theta_hinf = theta_ls;
    gamma_hinf = gamma_ls;
end

%% ====== 1.4 控制器频率响应 ======
b_K = [0, theta_hinf(:)'];  % FIR 系数 (含 z^0=0)
[Hk, fk] = freqz(b_K, 1, 1024, fs);

figure('Name','固定FIR控制器');
subplot(2,1,1);
semilogx(fk, 20*log10(abs(Hk)), 'b-', 'LineWidth',1.5); hold on;
xline(omega0*fs/(2*pi), 'r--', sprintf('%.1f Hz', omega0*fs/(2*pi)));
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on;
title(sprintf('K(z) 幅频响应 (N=%d, γ=%.3f)', N_fir, gamma_hinf));

subplot(2,1,2);
[H_loop, f_loop] = freqz(conv(G0.num{1}, b_K), G0.den{1}, 1024, fs);
semilogx(f_loop, 20*log10(abs(H_loop)), 'r-', 'LineWidth',1.5); hold on;
yline(0, 'k--', '0 dB');
xline(omega0*fs/(2*pi), 'r--');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on;
title(sprintf('|G_0(z)K(z)| 频响 (γ_{inf}=%.3f)', gamma_hinf));

%% ====== 1.5 固定控制器闭环仿真 ======
b = G0.num{1};  a = G0.den{1};
u_fixed = zeros(1, Nsim);
y_fixed = zeros(1, Nsim);
x_fixed = zeros(1, Nsim);
antinoise = zeros(1, Nsim);

for k = N_fir+1:Nsim
    if t(k) >= 10  % 10秒后开启
        u_fixed(k) = -sum(theta_hinf'.*x_fixed(k-(1:N_fir)));
        if isnan(u_fixed(k)), u_fixed(k) = 0; end
    end
    antinoise(k) = (-a(2:end)*antinoise(k-1:-1:k-length(a)+1)' + ...
                     b*u_fixed(k-1:-1:k-length(b))')/a(1);
    y_fixed(k) = antinoise(k) + d(k);
    x_fixed(k) = y_fixed(k) - antinoise(k);
end

figure('Name','固定FIR控制器仿真');
subplot(2,1,1);
plot(t, u_fixed, 'r-', 'LineWidth',1); xline(10, 'k:', 'Controller ON');
xlabel('Time (s)'); ylabel('Control Input u'); title('固定FIR控制输入'); grid on;
subplot(2,1,2);
plot(t, d, 'k:', 'LineWidth',1); hold on;
plot(t, y_fixed, 'b-', 'LineWidth',1); xline(10, 'k--');
xlabel('Time (s)'); ylabel('Output y'); title('闭环输出'); grid on;

rms_fixed = rms(y_fixed(15*fs:end));
rms_d = rms(d(15*fs:end));
fprintf('固定控制器: 扰动RMS=%.4f, 输出RMS=%.4f, 抑制=%.1f dB\n', ...
    rms_d, rms_fixed, 20*log10(rms_d/rms_fixed));

%%
%[text] ### （1）阶次对控制效果的影响分析
%[text] 考察 FIR 阶数 $N$ 对 $H_\infty$ 范数和噪声抑制性能的影响。
N_list = [5, 10, 15, 20, 30, 50];
gamma_N = zeros(size(N_list));
supp_N = zeros(size(N_list));

for ni = 1:length(N_list)
    Nn = N_list(ni);
    Phi_n = zeros(length(w_grid), Nn);
    for k = 1:Nn, Phi_n(:,k) = z_grid(:) .^ (-k); end
    Phi_i = z_interp .^ (-(1:Nn));
    Gi = squeeze(freqresp(G0, w_interp));
    Ae = [real(Gi * Phi_i); imag(Gi * Phi_i)];
    Be = [1; 0];
    th_n = Ae \ Be;
    gamma_N(ni) = max(abs(G_grid(:) .* (Phi_n * th_n)));

    % 快速仿真
    y_t = zeros(1, Nsim);
    anti_t = zeros(1, Nsim);
    u_hist = zeros(1, Nsim);  % 存储历史控制输入
    for k = Nn+1:Nsim
        if t(k) >= 10
            u_t = -sum(th_n'.*(-anti_t(k-(1:Nn)) + d(k-(1:Nn))));
        else
            u_t = 0;
        end
        u_hist(k) = u_t;
        anti_t(k) = (-a(2:end)*anti_t(k-1:-1:k-length(a)+1)' + ...
                      b*u_hist(k-1:-1:k-length(b))')/a(1);
        y_t(k) = anti_t(k) + d(k);
    end
    supp_N(ni) = 20*log10(rms(y_t(15*fs:end))/rms_d);
end

figure('Name','阶次影响分析');
yyaxis left;
plot(N_list, gamma_N, 'b-o', 'LineWidth',1.5); ylabel('H_\infty 范数 \gamma');
yyaxis right;
plot(N_list, supp_N, 'r-s', 'LineWidth',1.5); ylabel('抑制 (dB)');
xlabel('FIR 阶数 N'); grid on;
title('FIR 阶数对 H_\infty 范数和抑制效果的影响');
legend('\gamma = ||G_0K||_\infty', '噪声抑制', 'Location','best');

%%
%[text] ### （2）控制器的频率范围鲁棒性分析
%[text] 当扰动频率偏离设计频率 $\omega_0$ 时，检查控制器性能的退化情况。
freq_offset = (-0.3:0.01:0.3) * omega0;  % 频率偏移
supp_robust = zeros(size(freq_offset));

for fi = 1:length(freq_offset)
    w_test = omega0 + freq_offset(fi);
    if w_test <= 0 || w_test >= pi, supp_robust(fi) = NaN; continue; end
    z_test = exp(1j*w_test);
    G_test = squeeze(freqresp(G0, w_test));
    K_test = sum(theta_hinf(:)' .* (z_test .^ (-(1:N_fir))'));
    loop_test = abs(G_test * K_test);
    % 理想抑制取决于 |1 + loop| 的大小
    if loop_test > 0.1
        supp_robust(fi) = -20*log10(abs(1 + loop_test) / loop_test);
    else
        supp_robust(fi) = NaN;
    end
end

figure('Name','频率鲁棒性');
plot((omega0 + freq_offset)*fs/(2*pi), supp_robust, 'b-', 'LineWidth',1.5);
xline(omega0*fs/(2*pi), 'r--', 'Design freq.');
xlabel('Frequency (Hz)'); ylabel('理论抑制 (dB)'); grid on;
title(sprintf('控制器抗频率偏移鲁棒性 (N=%d)', N_fir));

%%
%[text] ## Example 2 自适应滤波器设计
%[text] 基于离散时间递推最小二乘 (RLS) 的自适应 FIR 控制器：
%[text] $K(z,\theta(k)) = \sum_{i=1}^N \theta_i(k) z^{-i}$，
%[text] 回归向量 $\phi(k) = G_0(z)F(z)\alpha(z)[z(k)]$，
%[text] 自适应律采用正规化 RLS + 参数投影。

%% ====== 2.1 自适应参数设置 ======
N_adapt = 30;                     % 自适应FIR阶数
theta_adapt = zeros(N_adapt, 1);  % 参数初值
% 注：G0 在 ω₀ 处增益 ≈ 2.3×10⁻⁴，控制器需 |K(e^{jω₀})| ≈ 4300，
%     对应 ||θ|| ≈ 1600+，故 θ_max 和 P(0) 需匹配此量级
gamma0 = 1.0;                     % 归一化增益
theta_max = 1e7;                  % 参数投影界（匹配固定控制器 ||θ||≈7e5）
P = eye(N_adapt) * 1e6;           % 协方差矩阵 P(0)，匹配参数尺度

%% ====== 2.2 滤波器离散化与状态初始化 ======
Fd = F;  % F 已是离散模型，无需 c2d
[numFd, denFd] = tfdata(Fd, 'v');
Gd = G0;  % G0 已是离散模型，无需 c2d
[numGd, denGd] = tfdata(Gd, 'v');

nF = max(length(numFd), length(denFd)) - 1;
nG = max(length(numGd), length(denGd)) - 1;
G_zf = zeros(nG, 1);    G0_zf = zeros(nG, 1);
F_ctrl_zf = zeros(nF, 1);
% 回归量中每个延迟分量的独立 F 和 G0 滤波器状态
phi_F_zf = zeros(N_adapt, nF);
phi_G0_zf = zeros(N_adapt, nG);

%% ====== 2.3 信号存储 ======
y_adapt = zeros(1, Nsim);  u_adapt = zeros(1, Nsim);
plant_out = zeros(1, Nsim);  z_hist = zeros(1, Nsim);
theta_hist = zeros(N_adapt, Nsim);  P_tr_hist = zeros(1, Nsim);
phi_norm_hist = zeros(1, Nsim);  theta_norm_hist = zeros(1, Nsim);
pred_err_hist = zeros(1, Nsim);  epsilon_hist = zeros(1, Nsim);

fprintf('\n===== 离散时间自适应 FIR 控制器 =====\n');
fprintf('N=%d, P(0)=%.0e I, θ_{max}=%.0e, γ_0=%.0f\n', N_adapt, P(1,1), theta_max, gamma0);
fprintf('回归量: φ_i(k) = G_0(z)F(z)[z(k-i)]\n');

%% ====== 2.4 仿真循环 ======
% 因果顺序: y(k)=G0[u(k-1)]+d(k) → z(k)=y(k)-G0[u(k-1)] → φ_i=G0·F[z(k-i)]
%           → RLS更新θ → u(k)=-F[θ^T z(k-i)]  (控制律用原始 z 延迟，非 φ)

phi_reg = zeros(N_adapt, 1);
z_delays = zeros(N_adapt, 1);  % 原始 z 延迟向量，供控制律使用

for k = 2:Nsim
    % 1. 系统输出: y(k) = G0[u(k-1)] + d(k)
    [y_G_k, G_zf] = filter(numGd, denGd, u_adapt(k-1), G_zf);
    plant_out(k) = y_G_k;
    y_adapt(k) = y_G_k + d(k);

    % 2. 观测信号: z(k) = y(k) - G0[u(k-1)]
    [G0_u_k, G0_zf] = filter(numGd, denGd, u_adapt(k-1), G0_zf);
    z_k = y_adapt(k) - G0_u_k;
    z_hist(k) = z_k;

    % 3. 回归向量: φ_i(k) = G0(z)F(z)[z(k-i)]，同时收集原始 z 延迟
    for i = 1:N_adapt
        z_delay = 0;
        if k > i, z_delay = z_hist(k - i); end
        z_delays(i) = z_delay;  % 保存原始延迟，供控制律 θ^T z(k-i) 使用
        [f_out, phi_F_zf(i,:)] = filter(numFd, denFd, z_delay, phi_F_zf(i,:)');
        [phi_reg(i), phi_G0_zf(i,:)] = filter(numGd, denGd, f_out, phi_G0_zf(i,:)');
    end

    phi_norm_hist(k) = norm(phi_reg);

    % 4. 离散时间正规化 RLS + 参数投影
    m_s_sq = 1 + gamma0 * (phi_reg' * phi_reg);
    pred_err = z_k - theta_adapt' * phi_reg;
    epsilon = pred_err / m_s_sq;
    pred_err_hist(k) = pred_err;
    epsilon_hist(k) = epsilon;

    P_phi = P * phi_reg;
    denom = m_s_sq + phi_reg' * P_phi;
    if denom > 1e-12
        K_gain = P_phi / denom;
        P = P - (K_gain * P_phi');
        theta_new = theta_adapt + K_gain * epsilon;  % 增益向量×归一化误差
        th_norm = norm(theta_new);
        theta_adapt = theta_new * min(1, theta_max / th_norm);  % 投影
    end

    theta_hist(:, k) = theta_adapt;
    P_tr_hist(k) = trace(P);
    theta_norm_hist(k) = norm(theta_adapt);

    % 5. 控制律: u(k) = -F(z)[θ^T(k) z(k-i)]
    %    注意：这里用原始 z 延迟，而非经 G0·F 滤波后的 φ
    Kz_k = theta_adapt' * z_delays;
    [u_adapt(k), F_ctrl_zf] = filter(numFd, denFd, -Kz_k, F_ctrl_zf);

    if mod(k, round(Nsim/8)) == 0
        fprintf('  t=%.1fs: ||θ||=%.3f, tr(P)=%.2e, z=%.4f\n', ...
            t(k), th_norm, trace(P), z_k);
    end
end

%% ====== 2.5 性能分析 ======
final_idx = round(0.5*Nsim):Nsim;
rms_adapt = sqrt(mean(y_adapt(final_idx).^2));
rms_yG = sqrt(mean(plant_out(final_idx).^2));
supp_adapt = 20*log10(rms_adapt/rms_d);
rms_u_adapt = sqrt(mean(u_adapt(final_idx).^2));

fprintf('\n=== 自适应控制性能 ===\n');
fprintf('扰动RMS: %.4f\n', rms_d);
fprintf('总输出RMS: %.4f (抑制 %.1f dB)\n', rms_adapt, supp_adapt);
fprintf('控制RMS: %.1f\n', rms_u_adapt);
fprintf('最终 ||θ||=%.4f, tr(P)=%.4e\n', norm(theta_adapt), trace(P));

% ---- 图5: 时域对比 ----
figure('Name','自适应控制器仿真');
subplot(3,1,1);
plot(t, d, 'Color',[0.7 0.7 0.7]); hold on;
plot(t, plant_out, 'r--', 'LineWidth',1); plot(t, y_adapt, 'b-', 'LineWidth',1);
legend('d(t)', 'G输出', '总输出 y(t)'); title('输出信号对比'); grid on;

subplot(3,1,2);
plot(t, u_adapt, 'g-', 'LineWidth',1); title('控制信号 u(t)'); grid on;

subplot(3,1,3);
plot(t, z_hist, 'k-', 'LineWidth',1); title('观测信号 z(t)'); grid on;

% ---- 图6: 参数演化与收敛 ----
figure('Name','参数演化与收敛');
subplot(2,2,1);
plot(t, theta_hist(1:5:end, :)'); title('θ(t) (每5个)'); grid on;
subplot(2,2,2);
semilogy(t, max(P_tr_hist, 1e-10), 'LineWidth',1); title('trace(P)'); grid on;
subplot(2,2,3);
plot(t, theta_norm_hist, 'LineWidth',1); title('||θ(t)||'); grid on;
subplot(2,2,4);
plot(t, phi_norm_hist, 'LineWidth',1); title('||φ(t)||'); grid on;

% ---- 图7: PSD 频谱对比 ----
figure('Name','PSD频谱');
Nfft = 2048;
[Pyy, f_p] = pwelch(y_adapt(end-Nfft+1:end).*hann(Nfft)', hann(512), 256, 512, fs);
[Pdd, ~]  = pwelch(d(end-Nfft+1:end).*hann(Nfft)', hann(512), 256, 512, fs);
semilogy(f_p(f_p<=50), Pdd(f_p<=50), 'Color',[0.7 0.7 0.7]); hold on;
semilogy(f_p(f_p<=50), Pyy(f_p<=50), 'b-', 'LineWidth',1.2);
xline(omega0*fs/(2*pi), 'r--', sprintf('%.1f Hz', omega0*fs/(2*pi)));
xlabel('Hz'); title('稳态PSD频谱'); grid on;
legend('扰动', '输出'); xlim([0 20]);

% ---- 自适应 vs 固定 ----
fprintf('\n===== 固定 vs 自适应 对比 =====\n');
fprintf('固定FIR: 抑制=%.1f dB\n', 20*log10(rms_d/rms_fixed));
fprintf('自适应:  抑制=%.1f dB\n', supp_adapt);
fprintf('调试完成。\n');
