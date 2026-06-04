%% 参数调优：θ_max 与 κ₀ 对抑制效果和控制能量的影响
clear; close all; clc;
addpath('../../functions');

%% ====== 系统定义 ======
s = tf('s');
fs = 10000; Ts = 1/fs;
G0 = 0.5*(s-0.2)/(s^2+s+1.25);
F_bar = 2*500^2*(s^2+s+1.25)/((s+500)^2*(s+0.2));  % F = κ₀·F̄

% 扰动
Tsim = 2; Nsim = Tsim*fs; t = (0:Nsim-1)*Ts;
omega_d = [70, 187]; A = [0.6, 0.7]; phi0 = [pi/4, pi/2];
d = zeros(1, Nsim);
for k = 1:length(omega_d)
    d = d + A(k)*sin(omega_d(k)*t + phi0(k));
end
d = d + 0.02*randn(1, Nsim);

%% ====== 参数扫描 ======
kappa_list = [0.05, 0.1, 0.2, 0.5];
theta_max_list = [3, 5, 10, 20];
Nparam = 20; lambda_val = 500;

results = struct('kappa',{}, 'theta_max',{}, 'suppression',{}, 'rms_u',{}, 'rms_y',{}, 'theta_norm',{}, 'trP',{}, 'robust_LHS',{});
result_idx = 1;

for ki = 1:length(kappa_list)
for ti = 1:length(theta_max_list)
    kappa0 = kappa_list(ki);
    theta_max_val = theta_max_list(ti);

    F = kappa0 * F_bar;
    Fd = c2d(F, Ts, 'tustin');
    [numFd, denFd] = tfdata(Fd, 'v');
    Gd = c2d(G0, Ts, 'tustin');
    [numGd, denGd] = tfdata(Gd, 'v');
    G0d = c2d(G0, Ts, 'tustin');
    [numG0d, denG0d] = tfdata(G0d, 'v');
    G1_s = tf(lambda_val, [1 lambda_val]);
    G1_d = c2d(G1_s, Ts, 'tustin');
    [numG1, denG1] = tfdata(G1_d, 'v');

    Lambda_order = Nparam:-1:1;

    hist_len = max([length(denFd), length(numFd), length(denGd), length(numGd)]) + Nparam + 5;
    G_u_hist = zeros(1, hist_len); G_y_hist = zeros(1, hist_len);
    G0_u_hist = zeros(1, hist_len); G0_y_hist = zeros(1, hist_len);
    F_ctrl_u_hist = zeros(1, hist_len); F_ctrl_y_hist = zeros(1, hist_len);
    phi_Lambda_u = zeros(Nparam, Nparam, hist_len);
    phi_Lambda_y = zeros(Nparam, Nparam, hist_len);
    phi_F_u = zeros(Nparam, hist_len); phi_F_y = zeros(Nparam, hist_len);
    phi_G0_u = zeros(Nparam, hist_len); phi_G0_y = zeros(Nparam, hist_len);

    y = zeros(1, Nsim); u = zeros(1, Nsim);
    theta = zeros(Nparam, 1); P = eye(Nparam)*500;
    gamma0 = 1.0;

    for k = hist_len+1:Nsim
        [y_G_k, G_u_hist, G_y_hist] = IIR_filter(numGd, denGd, u(k-1), G_u_hist, G_y_hist);
        y(k) = y_G_k + d(k);
        [G0_u_k, G0_u_hist, G0_y_hist] = IIR_filter(numG0d, denG0d, u(k-1), G0_u_hist, G0_y_hist);
        z_k = y(k) - G0_u_k;

        phi_reg = zeros(Nparam, 1);
        lambda_vals = zeros(Nparam, 1);
        for i = 1:Nparam
            casc_out = z_k;
            for stage = 1:Lambda_order(i)
                uu = reshape(phi_Lambda_u(i,stage,:), 1, []);
                yy = reshape(phi_Lambda_y(i,stage,:), 1, []);
                [casc_out, uu_n, yy_n] = IIR_filter(numG1, denG1, casc_out, uu, yy);
                phi_Lambda_u(i,stage,:) = uu_n; phi_Lambda_y(i,stage,:) = yy_n;
            end
            lambda_vals(i) = casc_out;
            [f_out, phi_F_u(i,:), phi_F_y(i,:)] = IIR_filter(numFd, denFd, casc_out, phi_F_u(i,:), phi_F_y(i,:));
            [phi_reg(i), phi_G0_u(i,:), phi_G0_y(i,:)] = IIR_filter(numG0d, denG0d, f_out, phi_G0_u(i,:), phi_G0_y(i,:));
        end

        m_s_sq = 1 + gamma0 * (phi_reg.' * phi_reg);
        pred_err = z_k - theta.' * phi_reg;
        epsilon = pred_err / m_s_sq;
        P = P - Ts * (P * (phi_reg * phi_reg.') * P) / m_s_sq;
        theta_new = theta + Ts * P * epsilon * phi_reg;
        if norm(theta_new) > theta_max_val
            theta = theta_new * (theta_max_val / norm(theta_new));
        else
            theta = theta_new;
        end
        K_z_k = theta.' * lambda_vals;
        [u(k), F_ctrl_u_hist, F_ctrl_y_hist] = IIR_filter(numFd, denFd, -K_z_k, F_ctrl_u_hist, F_ctrl_y_hist);
    end

    final_idx = round(0.5*Nsim):Nsim;
    rms_d = sqrt(mean(d(final_idx).^2));
    rms_y = sqrt(mean(y(final_idx).^2));
    rms_u = sqrt(mean(u(final_idx).^2));
    supp_db = 20*log10(rms_y/rms_d);
    theta_final_norm = norm(theta);
    trP_final = trace(P);

    % 鲁棒裕度
    T_rob = F_bar * G0 * (-0.001*s);
    [mag_T, ~] = bode(T_rob, logspace(-1, 4, 500));
    Hinf_T = max(squeeze(mag_T));
    robust_LHS = kappa0 * theta_max_val * Hinf_T;

    results(result_idx) = struct('kappa', kappa0, 'theta_max', theta_max_val, ...
        'suppression', supp_db, 'rms_u', rms_u, 'rms_y', rms_y, ...
        'theta_norm', theta_final_norm, 'trP', trP_final, 'robust_LHS', robust_LHS);

    fprintf('κ₀=%.2f, θ_max=%.0f: supp=%.1f dB, u_rms=%.0f, ||θ||=%.2f, rob_LHS=%.2f\n', ...
        kappa0, theta_max_val, supp_db, rms_u, theta_final_norm, robust_LHS);

    result_idx = result_idx + 1;
end
end

%% ====== 汇总表 ======
fprintf('\n========== 参数扫描汇总 ==========\n');
fprintf('%-8s %-10s %-10s %-10s %-10s %-12s\n', 'κ₀', 'θ_max', '抑制(dB)', 'u_RMS', '||θ||', '鲁棒LHS');
fprintf('%s\n', repmat('-', 1, 65));
for i = 1:length(results)
    r = results(i);
    flag = '';
    if r.robust_LHS < 1, flag = ' ✅'; else, flag = ' ❌'; end
    fprintf('%-8.2f %-10.0f %-10.1f %-10.0f %-10.2f %-12.2f%s\n', ...
        r.kappa, r.theta_max, r.suppression, r.rms_u, r.theta_norm, r.robust_LHS, flag);
end

%% ====== 二维可视化 ======
out_dir = 'E:/Code/MATLAB_ANC/ANC-Classic-Control-Simulations/demo3_Robust/demos/';
nk = length(kappa_list); nt = length(theta_max_list);
supp_mat = reshape([results.suppression], nt, nk)';
u_mat = reshape([results.rms_u], nt, nk)';
rob_mat = reshape([results.robust_LHS], nt, nk)';

fig = figure('Position', [100 100 1400 450]);

subplot(1,3,1);
imagesc(supp_mat); colorbar;
set(gca, 'XTick', 1:nt, 'XTickLabel', theta_max_list);
set(gca, 'YTick', 1:nk, 'YTickLabel', kappa_list);
xlabel('θ_{max}'); ylabel('κ₀');
title('抑制效果 (dB)');
for i = 1:nk, for j = 1:nt
    text(j, i, sprintf('%.0f', supp_mat(i,j)), 'HorizontalAlign','center', 'FontSize', 9, 'Color', 'w');
end, end

subplot(1,3,2);
imagesc(u_mat); colorbar;
set(gca, 'XTick', 1:nt, 'XTickLabel', theta_max_list);
set(gca, 'YTick', 1:nk, 'YTickLabel', kappa_list);
xlabel('θ_{max}'); ylabel('κ₀');
title('控制 RMS');
for i = 1:nk, for j = 1:nt
    text(j, i, sprintf('%.0f', u_mat(i,j)), 'HorizontalAlign','center', 'FontSize', 9, 'Color', 'w');
end, end

subplot(1,3,3);
imagesc(min(rob_mat, 10)); colorbar;
set(gca, 'XTick', 1:nt, 'XTickLabel', theta_max_list);
set(gca, 'YTick', 1:nk, 'YTickLabel', kappa_list);
xlabel('θ_{max}'); ylabel('κ₀');
title('鲁棒 LHS (<1=安全)');
for i = 1:nk, for j = 1:nt
    v = rob_mat(i,j);
    lbl = sprintf('%.1f', v);
    if v < 1, lbl = [lbl '✓']; end
    text(j, i, lbl, 'HorizontalAlign','center', 'FontSize', 8, 'Color', 'w');
end, end

saveas(fig, [out_dir 'param_sweep.png']);

% 推荐配置
fprintf('\n===== 推荐配置 =====\n');
% 找鲁棒安全且抑制最好的 (suppression越负越好 → 用min)
safe_idx = find([results.robust_LHS] < 1);
if ~isempty(safe_idx)
    [best_supp, best_i] = min([results(safe_idx).suppression]);
    best = results(safe_idx(best_i));
    fprintf('✓ 最优安全配置: κ₀=%.2f, θ_max=%.0f → 抑制 %.1f dB, u_RMS=%.0f, rob_LHS=%.2f\n', ...
        best.kappa, best.theta_max, best.suppression, best.rms_u, best.robust_LHS);

    % 找控制能量最低但仍>10dB的安全配置
    mid_safe = safe_idx([results(safe_idx).suppression] < -10);
    if ~isempty(mid_safe)
        [~, mid_i] = min([results(mid_safe).rms_u]);
        mid = results(mid_safe(mid_i));
        fprintf('○ 低控制安全配置:  κ₀=%.2f, θ_max=%.0f → 抑制 %.1f dB, u_RMS=%.0f, rob_LHS=%.2f\n', ...
            mid.kappa, mid.theta_max, mid.suppression, mid.rms_u, mid.robust_LHS);
    end
end

% 对比原始配置
orig = results([results.kappa] == 0.5 & [results.theta_max] == 1000);
fprintf('原始配置:     κ₀=%.2f, θ_max=%.0f → 抑制 %.1f dB, u_RMS=%.0f\n', ...
    orig.kappa, orig.theta_max, orig.suppression, orig.rms_u);

fprintf('\n图片已保存到 %s\n', out_dir);
