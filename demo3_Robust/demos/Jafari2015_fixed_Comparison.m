% Direct Hinf minimization with interpolation equalities (CVX)
% Inputs: G0 (tf with z variable), N_theta (N), w_interp (vector), r0 (optional)
clear; close all;clc;
% ------- User params --------
z = tf('z', -1);
G0 = (-0.00146*(z-0.1438)*(z-1)) / ((z - 0.7096)*(z^2 - 0.04369*z + 0.01392));
N = 30;                        % number of FIR coeffs

w_interp = [0.0521]; % interpolation frequencies (rad)
nf = length(w_interp);
z_interp = exp(1j * w_interp);


% frequency grid for Hinf constraints
Nfreq = 2000;                     
w_grid = linspace(0, pi, Nfreq);
z_grid = exp(1j*w_grid);
% precompute frequency responses
G_grid = squeeze(freqresp(G0, w_grid)).';        % 1 x Nfreq complex -> transpose to Nfreq x 1
G_interp = squeeze(freqresp(G0, w_interp)).';   % N x 1 complex
% build Phi matrices
Phi_grid = zeros(Nfreq, N);
Phi_interp = zeros(nf, N);
for k = 1:N
    Phi_grid(:,k)   = z_grid .^ (k-N-1);
    Phi_interp(:,k) = z_interp .^ (k-N-1); % equivalent
end

% build interpolation A * theta = b (complex)
A = diag(G_interp) * Phi_interp; % size N x N (complex)
b = ones(size(G_interp));        % complex ones

% if theta real, convert to real linear eqns
Aeq = [ real(A); imag(A) ];
beq = [ real(b)'; imag(b)' ];

% ---------- 可行性检查（先看是否存在满足等式的 theta） ----------
feasible_flag = true;
% 线性可行性：是否有 theta s.t. Aeq * theta = beq
% 简单尝试解方程并测残差
if size(Aeq,1) == size(Aeq,2)
    theta_test = Aeq \ beq;
    resnorm = norm(Aeq*theta_test - beq);
else
    theta_test = Aeq \ beq; % least squares if over/under
    resnorm = norm(Aeq*theta_test - beq);
end
fprintf('线性可行性检验: 等式最小残差 = %.3e\n', resnorm);
if resnorm > 1e-8
    warning('插值等式可能无法严格满足 (残差 %.3e). 若要强制满足，可考虑增大 N 或放宽约束。', resnorm);
    % 依使用场景，可视为不可行 or 允许松弛
end

% ---------- 最小二乘基准解（仅供对比） ----------
fprintf('\n--- 最小二乘基准解分析 ---\n');
theta_ls = theta_test;  % 使用前面计算的最小二乘解
fprintf('LS解的插值残差: %.3e\n', resnorm);

% 计算最小二乘解的H∞性能
mag_ls = abs(G_grid(:) .* (Phi_grid * theta_ls));
gamma_ls = max(mag_ls);
fprintf('LS解的H∞范数: %.6f\n', gamma_ls);
fprintf('LS解的参数范数: %.6f\n', norm(theta_ls));

% ---------- 多种求解器对比 ----------
solvers = {}; % 存储所有求解结果
solver_names = {}; % 存储求解器名称

% 1. 最小二乘解（基准）
solvers{end+1} = struct('name', 'LS', 'theta', theta_ls, 'gamma', gamma_ls, ...
                       'time', 0, 'status', 'success');
solver_names{end+1} = 'LS';

% 2. CVX求解
if exist('cvx_begin','file') == 2
    fprintf('\n--- CVX 求解 ---\n');
    tic;
    r0 = 1e4;
    cvx_begin quiet
        cvx_precision high
        variable theta_cvx(N,1)
        variable gamma_cvx nonnegative
        minimize( gamma_cvx )
        subject to
            for i = 1:Nfreq
                abs( G_grid(i) * (Phi_grid(i,:) * theta_cvx) ) <= gamma_cvx;
            end
            Aeq * theta_cvx == beq;
            norm(theta_cvx,2) <= r0;
    cvx_end
    time_cvx = toc;

    if strcmp(cvx_status, 'Solved')
        solvers{end+1} = struct('name', 'CVX', 'theta', theta_cvx, 'gamma', gamma_cvx, ...
                               'time', time_cvx, 'status', 'success');
        solver_names{end+1} = 'CVX';
        fprintf('CVX 求解成功: γ = %.6f, 用时 %.2fs\n', gamma_cvx, time_cvx);
    else
        fprintf('CVX 求解失败: %s\n', cvx_status);
    end
else
    fprintf('CVX 未安装，跳过CVX求解\n');
end

% 3. YALMIP + MOSEK 求解
if exist('yalmip','file') == 2
    fprintf('\n--- YALMIP + MOSEK 求解 ---\n');
    yalmip('clear');
    tic;
    
    % YALMIP 变量
    theta_mosek = sdpvar(N,1,'full');
    gamma_mosek = sdpvar(1,1);
    
    % 约束
    Constraints = [];
    Constraints = [Constraints, Aeq*theta_mosek == beq];
    for k = 1:Nfreq
        Fk = Phi_grid(k,:) * theta_mosek;
        Constraints = [Constraints, norm([real(G_grid(k)*Fk), imag(G_grid(k)*Fk)],2) <= gamma_mosek];
    end
    
    % 求解
    options = sdpsettings('solver','mosek','verbose',0);
    sol = optimize(Constraints, gamma_mosek, options);
    time_mosek = toc;
    
    if sol.problem == 0
        theta_mosek_val = value(theta_mosek);
        gamma_mosek_val = value(gamma_mosek);
        solvers{end+1} = struct('name', 'MOSEK', 'theta', theta_mosek_val, 'gamma', gamma_mosek_val, ...
                               'time', time_mosek, 'status', 'success');
        solver_names{end+1} = 'MOSEK';
        fprintf('MOSEK 求解成功: γ = %.6f, 用时 %.2fs\n', gamma_mosek_val, time_mosek);
    else
        fprintf('MOSEK 求解失败: problem = %d\n', sol.problem);
    end
else
    fprintf('YALMIP 未安装，跳过MOSEK求解\n');
end

% 4. YALMIP + SeDuMi 求解（备用）
if exist('yalmip','file') == 2 && length(solvers) == 1 % 如果只有LS解，尝试SeDuMi
    fprintf('\n--- YALMIP + SeDuMi 求解 ---\n');
    yalmip('clear');
    tic;
    
    theta_sedumi = sdpvar(N,1,'full');
    gamma_sedumi = sdpvar(1,1);
    
    Constraints = [];
    Constraints = [Constraints, Aeq*theta_sedumi == beq];
    % 减少约束数量以提高稳定性
    Nfreq_reduced = 200;
    w_reduced = linspace(0, pi, Nfreq_reduced);
    G_reduced = squeeze(freqresp(G0, w_reduced)).';
    Phi_reduced = zeros(Nfreq_reduced, N);
    for k = 1:N
        Phi_reduced(:,k) = exp(-1j*w_reduced*(k-1));
    end
    
    for k = 1:Nfreq_reduced
        Fk = Phi_reduced(k,:) * theta_sedumi;
        Constraints = [Constraints, norm([real(G_reduced(k)*Fk), imag(G_reduced(k)*Fk)],2) <= gamma_sedumi];
    end
    
    options = sdpsettings('solver','sedumi','verbose',0);
    sol = optimize(Constraints, gamma_sedumi, options);
    time_sedumi = toc;
    
    if sol.problem == 0
        theta_sedumi_val = value(theta_sedumi);
        gamma_sedumi_val = value(gamma_sedumi);
        solvers{end+1} = struct('name', 'SeDuMi', 'theta', theta_sedumi_val, 'gamma', gamma_sedumi_val, ...
                               'time', time_sedumi, 'status', 'success');
        solver_names{end+1} = 'SeDuMi';
        fprintf('SeDuMi 求解成功: γ = %.6f, 用时 %.2fs\n', gamma_sedumi_val, time_sedumi);
    else
        fprintf('SeDuMi 求解失败: problem = %d\n', sol.problem);
    end
end

% ---------- 多求解器对比分析 ----------
compare_multiple_solvers(solvers, Aeq, beq);

function compare_multiple_solvers(solvers, Aeq, beq)
    fprintf('\n--- 多求解器性能对比 ---\n');
    fprintf('%-10s | %-10s | %-10s | %-10s | %-8s\n', '求解器', 'H∞范数', '参数范数', '插值残差', '用时(s)');
    fprintf('-----------|------------|------------|------------|----------\n');
    
    % 找到最优解
    gammas = [];
    for i = 1:length(solvers)
        if strcmp(solvers{i}.status, 'success')
            gammas(end+1) = solvers{i}.gamma;
        end
    end
    best_gamma = min(gammas);
    
    % 显示结果
    for i = 1:length(solvers)
        solver = solvers{i};
        residual = norm(Aeq*solver.theta - beq);
        
        % 标记最优解
        if abs(solver.gamma - best_gamma) < 1e-6
            marker = ' ★';
        else
            marker = '  ';
        end
        
        fprintf('%-10s | %-10.6f | %-10.6f | %-10.3e | %-8.2f%s\n', ...
                solver.name, solver.gamma, norm(solver.theta), residual, solver.time, marker);
    end
    
    % 性能改善分析
    if length(solvers) > 1
        ls_gamma = solvers{1}.gamma; % 假设第一个是LS
        fprintf('\n--- 性能改善分析 ---\n');
        for i = 2:length(solvers)
            improvement = ls_gamma / solvers{i}.gamma;
            fprintf('%s 相对LS改善: %.2f倍', solvers{i}.name, improvement);
            if improvement > 1.2
                fprintf(' ✓ 显著改善\n');
            elseif improvement > 1.05
                fprintf(' △ 略有改善\n');
            else
                fprintf(' ○ 性能接近\n');
            end
        end
    end
end

% ---------- 验证所有求解器结果 ----------
fprintf('\n--- 密集网格验证 ---\n');
Mchk = 3000;
w_chk = linspace(0, pi, Mchk);
G_chk = squeeze(freqresp(G0, w_chk)).';
Phi_chk = zeros(Mchk, N);
for k=1:N, Phi_chk(:,k) = exp(-1j*w_chk*(k-1)); end

for i = 1:length(solvers)
    solver = solvers{i};
    mag = abs( G_chk(:) .* (Phi_chk * solver.theta) );
    gamma_check = max(mag);
    fprintf('%s: 理论γ=%.6f, 验证γ=%.6f, 误差=%.2e\n', ...
            solver.name, solver.gamma, gamma_check, abs(solver.gamma - gamma_check));
end

% 选择最优解用于后续仿真
[~, best_idx] = min([solvers{1}.gamma,solvers{2}.gamma,solvers{3}.gamma]);
theta_opt = solvers{best_idx}.theta;
gamma_opt = solvers{best_idx}.gamma;
best_solver_name = solvers{best_idx}.name;
fprintf('\n选择 %s 解进行后续仿真 (γ = %.6f)\n', best_solver_name, gamma_opt);


%%
% 扰动参数
Ts = 1/480;
fs = 1/Ts;
omega0 = 0.0521; % rad/sample, 25rad/sec, 3.98Hz
Nsim = 480*50;   % 50秒
t = (0:Nsim-1)*Ts;
d = sin(omega0*t*fs) + 0.02*randn(1,Nsim);
figure('Name', '扰动信号');
plot(t, d); grid on;
xlabel('Time(s)'); ylabel('Magnitude');
xlim([0,5]);

% ---------- 多求解器频率响应对比 ----------
nfft = 1024;
colors = {'r-', 'b-', 'g-', 'm-', 'c-', 'k-'}; % 支持最多6个求解器
styles = {'-', '.-', '--', ':', '-', '--'}; % 不同线型

figure('Name', '多求解器控制器频率响应对比');
subplot(2,1,1);
for i = 1:length(solvers)
    solver = solvers{i};
    b_K = [0; flip(solver.theta(:))];
    [Hk, fk] = freqz(b_K, 1, nfft, fs);
    
    color_idx = mod(i-1, length(colors)) + 1;
    style_idx = mod(i-1, length(styles)) + 1;
    plot_style = [colors{color_idx}(1) styles{style_idx}];
    
    semilogx(fk, 20*log10(abs(Hk)), plot_style, 'LineWidth', 2, ...
             'DisplayName', sprintf('%s (γ=%.3f)', solver.name, solver.gamma));
    hold on;
end
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)'); grid on;
title('控制器幅频响应对比');
legend('Location', 'best');

subplot(2,1,2);
for i = 1:length(solvers)
    solver = solvers{i};
    b_K = [0; flip(solver.theta(:))];
    [Hk, ~] = freqz(b_K, 1, nfft, fs);
    
    color_idx = mod(i-1, length(colors)) + 1;
    style_idx = mod(i-1, length(styles)) + 1;
    plot_style = [colors{color_idx}(1) styles{style_idx}];
    
    semilogx(fk, angle(Hk)*180/pi, plot_style, 'LineWidth', 2, ...
             'DisplayName', sprintf('%s', solver.name));
    hold on;
end
xlabel('Frequency (Hz)'); ylabel('Phase (deg)'); grid on;
title('控制器相频响应对比');
legend('Location', 'best');

% ---------- 多求解器闭环仿真对比 ----------
% 获取G0(z)分子分母系数
b_plant = G0.num{1};
a_plant = G0.den{1};

% 为每个求解器准备仿真变量
sim_results = {};
for i = 1:length(solvers)
    solver = solvers{i};
    
    % 初始化仿真变量
    u = zeros(1,Nsim);
    y = zeros(1,Nsim);
    x = zeros(1,Nsim);
    antinoise = zeros(1,Nsim);
    
    % 仿真循环
    for k = max(N,length(a_plant))+1:Nsim
        if t(k) >= 20  % 20秒后开启控制器
            u(k) = -sum(solver.theta'.*x(k-(N:-1:1)));
            if isnan(u(k))
                u(k) = 0; % 防止NaN
            end
        end
        % 离散系统差分方程递推
        antinoise(k) = (-a_plant(2:end)*antinoise(k-1:-1:k-length(a_plant)+1)' + ...
                        b_plant*u(k-1:-1:k-length(b_plant))')/a_plant(1);
        y(k) = antinoise(k) + d(k);
        x(k) = y(k) - antinoise(k);
    end
    
    % 存储结果
    sim_results{i} = struct('name', solver.name, 'gamma', solver.gamma, ...
                           'u', u, 'y', y, 'x', x);
end


% ---------- 多求解器时域仿真结果对比 ----------
% 控制输入对比
figure('Name', '多求解器控制输入信号对比');
for i = 1:length(sim_results)
    result = sim_results{i};
    color_idx = mod(i-1, length(colors)) + 1;
    style_idx = mod(i-1, length(styles)) + 1;
    plot_style = [colors{color_idx}(1) styles{style_idx}];
    
    plot(t, result.u, plot_style, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s控制输入', result.name));
    hold on;
end
xline(20, 'k:', 'Controller ON');
xlabel('Time (s)'); ylabel('Control Input u');
title('多求解器控制输入信号对比');
legend show; grid on;

% 闭环输出对比
figure('Name', '多求解器闭环输出性能对比');
plot(t, d, 'k:', 'LineWidth', 1.0, 'DisplayName', 'Disturbance');
hold on;
for i = 1:length(sim_results)
    result = sim_results{i};
    color_idx = mod(i-1, length(colors)) + 1;
    style_idx = mod(i-1, length(styles)) + 1;
    plot_style = [colors{color_idx}(1) styles{style_idx}];
    
    plot(t, result.y, plot_style, 'LineWidth', 1.5, ...
         'DisplayName', sprintf('%s (γ=%.3f)', result.name, result.gamma));
end
xline(20, 'k--', 'Controller ON');
xlabel('Time (s)'); ylabel('Plant output y');
title('多求解器闭环输出性能对比');
legend show; grid on;

% 性能统计（可扩展）
fprintf('\n--- 多求解器时域仿真性能统计 ---\n');
rms_d = rms(d(20*fs:end));  % 扰动RMS
fprintf('扰动信号RMS: %.4f\n', rms_d);
fprintf('%-10s | %-12s | %-12s\n', '求解器', '控制后RMS', '噪声抑制(dB)');
fprintf('-----------|--------------|-------------\n');

for i = 1:length(sim_results)
    result = sim_results{i};
    rms_controlled = rms(result.y(20*fs:end));
    suppression_db = 20*log10(rms_d/rms_controlled);
    fprintf('%-10s | %-12.4f | %-12.1f\n', result.name, rms_controlled, suppression_db);
end
