clear; close all; clc;
yalmip('clear');

% ------- User params --------
z = tf('z', -1);
G0 = (-0.00146*(z-0.1438)*(z-1)) / ((z - 0.7096)*(z^2 - 0.04369*z + 0.01392));
N = 7;                        % FIR 阶数

w_interp = [0.0521]; % 干扰频率 (rad)   (25 rad/s * Ts, Ts=1/480 -> 0.0521)
nf = length(w_interp);
z_interp = exp(1j * w_interp);

% ---------- 频率网格 ----------
Nfreq = 2000;  % 减少频率点，提高数值稳定性
w_grid = linspace(0, pi, Nfreq);
z_grid = exp(1j*w_grid);

% ---------- G0 的频率响应 ----------
G_grid = squeeze(freqresp(G0, w_grid)).';      % Nfreq x 1
G_interp = squeeze(freqresp(G0, w_interp)).'; % nf x 1

% ---------- Phi 矩阵 ----------
Phi_grid = zeros(Nfreq, N);
Phi_interp = zeros(nf, N);
for k = 1:N
    Phi_grid(:,k)   = z_grid .^ (-(k-1));
    Phi_interp(:,k) = z_interp .^ (-(k-1));
end

% ---------- 插值约束 A * theta = b ----------
A = diag(G_interp) * Phi_interp; % nf x N (complex)
b = ones(size(G_interp));        % 1's (complex)

Aeq = [ real(A); imag(A) ];
beq = [ real(b); imag(b) ];

% ---------- YALMIP 变量 ----------
theta = sdpvar(N,1,'full');   % FIR 系数 (实数)
gamma = sdpvar(1,1);          % H∞ upper bound

% ---------- 约束 ----------
Constraints = [];
% 插值约束
Constraints = [Constraints, Aeq*theta == beq];

% H∞ 约束: |G0(z_k)F(z_k)| <= gamma
for k = 1:Nfreq
    % F(z_k,theta) = Phi_grid(k,:)*theta
    Fk = Phi_grid(k,:) * theta;
    Constraints = [Constraints, norm([real(G_grid(k)*Fk), imag(G_grid(k)*Fk)],2) <= gamma];
end

% ---------- 目标函数 ----------
Objective = gamma;

% ---------- 求解 ----------
options = sdpsettings('solver','mosek','verbose',1);  
sol = optimize(Constraints, Objective, options);

if sol.problem == 0
    fprintf('优化成功!\n');
    theta_opt = value(theta);
    gamma_opt = value(gamma);
    fprintf('最优 gamma = %.4f\n', gamma_opt);
else
    disp('优化失败:');
    sol.info
    % 如果优化失败，使用最小二乘解作为备用方案
    theta_opt = Aeq \ beq;  % 最小二乘解
    gamma_opt = NaN;
    fprintf('使用最小二乘备用解\n');
end

% ---------- 验证 ----------
Fz = zeros(Nfreq,1);
for k = 1:Nfreq
    Fz(k) = Phi_grid(k,:) * theta_opt;
end
H = G_grid(:) .* Fz;

figure;
plot(w_grid, 20*log10(abs(H)), 'LineWidth',1.5);
xlabel('\omega (rad/sample)'); ylabel('|G0F| (dB)');
title('频率响应 |G0F|');
grid on;

