function [u, u_demand, clipped, diag_out] = anc_demo3_imc_fxnlms_step( ...
        e, ref, enable, theta_init, mu_nlms, Ahat_poly, Bhat_star, ...
        k_start, ramp_samples, control_scale, theta_norm_limit, ...
        u_limit) %#codegen
%ANC_DEMO3_IMC_FXNLMS_STEP 源参考 filtered-x NLMS 步函数（Demo3 自适应）
%
%   移植自 tests/internal/controller_demo3.m::sim_imc_fxnlms。
%   前馈结构 —— 需要参考信号 ref（仿真中 = 扰动源 d 本身；台架上为
%   参考传感器，需先验证与误差通道的相干性）。每样本 (k ≥ k_start):
%     x_f(k)  = Ŝ(z)·ref             滤波参考（次级路径模型递推）
%     若上一样本未削顶: θ ← θ + μ·e·x_f_hist / (x_f_hist'x_f_hist + 0.01)
%                        + 范数投影
%     u = −control_scale·ramp·(θ'·ref_hist)
%
%   与脚本一致: 更新/输出只用 ref 与 x_f 的历史（不含当前样本）；
%   ref 历史从 k=1 就开始累积，x_f 递推从 k_start 开始（此前保持 0）。
%
%   诊断输出 diag_out = [‖θ‖; x_f(k)]。

N_adapt = numel(theta_init);
nA = numel(Ahat_poly) - 1;
Lb = numel(Bhat_star);

persistent ref_hist xf_hist theta clipped_prev k_count
if isempty(ref_hist)
    ref_hist = zeros(1, max(N_adapt, Lb));
    xf_hist = zeros(1, max(N_adapt, nA));
    theta = theta_init(:);
    clipped_prev = false;
    k_count = 0;
end

if enable <= 0.5
    ref_hist(:) = 0;
    xf_hist(:) = 0;
    theta = theta_init(:);
    clipped_prev = false;
    k_count = 0;
    u = 0;  u_demand = 0;  clipped = false;  diag_out = [0; 0];
    return;
end

k_count = k_count + 1;

if k_count < k_start
    % 脚本循环从 kStart 才开始；此前仅累积参考历史
    ref_hist = [ref, ref_hist(1:end-1)];
    u = 0;  u_demand = 0;  clipped = false;
    diag_out = [norm(theta); 0];
    return;
end

% ---- 滤波参考: x_f = Ŝ(z)·ref（只用历史，与脚本索引一致）----
xf_new = -fir_dot(Ahat_poly(2:end), xf_hist) + fir_dot(Bhat_star, ref_hist);

% ---- NLMS 更新（用历史向量，不含当前样本）----
referenceVector = ref_hist(1:N_adapt)';
filteredVector = xf_hist(1:N_adapt)';
if ~clipped_prev
    theta = theta + mu_nlms * e * filteredVector ...
        / (filteredVector' * filteredVector + 0.01);
    theta_norm = norm(theta);
    if isfinite(theta_norm_limit) && theta_norm > theta_norm_limit
        theta = theta * (theta_norm_limit / theta_norm);
    end
end

% ---- 控制律 ----
ramp = min(1, max(0, (k_count - k_start) / max(1, ramp_samples)));
u_request = -control_scale * ramp * (theta' * referenceVector);
if ~isfinite(u_request)
    u_request = 0;
end
[u_new, was_clipped] = apply_limit(u_request, u_limit);
clipped_prev = was_clipped;

% ---- 历史推进 ----
ref_hist = [ref, ref_hist(1:end-1)];
xf_hist = [xf_new, xf_hist(1:end-1)];

u = u_new;
u_demand = u_request;
clipped = was_clipped;
diag_out = [norm(theta); xf_new];
end

function y = fir_dot(b, x)
% FIR 点积滤波: y = Σ b(i)·x(i)（与 tests/internal 的 FIR() 一致）
L = min(length(b), length(x));
if L < 1
    y = 0;
else
    y = b(1:L) * x(1:L)';
end
end

function [u_actual, clipped] = apply_limit(u_demand, limit)
% 归一化执行器硬限幅（与 tests/internal 的 apply_actuator_limit 一致）
u_actual = min(max(u_demand, -limit), limit);
clipped = abs(u_actual - u_demand) > 1e-10;
end
