function [u, u_demand, clipped, diag_out] = anc_demo2_lqg_adaptive_step( ...
        e, enable, A_design, B_design, C_design, K_lqr, L_kal, ...
        A_t12, B_t12, C_t12, bp_b, bp_a, nQ, lambda, F_init, ...
        adaptation_gain, method_is_lms, lms_epsilon, lms_leakage, ...
        rls_regularization, theta_norm_limit, warmup_samples, ...
        ramp_samples, u_limit) %#codegen
%ANC_DEMO2_LQG_ADAPTIVE_STEP 增广 LQG + filtered-x Youla-Q (RLS/NLMS) 步函数
%
%   移植自 tests/internal/controller_demo2.m::sim_lqg_augmented_adaptive
%   （被控对象状态递推留在 Simulink）。每样本:
%     r      = y − C_design·x̂            残差
%     r_filt = butter 带通 (DF2T，复现 MATLAB filter() 状态式)
%     x_f    = (C_t12 · X_t12)'           滤波回归向量
%     y_Q    = θ'·q_basis                 Q 附加输出
%     u      = ramp·(−K·x̂ + y_Q)
%     θ 更新: NLMS 或 RLS（k ≥ warmup 且未削顶），范数投影
%     x̂/X_t12 状态推进
%
%   诊断输出 diag_out = [‖θ‖; y_Q]。

n_design = size(A_design, 1);
n_bp = max(numel(bp_a), numel(bp_b)) - 1;

persistent x_hat theta P_mat q_basis t12_states filter_state k_count
if isempty(x_hat)
    x_hat = zeros(n_design, 1);
    theta = zeros(nQ, 1);
    P_mat = F_init * eye(nQ);
    q_basis = zeros(nQ, 1);
    t12_states = zeros(size(A_t12, 1), nQ);
    filter_state = zeros(n_bp, 1);
    k_count = 0;
end

if enable <= 0.5
    x_hat = zeros(n_design, 1);
    theta = zeros(nQ, 1);
    P_mat = F_init * eye(nQ);
    q_basis = zeros(nQ, 1);
    t12_states = zeros(size(A_t12, 1), nQ);
    filter_state = zeros(n_bp, 1);
    k_count = 0;
    u = 0;  u_demand = 0;  clipped = false;  diag_out = [0; 0];
    return;
end

k_count = k_count + 1;
y_k = e;

% ---- 残差带通（DF2T，与 MATLAB filter(b,a,x,zi) 状态式逐位一致）----
r_k = y_k - C_design * x_hat;
[r_filt, filter_state] = df2t_step(bp_b, bp_a, r_k, filter_state);
q_basis = [r_filt; q_basis(1:end-1)];

% ---- 控制律 ----
filtered_x = (C_t12 * t12_states).';
y_Q_k = theta' * q_basis;
ramp = min(1, max(0, (k_count - 1) / max(1, ramp_samples)));
u_request = ramp * (-K_lqr * x_hat + y_Q_k);
[u_k, was_clipped] = apply_limit(u_request, u_limit);

% ---- 自适应更新 ----
if k_count >= warmup_samples && ~was_clipped
    if method_is_lms > 0.5
        phi_energy = filtered_x' * filtered_x;
        lms_step = adaptation_gain / max(lms_epsilon + phi_energy, 1e-12);
        theta = (1 - lms_leakage) * theta - lms_step * filtered_x * y_k;
    else
        denom = lambda + filtered_x' * P_mat * filtered_x;
        gain = P_mat * filtered_x / max(denom, 1e-10);
        theta = theta - adaptation_gain * gain * y_k;
        P_mat = (P_mat - gain * filtered_x' * P_mat) / lambda;
    end

    theta_norm = norm(theta);
    if theta_norm > theta_norm_limit
        theta = theta * (theta_norm_limit / theta_norm);
    end
end
if method_is_lms <= 0.5
    P_mat = 0.5 * (P_mat + P_mat') + rls_regularization * eye(nQ);
end
if any(isnan(theta)) || any(isinf(theta))
    theta = zeros(nQ, 1);
    P_mat = F_init * eye(nQ);
end

% ---- 状态推进 ----
x_hat = (A_design - L_kal * C_design) * x_hat + B_design * u_k + L_kal * y_k;
t12_states = A_t12 * t12_states + B_t12 * q_basis.';

u = u_k;
u_demand = u_request;
clipped = was_clipped;
diag_out = [norm(theta); y_Q_k];
end

function [y, z] = df2t_step(b, a, x, z)
% Direct Form II Transposed 单样本步进，复现 MATLAB filter(b,a,x,zi):
%   y    = b(1)·x + z(1)
%   z(i) = b(i+1)·x + z(i+1) − a(i+1)·y   (a(1)=1，butter 输出已归一)
n = numel(z);
y = b(1) * x + z(1);
for i = 1:n-1
    z(i) = b(i+1) * x + z(i+1) - a(i+1) * y;
end
z(n) = b(n+1) * x - a(n+1) * y;
end

function [u_actual, clipped] = apply_limit(u_demand, limit)
% 归一化执行器硬限幅（与 tests/internal 的 apply_actuator_limit 一致）
u_actual = min(max(u_demand, -limit), limit);
clipped = abs(u_actual - u_demand) > 1e-10;
end
