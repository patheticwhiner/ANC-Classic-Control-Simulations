function [u, u_demand, clipped, diag_out] = anc_demo2_lqg_step( ...
        e, enable, A_design, B_design, C_design, K_lqr, L_kal, ...
        ramp_samples, u_limit) %#codegen
%ANC_DEMO2_LQG_STEP 增广 LQG 固定控制器步函数
%
%   移植自 tests/internal/controller_demo2.m::sim_lqg_augmented
%   （被控对象状态递推留在 Simulink 的 plant_secondary_ss 块中）:
%     ramp      = min(1, max(0, (k-1)/max(1, ramp_samples)))
%     u_request = -ramp * K * x_hat
%     x_hat     = (A_design - L*C_design)*x_hat + B_design*u + L*y
%
%   诊断输出 diag_out = [‖x_hat‖; 0]。

persistent x_hat k_count
if isempty(x_hat)
    x_hat = zeros(size(A_design, 1), 1);
    k_count = 0;
end

if enable <= 0.5
    x_hat = zeros(size(A_design, 1), 1);
    k_count = 0;
    u = 0;  u_demand = 0;  clipped = false;  diag_out = [0; 0];
    return;
end

k_count = k_count + 1;

% ---- 控制律（与脚本逐字一致，含求值顺序）----
ramp = min(1, max(0, (k_count - 1) / max(1, ramp_samples)));
u_request = -ramp * K_lqr * x_hat;
[u_new, was_clipped] = apply_limit(u_request, u_limit);

% ---- Kalman 更新 ----
x_hat = (A_design - L_kal * C_design) * x_hat ...
    + B_design * u_new + L_kal * e;

u = u_new;
u_demand = u_request;
clipped = was_clipped;
diag_out = [norm(x_hat); 0];
end

function [u_actual, clipped] = apply_limit(u_demand, limit)
% 归一化执行器硬限幅（与 tests/internal 的 apply_actuator_limit 一致）
u_actual = min(max(u_demand, -limit), limit);
clipped = abs(u_actual - u_demand) > 1e-10;
end
