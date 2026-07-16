function [u, u_demand, clipped, diag_out] = anc_demo3_fixed_step( ...
        e, enable, theta_fixed, F_num, F_den, Ahat_poly, Bhat_star, ...
        k_start, ramp_samples, control_scale, u_limit) %#codegen
%ANC_DEMO3_FIXED_STEP F(z) 频谱展平 + 固定 FIR 控制器步函数（Jafari）
%
%   移植自 tests/internal/controller_demo3.m::sim_fixed_jafari。
%   脚本中 z(k) = y(k) − g(k)（g 为真实次级路径输出）；部署时控制器
%   用内部次级路径模型 ĝ 重构 z = e − ĝ。离线模型中模型 = 被控对象且
%   两侧算术一致，故 ĝ ≡ g，逐位等效。
%
%   时序与脚本一致:
%     k < k_start          → 输出 0，不推进任何内部状态（z 历史保持 0）
%     k ≥ k_start          → ĝ 递推、z 重构、FIR→F(z)→ramp→限幅
%     （k ≤ N_fir 的 v/F(z) 支路不执行 —— 冻结参数下 k_start > N_fir，
%       保护性保留该守卫）
%
%   诊断输出 diag_out = [‖θ‖(常量); F(z) 输出 uk]。

N_fir = numel(theta_fixed);
nA = numel(Ahat_poly) - 1;
Lb = numel(Bhat_star);
nF_num = numel(F_num);
nF_den = numel(F_den) - 1;

persistent u_hist g_hist z_hist xF yF k_count
if isempty(u_hist)
    u_hist = zeros(1, max(Lb, 1));
    g_hist = zeros(1, max(nA, 1));
    z_hist = zeros(1, max(N_fir, 1));
    xF = zeros(nF_num, 1);
    yF = zeros(max(1, nF_den), 1);
    k_count = 0;
end

if enable <= 0.5
    u_hist(:) = 0;  g_hist(:) = 0;  z_hist(:) = 0;
    xF(:) = 0;  yF(:) = 0;  k_count = 0;
    u = 0;  u_demand = 0;  clipped = false;  diag_out = [0; 0];
    return;
end

k_count = k_count + 1;

if k_count < k_start
    u = 0;  u_demand = 0;  clipped = false;
    diag_out = [norm(theta_fixed); 0];
    return;
end

% ---- 内部次级路径估计: ĝ(k) = b̂*u_hist − â(2:end)*ĝ_hist ----
g_new = fir_dot(Bhat_star, u_hist) - fir_dot(Ahat_poly(2:end), g_hist);
z_new = e - g_new;

u_request = 0;
uk = 0;
was_clipped = false;
if k_count > N_fir
    % FIR 输出（只用 z 的历史，不含当前样本）→ F(z) → ramp
    v = fir_dot(theta_fixed(:).', z_hist);
    xF = [v; xF(1:nF_num-1)];
    uk = F_num * xF - F_den(2:end) * yF(1:nF_den);
    yF = [uk; yF(1:end-1)];
    ramp = min(1, max(0, (k_count - k_start) / max(1, ramp_samples)));
    u_request = -control_scale * ramp * uk;
    if ~isfinite(u_request)
        u_request = 0;
    end
    [u_new, was_clipped] = apply_limit(u_request, u_limit);
else
    u_new = 0;
end

% ---- 历史推进 ----
u_hist = [u_new, u_hist(1:end-1)];
g_hist = [g_new, g_hist(1:end-1)];
z_hist = [z_new, z_hist(1:end-1)];

u = u_new;
u_demand = u_request;
clipped = was_clipped;
diag_out = [norm(theta_fixed); uk];
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
