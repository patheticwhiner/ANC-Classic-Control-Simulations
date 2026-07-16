function [u, u_demand, clipped, diag_out] = anc_demo14_youla_rls_step( ...
        e, enable, A_poly, B_star, R0, S0, P0, nQ, lambda1, lambda2, ...
        F_diag, s0_norm, control_scale, u_limit) %#codegen
%ANC_DEMO14_YOULA_RLS_STEP RST + Youla-Q RLS 自适应步函数（Demo1/4 共用）
%
%   移植自 tests/internal/controller_demo1.m::sim_adaptive_landau
%   （controller_demo4 同名引擎差异仅两行:
%      demo4: u_request = u_request / So(1); u_request = control_scale * u_request;
%    demo1 传 s0_norm=1, control_scale=1 —— IEEE754 除以/乘以 1.0 为恒等，
%    两者位级一致）。运行时 HR=HS=1（demo1 低阶设计与 demo4 eMOPSO 设计
%    均把固定因子并入 R0/S0），故 conv(HR_HS, Q) ≡ Q, conv(HR_HS, B*) ≡ B*。
%
%   被控对象递推留在 Simulink；本函数为脚本引擎的第 2–7 步:
%     2. 扰动观测   w  = A·y − B*·u
%     3. 滤波扰动   w1 = (S0/P0)·w
%     4/5. RLS      e_post 后验误差, Q ← Q + F·Φ'·e_post,
%                   F ← inv(λ1·inv(F) + λ2·Φ'Φ + 1e-8·I)
%     6. 控制律     u = (−R0·Y − Q·W − S0(2:end)·U) / s0_norm · control_scale
%     7. 回归向量   w2 = (B*/P0)·w, Φ 滚动
%
%   参数（Parameter 作用域，编译期常量）: A_poly, B_star, R0, S0, P0,
%     nQ, lambda1, lambda2, F_diag, s0_norm, control_scale, u_limit
%   诊断输出 diag_out = [‖Q‖; e_post]。

% ---- 缓冲区长度（与脚本 buf_len 一致，编译期常量）----
buf_len = max([numel(A_poly), numel(B_star), numel(R0), numel(S0), ...
    numel(P0), nQ + 3]) + 5;

persistent U Y W W1 W2 Q_vec Phi F_mat
if isempty(U)
    U  = zeros(1, buf_len);
    Y  = zeros(1, buf_len);
    W  = zeros(1, buf_len);
    W1 = zeros(1, buf_len);
    W2 = zeros(1, buf_len);
    Q_vec = zeros(1, nQ + 1);
    Phi   = zeros(1, nQ + 1);
    F_mat = F_diag * eye(nQ + 1);
end

if enable <= 0.5
    U(:) = 0;  Y(:) = 0;  W(:) = 0;  W1(:) = 0;  W2(:) = 0;
    Q_vec = zeros(1, nQ + 1);
    Phi   = zeros(1, nQ + 1);
    F_mat = F_diag * eye(nQ + 1);
    u = 0;  u_demand = 0;  clipped = false;  diag_out = [0; 0];
    return;
end

% ---- 2. 扰动观测: w = A*y - B_star*u ----
Y = [e, Y(1:end-1)];
w_new = fir_dot(A_poly, Y) - fir_dot(B_star, U);
W = [w_new, W(1:end-1)];

% ---- 3. 滤波扰动: w1 = (S0/P0)*w ----
w1_new = (fir_dot(S0, W) - fir_dot(P0(2:end), W1)) / P0(1);
W1 = [w1_new, W1(1:end-1)];

% ---- 4. RLS 先验/后验误差 ----
e0 = w1_new - fir_dot(Q_vec, Phi);
e_post = e0 / (1 + Phi * F_mat * Phi');

% ---- 5. RLS 更新 ----
Q_vec = Q_vec + (F_mat * Phi' * e_post)';

F_new = lambda1 * inv(F_mat) + lambda2 * (Phi' * Phi); %#ok<MINV>
F_new = F_new + 1e-8 * eye(nQ + 1);
F_mat = inv(F_new); %#ok<MINV>

if any(isnan(F_mat(:))) || any(isinf(F_mat(:))) || rcond(F_mat) < 1e-12
    Q_vec = zeros(1, nQ + 1);
    F_mat = F_diag * eye(nQ + 1);
end

% ---- 6. 控制律（HR_HS_Q ≡ Q_vec，见头注）----
u_request = -fir_dot(R0, Y) - fir_dot(Q_vec, W) - fir_dot(S0(2:end), U);
u_request = u_request / s0_norm;
u_request = control_scale * u_request;
[u_new, was_clipped] = apply_limit(u_request, u_limit);
U = [u_new, U(1:end-1)];

% ---- 7. 回归向量: w2 = (B_star/P0)*w ----
w2_new = (fir_dot(B_star, W) - fir_dot(P0(2:end), W2)) / P0(1);
W2 = [w2_new, W2(1:end-1)];
Phi = [w2_new, Phi(1:end-1)];

u = u_new;
u_demand = u_request;
clipped = was_clipped;
diag_out = [norm(Q_vec); e_post];
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
