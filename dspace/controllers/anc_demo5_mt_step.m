function [u, u_demand, clipped, diag_out] = anc_demo5_mt_step( ...
        e, enable, theta_init, theta_min, theta_max, sign_S, sign_S0, ...
        k_gain, epsilon, dc_cancel, k_start, ramp_samples, ...
        u_limit) %#codegen
%ANC_DEMO5_MT_STEP Marino-Tomei (2016) 自适应内模 ANC 步函数
%
%   移植自 tests/internal/controller_demo5.m::sim_marino_tomei，
%   固定为冻结配置 method='exact'、output_timing='updated'
%   （export_frozen_params 对其他配置直接报错）。每样本 (k ≥ k_start):
%     α = √θᵢ, 内模旋转 + 增益注入:
%       w1ᵢ ←  cosα·w1ᵢ + (sinα/α)·w2ᵢ + sgnᵢ·k·(sinα/α)·e
%       w2ᵢ ← −α·sinα·w1ᵢ + cosα·w2ᵢ + sgnᵢ·k·(cosα−1)·e
%     频率自适应 (ε>0): θᵢ ← θᵢ + ε·sgnᵢ·w2ᵢ_old·e, 夹紧 [θmin, θmax]
%     u = −ramp·(w0 + Σ w1ᵢ)
%
%   fixed 变体 = 冻结参数 epsilon=0，共用本函数。
%   限幅为 demo5 变体（NaN → 0 且记削顶）。
%
%   诊断输出 diag_out = [√θ₁ (= ω̂₁·Ts 数字频率); w1₁]。
%   （换算 Hz: f = diag_out(1)·fs/(2π)，可在模型中用 Gain 块实现）

q = numel(theta_init);

persistent w w0 theta k_count
if isempty(w)
    w = zeros(2, q);
    w0 = 0;
    theta = theta_init(:);
    k_count = 0;
end

if enable <= 0.5
    w = zeros(2, q);
    w0 = 0;
    theta = theta_init(:);
    k_count = 0;
    u = 0;  u_demand = 0;  clipped = false;  diag_out = [0; 0];
    return;
end

k_count = k_count + 1;

if k_count >= k_start
    w_old = w;
    theta_old = theta;

    if dc_cancel > 0.5
        w0 = w0 + sign_S0 * k_gain * e;
    end

    for i = 1:q
        % method='exact' 离散化（与脚本逐字一致）
        alpha = sqrt(max(theta_old(i), 1e-12));
        ca = cos(alpha);
        sa = sin(alpha);

        w(1,i) = ca * w_old(1,i) + (sa / alpha) * w_old(2,i) ...
               + sign_S(i) * k_gain * (sa / alpha) * e;
        w(2,i) = -alpha * sa * w_old(1,i) + ca * w_old(2,i) ...
               + sign_S(i) * k_gain * (ca - 1) * e;

        if epsilon > 0
            theta(i) = theta_old(i) + epsilon * sign_S(i) * w_old(2,i) * e;
            theta(i) = min(max(theta(i), theta_min), theta_max);
        end
    end

    % output_timing='updated': 输出更新后状态
    u_request = -(w0 + sum(w(1,:)));

    ramp = min(1, max(0, (k_count - k_start) / max(1, ramp_samples)));
    u_request = ramp * u_request;
else
    u_request = 0;
end

u_demand = u_request;
[u, clipped] = apply_limit_demo5(u_request, u_limit);
diag_out = [sqrt(theta(1)); w(1,1)];
end

function [u_actual, was_clipped] = apply_limit_demo5(u_demand, limit)
% demo5 变体限幅（与 tests/internal/controller_demo5.m::apply_limit 一致）
if isnan(u_demand)
    u_actual = 0;
    was_clipped = true;
else
    u_actual = min(max(u_demand, -limit), limit);
    was_clipped = ~isfinite(u_demand) || abs(u_actual - u_demand) > 1e-10;
end
end
