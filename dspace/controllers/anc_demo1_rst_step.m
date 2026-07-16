function [u, u_demand, clipped, diag_out] = anc_demo1_rst_step( ...
        e, enable, R0, S0, control_scale, u_limit) %#codegen
%ANC_DEMO1_RST_STEP 固定 RST 控制器逐样本步函数（Demo1/Demo4 固定版共用）
%
%   移植自 tests/internal/controller_demo1.m::sim_fixed_rst（controller_demo4
%   同名引擎仅多一个 control_scale 乘法；demo1 传 control_scale=1，
%   IEEE754 乘以 1.0 为恒等，两者位级一致）。
%
%   被控对象递推留在 Simulink（Unit Delay + Discrete Filter），本函数只含
%   控制律:  u = control_scale * ( -(R0*[y] + S0(2:end)*[u]) / S0(1) )
%
%   输入:
%     e       — 误差麦克风采样 y(k)
%     enable  — >0.5 使能；<=0.5 时输出 0 并复位全部内部状态（重新武装）
%   参数（Parameter 作用域，编译期常量，模型工作区解析）:
%     R0, S0          — 中央 RST 多项式（行向量）
%     control_scale   — 输出缩放（demo1: 1）
%     u_limit         — 执行器硬限幅
%   输出:
%     u        — 限幅后控制量
%     u_demand — 限幅前需求
%     clipped  — 本样本是否触发限幅
%     diag_out — [0; 0]（与自适应控制器接口对齐的诊断占位）

persistent y_hist u_hist
if isempty(y_hist)
    y_hist = zeros(1, numel(R0));
    u_hist = zeros(1, max(numel(S0) - 1, 1));
end

if enable <= 0.5
    y_hist(:) = 0;
    u_hist(:) = 0;
    u = 0;
    u_demand = 0;
    clipped = false;
    diag_out = [0; 0];
    return;
end

% ---- 控制律（与 sim_fixed_rst 运算顺序一致）----
y_vec = [e, y_hist(1:end-1)];
u_request = -(fir_dot(R0, y_vec) + fir_dot(S0(2:end), u_hist)) / S0(1);
u_request = control_scale * u_request;

[u_new, was_clipped] = apply_limit(u_request, u_limit);

% ---- 滚动缓冲区更新 ----
u_hist = [u_new, u_hist(1:end-1)];
y_hist = y_vec;

u = u_new;
u_demand = u_request;
clipped = was_clipped;
diag_out = [0; 0];
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
