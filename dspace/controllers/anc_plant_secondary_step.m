function anti = anc_plant_secondary_step(u_prev, PLANT_A, PLANT_B_STAR) %#codegen
%ANC_PLANT_SECONDARY_STEP 离线仿真被控对象（次级路径）逐样本步函数
%
%   与 tests/internal 各引擎中的被控对象递推逐位一致:
%     anti(k) = FIR(B(2:end), [u(k-1), u(k-2), ...]) - FIR(A(2:end), anti_buf)
%
%   仅用于离线闭环模型（RT 模型中被控对象是真实声学路径）。
%   上游必须串一个 Unit Delay: 本函数输入 u_prev = u(k-1)，
%   从而 anti(k) 只依赖 u 的历史，与脚本先算被控对象、后算控制律的
%   样本内顺序一致，且无代数环。
%
%   参数:
%     PLANT_A      — 分母多项式 A（FIR 模型时为 1）
%     PLANT_B_STAR — B_full(2:end)（B_full 含前导零延迟约定）

persistent u_hist anti_hist
if isempty(u_hist)
    u_hist = zeros(1, max(numel(PLANT_B_STAR), 1));
    anti_hist = zeros(1, max(numel(PLANT_A) - 1, 1));
end

u_vec = [u_prev, u_hist(1:end-1)];
anti_new = fir_dot(PLANT_B_STAR, u_vec) - fir_dot(PLANT_A(2:end), anti_hist);

anti_hist = [anti_new, anti_hist(1:end-1)];
u_hist = u_vec;

anti = anti_new;
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
