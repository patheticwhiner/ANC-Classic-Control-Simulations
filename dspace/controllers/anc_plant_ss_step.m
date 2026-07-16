function anti = anc_plant_ss_step(u_prev, PLANT_AF, PLANT_BF, PLANT_CF) %#codegen
%ANC_PLANT_SS_STEP 离线仿真被控对象（状态空间形式，demo2 专用）
%
%   与 tests/internal/controller_demo2.m 仿真引擎的被控对象逐位一致:
%     y_k = Cf * x_plant + d(k)         (取 x_k)
%     x_plant = Af * x_plant + Bf * u_k (控制量作用后更新)
%
%   上游必须串一个 Unit Delay（u_prev = u(k-1)），本步先用 u(k-1) 推进
%   状态到 x_k 再输出 Cf*x_k —— 与脚本"先输出后更新"的时序等价:
%     x_k = Af*x_{k-1} + Bf*u(k-1);  anti_k = Cf*x_k
%
%   仅用于离线闭环模型。

persistent x_plant
if isempty(x_plant)
    x_plant = zeros(size(PLANT_AF, 1), 1);
end

x_plant = PLANT_AF * x_plant + PLANT_BF * u_prev;
anti = PLANT_CF * x_plant;
end
