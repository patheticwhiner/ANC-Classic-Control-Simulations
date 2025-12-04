%% 1 系统建模
% 物理参数定义
m = 1;   % 质量 (kg)
k = 1;   % 弹簧刚度 (N/m)
b = 0.5; % 阻尼系数 (N·s/m)

% 状态空间方程
A = [0 1; -k/m -b/m];  % 系统矩阵 [3](@ref)
B = [0; 1/m];           % 输入矩阵 [3](@ref)
C = eye(2);              % 输出矩阵（全状态观测）
D = 0;                  % 直接传递项
sys = ss(A, B, C, D);   % 构建状态空间模型

%% 2 权重矩阵设置
Q = diag([10, 1]);  % 状态权重：更重视位置误差[2](@ref)
R = 0.1;            % 控制输入权重：允许适度控制能耗[3](@ref)

%% 3 LQR反馈增益计算
[K, S, E] = lqr(A, B, Q, R);  % 求解LQR最优反馈矩阵[1,3](@ref)
fprintf('反馈矩阵 K = [%.4f, %.4f]\n', K(1), K(2));

%% 4 闭环系统仿真
% 构建闭环系统
A_closed = A - B*K;           % 闭环系统矩阵[3](@ref)
sys_closed = ss(A_closed, B, C, D);

% 初始条件与时间范围
x0 = [1; 0];                  % 初始位置1m，速度0m/s
t = 0:0.01:5;                 % 仿真时间5秒

% 仿真动态响应
[y, t, x] = initial(sys_closed, x0, t);
u = -K * x';                  % 计算控制输入序列

% 绘图
figure;
subplot(3,1,1); plot(t, x(:,1)); title('位置响应'); ylabel('x1 (m)');
subplot(3,1,2); plot(t, x(:,2)); title('速度响应'); ylabel('x2 (m/s)');
subplot(3,1,3); plot(t, u);     title('控制输入');  ylabel('u (N)'); xlabel('时间 (s)');

%% 5 结果验证
% 验证1：闭环系统稳定性
fprintf('闭环系统特征值实部：%.3f, %.3f\n', real(E(1)), real(E(2))); 
assert(all(real(E) < 0), '闭环系统不稳定！');

% 验证2：稳态误差分析
steady_state_error = abs(x(end,1)); 
fprintf('稳态位置误差：%.4f m\n', steady_state_error);
assert(steady_state_error < 0.05, '稳态误差超出阈值！');

% 验证3：控制能量计算
control_energy = trapz(t, u.^2); 
fprintf('控制能量消耗：%.3f J\n', control_energy);

%% 功能扩展
% 频域分析：绘制闭环系统Bode图
figure; bode(sys_closed); title('闭环系统频域特性');

% 数据导出：保存仿真结果供进一步分析
results = struct('time',t, 'position',x(:,1), 'control',u);
save('lqr_test_results.mat', 'results');