function y_current = realtimeIIR(x_new, b, a, reset)
% REALTIMEIIR 实时IIR滤波器函数
% 输入参数:
%   x_new  - 当前输入样本
%   b      - 分子系数向量
%   a      - 分母系数向量
%   reset  - 重置滤波器状态（可选，默认0）
% 输出参数:
%   y_current - 当前输出样本

persistent x_buf y_buf
if isempty(x_buf) || (nargin>3 && reset)
    % 初始化缓冲区
    x_buf = zeros(1, length(b));    % 输入缓冲区
    y_buf = zeros(1, length(a)-1);  % 输出缓冲区
end

% 更新输入缓冲区
x_buf = [x_new, x_buf(1:end-1)];

% 计算当前输出
y_current = sum(b .* x_buf) - sum(a(2:end) .* y_buf);

% 更新输出缓冲区
y_buf = [y_current, y_buf(1:end-1)];
end