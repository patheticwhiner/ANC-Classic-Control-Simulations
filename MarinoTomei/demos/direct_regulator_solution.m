function [Pi, Gamma] = direct_regulator_solution(A, b, C, R, P, q, beta)
% DIRECT_REGULATOR_SOLUTION 直接求解线性调节器方程
%   该函数使用MATLAB内置函数求解线性调节器方程
%   ΠR = AΠ + 1/β*b*Γ + P, q = CΠ
%
% 输入参数:
%   A - 系统状态矩阵
%   b - 输入矩阵
%   C - 输出矩阵
%   R - 扰动系统状态矩阵
%   P - 扰动到状态的映射矩阵
%   q - 扰动到输出的映射矩阵
%   beta - 系统参数
%
% 输出参数:
%   Pi - 调节器方程的解Π
%   Gamma - 调节器方程的解Γ (控制输入u = Γw)

% 检查输入参数维度
[n_state, ~] = size(A);
[~, n_input] = size(b);
[n_output, ~] = size(C);
[n_dist, ~] = size(R);

assert(n_input == 1, '当前实现仅支持单输入系统');
assert(n_output == 1, '当前实现仅支持单输出系统');

% 方法1: 构建线性方程组求解

% 创建系数矩阵
coef_matrix = zeros(n_state*n_dist + n_dist, n_state*n_dist + n_input*n_dist);

% 构造Sylvester方程部分
for i = 1:n_state
    for j = 1:n_dist
        row_idx = (i-1)*n_dist + j;
        
        % Π*R项
        for k = 1:n_dist
            col_idx = (i-1)*n_dist + k;
            coef_matrix(row_idx, col_idx) = R(k, j);
        end
        
        % -A*Π项
        for k = 1:n_state
            for l = 1:n_dist
                col_idx = (k-1)*n_dist + l;
                coef_matrix(row_idx, col_idx) = coef_matrix(row_idx, col_idx) - A(i, k) * (l == j);
            end
        end
    end
end

% 填充输入项 -(1/β)*b*Γ
for i = 1:n_state
    for j = 1:n_dist
        row_idx = (i-1)*n_dist + j;
        col_idx = n_state*n_dist + j;
        coef_matrix(row_idx, col_idx) = -(1/beta)*b(i);
    end
end

% 填充输出约束 C*Π = q
for j = 1:n_dist
    row_idx = n_state*n_dist + j;
    for i = 1:n_state
        col_idx = (i-1)*n_dist + j;
        coef_matrix(row_idx, col_idx) = C(i);
    end
end

% 构造右侧向量
rhs = zeros(n_state*n_dist + n_dist, 1);

% 填充Sylvester方程右侧 (P)
for i = 1:n_state
    for j = 1:n_dist
        row_idx = (i-1)*n_dist + j;
        if i == 1
            rhs(row_idx) = P(j);  % P行对应的值
        else
            rhs(row_idx) = 0;     % 其他行为0
        end
    end
end

% 填充输出约束右侧 (q)
for j = 1:n_dist
    row_idx = n_state*n_dist + j;
    rhs(row_idx) = q(j);
end

% 求解线性方程组
sol = coef_matrix \ rhs;

% 提取解：Π和Γ
Pi = zeros(n_state, n_dist);
for i = 1:n_state
    for j = 1:n_dist
        idx = (i-1)*n_dist + j;
        Pi(i, j) = sol(idx);
    end
end

Gamma = zeros(1, n_dist);
for j = 1:n_dist
    idx = n_state*n_dist + j;
    Gamma(j) = sol(idx);
end

% 验证解的正确性
residual1 = Pi*R - A*Pi - (1/beta)*b*Gamma - [P; zeros(size(A,1)-1, size(P,2))];
residual2 = C*Pi - q;

fprintf('--- 直接求解结果验证 ---\n');
fprintf('方程1残差范数: %.6e\n', norm(residual1, 'fro'));
fprintf('方程2残差范数: %.6e\n', norm(residual2, 'fro'));

% 方法2: 可以尝试使用其他MATLAB内置函数求解
% 例如，对于特定结构的问题，可以尝试sylvester函数
% 但由于本问题有附加约束C*Π = q，一般需要联立求解

end
