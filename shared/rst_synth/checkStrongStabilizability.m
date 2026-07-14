function [isStrong, unstable_zeros, multiplicities, Bminus, Bplus] = checkStrongStabilizability(B, A)
% CHECKSTRONGSTABILIZABILITY 判定 P(z) = B(z)/A(z) 是否满足 strong stabilizability
% 输入:
%   B, A - z⁻¹ 表达下的多项式系数
% 输出:
%   isStrong - 布尔值，表示系统是否满足 strong stabilizability
%   unstable_zeros - 不稳定零点列表 (cell数组)
%                   实根: {z1}, {z2}, ...
%                   复根: {[z3, conj(z3)]}, {[z4, conj(z4)]}, ...
%   multiplicities - 对应的重数列表 (cell数组)
%   z_B - 所有零点列表
%   Bminus - 包含所有不稳定零点的多项式
%   Bplus - 包含所有稳定零点的多项式
%
% 注意: 传递函数满足 strong stabilizability 当且仅当在不稳定区域（|z|<1）
%       任意两个相邻的实零点之间，存在偶数个（可以为0个）实极点

    % 先将 z⁻¹ 表达的多项式变为 z 表达
    B_rev = fliplr(B);  % B(z)
    A_rev = fliplr(A);  % A(z)

    % 计算零点与极点
    z_B = roots(B_rev);
    z_A = roots(A_rev);
    
    % ===== 第一部分：打印系统信息和绘图 =====
    fprintf('============= 系统强可稳定性分析 =============\n');
    
    % 打印零点信息
    fprintf('系统零点:\n');
    for i = 1:length(z_B)
        if imag(z_B(i)) == 0
            fprintf('零点 %d: %.4f (|λ| = %.4f)\n', i, real(z_B(i)), abs(z_B(i)));
        else
            fprintf('零点 %d: %.4f%+.4fi (|λ| = %.4f)\n', i, real(z_B(i)), imag(z_B(i)), abs(z_B(i)));
        end
    end
    
    % 打印极点信息
    fprintf('系统极点:\n');
    for i = 1:length(z_A)
        if imag(z_A(i)) == 0
            fprintf('极点 %d: %.4f (|λ| = %.4f)\n', i, real(z_A(i)), abs(z_A(i)));
        else
            fprintf('极点 %d: %.4f%+.4fi (|λ| = %.4f)\n', i, real(z_A(i)), imag(z_A(i)), abs(z_A(i)));
        end
    end
    
    % 识别不稳定实零点和不稳定实极点
    unstable_real_zeros = z_B(abs(z_B) < 1 & imag(z_B) == 0);
    unstable_real_zeros = sort(unstable_real_zeros);
    real_poles = z_A(imag(z_A) == 0);
    unstable_real_poles = z_A(abs(z_A) < 1 & imag(z_A) == 0);
    unstable_real_poles = sort(unstable_real_poles);
    
    % 绘图：实极点和不稳定实零点分布
    figure; hold on;

    % 绘制单位圆
    theta = linspace(0, 2*pi, 100);
    plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1, 'DisplayName', '单位圆');

    % 绘制实极点
    if ~isempty(unstable_real_poles)
        plot(unstable_real_poles, zeros(size(unstable_real_poles)), 'bx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '不稳定实极点');
    end

    % 绘制不稳定实零点
    if ~isempty(unstable_real_zeros)
        plot(unstable_real_zeros, zeros(size(unstable_real_zeros)), 'ro', 'MarkerSize', 8, 'LineWidth', 2, 'DisplayName', '不稳定实零点');
    end

    % 添加标签和图例
    title('系统实极点和不稳定实零点分布');
    xlabel('实部'); ylabel('虚部');
    axis equal;
    axis([-1.5 1.5 -1.5 1.5]);    
    legend('Location', 'Best');
    grid on;
    
    % ===== 第二部分：判断强可稳定性 =====
    % 如果不稳定实零点数 <= 1，一定满足
    if length(unstable_real_zeros) <= 1
        isStrong = true;
        disp('系统 Strong Stabilizable');
    else
        % 对每对相邻零点，统计其间实极点数量
        isStrong = true;
        for i = 1:length(unstable_real_zeros)-1
            z1 = unstable_real_zeros(i);
            z2 = unstable_real_zeros(i+1);

            count = sum(real_poles > z1 & real_poles < z2);
            if mod(count, 2) ~= 0
                isStrong = false;
                disp('系统 Strong UNStabilizable!');
                break;
            end
        end
        
        if isStrong
            disp('系统 Strong Stabilizable');
        end
    end
    
    % ===== 第三部分：分析不稳定零点 =====
    % 将零点分类为稳定和不稳定
    unstable_indices = abs(z_B) < 1;
    unstable_zeros_array = z_B(unstable_indices);
    stable_zeros_array = z_B(~unstable_indices);
    
    if isempty(unstable_zeros_array)
        unstable_zeros = {};
        multiplicities = {};
        Bminus = 1;
        Bplus = B;
    else
        % 分类为实零点和复零点
        real_zeros = unstable_zeros_array(imag(unstable_zeros_array) == 0);
        complex_zeros = unstable_zeros_array(imag(unstable_zeros_array) ~= 0);
        
        unstable_zeros = {};
        multiplicities = {};
        
        % 处理实零点
        if ~isempty(real_zeros)
            % 分组相似的实零点并计数
            unique_real = [];
            counts_real = [];
            
            temp_real = real_zeros;
            while ~isempty(temp_real)
                current = temp_real(1);
                similar = abs(temp_real - current) < 1e-4;
                
                unique_real = [unique_real; current];
                counts_real = [counts_real; sum(similar)];
                
                temp_real = temp_real(~similar);
            end
            
            for i = 1:length(unique_real)
                unstable_zeros{end+1} = unique_real(i);
                multiplicities{end+1} = counts_real(i);
            end
        end
        
        % 处理复零点
        if ~isempty(complex_zeros)
            % 只处理上半平面的复零点
            upper_half = complex_zeros(imag(complex_zeros) > 0);
            
            unique_complex = [];
            counts_complex = [];
            
            for i = 1:length(upper_half)
                z = upper_half(i);
                
                % 检查是否已处理过这个值
                if any(abs(unique_complex - z) < 1e-10)
                    continue;
                end
                
                % 计算重数
                similar_to_z = sum(abs(complex_zeros - z) < 1e-10);
                similar_to_conj_z = sum(abs(complex_zeros - conj(z)) < 1e-10);
                
                % 取较小的数作为共轭对的重数
                pair_count = min(similar_to_z, similar_to_conj_z);
                
                unique_complex = [unique_complex; z];
                counts_complex = [counts_complex; pair_count];
                
                % 将共轭复根对放在同一个cell中
                unstable_zeros{end+1} = [z, conj(z)];
                multiplicities{end+1} = pair_count;
            end
        end
        
        % ===== 第四部分：构造多项式因式分解 =====
        % 首先翻转B为z表示法
        B_z = fliplr(B);
        
        % 构造包含不稳定/稳定零点的多项式
        if ~isempty(unstable_zeros_array)
            Bminus_z = poly(unstable_zeros_array);
        else
            Bminus_z = 1;
        end
        
        if ~isempty(stable_zeros_array)
            Bplus_z = poly(stable_zeros_array);
        else
            Bplus_z = 1;
        end
        
        % 计算缩放因子使得 B = Bminus * Bplus
        if ~isempty(unstable_zeros_array) && ~isempty(stable_zeros_array)
            scaling = B_z(1) / (Bminus_z(1) * Bplus_z(1));
            Bminus_z = scaling * Bminus_z;
        elseif ~isempty(unstable_zeros_array) % 只有不稳定零点
            scaling = B_z(1) / Bminus_z(1);
            Bminus_z = scaling * Bminus_z;
            Bplus_z = 1;
        else % 只有稳定零点
            scaling = B_z(1) / Bplus_z(1);
            Bplus_z = scaling * Bplus_z;
            Bminus_z = 1;
        end
        
        % 验证因式分解是否正确
        B_reconstructed = conv(Bminus_z, Bplus_z);
        error_margin = max(abs(B_z - B_reconstructed));
        
        if error_margin > 1e-10
            warning('多项式因式分解可能不准确，最大误差: %e', error_margin);
        end
        
        % 转换回z^(-1)表示法
        Bminus = fliplr(Bminus_z);
        Bplus = fliplr(Bplus_z);
    end
end