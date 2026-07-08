% % 定义单频率信号
% function d_t = single_frequency_signal(t)
%     d_t = 0.3 + sin(0.5 * t);
% end
% 
% % 修改 estimator_dynamics 函数中的输入信号
% function dydt = estimator_dynamics(t, y, k, epsilon, num_frequencies)
%     % 计算输入信号 d(t)
%     d_t = single_frequency_signal(t);
%     
%     % ...（其他代码保持不变）
% end

function frequency_estimator()
    % 定义参数
    k = 0.5;
    epsilon = 0.05;
    num_frequencies = 2; % 干扰频率数
    initial_conditions = [0; 0; 0; 0; 0; 0.3; 4.2]; % 初始条件
    
    % 定义时间向量
    t_span = [0, 200];
    [t, y] = ode45(@(t, y) estimator_dynamics(t, y, k, epsilon, num_frequencies), t_span, initial_conditions);
    
    % 提取结果
    w0_hat = y(:, 1);
    w1_hat = y(:, 2);
    w2_hat = y(:, 3);
    w3_hat = y(:, 4);
    w4_hat = y(:, 5);
    theta1_hat = y(:, 6);
    theta2_hat = y(:, 7);
    
    % 绘制结果
    figure;
    subplot(2, 2, 1);
    plot(t, w0_hat);
    title('估计的偏置项 \hat{w}_0');
    xlabel('时间');
    grid on;
    
    subplot(2, 2, 2);
    plot(t, w1_hat);
    title('估计的正弦项 \hat{w}_1');
    xlabel('时间');
    grid on;
    
    subplot(2, 2, 3);
    plot(t, theta1_hat);
    title('频率估计 \hat{\theta}_1');
    xlabel('时间');
    grid on;
    
    subplot(2, 2, 4);
    plot(t, theta2_hat);
    title('频率估计 \hat{\theta}_2');
    xlabel('时间');
    grid on;
end

function dydt = estimator_dynamics(t, y, k, epsilon, num_frequencies)
    % 计算输入信号 d(t)
    d_t = 0.3 + sin(0.5 * t) + sin(2 * t);
    
    % 计算估计的 d_hat(t)
    d_hat = y(1) + y(2) + y(4);
    
    % 计算误差
    error = d_t - d_hat;
    
    % 更新估计器状态
    dydt(1) = k * error;
    dydt(2) = y(3) + k * error;
    dydt(3) = -y(6) * y(2);
    dydt(4) = y(5) - k * error;
    dydt(5) = -y(7) * y(4);
    dydt(6) = epsilon * y(3) * error;
    dydt(7) = -epsilon * y(5) * error;
    
    % 转换为列向量
    dydt = dydt';
end