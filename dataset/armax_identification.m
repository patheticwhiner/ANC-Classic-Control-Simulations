%% ARMAX 模型辨识与验证 (使用系统辨识工具箱)
%
% 该脚本加载由 'id_validation_workflow.m' 准备的数据，
% 并使用 MATLAB 系统辨识工具箱来辨识和验证 ARMAX 模型。

%% 1. 初始化
clear; close all; clc;

%% 2. 加载和预处理数据
% 加载仿真数据
info = DataManager('raw_dspace');
x = info.signal.x; y = info.signal.y; t = info.signal.t; Fs = info.fs;
fs = Fs;

% 预处理 (与原脚本保持一致)
% 去趋势处理
data_x = detrend(x);
data_y = detrend(y);

% 陷波滤波 - 滤除50Hz及谐波
f_notch = [50, 100, 150];
for i = 1:length(f_notch)
    wo = f_notch(i)/(fs/2);
    bw = wo/35;
    [b, a] = iirnotch(wo, bw);
    data_y = filtfilt(b, a, data_y);
end

%% 3. 数据集划分
% 将数据分为辨识集和验证集（例如，各占50%）
N_total = length(data_x);
split_point = round(N_total / 2);

% 辨识集
id_data_x = data_x(1:split_point);
id_data_y = data_y(1:split_point);

% 验证集
val_data_x = data_x(split_point+1:end);
val_data_y = data_y(split_point+1:end);
val_t = t(split_point+1:end); % 验证集对应的时间轴

fprintf('数据总量: %d 点\n', N_total);
fprintf('辨识集大小: %d 点\n', length(id_data_x));
fprintf('验证集大小: %d 点\n', length(val_data_x));

%% 4. ARMAX 模型辨识与验证
% 准备数据
% 将辨识数据和验证数据打包成系统辨识工具箱所需的 iddata 格式
identification_data = iddata(id_data_y, id_data_x, 1/fs);
validation_data = iddata(val_data_y, val_data_x, 1/fs);

% 设置ARMAX模型阶数 [na nb nc nk]
% na: 自回归部分的阶数
% nb: 外部输入的阶数
% nc: 移动平均噪声模型的阶数
% nk: 输入延迟
% 注意: 这些阶数需要根据实际系统进行调整以获得最佳效果
armax_orders = [30 30 30 22]; 

% 使用辨识集进行ARMAX模型辨识
fprintf('\n正在辨识ARMAX模型...\n');
armax_model = armax(identification_data, armax_orders);

% 在辨识集上评估模型
[~, fit_id_armax, ~] = compare(identification_data, armax_model);
fprintf('ARMAX模型在辨识集上的Best Fit为: %.2f%%\n', fit_id_armax);

% 使用验证集进行模型验证
fprintf('正在使用验证集验证ARMAX模型...\n');
[y_val_pred_armax, fit_val_armax, ~] = compare(validation_data, armax_model);

% --- 绘图 ---
figure;
plot(val_t, val_data_y, 'b', 'LineWidth', 1);
hold on;
plot(val_t, y_val_pred_armax.OutputData, 'r--', 'LineWidth', 1.5);
title(sprintf('ARMAX 模型验证结果 (拟合度: %.2f%%)', fit_val_armax(1)));
legend('实际输出 (验证集)', 'ARMAX模型预测输出');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

fprintf('ARMAX模型在验证集上的Best Fit为: %.2f%%\n', fit_val_armax(1));

% 显示模型参数
disp('辨识出的ARMAX模型参数:');
disp(armax_model);

% 绘制ARMAX模型的频率响应 (传递函数 G = B/A)
figure;
freqz(armax_model.B, armax_model.A, 2048, fs);
title('辨识得到的ARMAX模型频率响应 (G = B/A)');
grid on;

%% 5. 手动计算拟合度
% 为了与LMS脚本的验证方法保持一致，此处额外手动计算拟合度

% 辨识集
y_id_pred_manual = predict(armax_model, identification_data);
figure;
plot(val_t, id_data_y); hold on;
plot(val_t, y_id_pred_manual.OutputData);
id_error_manual = id_data_y - y_id_pred_manual.OutputData;
plot(val_t, id_error_manual);
fit_percent_id_manual = (1 - norm(id_error_manual)^2 / norm(id_data_y)^2) * 100;
fprintf('\n手动计算 --> ARMAX模型在辨识集上的拟合度为: %.2f%%\n', fit_percent_id_manual);

% 验证集
y_val_pred_manual = predict(armax_model, validation_data);
val_error_manual = val_data_y - y_val_pred_manual.OutputData;
figure;
plot(val_t, val_data_y); hold on;
plot(val_t, y_val_pred_manual.OutputData); 
plot(val_t, val_error_manual);
fit_percent_val_manual = (1 - norm(val_error_manual)^2 / norm(val_data_y)^2) * 100;
fprintf('手动计算 --> ARMAX模型在验证集上的拟合度为: %.2f%%\n', fit_percent_val_manual);

%% 6. 仿真验证 (纯模型预测)
% 使用 sim 函数进行纯仿真，这更能反映模型在开环应用中的表现
fprintf('\n正在进行纯仿真验证 (sim)...\n');

% sim 函数只使用输入信号 u 来模拟输出，不使用真实的 y 进行校正
y_val_sim = sim(armax_model, validation_data.u);
sim_error = val_data_y - y_val_sim;
fit_percent_sim = (1 - norm(sim_error)^2 / norm(val_data_y)^2) * 100;

figure('Name', '纯仿真验证(验证集)');
plot(val_t, val_data_y, 'b'); hold on;
plot(val_t, y_val_sim, 'm-.');
plot(val_t, sim_error, 'g:');
title(sprintf('纯仿真验证 (sim) 结果 (拟合度: %.2f%%)', fit_percent_sim));
legend('实际输出', '仿真输出 (sim)', '仿真误差');
xlabel('时间 (s)');
ylabel('幅值');
grid on;

fprintf('纯仿真 (sim) --> ARMAX模型在验证集上的拟合度为: %.2f%%\n', fit_percent_sim);

%% 7. 保存辨识结果
% 创建一个结构体来保存模型和相关信息
ARMAXmodel.model = armax_model;
ARMAXmodel.fs = fs;
ARMAXmodel.orders = armax_orders;
ARMAXmodel.description = 'ARMAX model identified from experimental data. G(q) = B(q)/A(q), H(q) = C(q)/A(q)';

% 生成符合命名规范的文件名: armax_<orders>_<date>.mat
dateStr = datestr(now, 'yyyy-mm-dd');
ordersStr = sprintf('%d%d%d%d', armax_orders(1), armax_orders(2), ...
                       armax_orders(3), armax_orders(4));
filename = fullfile('dataset', sprintf('armax_%s_%s.mat', ordersStr, dateStr));

% 保存到 .mat 文件
save(filename, 'ARMAXmodel');

fprintf('\nARMAX模型已成功保存到 %s 文件中。\n', filename);;