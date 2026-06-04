% 图像Bode曲线数据提取与传递函数拟合脚本
% 需要：Image Processing Toolbox, System Identification Toolbox

clear; clc; close all;

%% 1. 读取并显示Bode图图片
img = imread('BodeCurve_3.png');
if size(img,3) == 3
    img_gray = rgb2gray(img);
else
    img_gray = img;
end

% 1. 二值化（黑色为0，白色为1）
bw = imbinarize(img_gray, 'adaptive', 'ForegroundPolarity', 'dark', 'Sensitivity', 0.4);
figure; imshow(bw); title('二值化图像');

% 2. 取反，黑色曲线为1
bw_curve = ~bw;
figure; imshow(bw_curve); title('黑色曲线掩码');

% 2.5 形态学膨胀，连接虚线
se = strel('disk', 4); % 盘型结构元，半径可调整
bw_dilated = imdilate(bw_curve, se);
figure; imshow(bw_dilated); title('盘型膨胀后曲线掩码');

% 3. 去除小噪点
bw_clean = bwareaopen(bw_dilated, 10);
figure; imshow(bw_clean); title('去噪后曲线掩码');

% 3. 骨架化，细化曲线为单像素宽度
bw_skel = bwskel(bw_clean);
figure; imshow(bw_skel); title('骨架化后曲线掩码');

[y_pix, x_pix] = find(bw_skel);
margin = 15; % 可根据实际图片调整
img_h = size(bw_skel,1);
img_w = size(bw_skel,2);
valid_idx = (x_pix > margin) & (x_pix < img_w - margin) & ...
            (y_pix > margin) & (y_pix < img_h - margin);
x_pix = x_pix(valid_idx);
y_pix = y_pix(valid_idx);


% 4. 每列kmeans聚类分离两条曲线
x_unique = unique(x_pix);
curve1 = nan(length(x_unique),2); % 预分配，便于插值
curve2 = nan(length(x_unique),2);
for i = 1:length(x_unique)
    idx = find(x_pix == x_unique(i));
    y_vals = y_pix(idx);
    if length(y_vals) >= 2
        [~, C] = kmeans(double(y_vals), 2, 'Replicates', 3);
        C = sort(C);
        [~, idx1] = min(abs(y_vals - C(1)));
        [~, idx2] = min(abs(y_vals - C(2)));
        curve1(i,:) = [x_unique(i), y_vals(idx1)];
        curve2(i,:) = [x_unique(i), y_vals(idx2)];
    elseif length(y_vals) == 1
        % 只有一个点时，同时赋值给两条曲线
        curve1(i,:) = [x_unique(i), y_vals];
        curve2(i,:) = [x_unique(i), y_vals];
    else
        curve1(i,:) = [x_unique(i), nan];
        curve2(i,:) = [x_unique(i), nan];
    end
end

% 5. 对nan断点插值补全
curve1(:,2) = fillmissing(curve1(:,2), 'linear');
curve2(:,2) = fillmissing(curve2(:,2), 'linear');

% 6. 映射到物理坐标
img_size = size(img_gray);
x_min = 0; x_max = 1200; % 频率范围
y_min = -30; y_max = 10; % 幅值范围

f_curve1 = curve1(:,1) / img_size(2) * (x_max - x_min) + x_min;
mag_db_curve1 = (1 - curve1(:,2) / img_size(1)) * (y_max - y_min) + y_min;

f_curve2 = curve2(:,1) / img_size(2) * (x_max - x_min) + x_min;
mag_db_curve2 = (1 - curve2(:,2) / img_size(1)) * (y_max - y_min) + y_min;

figure;
plot(f_curve1, mag_db_curve1, 'b.-', 'DisplayName', '曲线1');
hold on;
plot(f_curve2, mag_db_curve2, 'r.-', 'DisplayName', '曲线2');
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
legend; grid on; title('骨架化+聚类+插值分离的两条曲线');
axis([0, 1200, -15, 10]);

%% 3. 数据处理
w1 = 2 * pi * f_curve1; % 曲线1角频率
mag1 = 10.^(mag_db_curve1/20); % 曲线1线性幅值
w2 = 2 * pi * f_curve2; % 曲线2角频率
mag2 = 10.^(mag_db_curve2/20); % 曲线2线性幅值

%% 4. 传递函数拟合（以2阶为例，可调整阶数）
order = 2;
Ts = 0; % 连续系统

data_idfrd1 = idfrd(mag1 .* exp(1j*zeros(size(mag1))), w1, Ts);
data_idfrd2 = idfrd(mag2 .* exp(1j*zeros(size(mag2))), w2, Ts);

sys1 = tfest(data_idfrd1, order);
sys2 = tfest(data_idfrd2, order, order);

figure;
bode(sys1); hold on; bode(sys2); title('两条曲线的tfest拟合Bode图'); legend('曲线1拟合','曲线2拟合');

%% 5. 与原始数据对比（全部用plot，频率轴一致）
mag_fit1 = squeeze(abs(freqresp(sys1, w1)));
mag_fit2 = squeeze(abs(freqresp(sys2, w2)));
figure;
plot(f_curve1, mag_db_curve1, 'bo-', 'DisplayName', '曲线1原始数据'); hold on;
plot(f_curve1, 20*log10(mag_fit1), 'r--', 'LineWidth', 1.5, 'DisplayName', '曲线1 tfest拟合');
plot(f_curve2, mag_db_curve2, 'go-', 'DisplayName', '曲线2原始数据');
plot(f_curve2, 20*log10(mag_fit2), 'k--', 'LineWidth', 1.5, 'DisplayName', '曲线2 tfest拟合');
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('两条曲线原始数据与tfest拟合对比');
grid on;
legend('show');

%% 6. 显示传递函数
fprintf('曲线1 tfest拟合得到的传递函数为：\n');
sys1
fprintf('曲线2 tfest拟合得到的传递函数为：\n');
sys2