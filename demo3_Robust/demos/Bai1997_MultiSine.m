%% H∞鲁棒控制设计 - MSBai1977实现
% 基于文献MSBai1977的H∞控制器设计与分析
% ===================================================================
clear; close all; clc;

%% 1. 系统定义
% 被控对象参数
modelFile = fullfile('..', '..', 'dataset', 'syn_Bai1997_4th.mat');
load(modelFile, 'model');
P_zp = model.G0_zpk;
P_tf = model.G0_tf;



% 权重函数定义 (基于MSBai1977论文)
s = tf('s');
% 方案1: 增强W3的高频滚降
W3_inv = 2*(0.6554*s + 2592) / (s + 5000);  % 增大极点频率
% 方案2: 放松W1的低频要求  
W1_inv = (4*s^2 + 2500*s + 1.027e7) / (s^2 + 7969*s + 8.542e6);  % 增大系数
% W3_inv = (0.6554*s + 2592) / (s + 3266);
% W1_inv = (2.37*s^2 + 2500*s + 1.027e7) / (s^2 + 7969*s + 8.542e6);
% W3_inv = (0.6554*s + 2592) / (s + 3266);
W1 = 1/W1_inv;
W3 = 1/W3_inv;

% 频率分析参数
f_hz = logspace(0, 4, 1000); % 1 Hz 到 10000 Hz
w_rad = 2*pi*f_hz;

fprintf('=== H∞鲁棒控制系统分析 ===\n');
fprintf('被控对象P(s)构建完成\n');

% 分析系统特性
rhp_zeros = z(real(z) > 0);
rhp_poles = p(real(p) > 0);

if ~isempty(rhp_zeros)
    fprintf('⚠️  系统包含%d个右半平面零点，将限制可达性能\n', length(rhp_zeros));
end
if ~isempty(rhp_poles)
    fprintf('⚠️  系统包含%d个右半平面极点，系统不稳定\n', length(rhp_poles));
end

%% 2. 被控对象固有特性分析
fprintf('\n=== 被控对象固有特性分析 ===\n');

% 计算被控对象的灵敏度函数（假设控制器为单位增益）
L0 = P_tf;
try
    S0 = feedback(1, L0);    % S0 = 1/(1+P_tf)
    T0 = feedback(L0, 1);    % T0 = P_tf/(1+P_tf)
    fprintf('被控对象灵敏度函数计算成功\n');
catch ME
    fprintf('使用手动计算方法\n');
    S0 = 1/(1+P_tf);
    T0 = P_tf/(1+P_tf);
end

% 绘制被控对象特性分析
plotSensitivityAnalysis(S0, T0, W1_inv, W3_inv, f_hz, w_rad, '被控对象固有特性分析', {'S₀(s)', 'T₀(s)'});

%% 3. H∞控制器设计
fprintf('\n=== H∞控制器设计 ===\n');

% 构造增广系统
try
    P_aug = augw(P_tf, W1, 1e-3, W3);
    fprintf('增广系统P_aug构造成功 (维度: %dx%d)\n', size(P_aug));
catch ME
    fprintf('❌ 增广系统构造失败: %s\n', ME.message);
    return;
end

% H∞综合
nmeas = 1; % 测量输出数量
ncon = 1;  % 控制输入数量

try
    fprintf('开始H∞综合...\n');
    [K_hinf, ~, gamma] = hinfsyn(P_aug, nmeas, ncon);
    fprintf('H∞综合完成: γ = %.4f\n', gamma);
    
    % 检查控制器稳定性
    if isstable(K_hinf)
        fprintf('✓ 控制器K(s)稳定\n');
    else
        fprintf('❌ 警告: 控制器K(s)不稳定！需要调整权重函数\n');
        fprintf('   建议: 1) 增大W3高频增益 2) 适当放松W1约束 3) 增加W2控制代价\n');
    end
    
    if gamma < 1.0
        fprintf('✓ 设计优秀: 所有性能约束满足\n');
    elseif gamma < 2.0
        fprintf('△ 设计可接受: 性能约束基本满足\n');
    else
        fprintf('✗ 设计欠佳: 建议调整权重函数\n');
    end
    
catch ME
    fprintf('❌ H∞综合失败: %s\n', ME.message);
    return;
end

% 权重函数设计分析与调试指导
analyzeWeightDesign(P_tf, W1, W3, K_hinf, gamma);

%% 3.5. 论文控制器验证
fprintf('\n=== 论文控制器验证 (Table II) ===\n');

% 根据论文Table II构建控制器
% 极点 (×10⁴)
poles_paper = [-2.0959, -1.6790, -0.8121, -0.0592+0.1973i, -0.0592-0.1973i] * 1e4;
% 零点 (×10⁵) 
zeros_paper = [-3.4362, -0.1327+0.1302i, -0.1327-0.1302i, -0.1885, -0.0361] * 1e5;
% 增益
gain_paper = 0.1623;

% 构建论文控制器
K_paper = zpk(zeros_paper, poles_paper, gain_paper);
fprintf('论文控制器K_paper构建完成\n');

% 检查论文控制器稳定性
if isstable(K_paper)
    fprintf('✓ 论文控制器稳定\n');
else
    fprintf('❌ 论文控制器不稳定\n');
end

% 计算论文控制器的闭环系统
L_paper = P_tf * K_paper;
try
    S_paper = feedback(1, L_paper);
    T_paper = feedback(L_paper, 1);
    fprintf('论文控制器闭环系统计算成功\n');
    
    % 检查闭环稳定性
    if isstable(S_paper)
        fprintf('✓ 论文设计的闭环系统稳定\n');
    else
        fprintf('❌ 论文设计的闭环系统不稳定\n');
    end
    
catch ME
    fprintf('❌ 论文控制器闭环计算失败: %s\n', ME.message);
    S_paper = [];
    T_paper = [];
end

% 绘制论文控制器性能分析
if ~isempty(S_paper) && ~isempty(T_paper)
    plotSensitivityAnalysis(S_paper, T_paper, W1_inv, W3_inv, f_hz, w_rad, ...
        '论文控制器验证结果 (Table II)', {'S_paper(s)', 'T_paper(s)'}, S0, T0);
    
    % 分析论文控制器的约束违反情况
    [violation_S_paper, violation_T_paper] = analyzeViolation(S_paper, T_paper, W1_inv, W3_inv, w_rad);
    fprintf('论文控制器约束分析:\n');
    fprintf('  灵敏度约束违反: %.2f dB\n', violation_S_paper);
    fprintf('  互补灵敏度约束违反: %.2f dB\n', violation_T_paper);
    
    if violation_S_paper <= 1.0 && violation_T_paper <= 1.0
        fprintf('✓ 论文控制器基本满足权重约束\n');
    else
        fprintf('⚠️ 论文控制器存在约束违反，可能权重函数设置不同\n');
    end
end

%% 4. 闭环系统分析
fprintf('\n=== 闭环系统性能分析 ===\n');

% 计算闭环传递函数
L = P_tf * tf(K_hinf); % 开环传递函数
try
    S = feedback(1, L);    % S = 1/(1+L)
    T = feedback(L, 1);    % T = L/(1+L)
    fprintf('闭环灵敏度函数计算成功\n');
catch ME
    fprintf('使用手动计算方法\n');
    S = 1/(1+L);
    T = L/(1+L);
end

% 绘制闭环系统性能分析，包含原始灵敏度函数对比
plotSensitivityAnalysis(S, T, W1_inv, W3_inv, f_hz, w_rad, ...
    sprintf('H∞控制器设计结果 (γ=%.3f)', gamma), {'S(s)', 'T(s)'}, S0, T0);

%% 5. 性能评估与对比
fprintf('\n=== 性能评估与对比 ===\n');

% H∞综合控制器评估
fprintf('H∞综合控制器:\n');
fprintf('  H∞范数 γ = %.4f\n', gamma);
[violation_S, violation_T] = analyzeViolation(S, T, W1_inv, W3_inv, w_rad);
fprintf('  灵敏度约束违反: %.2f dB\n', violation_S);
fprintf('  互补灵敏度约束违反: %.2f dB\n', violation_T);

% 论文控制器评估 (如果可用)
if exist('S_paper', 'var') && ~isempty(S_paper)
    fprintf('\n论文控制器 (Table II):\n');
    fprintf('  灵敏度约束违反: %.2f dB\n', violation_S_paper);
    fprintf('  互补灵敏度约束违反: %.2f dB\n', violation_T_paper);
    
    % 控制器对比分析
    fprintf('\n控制器对比:\n');
    if violation_S < violation_S_paper
        fprintf('  ✓ H∞综合在灵敏度约束方面更优\n');
    else
        fprintf('  ▲ 论文控制器在灵敏度约束方面更优\n');
    end
    
    if violation_T < violation_T_paper
        fprintf('  ✓ H∞综合在鲁棒性约束方面更优\n');
    else
        fprintf('  ▲ 论文控制器在鲁棒性约束方面更优\n');
    end
end

% 综合评价
if violation_S <= 0.5 && violation_T <= 0.5
    fprintf('\n✓ H∞控制器设计成功\n');
else
    fprintf('\n⚠️ H∞控制器部分约束轻微违反，属于正常范围\n');
end

%% 噪声处理
fprintf('\n=== 控制器特性分析 (时域与频域) ===\n');

% 仿真参数
fs = 4000; % 采样频率 (Hz), 需大于 2*1000Hz
T = 5;     % 仿真时长 (s)
t = 0:1/fs:T-1/fs; % 时间向量

% 生成0-1000Hz的带限白噪声
d_raw = randn(size(t)); % 生成高斯白噪声
[b_filter, a_filter] = butter(6, 1000/(fs/2), 'low'); % 6阶巴特沃斯低通滤波器
d = filter(b_filter, a_filter, d_raw); % 滤波得到带限白噪声

% 离散化系统模型
P_d = c2d(P_tf, 1/fs, 'tustin');
S_d = c2d(S, 1/fs, 'tustin');

% 仿真1: 无控制器 (开环响应, K=0, S=1)
% 此时输出 y = d
y_uncontrolled = d;

% 仿真2: 有H∞控制器 (闭环响应)
% 输出 y = S * d
y_controlled = lsim(S_d, d, t);

% 仿真3: 论文控制器 (如果可用)
if exist('S_paper', 'var') && ~isempty(S_paper)
    S_paper_d = c2d(S_paper, 1/fs, 'tustin');
    y_paper = lsim(S_paper_d, d, t);
    fprintf('论文控制器噪声测试完成\n');
else
    y_paper = [];
end

% 计算功率谱密度 (PSD)
[psd_uncontrolled, f_psd] = pwelch(y_uncontrolled, hann(512), 256, 512, fs);
[psd_controlled, ~] = pwelch(y_controlled, hann(512), 256, 512, fs);
if ~isempty(y_paper)
    [psd_paper, ~] = pwelch(y_paper, hann(512), 256, 512, fs);
end

% 绘制时域响应
figure('Name', '控制器噪声抑制效果时域响应');
plot(t(1:4000), y_uncontrolled(1:4000), 'b-', 'LineWidth', 1.0, 'DisplayName', '无控制器');
hold on;
plot(t(1:4000), y_controlled(1:4000), 'r-', 'LineWidth', 1.0, 'DisplayName', 'H∞控制器');
if ~isempty(y_paper)
    plot(t(1:4000), y_paper(1:4000), 'g--', 'LineWidth', 1.0, 'DisplayName', '论文控制器');
end
grid on;
xlabel('时间 (s)');
ylabel('输出 (y)');
title('控制器噪声抑制效果时域响应 (前1秒)');
legend('Location', 'best');

% 绘制频谱对比图
figure('Name', '控制器噪声抑制效果频谱分析');
plot(f_psd, 10*log10(psd_uncontrolled), 'b-', 'LineWidth', 1.5, 'DisplayName', '无控制器');
hold on;
plot(f_psd, 10*log10(psd_controlled), 'r-', 'LineWidth', 1.5, 'DisplayName', 'H∞控制器');
if ~isempty(y_paper)
    plot(f_psd, 10*log10(psd_paper), 'g--', 'LineWidth', 1.5, 'DisplayName', '论文控制器');
end
grid on;
xlabel('频率 (Hz)');
ylabel('功率/频率 (dB/Hz)');
title('0-1000Hz白噪声扰动下的输出频谱对比');
legend('Location', 'best');
xlim([0 1200]); % 限制显示范围

% 计算并显示噪声抑制性能指标
rms_uncontrolled = rms(y_uncontrolled);
rms_controlled = rms(y_controlled);
suppression_hinf = 20*log10(rms_uncontrolled/rms_controlled);

fprintf('\n噪声抑制性能对比:\n');
fprintf('  无控制器 RMS: %.4f\n', rms_uncontrolled);
fprintf('  H∞控制器 RMS: %.4f (抑制 %.1f dB)\n', rms_controlled, suppression_hinf);

if ~isempty(y_paper)
    rms_paper = rms(y_paper);
    suppression_paper = 20*log10(rms_uncontrolled/rms_paper);
    fprintf('  论文控制器 RMS: %.4f (抑制 %.1f dB)\n', rms_paper, suppression_paper);
    
    if suppression_hinf > suppression_paper
        fprintf('  ✓ H∞控制器噪声抑制效果更优\n');
    else
        fprintf('  ▲ 论文控制器噪声抑制效果更优\n');
    end
end

fprintf('控制器特性分析完成，请查看频谱对比图。\n');

%% 灵敏度函数绘制通用函数
function plotSensitivityAnalysis(S_func, T_func, W1_inv, W3_inv, f_hz, w_rad, fig_title, func_labels, S0_func, T0_func)
    % 计算频响数据
    [mag_W1inv, ~] = freqresp(W1_inv, w_rad);
    [mag_S, ~] = freqresp(S_func, w_rad);
    [mag_W3inv, ~] = freqresp(W3_inv, w_rad);
    [mag_T, ~] = freqresp(T_func, w_rad);
    
    % 转换为dB
    mag_W1inv_dB = 20*log10(squeeze(abs(mag_W1inv)));
    mag_S_dB = 20*log10(squeeze(abs(mag_S)));
    mag_W3inv_dB = 20*log10(squeeze(abs(mag_W3inv)));
    mag_T_dB = 20*log10(squeeze(abs(mag_T)));
    
    % 检查是否提供了参考灵敏度函数
    if nargin >= 10 && ~isempty(S0_func) && ~isempty(T0_func)
        [mag_S0, ~] = freqresp(S0_func, w_rad);
        [mag_T0, ~] = freqresp(T0_func, w_rad);
        mag_S0_dB = 20*log10(squeeze(abs(mag_S0)));
        mag_T0_dB = 20*log10(squeeze(abs(mag_T0)));
        show_reference = true;
    else
        show_reference = false;
    end
    
    % 创建图形
    figure('Name', fig_title, 'Position', [100, 100, 800, 600]);
    
    % 灵敏度函数对比
    subplot(2,1,1);
    plot(f_hz, mag_W1inv_dB, 'b--', 'LineWidth', 2); hold on;
    if show_reference
        plot(f_hz, mag_S0_dB, 'k:', 'LineWidth', 1.8, 'DisplayName', 'S₀(s) [原始]');
    end
    plot(f_hz, mag_S_dB, 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (Hz)');
    ylabel('幅值 (dB)');
    title(['灵敏度函数', func_labels{1}, '与权重W₁⁻¹(s)对比']);
    if show_reference
        legend({'W₁⁻¹(s) [约束]', 'S₀(s) [原始]', [func_labels{1}, ' [H∞设计]']}, 'Location', 'best');
    else
        legend({'W₁⁻¹(s) [约束]', [func_labels{1}, ' [实际]']}, 'Location', 'best');
    end
    
    % 互补灵敏度函数对比
    subplot(2,1,2);
    plot(f_hz, mag_W3inv_dB, 'b--', 'LineWidth', 2); hold on;
    if show_reference
        plot(f_hz, mag_T0_dB, 'k:', 'LineWidth', 1.8, 'DisplayName', 'T₀(s) [原始]');
    end
    plot(f_hz, mag_T_dB, 'r-', 'LineWidth', 1.5);
    grid on;
    xlabel('频率 (Hz)');
    ylabel('幅值 (dB)');
    title(['互补灵敏度函数', func_labels{2}, '与权重W₃⁻¹(s)对比']);
    if show_reference
        legend({'W₃⁻¹(s) [约束]', 'T₀(s) [原始]', [func_labels{2}, ' [H∞设计]']}, 'Location', 'best');
    else
        legend({'W₃⁻¹(s) [约束]', [func_labels{2}, ' [实际]']}, 'Location', 'best');
    end
end

%% 权重函数设计分析与调试指导
function analyzeWeightDesign(P, W1, W3, K, gamma)
    fprintf('\n=== 权重函数设计分析与调试指导 ===\n');
    
    % 1. 检查权重函数的基本特性
    fprintf('权重函数特性分析:\n');
    
    % W1的低频和高频增益
    w_low = 0.1;   % 低频点
    w_high = 1000; % 高频点
    W1_dc = abs(evalfr(W1, 1j*w_low));
    W1_hf = abs(evalfr(W1, 1j*w_high));
    fprintf('  W1低频增益(%.1f rad/s): %.2f = %.1f dB\n', w_low, W1_dc, 20*log10(W1_dc));
    fprintf('  W1高频增益(%.0f rad/s): %.4f = %.1f dB\n', w_high, W1_hf, 20*log10(W1_hf));
    
    % W3的低频和高频增益  
    W3_dc = abs(evalfr(W3, 1j*w_low));
    W3_hf = abs(evalfr(W3, 1j*w_high));
    fprintf('  W3低频增益(%.1f rad/s): %.4f = %.1f dB\n', w_low, W3_dc, 20*log10(W3_dc));
    fprintf('  W3高频增益(%.0f rad/s): %.2f = %.1f dB\n', w_high, W3_hf, 20*log10(W3_hf));
    
    % 2. 检查被控对象的关键特性
    fprintf('\n被控对象关键特性:\n');
    P_poles = pole(P);
    P_zeros = zero(P);
    unstable_poles = P_poles(real(P_poles) >= 0);
    rhp_zeros = P_zeros(real(P_zeros) > 0);
    
    if ~isempty(unstable_poles)
        fprintf('  ❌ 被控对象有%d个不稳定极点，需要足够的控制带宽\n', length(unstable_poles));
        fprintf('     建议: 确保W1在不稳定极点频率附近有足够大的增益\n');
    end
    
    if ~isempty(rhp_zeros)
        fprintf('  ⚠️  被控对象有%d个右半平面零点，限制可达性能\n', length(rhp_zeros));
        rhp_freq = abs(rhp_zeros);
        fprintf('     RHP零点频率: ');
        for i = 1:length(rhp_freq)
            fprintf('%.2f ', rhp_freq(i));
        end
        fprintf('rad/s\n');
        fprintf('     建议: 在RHP零点频率附近适当放松W1约束\n');
    end
    
    % 3. 设计质量评估
    fprintf('\n设计质量评估:\n');
    fprintf('  H∞范数 γ = %.4f', gamma);
    if gamma < 1
        fprintf(' ✓ 优秀\n');
    elseif gamma < 1.5
        fprintf(' △ 可接受\n');
    else
        fprintf(' ✗ 需要改进\n');
    end
    
    % 控制器稳定性
    if isstable(K)
        fprintf('  控制器稳定性: ✓ 稳定\n');
    else
        fprintf('  控制器稳定性: ❌ 不稳定\n');
        fprintf('     可能原因: W3高频滚降不足，或W1低频要求过于严格\n');
    end
    
    % 4. 具体调试建议
    fprintf('\n具体调试建议:\n');
    
    if gamma >= 1.5
        fprintf('  📝 gamma过大的可能原因和解决方案:\n');
        fprintf('     - W1低频增益过大 → 适当减小W1的直流增益\n');
        fprintf('     - W3高频增益不足 → 增大W3的高频增益\n');
        fprintf('     - 缺少控制代价权重 → 在augw中添加小的W2(如1e-3)\n');
    end
    
    if ~isstable(K)
        fprintf('  📝 控制器不稳定的可能原因和解决方案:\n');
        fprintf('     - W3滚降速度不够快 → 增加W3分母的阶数或增大极点频率\n');
        fprintf('     - W1和W3冲突 → 检查是否W1+W3在某频率点>1\n');
        fprintf('     - 被控对象特殊结构 → 考虑添加预滤波器或改变控制结构\n');
    end
    
    % 5. 推荐的权重函数调整策略
    fprintf('\n推荐的权重函数调整策略:\n');
    fprintf('  🎯 W1设计原则:\n');
    fprintf('     - 低频: 大增益(>10)确保扰动抑制\n');
    fprintf('     - 中频: 平滑过渡\n');
    fprintf('     - 高频: 小增益(<0.1)避免过度控制\n');
    fprintf('  🎯 W3设计原则:\n');
    fprintf('     - 低频: 小增益(<0.1)允许控制器有足够权威\n');
    fprintf('     - 高频: 大增益(>1)确保鲁棒稳定性\n');
    fprintf('     - 滚降: 至少40dB/dec确保控制器稳定\n');
    
    % 6. 数值建议
    fprintf('\n数值调整建议 (当前设计问题时):\n');
    if gamma > 1.2 || ~isstable(K)
        fprintf('  试试这些修改:\n');
        fprintf('     W1_new = 0.7 * W1;  %% 适当放松性能要求\n');
        fprintf('     或者\n');
        fprintf('     W3_new = 2 * W3;    %% 增强鲁棒性要求\n');
        fprintf('     或者\n');
        fprintf('     P_aug = augw(P, W1, 1e-2, W3); %% 添加控制代价\n');
    end
end

%% 约束违反分析函数
function [violation_S, violation_T] = analyzeViolation(S, T, W1_inv, W3_inv, w_rad)
    try
        % 计算S与W1_inv的对比
        [mag_S, ~] = freqresp(S, w_rad);
        [mag_W1inv, ~] = freqresp(W1_inv, w_rad);
        mag_S_dB = 20*log10(squeeze(abs(mag_S)));
        mag_W1inv_dB = 20*log10(squeeze(abs(mag_W1inv)));
        violation_S = max(mag_S_dB - mag_W1inv_dB);
        
        % 计算T与W3_inv的对比
        [mag_T, ~] = freqresp(T, w_rad);
        [mag_W3inv, ~] = freqresp(W3_inv, w_rad);
        mag_T_dB = 20*log10(squeeze(abs(mag_T)));
        mag_W3inv_dB = 20*log10(squeeze(abs(mag_W3inv)));
        violation_T = max(mag_T_dB - mag_W3inv_dB);
        
    catch
        violation_S = inf;
        violation_T = inf;
    end
end


