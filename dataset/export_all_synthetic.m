% dataset/export_all_synthetic.m
% 一次性导出仓库中所有硬编码理论模型到 dataset/ 目录
% 命名规范: syn_{FEATURE}.mat  (与现有 syn_whitenoise_ssmodel.mat 一致)
% 结构: 每个 .mat 文件含 model struct, 与 DataManager 兼容
% 注意: 论文年份(如 TAC2015)是来源标识, 不是数据采集日期

clear; close all; clc;
script_dir = fileparts(mfilename('fullpath'));
if isempty(script_dir), script_dir = pwd; end

fprintf('========== 导出所有合成模型到 dataset/ ==========\n');

%% Model 1: Jafari TAC 2015 — 离散 AVC (3阶, Ts=1/480)
Ts = 1/480;  z = tf('z', Ts);
G0 = (-0.00146*(z-0.1438)*(z-1)) / ((z - 0.7096)*(z^2 - 0.04369*z + 0.01392));
model = struct('name','JafariTAC2015_AVC','type','syn_tf', ...
    'domain','discrete','G0',G0,'Ts',Ts,'fs',1/Ts,'orders',[3], ...
    'source','Jafari, Ioannou & Rudd (2015) IEEE TAC', ...
    'desc','Discrete 3rd-order AVC plant: G0(z)=-0.00146(z-0.1438)(z-1)/((z-0.7096)(z^2-0.04369z+0.01392))');
save(fullfile(script_dir, 'syn_TAC2015_3rd.mat'), 'model');
fprintf('  [1/7] syn_TAC2015_3rd.mat\n');

%% Model 2: Jafari JVC 2017 — 连续 CC (3阶)
s = tf('s');
G0 = 0.5*(s-0.2)/(s^2+s+1.25);
model = struct('name','JafariJVC2017_CC','type','syn_tf', ...
    'domain','continuous','G0',G0,'fs_nominal',10000,'orders',[3], ...
    'source','Jafari & Ioannou (2016) JVC', ...
    'desc','Continuous 3rd-order CC plant: G0(s)=0.5(s-0.2)/(s^2+s+1.25), NMP zero at s=0.2');
save(fullfile(script_dir, 'syn_JVC2017_3rd.mat'), 'model');
fprintf('  [2/7] syn_JVC2017_3rd.mat\n');

%% Model 3: Bai 1997 — 耳机 H∞ (4阶, 连续)
zp = [-3.0841, 1.0320, -0.4387, 0.0034];
pp = [0.6612+0.3483i, 0.6612-0.3483i, -0.4426+0.3324i, -0.4426-0.3324i];
P_zpk = zpk(zp, pp, 0.3921);
P_tf = tf(P_zpk);
model = struct('name','Bai1997_Headphone','type','syn_zpk', ...
    'domain','continuous','G0_zpk',P_zpk,'G0_tf',P_tf,'fs_nominal',4000,'orders',[4,4], ...
    'source','Bai & Lee (1997) IEEE Trans. SAP', ...
    'desc','4th-order headphone plant: 4 zeros (2 RHP), 4 poles (2 RHP), k=0.3921');
save(fullfile(script_dir, 'syn_Bai1997_4th.mat'), 'model');
fprintf('  [3/7] syn_Bai1997_4th.mat\n');

%% Model 4: Jafari JVC 2017 — 高阶示例 (6阶, 连续)
G0 = (s*(s^2+4)*(s-0.8)*(s+1.4)) / ((s+0.5)^3*(s+2)^2*(s+3));
model = struct('name','JafariJVC2017_Example','type','syn_tf', ...
    'domain','continuous','G0',G0,'fs_nominal',10000,'orders',[5,6], ...
    'source','Jafari & Ioannou (2016) JVC', ...
    'desc','6th-order continuous example: 5 zeros (RHP at s=0.8), 6 poles');
save(fullfile(script_dir, 'syn_JVC2017_6th.mat'), 'model');
fprintf('  [4/7] syn_JVC2017_6th.mat\n');

%% Model 5: Carmona 2000 — 管道 ANC (7阶, 离散, Fs=2000)
A_coeffs = [1 -1.3941 -0.0389 1.2131 -1.1895 0.0430 1.0517 -0.6267];
B_coeffs = [0.0304 0.0709 -0.0947 -0.0170 -0.0104 -0.0787 0.0414 0.0380 0.0250 0.0366 -0.0584 0.0540 -0.0862 -0.6267];
d_delay = 6;
Fs_c = 2000;
G0_c = tf(B_coeffs, A_coeffs, 1/Fs_c);
model = struct('name','Carmona2000_Pipeline','type','syn_tf', ...
    'domain','discrete','G0',G0_c,'Ts',1/Fs_c,'fs',Fs_c, ...
    'A_coeffs',A_coeffs,'B_coeffs',B_coeffs,'d_delay',d_delay,'orders',[7,13], ...
    'source','Carmona & Alvarado (2000) ASME', ...
    'desc','7th-order discrete pipeline ANC model, d=6, Fs=2000Hz');
save(fullfile(script_dir, 'syn_Carmona2000_7th.mat'), 'model');
fprintf('  [5/7] syn_Carmona2000_7th.mat\n');

%% Model 6: 质量-弹簧-阻尼器 (2阶 SS, 连续)
m=1; k_s=1; b_s=0.5;
A_ss = [0 1; -k_s/m -b_s/m];
B_ss = [0; 1/m];
C_ss = eye(2);
D_ss = 0;
sys_ss = ss(A_ss, B_ss, C_ss, D_ss);
model = struct('name','MassSpringDamper','type','syn_ss', ...
    'domain','continuous','sys',sys_ss, ...
    'A',A_ss,'B',B_ss,'C',C_ss,'D',D_ss,'orders',[2], ...
    'source','Standard textbook model', ...
    'desc','2nd-order mass-spring-damper: m=1, k=1, b=0.5');
save(fullfile(script_dir, 'syn_MassSpringDamper_2nd.mat'), 'model');
fprintf('  [6/7] syn_MassSpringDamper_2nd.mat\n');

%% Model 7: RST 教学模型 (2阶, 离散, Ts=1)
B_rst = [0.2, 0.15];
A_rst = [1, -1.2, 0.45];
d_rst = 1;  Ts_rst = 1;
G0_rst = tf(B_rst, A_rst, Ts_rst);
model = struct('name','RST_Toy','type','syn_tf', ...
    'domain','discrete','G0',G0_rst,'Ts',Ts_rst,'fs',1/Ts_rst, ...
    'B_poly',B_rst,'A_poly',A_rst,'d_delay',d_rst,'orders',[2,1], ...
    'source','MOEA benchmark (demo4_Robust)', ...
    'desc','2nd-order discrete toy model: B=[0.2,0.15], A=[1,-1.2,0.45], d=1, Ts=1');
save(fullfile(script_dir, 'syn_RSTtoy_2nd.mat'), 'model');
fprintf('  [7/7] syn_RSTtoy_2nd.mat\n');

fprintf('\n========== 全部导出完成 ==========\n');
fprintf('文件位置: %s\n', script_dir);
