function cal = rt_calibration()
%RT_CALIBRATION 台架标定增益（V ↔ 模型单位）—— 上板 Stage 0 后人工填写
%
%   ⚠ 下列值为占位默认（增益 1、无失调）。首次上板必须按
%   docs/MICROLABBOX_EXPERIMENT_PROTOCOL.md Stage 0/1 的标定流程实测后
%   更新，并在提交信息中注明标定日期与测量方法。
%
%   MicroLabBox RTI 块约定: ADC 输出 = V_in/10（±10 V → ±1），
%   DAC 输入 = V_out/10。模型内在 RTI 块外补 ×10 / ×0.1，使
%   IO 子系统边界统一为“物理伏特”。

cal = struct();

% --- 输入链: 误差麦克风电压 → 模型单位 e ---
% e_model = cal.e_per_volt * V_mic（含前置放大器增益的总折算）
cal.e_per_volt = 1.0;          % [模型单位/V] 占位，Stage 1 标定

% --- 输入链: 参考通道（demo3）---
cal.ref_per_volt = 1.0;        % [模型单位/V] 占位

% --- 输出链: 模型单位 u → DAC 电压 ---
% V_dac = cal.volt_per_u * u（功放增益在链路之外，见规程接线表）
cal.volt_per_u = 1.0;          % [V/模型单位] 占位，从低增益阶梯上探

% --- 安全 ---
cal.dac_volt_limit = 2.0;      % [V] DAC 输出绝对限幅（保守起步，逐步放开）
cal.watchdog_window_s = 0.25;  % [s] e 滑动 RMS 窗口
cal.watchdog_rms_limit = 2.0;  % [模型单位] e_rms 超此值锁零（Stage 2 校准）

% --- 通道分配（MicroLabBox ADC/DAC Class 1，通道号按接线表）---
cal.adc_ch_error = 1;
cal.adc_ch_ref = 2;
cal.dac_ch_control = 1;
end
