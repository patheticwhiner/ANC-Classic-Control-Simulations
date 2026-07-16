function [gate, e_rms, latched_out] = anc_safety_gate_step( ...
        e, panic, latch_reset, window_samples, rms_limit) %#codegen
%ANC_SAFETY_GATE_STEP RT 安全门: e 滑动 RMS 看门狗 + 软件急停
%
%   gate ∈ {0,1} 乘在控制输出上:
%     - panic > 0.5（ControlDesk 可调）→ gate = 0（软件急停，非物理急停！
%       物理急停 = 功放静音，见实验规程安全附录）
%     - e 的 window_samples 滑动 RMS 超过 rms_limit → 锁零（latched），
%       直到 latch_reset 上升沿才解除
%
%   参数: window_samples（编译期常量，决定缓冲区尺寸）、rms_limit。
%   rms_limit 在上板 Stage 2 用开环基线标定（见规程）。

persistent buf idx acc latched prev_reset
if isempty(buf)
    buf = zeros(1, window_samples);
    idx = 1;
    acc = 0;
    latched = false;
    prev_reset = 0;
end

% ---- 滑动均方（环形缓冲 + 运行和）----
e2 = e * e;
acc = acc - buf(idx) + e2;
buf(idx) = e2;
idx = idx + 1;
if idx > window_samples
    idx = 1;
end
e_rms = sqrt(max(acc, 0) / window_samples);

% ---- 锁零/复位 ----
if e_rms > rms_limit
    latched = true;
end
if latch_reset > 0.5 && prev_reset <= 0.5
    latched = false;
    acc = 0;
    buf(:) = 0;
end
prev_reset = latch_reset;

if latched || panic > 0.5
    gate = 0;
else
    gate = 1;
end
latched_out = double(latched);
end
