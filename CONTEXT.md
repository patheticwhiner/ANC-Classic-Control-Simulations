# ANC Classic Control Simulations

主动噪声控制（ANC）经典控制方法的对比实验平台——用统一被控对象和测试基准评估四种控制器设计方法（Youla 极点配置、LQG、Jafari F(z)+FIR+NLMS、ε-MOPSO RST）的性能。

## Language

**被控对象 (Plant)**:
由 ARMAX(30,30,30,22) 辨识模型描述的实测声学管道，48kHz 采样，d=22 纯延迟。
_Avoid_: 系统、传递函数（仅指具体数学模型时不用）

**抑制 (Suppression)**:
闭环输出相对开环输出的 RMS 比，以 dB 表示：`20·log₁₀(RMS(y_open) / RMS(y_closed))`。正值表示噪声降低。
_Avoid_: 降噪量、衰减（attenuation 指声压级差，不同概念）

**稳态抑制 (Steady-state Suppression)**:
仿真时间 `[3, Tsim]` 的后 80% 窗口内计算的抑制 dB。固定控制器从第 0.5 秒开始取窗。
_Avoid_: 最终抑制、收敛后抑制

**收敛时间 (Convergence Time)**:
自适应控制器从启动到抑制达到稳态抑制 50% 的时间。固定控制器填 0。
_Avoid_: 响应时间、建立时间

**测试信号 (Test Signal)**:
Phase 2 预生成并存储的三组标准扰动信号，RMS 归一化到 0.8。
_Avoid_: 扰动、干扰（interference 含义不同）

**固定控制器 (Fixed Controller)**:
参数不随时间变化的控制器。包括 demo1 的 R₀/S₀、demo2 的 LQG、demo3 的 FIR K(z)。
_Avoid_: 静态控制器、非自适应控制器

**自适应控制器 (Adaptive Controller)**:
参数在线更新的控制器。包括 demo1 的 Youla-Q RLS、demo3 的 NLMS。
_Avoid_: 在线控制器、学习控制器

**控制器设计频率 (Design Frequency)**:
固定控制器优化时目标抑制的频率点，取被控对象第一共振峰（~334 Hz）。
_Avoid_: 调谐频率、中心频率

**频率鲁棒性 (Frequency Robustness)**:
固定控制器在干扰频率偏离设计频率时维持抑制的能力。由 T2 扫频测试评估。
_Avoid_: 频率敏感度

**有效抑制带宽 (Effective Suppression Bandwidth)**:
灵敏度 |S(f)| < −10 dB 的频率范围。受限于延迟相位滞后（>500 Hz 不可控）。
_Avoid_: 控制带宽

**开环响应 (Open-loop Response)**:
无控制时被控对象对测试信号的输出 y_open = P(z)·d。预计算并存储在测试信号文件中。
_Avoid_: 未控响应、被动响应

**测试基准 (Test Suite)**:
`tests/` 目录下的标准化评估框架：信号生成脚本、指标计算函数、控制器接口规范。
_Avoid_: 测试框架、验证平台（validator 指 GUI 工具，不同概念）

**控制器验证器 (Controller Validator)**:
根目录的 `controller_validator.m`——RST 专用 GUI 验证工具（频域/时域/鲁棒性标签页）。与 tests/ 测试基准互补，不做批量评估。
_Avoid_: 测试工具、评估器
