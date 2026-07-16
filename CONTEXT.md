# ANC Classic Control Simulations

主动噪声控制（ANC）经典控制方法的对比实验平台——用统一被控对象和测试基准评估四种控制器设计方法（Youla 极点配置、LQG、Jafari F(z)+FIR+NLMS、ε-MOPSO RST）的性能。

## Language

**被控对象 (Plant)**:
由 ARMAX(30,30,30,22) 辨识模型描述的实测声学管道，48kHz 采样，d=22 纯延迟。
_Avoid_: 系统、传递函数（仅指具体数学模型时不用）

**抑制 (Suppression)**:
闭环输出相对开环输出的 RMS 比，以 dB 表示：`20·log₁₀(RMS(y_open) / RMS(y_closed))`。正值表示噪声降低。
_Avoid_: 降噪量、衰减（attenuation 指声压级差，不同概念）

**成功控制器设计 (Successful Controller Design)**:
控制器及闭环内部稳定、控制量不触发饱和，并在至少一个被控对象主导频点实现不低于 20 dB 的抑制。不满足时应报告为算法不适用或设计失败，而不是以完成仿真代替成功。
_Avoid_: 仿真跑通、结果有效、有效果

**稳态抑制 (Steady-state Suppression)**:
仿真时间 `[3, Tsim]` 的后 80% 窗口内计算的抑制 dB。固定控制器从第 0.5 秒开始取窗。
_Avoid_: 最终抑制、收敛后抑制

**收敛时间 (Convergence Time)**:
自适应控制器从启动到抑制达到稳态抑制 50% 的时间。固定控制器填 0。
_Avoid_: 响应时间、建立时间

**测试信号 (Test Signal)**:
Phase 2 预生成并存储的三组标准扰动信号，RMS 归一化到 0.8。
_Avoid_: 扰动、干扰（interference 含义不同）

**校准集 (Calibration Set)**:
仅用于为每种控制器独立选择参数的共享信号集合，其性能不得作为最终算法结论。
_Avoid_: 训练集、测试集、调参结果

**评价集 (Evaluation Set)**:
不参与参数选择、仅用于最终性能报告的共享信号集合，包括工程目标频带评价和近 Nyquist 压力场景。
_Avoid_: 校准集、验证集、调参集

**场景专用控制器设计 (Scenario-specific Controller Design)**:
针对 T1、T2、T3 等测试信号类别分别设计并冻结控制器；每套控制器只能使用对应类别的校准子集选择参数，随后在同类但未参与调参的评价信号上测试。
_Avoid_: 针对评价信号调参、每次运行重新调参、单一通用控制器

**固定控制器 (Fixed Controller)**:
参数不随时间变化的控制器。包括 demo1 的 R₀/S₀、demo2 的 LQG、demo3 的 FIR K(z)。
_Avoid_: 静态控制器、非自适应控制器

**自适应控制器 (Adaptive Controller)**:
参数在线更新的控制器。当前 Cylinder 1 dm 主实验包括 demo1 的 Youla-Q RLS、
demo2 的 normalized-LMS（filtered-x Q-RLS 保留为失败对照）、demo3 的源参考 IMC-FxNLMS；
宽带 T3 使用的旧结构在独立边界报告中单列。
_Avoid_: 在线控制器、学习控制器

**工程目标频带 (Engineering Target Band)**:
控制器性能结论主要适用的物理频率范围；Cylinder 1 dm 基准取 100–500 Hz。约 880 Hz 的模型峰值属于近 Nyquist 压力测试，不作为成功设计判据。
_Avoid_: 全频段、有效抑制带宽

**归一化执行器约束 (Normalized Actuator Constraint)**:
测试基准中的控制量采用归一化单位，统一硬限幅为 `u ∈ [-5, 5]`；控制器调参阶段要求未限幅控制峰值不超过 4，以保留 20% 裕量。最终评价出现任何限幅时，不认定该控制器为成功设计。
_Avoid_: 真实电压上限、控制量越大越好、仅报告限幅后的结果

**控制器设计频率 (Design Frequency)**:
每种控制器在工程目标频带内独立选择的主要优化频点，不固定继承其他被控对象的共振频率。
_Avoid_: 固定 334 Hz、调谐频率、中心频率

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
