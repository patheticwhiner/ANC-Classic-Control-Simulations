# MATLAB 仿真与 FPGA ANC 模块结合实验计划

## 1. 目标与边界

本计划把当前已经冻结的 MATLAB 结果，逐步映射到
`/home/dcol/Projects/FPGA/DCOL-ANC-Sources` 的 SystemVerilog RTL 和板级音频链路。
第一阶段只做只读分析、系数导出、bit-accurate 仿真和 controller-in-the-loop；不修改 FPGA 仓库，
也不把“RTL 能综合”当成“ANC 算法已在硬件有效”。

当前 MATLAB 主实验口径：

| 方法 | T1 fixed | T1 adaptive | T2 adaptive | 宽带 T3 |
|---|---:|---:|---:|---|
| Demo1 | 30.66 dB，成功 | 39.56 dB，成功 | 38.28 dB，成功 | 独立实验，无效 |
| Demo2 | 35.38 dB，成功 | 35.58 dB，成功 | 36.99 dB，成功 | 独立实验，无效 |
| Demo3 | 37.50 dB，成功 | 39.47 dB，成功 | 23.71 dB，成功 | 独立实验，无效 |
| Demo4 | 25.62 dB，需求超限，失败 | 39.67 dB，需求超限，失败 | 38.57 dB，需求超限，失败 | 独立实验，无效 |

T2 的“成功”仍需披露局部边缘失效：Demo2 最差局部约 4.53 dB，Demo3 最差局部约 13.55 dB。
硬件阶段必须保留同一套抑制、控制需求、饱和和稳定性指标，不能只比较单一 RMS。

## 2. FPGA 仓库现状

### 2.1 已有可复用模块

| 能力 | 现有实现 | 关键接口/假设 |
|---|---|---|
| 定点 FIR | `hdl/rtl/fir/fir_core.sv`、`fir_reload.sv`、`fir_update.sv` | `in/in_valid` 到 `out/out_valid`；多 lane；系数可重载或流水更新 |
| 定点 IIR | `hdl/rtl/iir/iir_core.sv`、`iir_reload.sv` | 分开存放 B/A taps；同样是 valid 驱动的串行样本接口 |
| FxLMS | `hdl/rtl/algorithms/FxLMS.sv` | 外部 `ref_in`、`error_in`；次级路径 FIR；自适应 FIR |
| FxNLMS | `hdl/rtl/algorithms/FxNLMS.sv` | 外部源参考和误差，结构上最接近 Demo3 的源参考 filtered-x 更新 |
| IMC-FxLMS | `hdl/rtl/algorithms/IMC_FxLMS.sv` | 用 `error - filtered_out` 形成内部参考，和 MATLAB Demo3 的源参考结构不同 |
| 采样率转换 | `hdl/rtl/sample_rate/fir_decimator.sv`、`fir_interpolator.sv` | 现有 48 kHz ↔ 2 kHz、24 倍、381-tap 低通路径 |
| 系数加载 | SCPI/UART、SD loader | FIR 支持 binary block；IIR 默认从 SD 扇区加载 |
| 行为级黄金模型 | `matlab/*_step.m`、`fxlms_float_sim.m` | 逐时钟模型，可用于 RTL 对拍，但当前参数需先校准 |

### 2.2 时钟和采样率事实

- 板级 I2S 默认 24-bit、48 kHz。
- `clk_main` 为 73.728 MHz，理论上每个音频样本为 **1536 个主时钟周期**，不是 MATLAB 脚本当前写的 1600。
- 已有 `low_rate_loopback_top.sv` 使用 `fir_decimator`/`fir_interpolator` 做 48 kHz ↔ 2 kHz，
  因此第一条落地路线应优先把控制器运行在 2 kHz，再通过这条链路回到音频速率。
- 381-tap 低通系数是 Q1.15 资源，不能自动视为 ANC 控制器系数；必须分别验证幅频、群延迟和通带边缘。

### 2.3 定点格式风险

当前 MATLAB 冻结控制器系数范围已经超过 Q1.23 的 `[-1,1)` 范围：

- Demo1 `R0/S0` 最大绝对值约 2.37。
- Demo4 `R0/S0` 最大绝对值约 3.22。
- Demo3 固定 FIR 最小系数约 -4.49，范数约 4.59。
- Demo2 的 LQG 增益最大约 4.02。

因此，第一版导出应采用显式的有符号格式，例如数据保持 Q1.23，控制器系数先用 24-bit、19 fractional bits
（可表示约 `[-16,16)`）作为候选，再按误差和资源重新选择；不能把所有系数都无条件量化成 Q1.23。
`FxLMS/FxNLMS` 内部自适应更新目前按 `MU_WIDTH + 2*DATA_WIDTH` 计算系数位宽，默认会达到 72 bit，
需要单独做 BRAM/DSP 资源评估和“内部高精度、写回缩放”的方案比较。

## 3. MATLAB 控制器到 FPGA 的映射

| MATLAB 方法 | 第一版 FPGA 映射 | 是否直接等价 | 需要补充的模块 |
|---|---|---|---|
| Demo3 fixed F(z)+FIR | `iir_reload` 实现 F(z)，`fir_reload` 实现 8-tap 控制 FIR；次级路径复用 FIR | 近似直接，需统一延迟和符号 | 系数导出、IIR B/A 符号验证、控制输出缩放 |
| Demo3 IMC-FxNLMS | 优先检查 `FxNLMS.sv` 的外部源参考端口；`IMC_FxLMS.sv` 只能作为结构对照 | **不直接等价** | 源参考/误差对齐、filtered-x 延迟、归一化和限幅感知更新 |
| Demo1 RST + Youla-Q | RST 的 R/S 多项式可由 IIR/FIR 级联实现；Youla-Q 先离线固定，不先做在线 RLS | 固定部分可实现，自适应部分暂不等价 | RST 递推封装、Q 更新资源和稳定性保护 |
| Demo2 LQG | 先把冻结 LQG 闭环化为等价 IIR/状态空间传递函数 | 固定版本可近似，在线版本不直接支持 | 状态空间 MAC/观测器，或离线 IIR realization |
| Demo4 eMOPSO RST | 仅作为失败对照，不进入首轮硬件验收 | 不建议首轮实现 | 先解决 MATLAB 中控制需求超限，再评估硬件价值 |

信息结构必须保持一致：Demo3 使用源参考；Demo1/2/4 是反馈或内部重构。不能把 `IMC_FxLMS` 的内部参考
结构描述成 MATLAB Demo3 的源参考算法，也不能把不同信息结构的硬件结果做绝对排行榜。

## 4. 分阶段实验

### P0：接口和时序契约冻结

1. 将 MATLAB 周期模型中的 `CLK_PER_SAMPLE` 改为由 `CLK_FRE/FS` 计算，默认得到 1536；同时保留一个参数化测试用于非整数比率。
2. 为每个模块记录：输入 valid 时刻、输出 valid 时刻、样本延迟、吞吐率、复位清空周期、系数加载完成条件。
3. 用一个冲激和一个单频正弦测量 FIR、IIR、decimator、interpolator 的实际延迟，并把延迟写入 MATLAB golden model。
4. 明确 ANC 闭环的样本对齐：`d(n)`、`u(n)`、次级路径输出、`e(n)` 和 filtered-x 参考必须使用同一时间索引。

**通过门槛**：所有模块的 valid 序列可由 MATLAB `*_step.m` 重现；对齐误差不得超过一个明确记录的 pipeline 延迟。

### P1：系数导出和定点误差基准

1. 从冻结 `demo1234_stage_results.mat` 导出 Demo3 fixed 的 `F_num/F_den` 和 FIR taps；从模型结果导出次级路径 FIR。
2. 对每组系数生成：浮点 CSV、定点整数 CSV、`readmemh` 文件和量化误差摘要。
3. 对每个候选 fractional bits（例如 15、19、23）计算：系数最大误差、频响幅相误差、闭环极点变化、控制输出峰值变化。
4. 对 IIR 分母系数单独检查符号约定。`iir_core` 的实现是直接累加 B 与 A 历史项，MATLAB 导出必须匹配该递推式，不能直接假设 `a` 是 Control System Toolbox 的分母符号。

**通过门槛**：量化后的单级频响在 100–500 Hz 的幅相误差、极点半径和控制需求均在预先登记的容差内，且无定点溢出。

### P2：FIR/IIR RTL 单元与 bit-accurate 对拍

1. 先复用 `sim/fir/tb_fir.sv`、`sim/iir/tb_iir.sv`，再增加 MATLAB 生成的随机、冲激、357 Hz 和 300→420 Hz 向量。
2. 覆盖系数加载、非法地址、最后一个 lane、复位后清空、`in_valid/out_valid` 间隔输入。
3. 对每个输出样本比较：RTL 原始整数、MATLAB 周期模型整数、浮点参考值；分离“时序错位”和“数值误差”。
4. 对 `fir_update` 增加系数更新地址连续/间断两种测试，确认三阶段 read-modify-write 不丢更新。

**通过门槛**：逐样本 bit-exact（允许已登记的 rounding/saturation 差异），所有 testbench 明确 PASS/FAIL，而不是只观察波形。

### P3：2 kHz controller-in-the-loop

1. 采用现有 48 kHz → 2 kHz decimator 和 2 kHz → 48 kHz interpolator 作为速率边界。
2. 先接入 Demo3 fixed：2 kHz 源参考、8-tap FIR、F(z) IIR、固定次级路径模型。
3. 在 MATLAB 中把 decimator/interpolator 的群延迟和缩放纳入同一闭环，再与 RTL 输出逐样本对拍。
4. 只在 T1 357 Hz 和 T2 300–420 Hz 做 controller-in-the-loop；宽带 T3 保持独立失败边界，不混入主验收。

**通过门槛**：T1 抑制、T2 总体抑制和最差局部抑制相对 MATLAB 冻结结果的差异有明确容差；控制需求不超过 4，硬限幅为 0。

### P4：自适应和资源评估

1. 在 P3 固定 FIR 通过后，再接 `FxNLMS.sv`；先使用短 tap（8/32）和低速仿真，确认 filtered-x 与 error 的 valid 对齐。
2. 比较三种实现：全精度 72-bit 更新、缩放后 24/32-bit 写回、块浮点/共享指数；报告 DSP、BRAM、Fmax 和收敛差异。
3. 在线更新必须加入 MATLAB 中已有的 norm limit、饱和感知更新和启动渐入，否则不能声称对应 Demo3 IMC-FxNLMS。
4. 只有当 `FxNLMS` 路径稳定后，才评估 `IMC_FxLMS`；两者分别报告源参考和内部参考信息结构。

### P5：Demo1/2/4 的硬件扩展

1. Demo1：先冻结 RST/Q 为固定系数 IIR/FIR，验证反馈闭环；在线 Youla-Q RLS 作为独立资源课题。
2. Demo2：先实现固定 LQG 的等价 IIR 或状态空间 MAC；观测器和矩阵乘法需要单独做吞吐率/资源预算。
3. Demo4：在 MATLAB 重新把控制需求压到 4 以内前，仅保留为失败对照，不进入板上成功验收。

### P6：hardware-in-the-loop

1. 通过现有 UDP/I2S 路径发送带序号的 `d/ref/e/u` 帧，记录 FPGA 侧 `out_valid`、丢包、FIFO 溢出和限幅计数。
2. 先用离线已知 secondary-path FIR，再接真实声学路径；真实路径辨识结果必须独立保存并重新验证。
3. 使用 T1、反向 T2、最后才是独立 T3；每次只冻结一组系数，禁止用硬件评价集重新调参。
4. 输出硬件版 CSV/PNG，沿用 MATLAB 报告的六项场景契约：目的、设置、比较、量化结果、解释、结论。

## 5. 首个可执行里程碑

推荐首个里程碑为 **Demo3 fixed + 2 kHz controller-in-the-loop**：

1. 生成 2 kHz 版 F(z)+8-tap FIR 和次级路径系数。
2. 选择并登记控制器系数格式（建议先试 24-bit/19 fractional，并保留 32-bit 参考）。
3. 完成 `fir_reload`、`iir_reload`、次级路径和采样率转换的 bit-accurate 对拍。
4. 在 T1 357 Hz 和反向 T2 上比较抑制、最差局部、控制需求、硬限幅和闭环稳定性。
5. 通过后再切换 `FxNLMS.sv` 自适应；不要从 `IMC_FxLMS.sv` 直接替代 MATLAB Demo3 的源参考结构。

## 6. 当前阻塞项和待确认问题

- 2 kHz 控制器究竟放在 FPGA 采样率转换之后，还是把全部控制器重离散化到 48 kHz；首轮按已有 24 倍转换链路执行。
- `CLK_PER_SAMPLE=1600` 的旧 MATLAB 参数需改为 1536 或由时钟/采样率自动计算。
- Demo3/Demo1/Demo2/Demo4 系数超出 Q1.23；必须先完成格式扫描，不能直接复用现有 Q1.23 假设。
- `FxNLMS` 的 72-bit 自适应系数和现有 `fir_update` 资源是否可接受，需综合报告确认。
- FPGA 仓库现有 testbench 已覆盖模块逻辑，但尚未证明当前 MATLAB Cylinder 1 dm 模型、采样率转换延迟和硬件定点系数三者已经统一；这正是 P0–P3 的验证目标。

## 7. 产物清单

- `fpga_plan/coefficients/<controller>_<format>.csv`
- `fpga_plan/coefficients/*.mem`
- `fpga_plan/vectors/{t1,t2,t3}/input_*.mem`
- `fpga_plan/rtl_compare/*.csv`
- `fpga_plan/reports/FPGA_BIT_ACCURATE_REPORT.md`
- `fpga_plan/reports/FPGA_CONTROLLER_IN_LOOP_REPORT.md`

计划依据：FPGA 仓库 `README.md`、`matlab/README.md`、FIR/IIR/算法 RTL、采样率转换 RTL 及其 testbench；
所有硬件结论仍需以实际仿真或板上冻结评价结果为准。
