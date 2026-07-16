# MicroLabBox (DS1202) ANC 上板实验规程

> 人工执行文档。每个步骤含 **目的 / 命令 / 预期结果 / 异常处理** 四要素。
> 开发机（Linux）只负责生成模型与分析数据；**构建、下载、运行全部在
> dSPACE 主机（Windows）上由实验人员手动完成**。

对应模型: `dspace/models/rt_<demo>_<test>_<variant>.slx`（由
`build_rt_model` 在 dSPACE 主机上重新生成，换入 RTI 块）。
成功判据与仿真报告同口径：稳态抑制 ≥ 20 dB、控制需求 ≤ 4（模型单位）、
零饱和、无失稳 —— 不得只比较单一 RMS。

---

## 0. 硬件映射与接线表

| MicroLabBox 通道 | 连接 | 途径 | 备注 |
|---|---|---|---|
| ADC Class-1 CH1 | 误差麦克风 | 前置放大 → **抗混叠 LPF（fc ≈ 800 Hz，模拟）** → BNC | 2 kHz 采样，混叠必须模拟域解决 |
| ADC Class-1 CH2 | 参考传感器（仅 demo3） | 同上调理链 | 上板前须验证与误差通道相干性 |
| DAC CH1 | 控制扬声器 | **重建滤波器（fc ≈ 800 Hz）** → 功放 → 扬声器 | 功放增益预置最小 |
| DAC CH2（可选） | 扰动源扬声器 | 功放 → 扬声器 | 噪声工况回放（T1 音、T2 扫频） |

接地: 信号地单点接 MicroLabBox 模拟地；功放与调理电源分开走线。
**物理急停 = 功放静音开关**（模型内 PanicSwitch 只是软件门，见安全附录）。

> ⚠ 16 kHz 过采样 + 模型内抽取是抗混叠的备选方案（一期不做）；若外部
> 模拟滤波器不可得，先按备选方案改造模型再上板，不允许裸 2 kHz 采样。

## 1. dSPACE 主机构建（每次参数/模型更新后）

**目的**：在有 RTI1202 的机器上生成含真实 ADC/DAC 块的 RT 模型并编译。

**命令**（dSPACE 主机 MATLAB）：
```matlab
cd <repo>
run('dspace/startup.m');
export_frozen_params;      % 若 tests 冻结结果有更新
build_all;                 % rti_available()==true → 换入 DS1202 块
% 打开目标模型, Ctrl+B 构建下载（例）
open_system('dspace/models/rt_demo1_T1_fixed.slx')
```

**预期结果**：构建日志显示 `rti1202.tlc` 目标；ControlDesk 中出现实验；
变量树可见 `EnableSwitch / PanicSwitch / LatchReset / MasterGain` 与监视
信号 `e / u / u_demand / clipped / diag / e_rms / watchdog_latched / frame_clock`。

**异常处理**：
- `RtiBlockNotFound`：打开 `rtilib1202` 库确认实际块名，更新
  `build_rt_model.m::add_rti_block` 的正则。
- 构建报 turnaround 超限：见 Stage 0 的负载检查。
- codegen 对 `inv()` 告警（demo1/4 自适应 RLS）：允许告警；若报错，
  改用 rank-1 协方差更新前**必须**重跑离线等效性并在报告披露偏差。

## 2. Stage 0 — I/O 回环与安全机制自检（不接功放）

**目的**：标定端到端增益/偏置/延迟，验证看门狗与软件急停，测量任务负载。

**命令**：DAC CH1 用 BNC 短接回 ADC CH1（不接功放！）。运行
`rt_demo1_T1_fixed`，`EnableSwitch=0`，在 ControlDesk 用 Test 面板向
DAC 写 0.5 V、1 kHz 以下正弦。

**预期结果**：
- ADC 读数 = DAC 写入 ×(1±2%)，偏置 < 5 mV；回环延迟 ≈ 1–2 个采样周期
  （记录实测值）。
- Task turnaround time（RTI 任务信息变量）≪ 500 µs（2 kHz 周期），
  典型应 < 50 µs。
- 手动把 `PanicSwitch` 置 1 → DAC 输出立即为 0。
- 向 ADC 注入大幅信号使 `e_rms` 超过 `rms_limit` → `watchdog_latched=1`
  且输出锁零；`LatchReset` 置 1 再回 0 → 解除。

**异常处理**：
- 增益偏差 > 2%：把实测折算系数写入 `dspace/params/rt_calibration.m`
  （`e_per_volt` / `volt_per_u`），重新生成模型。
- turnaround 接近 500 µs：检查是否误开信号记录全量高频采样；必要时
  减少 ControlDesk 采集信号数。
- 看门狗不触发：确认 `window_samples`/`rms_limit` 已随标定更新。

## 3. Stage 1 — 次级路径台架复辨识

**目的**：确认台架当前声学路径与 `cylinder1dm_2k_secondary_fir`
（2026-07-13 实测）仍一致；不一致则以新模型重新导出控制器。

**命令**：接好功放（低增益）。用 DAC CH1 播放 100–600 Hz 带限噪声
（`signal_excitation/` 生成，或 ControlDesk 回放），同步记录 DAC 激励与
ADC CH1 响应 ≥ 60 s，导出 `.mat`。在开发机用与
`sysid_models/20260713_cylinder1dm/` 相同的 LMS 管线拟合 FIR(16) @ 2 kHz。

**预期结果**：新旧 FIR 在 100–500 Hz 幅频偏差 < 3 dB、相位偏差 < 30°，
拟合误差比与原始档案同量级（≈ 6%）。

**异常处理**：
- 偏差超限（扬声器/麦克风位置变动、温度变化）：把新 FIR 注册进
  `dataset/model_registry.m`，重跑 `tests/run_cylinder1dm_stage()` 冻结新
  参数 → `export_frozen_params` → `build_all` → 重新过离线等效性门禁，
  然后才允许继续 Stage 2。**禁止**在旧控制器参数上"试试看"。
- 相干性差（< 0.9）：检查信噪比、接线、混叠滤波器。

## 4. Stage 2 — 固定 RST（demo1）低增益上板

**目的**：首个闭环实验；建立开环基线，标定看门狗阈值，验证抑制方向正确。

**命令**：噪声源播放 T1 工况（357 Hz 音 + 低幅噪声，声压从低起步）。
模型 `rt_demo1_T1_fixed`：
1. `EnableSwitch=0` 录 ≥ 10 s 开环基线（`rig_openloop_T1_<date>.mat`）。
2. 用开环 `e_rms` 设置看门狗：`rms_limit ≈ 1.5 × e_rms_open`（写回
   `rt_calibration.m` 并在 ControlDesk 同步修改）。
3. `MasterGain=0.2` → `EnableSwitch=1`，观察 30 s。
4. 增益阶梯 0.2 → 0.5 → 1.0，每档观察 ≥ 30 s 并录波。

**预期结果**：每档 `e_rms` 单调下降或持平，`clipped` 恒 0，`u_demand`
峰值 < 4×MasterGain；满增益时 357 Hz 处可听/可测衰减（仿真参考:
27.3 dB，台架允许打折但应 ≥ 20 dB）。

**异常处理**（任一发生 → `PanicSwitch=1`，回退上一档增益）：
- `e_rms` 增大（正反馈）：检查扬声器/麦克风极性；若极性正确而方向仍错，
  怀疑次级路径相位变化 → 回 Stage 1。
- 持续 `clipped=1`：`volt_per_u` 过大或声压过高。
- 看门狗频繁误触发：阈值过紧，用当前工况基线重标。

## 5. Stage 3 — 自适应控制器逐个上板

**目的**：在固定 RST 验证过的同一工况链上，逐个引入自适应结构。

**顺序与要点**（每个都先 `MasterGain=0.2` 阶梯上探，录波 ≥ 10 s）：

| # | 模型 | 观察量（diag） | 特有前置检查 |
|---|---|---|---|
| 1 | `rt_demo1_T1_adaptive` | ‖Q‖ 收敛、e_post → 0 | — |
| 2 | `rt_demo4_T1_adaptive` / `rt_demo4_T2_adaptive` | 同上 | demo4 需求余量小（仿真 3.73），增益阶梯更保守 |
| 3 | `rt_demo3_T1_adaptive`（需 ADC CH2） | ‖θ‖ 收敛 | **先测参考↔误差通道相干性 > 0.9**（100–500 Hz） |
| 4 | `rt_demo2_T1_adaptive` / `rt_demo2_T2_adaptive` | ‖θ‖、y_Q | demo2 T1 内稳定性指标在仿真未过——只做短时（≤60 s）观察性实验，全程人守 |
| 5 | `rt_demo5_T1_fixed` → `rt_demo5_T1_adaptive` | √θ₁（×fs/2π = f̂ Hz） | 仅 T1 单音；T2 扫频禁止（Case A 不适用） |

**预期结果**：自适应量在 3 s 内收敛后冻结波动；抑制优于对应固定版。

**异常处理**：
- 自适应发散（‖Q‖/‖θ‖ 持续增长）：`EnableSwitch=0`（会复位控制器状态
  ——enable 上升沿即重新武装），必要时降 `MasterGain`。
- demo5 频率估计撞界（f̂ 到达 20/500 Hz 界）：初值偏差过大或 ε 过大，
  按冻结参数复核 ControlDesk 中的实际值。

## 6. Stage 4 — 对比战役与数据归档

**目的**：采集可与仿真报告逐列对比的最终数据。

**命令**：每个 (demo, test, variant) 在同一工况下录 ≥ 10 s：
信号集 `{t, e, u, u_demand, clipped, diag, frame_clock}` @ 2 kHz，
导出命名 `rig_<demo>_<test>_<variant>_<yyyymmdd>.mat`，连同当日开环
基线 `rig_openloop_<test>_<yyyymmdd>.mat` 拷回开发机
`dspace/output/rig/`。分析：

```matlab
run('dspace/startup.m');
r = analyze_rig_capture('dspace/output/rig/rig_demo1_T1_fixed_20260801.mat', ...
                        'dspace/output/rig/rig_openloop_T1_20260801.mat');
```

**预期结果**：脚本输出 supp_db / 需求峰值 / 饱和计数 / 丢帧数；
丢帧 = 0；台架 supp_db 与仿真差 < 10 dB 量级（声学失配预期内），
结论按同一成功判据（≥20 dB、需求 ≤4、零饱和）给出。

**异常处理**：
- 丢帧 > 0：降低 ControlDesk 采集负载或用平台 recorder 而非在线绘图。
- 台架抑制显著劣于仿真：先查 Stage 1 路径失配，再查标定增益，最后才
  调控制器 —— 禁止直接在台架上重调冻结参数。

## 7. 安全附录

1. **物理急停**：功放静音/断电开关必须在实验人员手边；软件
   `PanicSwitch` 与看门狗只是第一层。
2. **电压上限**：`dac_volt_sat` 出厂 ±2 V（`rt_calibration.m`），只允许
   在增益阶梯全程无削顶后逐步放开，任何时候 ≤ ±10 V。
3. **声压**：扬声器近场佩戴听力防护；连续单音实验 SPL ≤ 85 dB(A)。
4. **首次上电双人制**：Stage 0–2 需两人，一人操作 ControlDesk，
   一人守物理急停。
5. **检查单**（每次上电前）：接线核对 → 功放最小增益 → `EnableSwitch=0`
   确认 → 看门狗阈值已按当日基线标定 → 录波已配置。
6. **本仓库红线**：任何自动化脚本（含 CI/Claude）不得执行下载、烧录、
   运行目标硬件的操作；上述步骤仅限实验人员手动执行。
