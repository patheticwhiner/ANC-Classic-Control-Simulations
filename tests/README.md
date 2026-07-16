# Cylinder 1 dm 控制器仿真

本目录以五个可直接阅读和逐段运行的 MATLAB 主实验脚本为主。每个脚本按照仓库外
原始 Demo 的风格使用 `%% Part` 组织模型加载、控制器参数、T1/T2 仿真和结果汇总。
100–500 Hz 宽带 T3 是独立失败边界实验，不混入 fixed/adaptive 主比较。

## 用户入口

```text
tests/
├── demo1_cylinder1dm.m       # 低阶 RST + Youla-Q RLS
├── demo2_cylinder1dm.m       # 增广 LQG + Q-RLS/LMS
├── demo3_cylinder1dm.m       # F(z)+FIR + Filtered-X LMS
├── demo4_cylinder1dm.m       # eMOPSO 灵敏度整形 RST
├── demo5_cylinder1dm.m       # Marino-Tomei 频估器 vs IMC-FxLMS (纯反馈)
├── run_cylinder1dm_stage.m   # T1/T2 校准、冻结评价和主报告
├── render_cylinder1dm.m      # 使用冻结结果重新生成主实验或宽带图像
├── test_display_options.m    # 图窗静默/显示开关
├── test_demo_overview.m      # 综合图自适应状态面板回归测试
├── test_demo5_diagnostics.m  # Marino-Tomei 符号与离散 Case A 边界回归测试
└── startup.m                 # 路径初始化
```

`internal/` 保存五个 Demo 控制器、信号加载/生成、指标计算和报告生成函数。
调参与对照实验是 `run_cylinder1dm_stage.m` 的局部函数，
单一调用方的辅助逻辑（RST 设计、限幅等）都内联在使用它们的文件里，
不存在额外的公共包装层。

## 直接运行 Demo

在 MATLAB 中执行：

```matlab
run('tests/startup.m');
demo1_cylinder1dm
demo2_cylinder1dm
demo3_cylinder1dm
demo4_cylinder1dm
demo5_cylinder1dm
```

也可以打开任一 Demo，只运行需要检查的 `%% Part`。控制器参数直接写在对应场景前，
例如 Demo2 T1 的 `R_lqr=1e-4`，不需要再查找另一层默认参数。

脚本使用冻结评价信号，主实验当前报告口径为：

| Demo | T1 fixed | T1 adaptive | T2 adaptive |
|---|---:|---:|---:|
| Demo1 | 27.33 dB | 38.97 dB | 36.58 dB |
| Demo2 | 36.18 dB | 36.59 dB | 36.16 dB（normalized-LMS） |
| Demo3 | 36.72 dB | 39.13 dB（IMC-FxNLMS） | 23.16 dB（IMC-FxNLMS） |
| Demo4 | 21.94 dB，需求 3.058 | 39.96 dB，需求 3.338 | 38.85 dB，需求 3.730 |
| Demo5 Marino-Tomei | 25.75 dB，需求 3.304 | 11.21 dB，需求 3.248 | -0.00 dB，Case A 不适用 |

## Demo5 — Marino-Tomei vs IMC-FxLMS (纯反馈)

Demo5 比较两种纯反馈（无参考麦克风）ANC 方法：
- **Marino-Tomei (2016)**：自适应内模在线估计窄带频率并抵消，仅适用固定频率正弦
- **IMC-FxLMS**：经典反馈 ANC 基线，IMC 结构重构内部参考 + filtered-x NLMS

主阶段默认统一使用 `LMS_20260713_*_N16_fs2k.mat` 的测量 FIR 主/次级路径。
冻结数字必须来自同一路径族的当前阶段结果，旧 ARMAX 结果不与 FIR 结果混用。

正式 FIR 实验使用 `output_timing='updated'`，与参考 Demo 一致：当前误差样本更新
内模状态后生成本拍控制命令，该命令经显式一拍次级路径延迟后作用于误差麦克风。
`previous` 仅保留为 ARMAX 离散相位边界诊断，不参与正式 FIR 候选选择。

T1 (357 Hz) 用于检验定频 Case A 是否真的满足抑制与执行器门限；成功状态始终从
冻结评价结果动态生成，不能由 Case A 适用性预设。
T2 (300→420 Hz 扫频) 的离散有效次级响应跨过零点，Marino-Tomei Case A 因固定符号条件不适用；
IMC-FxLMS 仍作为纯反馈基线运行。T3 (宽带) 不纳入 Demo5 主阶段。

### Demo5 实现核对结论

本次核对将原始 `demo5_MarinoTomei` 与 `tests/internal/controller_demo5` 放在同一条合成次级路径和双音设置下进行对照。内模状态更新、反馈符号和频率自适应律能够复现原始 Demo 的抑制效果，说明核心 Marino-Tomei 公式没有被实现成反号或失效版本。

此前 `tests/` 的正式候选使用 `output_timing='previous'`，而原始 Demo 使用更新后的内模状态生成本拍输出；这会引入额外一拍控制滞后。统一测量 FIR 后，正式候选改为 `output_timing='updated'`，`previous` 仅保留为 ARMAX 离散相位边界诊断。

此前可读脚本在没有兼容冻结 MAT 时会回放一个已知无效的 adaptive fallback，导致控制量接近零并产生“没有降噪”的表象；fallback 现已与当前校准候选一致。

当前可复现结论是：T1 fixed 达到 25.75 dB 且满足需求/限幅约束；T1 adaptive 只有 11.21 dB；T2 因 Case A 有效次级响应跨零而不适用。Demo5 因此是“固定窄带版本局部成功”，不是完整自适应扫频控制器成功。

```matlab
run('tests/startup.m');
demo5_cylinder1dm
```

Demo5 进入主阶段的机器可读冻结结果，但在解释中与 Demo1-4 分开：
Marino-Tomei 的 Case A 适用性和纯反馈 IMC-FxLMS 的信息结构不同，不能作为简单排行榜。
主阶段会保留 Demo5 T1/T2 的失败或不适用状态，并生成 `DEMO5_CYLINDER1DM_REPORT.md`。

Demo5 的离散诊断可单独运行：

```matlab
run('tests/test_demo5_diagnostics.m');
```

该回归测试同时检查论文反馈符号、反向符号、357 Hz 的离散有效 Case A 比值，
以及 T2 扫频的符号跨越边界。

## 完整阶段运行

需要重新执行校准、冻结评价和报告生成时运行：

```matlab
run('tests/startup.m');
stage = run_cylinder1dm_stage();
```

主要产物：

- `DEMO1234_STAGE_REPORT.md`
- `DEMO1_CYLINDER1DM_REPORT.md`
- `DEMO2_CYLINDER1DM_REPORT.md`
- `DEMO3_CYLINDER1DM_REPORT.md`
- `DEMO4_CYLINDER1DM_REPORT.md`
- `DEMO5_CYLINDER1DM_REPORT.md`
- `output/cylinder1dm_2k/demo1234_stage/demo1234_stage_results.mat`

每个独立报告还引用对应的 `figures/cylinder1dm_2k/demoN/demoN_analysis.png`，
用于联合解释 T1 fixed/adaptive、T2 局部分箱、自适应状态和控制需求。

Demo4 每个候选需要一次 eMOPSO 优化，完整阶段耗时明显长于以前的 Demo1-3 版本。

## 独立宽带边界运行

宽带设计不参与上述主阶段。单独执行：

```matlab
run('tests/startup.m');
broadband = run_cylinder1dm_stage('broadband');
```

该入口只运行 T3 校准与冻结评价，生成：

- `BROADBAND_CYLINDER1DM_REPORT.md`
- `output/cylinder1dm_2k/broadband_stage/broadband_stage_results.mat`
- `output/cylinder1dm_2k/broadband_stage/broadband_tuning.csv`
- `output/cylinder1dm_2k/broadband_stage/broadband_evaluation.csv`

当前四种宽带设计均无效：Demo1、Demo2、Demo3、Demo4 的冻结评价分别为
0.53、0.12、-4.32、3.41 dB，其中 Demo1/3/4 还违反控制需求或限幅约束。

## 图像显示

默认显示并保留图窗（同时保存 PNG）：

```matlab
render_cylinder1dm();
```

宽带图像独立重绘：

```matlab
render_cylinder1dm('broadband');
```

批量/无桌面运行时建议切换为静默保存：

```matlab
test_display_options('silent');
```

恢复默认显示模式执行 `test_display_options('reset')`。设置保存在当前
MATLAB 会话中，重启后回到默认。
