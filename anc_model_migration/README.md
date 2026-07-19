# ANC model migration designer

该目录把 `tests/internal/controller_demo*.m` 中的控制器设计方法变成可选择模型、噪声和控制器方案的迁移工具。工作流参考 `dSPACE_ANC/ANC_intro/controller_design_gui.m` 的分层方式，但增加了完整的控制器设计、校准选参和留出评价。

## GUI 入口

```matlab
run('project_init.m');
addpath('anc_model_migration');
anc_migration_designer_gui;
```

面向部署的上层工具可以限制 GUI 中可选的控制器，并注入经过审查的
部署导出函数。例如 dSPACE 固定反馈 IIR 页面只显示两个固定 RST 入口：

```matlab
anc_migration_designer_gui( ...
    'ControllerIds', {'demo1_rst_fixed','demo4_emopso_rst'}, ...
    'DeploymentExporter', @export_fixed_iir_controller);
```

部署按钮只在留出评价 `migration.summary.passed` 为真时启用。迁移后端
本身仍保存完整 `migration` 证据；具体的 RTI 系数格式、拓扑和采样率
检查由注入的部署导出函数负责，避免把反馈 `R/S` 与同名的前馈滤波器
系数混用。

界面提供：

- 自动发现 `sysid_models/20260718` 中匹配的主路径/次级路径模型对；
- FIR/ARMAX 与实际采样率（包括名为 `fs4k`、实际为 3999 Hz 的模型）选择；
- 定频、扫频、带限噪声、白噪声、多正弦和音频文件输入；
- Demo1–5 固定/自适应控制器及 IMC-FxLMS 方案；
- Quick/Standard 两级自动调参；
- 与实验设置并排的“控制器参数”分区：可关闭候选扫描，并按控制器逐项勾选、验证和覆盖手动参数；
- 校准集候选选择、独立随机种子或反向扫频的留出评价；
- MAT 结果、调参 CSV、残差/PSD/候选/控制需求图。

Marino–Tomei 的输入类型在 UI 和后端中受理论信息结构约束；不兼容组合不会运行。eMOPSO 方案需要离线优化，耗时明显长于其他方案。

## 非 GUI 入口

```matlab
addpath('anc_model_migration');
cfg = anc_migration_config();
catalog = scan_anc_migration_models(cfg.model_root);
pairs = pair_anc_migration_models(catalog);

cfg.models.primary_file = pairs(1).primary.file;
cfg.models.secondary_file = pairs(1).secondary.file;
cfg.noise.type = 'linear_chirp';
cfg.design.controller_id = 'demo2_lms';
cfg.output.plot = true;
migration = run_anc_migration(cfg);
```

脚本也可使用与 GUI 相同的手动参数后端：

```matlab
cfg.design.auto_tune = false;
cfg.design.manual_params = struct('nQ', 48, ...
    'adaptation_gain', 0.05, 'theta_norm_limit', 8);
```

开启自动调参时，手动值覆盖每个候选的对应字段；关闭自动调参时，仅保留一个采样率相关的默认候选，再施加手动覆盖。

`migration.calibration` 保存全部候选及选择证据，`migration.evaluation` 保存冻结超参数在留出信号上的结果。评价失败、饱和或控制需求超限不会被隐藏，`migration.summary.passed` 只在抑制、稳定性、需求和限幅门槛同时满足时为真。

## 算法实现文档

- TeX 源文件：[`docs/anc_model_migration_algorithm_guide.tex`](docs/anc_model_migration_algorithm_guide.tex)
- 已编译 PDF：[`docs/anc_model_migration_algorithm_guide.pdf`](docs/anc_model_migration_algorithm_guide.pdf)
- 控制系统框图 TeX：[`docs/anc_model_migration_block_diagrams.tex`](docs/anc_model_migration_block_diagrams.tex)
- 控制系统框图 PDF：[`docs/anc_model_migration_block_diagrams.pdf`](docs/anc_model_migration_block_diagrams.pdf)

算法指南按可执行源码说明 12 个控制器入口的公式、信息结构、实现步骤、自动候选缩放、手动参数语义以及迁移到其他辨识模型时的适用条件。控制系统框图文档说明框图绘制方法，并为每个入口给出独立框图。使用 XeLaTeX 编译：

```powershell
latexmk -xelatex docs/anc_model_migration_algorithm_guide.tex
latexmk -xelatex docs/anc_model_migration_block_diagrams.tex
```

## 延迟约定

ARMAX 使用模型 `orders(4)` 和分子中的前导零。LMS MAT 文件只保存 `s.w/fs/mu`，其文档中的 `nk` 作为辨识元数据读取；控制仿真另外显式增加 `cfg.models.control_delay_samples`（默认一拍）。两类延迟分别保存，避免把辨识延迟和控制计算延迟混为一谈。

## 回归测试

```matlab
addpath('anc_model_migration');
test_anc_model_migration();
```

常规测试检查 8 个模型、4 个匹配模型对、12 个控制器的实际调度、2 kHz FIR 端到端运行，以及不可见 GUI 中的 4 个模型对、12 个控制器和 6 类噪声选择项。

完整的 4 模型对 × 12 控制器组合矩阵可运行：

```matlab
matrix = test_anc_model_migration_matrix();
```

该矩阵要求所有 48 个组合都能执行并返回结果证据；它不把“所有方法都必须达到 20 dB”误当成迁移成功条件，性能是否通过仍由每个结果的 `summary.passed` 单独给出。
