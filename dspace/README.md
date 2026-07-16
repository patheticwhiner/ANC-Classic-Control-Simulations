# dspace/ — Simulink + dSPACE MicroLabBox 部署子项目

把 `tests/` 下已冻结的 ANC 控制器（Demo1/2/3/4/5，2 kHz 逐样本离散实现）迁移到
Simulink，最终部署到 dSPACE MicroLabBox (DS1202)。三层产物：

1. **离线闭环模型** `models/offline_*.slx` — 被控对象用 cylinder1dm 实测 FIR 路径，
   与 MATLAB 脚本仿真逐样本对拍（等效性门禁）。
2. **实时模型** `models/rt_*.slx` — 同一控制器库链接，I/O 换成 RTI DS1202 ADC/DAC
   （无 RTI 环境自动降级为 Inport/Outport 占位）。
3. **上板规程** `docs/MICROLABBOX_EXPERIMENT_PROTOCOL.md` — 接线、标定、安全、
   分阶段 bring-up（人工执行）。

## 原则

- **所有 .slx 由脚本生成**（`build/`），生成脚本是唯一提交物；`.slx`/`.mat` 均被
  `.gitignore` 忽略，任何人重跑脚本即可再生。
- **冻结参数不手抄**：`params/export_frozen_params.m` 从
  `tests/output/cylinder1dm_2k/demo1234_stage/demo1234_stage_results.mat` 出发，
  复用 `tests/internal/controller_demoN.m` 的设计代码路径收割系数。
- **等效性是硬门禁**：固定控制器 `max|Δe| ≤ 1e-12`，自适应 `≤ 1e-9`；
  超差按移植 bug 处理，不放宽容差。证据见 `SIMULINK_EQUIVALENCE_REPORT.md`。

## 再生成流程（本机 / 任何有 Simulink 的机器）

```bash
matlab -batch "run('dspace/startup.m'); export_frozen_params; build_all"
matlab -batch "run('dspace/startup.m'); verify_offline_equivalence"
```

## dSPACE 主机（Windows）构建

1. 拉取仓库，确认已安装 MATLAB + Simulink + Simulink Coder + dSPACE RTI1202。
2. 运行上面同样的两条命令（RTI 存在时 `build_rt_model` 自动换入 DS1202 ADC/DAC 块，
   并设置 `rti1202.tlc` 代码生成目标）。
3. 打开 `models/rt_demo1_fixed.slx` → Ctrl+B 生成并下载到 MicroLabBox。
4. 之后严格按 `docs/MICROLABBOX_EXPERIMENT_PROTOCOL.md` 分阶段执行。

## 目录

| 路径 | 内容 |
|---|---|
| `controllers/` | `%#codegen` MATLAB Function 源（逐样本步函数，persistent 状态） |
| `params/` | 冻结参数导出/加载 + 台架标定增益 |
| `build/` | 模型生成脚本（库 / 离线 / RT） |
| `models/` | 生成的 .slx（gitignored） |
| `verify/` | 等效性对拍 + 台架录波分析 |
| `output/` | 等效性 CSV 等运行产物（gitignored） |
| `docs/` | 上板实验规程 |
