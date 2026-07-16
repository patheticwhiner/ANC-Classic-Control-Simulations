# MATLAB 仿真到 Simulink/dSPACE MicroLabBox 的迁移计划

## 1. 目标与边界

把 `tests/` 已冻结的五个 ANC 控制器（2 kHz 逐样本离散实现）迁移到
Simulink，最终部署到 dSPACE MicroLabBox (DS1202)，输入/输出替换为
ADC/DAC。三层产物，全部在 `dspace/`：

1. **离线闭环模型** `models/offline_<demo>_<test>_<variant>.slx` —— 被控
   对象用 cylinder1dm 实测 FIR 路径，与脚本仿真逐样本对拍（硬门禁：
   fixed ≤ 1e-12，adaptive ≤ 1e-9，Δsupp < 0.01 dB）。
2. **实时模型** `models/rt_*.slx` —— 同一控制器库链接，I/O 换 RTI1202
   ADC/DAC（无 RTI 环境自动降级 Inport/Outport 占位），含
   MasterGain/软件急停/滑动 RMS 看门狗安全层。
3. **上板规程** `docs/MICROLABBOX_EXPERIMENT_PROTOCOL.md` —— Stage 0
   回环标定 → Stage 1 次级路径复辨识 → Stage 2 固定 RST 增益阶梯 →
   Stage 3 自适应逐个上板 → Stage 4 同口径对比战役。

边界：与 FPGA 计划同一原则 —— 硬件阶段保留同一套抑制/需求/饱和/稳定
指标；开发机绝不向硬件下载/运行；台架上禁止重调冻结参数。

## 2. 迁移矩阵

| 控制器 | 运行时结构 | 场景 | 依据（冻结评价） |
|---|---|---|---|
| demo1/demo4 | RST + Youla-Q RLS（共用步函数，s0_norm/control_scale 参数区分） | T1 fixed/adaptive, T2 adaptive | 全部成功，27–40 dB |
| demo2 | 增广 LQG (+ filtered-x Q, RLS/NLMS) | T1 fixed/adaptive, T2 adaptive | T2 成功；T1 内稳定性指标未过 → 台架仅短时观察实验 |
| demo3 | F(z)+固定 FIR / IMC-FxNLMS（需参考 ADC） | T1 fixed/adaptive, T2 adaptive | 成功；前馈结构，上板先验证相干性 |
| demo5 | Marino-Tomei 自适应内模 | 仅 T1（fixed/adaptive） | T2 Case A 不适用，禁止扫频上板 |
| demo5_imc | —— 不迁移 | —— | 评估需求超限（4.22/4.42 > 4） |

## 3. 设计决策

- **所有 .slx 由脚本生成**（`dspace/build/`），`.slx`/`.mat` 保持
  gitignore；生成脚本 + `dspace/controllers/*.m` 步函数源码是提交物。
- **冻结参数不手抄**：`export_frozen_params` 从 stage MAT 重跑
  `controller_demoN` 收割 `result.extra`（复用设计代码路径），并断言
  重跑 supp_db 与冻结值逐位一致（同时证明 demo2/demo5 的 extra 附加
  是纯输出）。产物 `anc_frozen_params.mat` 按需再生。
- **逐位一致移植**：步函数保持脚本引擎的运算顺序（含 RLS 的
  `inv()`、demo4 的 `/S0(1)` 与 `control_scale`——demo1 传 1，IEEE754
  恒等）；离线被控对象块复刻各引擎的递推算术（demo2 用状态空间块，
  其余多项式递推）。
- **参数经模型工作区**（Parameter 作用域、非可调 → 编译期常量，
  persistent 数组定尺寸，满足 dSPACE codegen 无动态分配约束）。
- **抗混叠是外部模拟滤波**（2 kHz 采样下模型内不可解）；16 kHz 过采样
  + 模型内抽取为备选风险项，一期不做。

## 4. 状态（2026-07-16）

- [x] Phase A 脚手架 + 参数导出管线（代码完成）
- [x] Phase B 五个控制器步函数 + 库/离线模型生成器 + 等效性验证脚本
- [x] Phase C RT 模型生成器（RTI/占位双分支 + 安全层）
- [x] Phase D 实验规程 + 台架录波分析
- [ ] **运行时验证未执行**：本机 `matlab -batch` 挂起（残留实例/license
  问题），需在可用 MATLAB 会话中执行 `dspace/README.md` 的再生成流程，
  过等效性门禁后提交 `SIMULINK_EQUIVALENCE_REPORT.md`
- [ ] dSPACE 主机构建 + 上板 Stage 0–4（人工，按规程）

## 5. 风险

| 风险 | 缓解 |
|---|---|
| MATLAB Function 脚本化注入（Stateflow API）版本行为差异 | demo1_fixed 探路案例先过 1e-12 门禁再信任全矩阵 |
| `filter()` 带通的 DF2T 复现在自适应回路中放大微差 | demo2 自适应容差 1e-9；超差即调查，不放宽 |
| RLS `inv()` 在 rti1202 codegen 报错 | 备选 rank-1 更新，须重过等效性并披露 |
| RTI 块名与假设不符 | `add_rti_block` 正则搜索 + 明确报错，主机上一次性修正 |
| 台架声学路径漂移 | Stage 1 复辨识门禁；超限走完整重冻结流程 |
