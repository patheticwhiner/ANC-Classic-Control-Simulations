# 模型来源与统一分析报告

> 生成日期: 2026-07-13  |  生成器: `tools/generate_model_report.m`  |  注册表: `dataset/model_registry.m`
>
> 本报告是模型资产、来源链路和适用范围的统一入口。公式级说明仍保留在 `dataset/ModelReference.md`；辨识过程细节以各 MAT 文件同名 Markdown 报告为准。

## 一、先建立模型地图

这些文件并不是一组可以直接按“抑制 dB”排序的同类对象。它们分别用来隔离不同控制难题：

| 模型族 | 代表模型 | 它刻意暴露的问题 | 它不能证明什么 |
|:---|:---|:---|:---|
| 实测反馈声学路径 | `armax_30303022` | 高阶、长延迟、NMP 零点、辨识不确定性同时存在 | 不能单独归因某一种算法机制 |
| 管道 RST 文献对象 | `syn_Carmona2000_7th` | 离散管道、明显延迟、非最小相位 | 不能代表当前 48 kHz 实验平台 |
| Jafari 理论对象 | `syn_TAC2015_3rd`, `syn_JVC2017_*` | 低增益、边界/RHP 零点、频谱展平和自适应律 | 不能给出真实声学工程性能 |
| 不稳定鲁棒对象 | `syn_Bai1997_4th` | 开环镇定与性能塑形必须同时完成 | 不能与稳定模型只比抑制量 |
| 前馈双路径对象 | `syn_Ho2020_ALE` | P/S 因果关系、次级通路 NMP、filtered-x 结构 | 不能直接验证反馈 ANC 控制器 |
| 教学对象 | `syn_RSTtoy_2nd`, `syn_MassSpringDamper_2nd` | 验证方程、接口和优化器是否基本正确 | 不能支撑 ANC 有效性结论 |
| 扰动/滤波模型 | `syn_whitenoise`, `syn_bpf` | 生成测试信号或早期 LQG 场景 | 不是与声学路径并列的被控对象 |
| 非线性输出调节 | demo5 脚本内嵌模型 | 未知频率、外系统、非线性参数自适应 | 不能用同一套 LTI 零极点指标概括 |

**最重要的分类轴不是算法名称，而是：单路径/双路径、反馈/前馈、稳定/不稳定、最小/非最小相位、解析模型/辨识模型。** 只有这些轴一致，性能数值才有横向意义。

## 二、动力学画像

下表直接从当前模型计算；标为“分析快照”的行来自既有分析文档。连续模型频率按 `f = omega/(2*pi)` 换算，仅表示模型内部时间尺度，不能当成声学 Hz 与离散模型横比。

| 模型/路径 | 域 | 极点与稳定性 | 零点约束 | 延迟 | 扫描内最大增益 | 分析依据 |
|:---|:---|:---|:---|:---|:---|:---|
| `armax_30303022` | 离散 / 48 kHz | 30 极点；稳定；max\|p\|=0.9958 | 29 零点；9 NMP | 22 点 / 0.458 ms | 2.50 dB @ 334.6 Hz | 分析快照（当前对象为空） |
| `syn_TAC2015_3rd` | 离散 / 480 Hz | 3 极点；稳定；max\|p\|=0.7096 | 2 零点；1 边界 | 0 | -54.67 dB @ 240 Hz | 当前模型计算 |
| `syn_JVC2017_3rd` | 连续 | 2 极点；稳定；max Re(p)=-0.5 | 1 零点；1 NMP | 0 | -5.88 dB @ 0.1769 Hz | 当前模型计算 |
| `syn_JVC2017_6th` | 连续 | 6 极点；稳定；max Re(p)=-0.5 | 5 零点；1 NMP + 3 边界 | 0 | -4.34 dB @ 0.06022 Hz | 当前模型计算 |
| `syn_Bai1997_4th` | 连续 | 4 极点；不稳定；max Re(p)=0.6612 | 4 零点；2 NMP | 0 | 3.55 dB @ 0.1002 Hz | 当前模型计算 |
| `syn_Carmona2000_7th` | 离散 / 2 kHz | 7 极点；稳定；max\|p\|=0.9534 | 13 零点；13 NMP | 6 点 / 3 ms | 21.24 dB at DC | 当前模型计算 |
| `syn_MassSpringDamper_2nd` | 连续 | 2 极点；稳定；max Re(p)=-0.25 | 0 零点；最小相位 | 0 | 6.30 dB @ 0.1489 Hz | 当前模型计算 |
| `syn_Ho2020_ALE:P` | 离散 / 4 kHz | 2 极点；稳定；max\|p\|=0 | 2 零点；最小相位 | 5 点 / 1.25 ms | 3.52 dB @ 2000 Hz | 当前模型计算 |
| `syn_Ho2020_ALE:S` | 离散 / 4 kHz | 2 极点；稳定；max\|p\|=0 | 2 零点；1 NMP | 2 点 / 0.5 ms | 7.96 dB @ 1000 Hz | 当前模型计算 |
| `syn_RSTtoy_2nd` | 离散 / 1 Hz | 2 极点；稳定；max\|p\|=0.6708 | 1 零点；最小相位 | 1 点 / 1 s | 3.00 dB @ 0.03564 Hz | 当前模型计算 |
| `syn_whitenoise` | 离散 / 1 kHz | 8 极点；稳定；max\|p\|=0.9734 | 7 零点；1 NMP | 0 | -0.00 dB @ 35.22 Hz | 当前模型计算 |
| `syn_bpf` | 离散 / 1 kHz | 8 极点；稳定；max\|p\|=0.9551 | 7 零点；5 NMP | 0 | 14.93 dB @ 40.59 Hz | 当前模型计算 |

读表时注意：

- “稳定”只说明自由响应不发散，不代表容易控制。ARMAX 与 Carmona 都稳定，但延迟和 NMP 零点仍限制可达带宽与稳定逆。
- `z=1` 或虚轴零点列为“边界零点”，与严格位于单位圆外/RHP 的 NMP 零点分开。
- FIR 路径的原点极点主要表示有限记忆，不应和接近单位圆的轻阻尼动态等量看待；纯延迟已单独列出。
- “扫描内最大增益”用于认识增益尺度，不自动等于共振峰；例如 Carmona 的最大值在 DC。

## 三、逐模型认识

### 3.1 `armax_30303022`：真实但不可重新计算的当前基准

- **怎么看它**：这是稳定的实测声学路径，不需要先解决开环镇定；主要动态集中在约 334.6 Hz，最大极点模约 0.9958，说明共振衰减慢、模型记忆长。
- **真正障碍**：30 阶、9 个 NMP 零点、22 点延迟。22 点在 48 kHz 下是 0.458 ms，经验反馈带宽上限约 727 Hz；全阶 Bezout 方程还会形成约 `1e16` 的病态矩阵。
- **研究含义**：它适合检验控制器能否同时承受延迟、NMP 和高阶辨识误差，不适合用来单独证明某个理论环节。
- **数据警告**：当前 MAT 中 `ARMAXmodel.model` 是空数组。报告里的 30/29 阶、9 个 NMP 零点和 334.6 Hz 来自 `About_Plant_Performance.md` 的既有快照；必须恢复 A/B/C 多项式后才能重新设计和独立复核。

### 3.2 `syn_Carmona2000_7th`：点数不大，但物理延迟更大

- 7 个动态极点、6 点纯延迟，采样率 2 kHz。6 点只看数字似乎小于 ARMAX 的 22 点，但对应 **3 ms**，约是 0.458 ms 的 6.5 倍。
- 当前系数计算得到 13 个单位圆外零点，说明稳定逆不是可用策略；RST/Youla 必须做闭环塑形，而不是直接消去对象。
- 它适合研究低阶多项式设计和延迟影响；不应因为同为“管道模型”就替代 48 kHz 实测对象。

### 3.3 Jafari 系列：把单一理论障碍逐步放大

- `syn_TAC2015_3rd` 全部极点稳定，`z=1` 是边界零点而非严格 NMP；DC 增益为零，整个频带增益也很低。控制难点是需要很大的补偿增益，同时避免在边界零点附近作不可靠求逆。
- `syn_JVC2017_3rd` 只有一对稳定极点和一个 RHP 零点 `s=0.2`，是理解“为什么要反射 NMP 零点、构造 F(s)”的最干净入口。
- `syn_JVC2017_6th` 保留同一 RHP 零点问题，同时加入多重极点、原点零点和虚轴零点，用于检查频谱展平在高阶/边界零点下是否仍成立。
- 三者是机制模型：先在这里理解算法为什么有效，再到 ARMAX 检查数值和工程边界。

### 3.4 `syn_Bai1997_4th`：唯一以“先镇定”为首要任务的对象

- 4 个极点中有 2 个 RHP 极点，4 个零点中有 2 个 RHP 零点。开环不稳定和非最小相位同时存在。
- 对它的第一项评价是闭环内部稳定性，其次才是灵敏度塑形和抑制；一个抑制较低但可靠镇定的控制器，不能简单判为差于稳定对象上的高抑制控制器。
- 它是 H-infinity 权重设计和水床效应的基准，不是管道 ANC 的缩小版。

### 3.5 `syn_Ho2020_ALE`：理解前馈 ANC 的因果窗口

- 主通路 P 延迟 5 点、次级通路 S 延迟 2 点，4 kHz 下分别为 1.25 ms 和 0.5 ms；参考信号因此提供约 **3 点/0.75 ms** 的因果余量。
- S 含一个单位圆外零点（最大模 2），所以 `S^-1` 不是稳定因果滤波器；filtered-x/ALE 的意义正是避免把控制问题简化成精确求逆。
- 这是双路径前馈对象。若拿它与反馈 ARMAX 模型比较，改变的不仅是模型，还包括信息结构，因此控制器抑制量不能直接横比。

### 3.6 教学与信号模型：知道它们证明到哪里为止

- `syn_RSTtoy_2nd` 低阶、稳定、最小相位，只能证明 RST 方程和优化器在良态问题上可运行。
- 质量-弹簧-阻尼器用于理解状态反馈、观测器和阻尼，不包含声传播延迟、NMP 次级通路或辨识误差。
- `syn_whitenoise` 是扰动生成模型，`syn_bpf` 是早期 LQG 场景对象；二者 MAT 中采样时间原为 `-1`（未指定），本报告按 DataManager 的名义 1 kHz 解释。其绝对频率结论因此低于明确采样率模型的可信度。

### 3.7 demo5：不是缺一张零极点表，而是模型范式不同

demo5 的核心是非线性输出调节、外系统和未知频率估计。它的“模型”包含被控对象、扰动生成器和自适应参数化，不能压缩成一个 SISO 传递函数后与前四个 demo 公平比较。正确对齐方式是统一扰动场景、可用测量和控制能量，而不是强行统一零极点阶数。

## 四、哪些结论可以比较

| 想回答的问题 | 应使用的模型 | 必须保持一致 | 不应混入 |
|:---|:---|:---|:---|
| 哪个控制器更适合当前实测反馈路径 | 恢复后的 `armax_30303022` 或同一新采集模型 | 对象、采样率、测试信号、窗口、控制限幅 | 论文低阶对象、前馈双路径对象 |
| 某篇论文方法是否复现正确 | 该论文对应 `syn_*` | 原始参数、扰动、初值、指标 | 其他论文模型上的更高 dB |
| NMP 零点处理是否有效 | JVC 3rd -> JVC 6th -> ARMAX | 相同 F/K 分解和稳定性判据 | 只在 RST toy 上成功 |
| 延迟如何限制反馈带宽 | Carmona 与 ARMAX 可作定性对照 | 比较毫秒而不是采样点数 | 直接比较阶数或抑制 dB |
| 前馈算法是否利用参考预览 | Ho P/S 或同一实测 P/S 对 | 主/次通路、参考噪声、因果余量 | 单路径反馈模型 |
| SYSID 表示哪一种更好 | 同一次采集的 FIR/ARMAX/OE/SS | 同一训练/验证段和采样率 | 跨日期、跨装置的 fit 数值 |

### 推荐认识顺序

1. 用 RST toy 和质量-弹簧模型确认基本方程、极点和状态反馈概念。
2. 用 JVC 3rd 单独理解一个 RHP 零点为何阻止稳定精确求逆。
3. 用 JVC 6th 与 TAC 观察边界零点、低增益和高阶频谱展平。
4. 用 Bai 区分“镇定问题”和“性能优化问题”。
5. 用 Carmona 理解离散延迟和多项式设计，再进入 ARMAX 的高阶与数值病态。
6. 最后用 Ho 与 demo5 分别扩展到前馈双路径和非线性外系统，不再强求单一 LTI 排名。

## 五、实验选型

- 跨 demo 横向对比暂时沿用 `armax_30303022`，因为当前测试基准已围绕它建立；但当前 MAT 的模型对象为空，修复前只能复用既有信号和分析快照，不能重新设计控制器。
- 论文方法复现应使用对应的 `syn_*` 文献基准，只能证明算法复现一致性，不能直接外推到真实声学平台。
- `sysid_models/` 中的新模型按采集批次动态盘点。它们适合候选模型筛选，但在残差白度、采集来源或验证报告不完整时，不应替换统一基准。
- demo5 仍有脚本内嵌对象。它们已登记来源与风险，但在抽取成独立模型并接入测试基准前，不宜与 `dataset/` 模型混为一类。

### 建议使用矩阵

| 任务 | 首选模型 | 理由 | 限制 |
|:---|:---|:---|:---|
| 四类控制器统一比较 | `armax_30303022` | tests 已统一信号、指标和加载入口 | 48 kHz、高阶、Bezout 数值病态 |
| RST/Youla 论文复现 | `syn_Carmona2000_7th` | 低阶管道 ANC 文献基准 | 不能代表当前实验平台 |
| Jafari 自适应/鲁棒复现 | `syn_TAC2015_3rd`, `syn_JVC2017_*` | 与推导和论文场景一致 | 合成理论模型 |
| H-infinity 复现 | `syn_Bai1997_4th` | 不稳定、非最小相位参考对象 | 耳机模型，不是管道模型 |
| 前馈 ALE/FXLMS 结构验证 | `syn_Ho2020_ALE` | 明确的 P(z)+S(z) 路径对 | 合成等价模型，尚无实测验证 |
| 最新硬件路径评估 | `sysid_models/<批次>/` | 保留采集批次和同名验证报告 | 需逐批审查，不能自动视为 canonical |

## 附录 A：资产概览

- 注册记录: **18**（DataManager 可加载 13，脚本内嵌/来源记录 5）
- 本机 `sysid_models/` MAT 文件: **39**（有可关联报告 19，无可关联报告 20）
- 注册路径缺失: **0**；`dataset/` 未登记 MAT 文件: **0**
- 注册资产内容问题: **1**
- 静态审计发现的其他脚本内嵌模型候选: **2**

### 来源分层

| 来源层级 | 含义 | 可支持的结论 |
|:---|:---|:---|
| `measured` / `measured-identification` | 来自硬件采集及其辨识结果 | 对指定装置、采样率和采集批次有效 |
| `literature` / `literature-reproduction` | 论文给定或按论文重建 | 支持方法复现，不自动支持工程泛化 |
| `synthetic` / `repository-synthetic` | 本工程构造的测试对象或扰动 | 支持功能与边界测试 |
| `educational` | 教材或低阶教学对象 | 支持概念演示，不支持 ANC 性能声明 |

## 附录 B：统一注册表

| ID | 角色 | 来源 | 表示 | 域/采样 | 成熟度 | 主资产 |
|:---|:---|:---|:---|:---|:---|:---|
| `armax_30303022` | plant | `measured-identification` | ARMAX(30,30,30,22) | discrete / 48000 Hz | `incomplete` | `dataset/armax_30303022_2026-01-20.mat` |
| `syn_TAC2015_3rd` | plant | `literature` | third-order denominator with one boundary zero | discrete / 480 Hz | `reference` | `dataset/syn_TAC2015_3rd.mat` |
| `syn_JVC2017_3rd` | plant | `literature` | second-order denominator with one RHP zero | continuous / nominal 10 kHz discretization | `reference` | `dataset/syn_JVC2017_3rd.mat` |
| `syn_JVC2017_6th` | plant | `literature` | sixth-order transfer function | continuous / not applicable | `reference` | `dataset/syn_JVC2017_6th.mat` |
| `syn_Bai1997_4th` | plant | `literature` | fourth-order unstable transfer function | continuous / nominal 4 kHz discretization | `reference` | `dataset/syn_Bai1997_4th.mat` |
| `syn_Carmona2000_7th` | plant | `literature` | seventh-order transfer function with delay | discrete / 2000 Hz | `reference` | `dataset/syn_Carmona2000_7th.mat` |
| `syn_MassSpringDamper_2nd` | plant | `educational` | second-order state space | continuous / not applicable | `teaching` | `dataset/syn_MassSpringDamper_2nd.mat` |
| `syn_Ho2020_ALE` | primary-secondary-path-pair | `literature` | FIR path pair | discrete / 4000 Hz | `reference` | `dataset/syn_Ho2020_ALE.mat` |
| `syn_RSTtoy_2nd` | plant | `educational` | second-order transfer function with delay | discrete / 1 Hz | `teaching` | `dataset/syn_RSTtoy_2nd.mat` |
| `syn_whitenoise` | disturbance | `synthetic` | state space | discrete / nominal 1000 Hz | `teaching` | `dataset/syn_whitenoise_ssmodel.mat` |
| `syn_bpf` | plant | `synthetic` | state space | discrete / nominal 1000 Hz | `teaching` | `dataset/syn_bpf_ssmodel.mat` |
| `lms_sysid` | identification-data | `measured` | four-channel LMS data bundle | discrete / 48000 Hz | `archive` | `dataset/lms_sysid_2026-01-20.mat` |
| `raw_dspace` | raw-identification-data | `measured` | input-output time series | discrete / 48000 Hz | `archive` | `dataset/raw_dspace_primpath.mat` |
| `marino_tomei_2011_inline` | plant-and-exosystem | `literature-reproduction` | nonlinear state equations | continuous / not applicable | `inline` | `demo5_MarinoTomei/MarinoTomei_2011_adaptive_regulation.m` |
| `marino_tomei_2016_inline` | plant-and-exosystem | `literature-reproduction` | second-order transfer function | continuous / not applicable | `inline` | `demo5_MarinoTomei/MarinoTomei_2016_adaptive_freq_est.m` |
| `marino_tomei_2023_inline` | plant-family | `literature-reproduction` | four selectable unstable transfer functions | continuous / not applicable | `review` | `demo5_MarinoTomei/MarinoTomei_2023_unstable.m` |
| `demo5_realistic_secondary_path_inline` | secondary-path | `repository-synthetic` | FIR path | discrete / 8000 Hz | `inline` | `demo5_MarinoTomei/MarinoTomei_ANC_IMC_FXLMS_compare.m` |
| `demo5_akhtar_paths_inline` | primary-secondary-path-pair | `literature-reproduction` | FIR path pair | discrete / 4000 Hz | `review` | `demo5_MarinoTomei/MarinoTomei_AkhtarBenchmark_compare.m` |

### 来源与证据

- **`armax_30303022`**: dSPACE acoustic-duct capture; identified by dataset/armax_identification.m  
  证据: `dataset/raw_dspace_primpath.mat`；使用方: demo2_LQG; demo3_Robust; demo4_Robust; tests；备注: Intended cross-demo benchmark; current MAT model object is empty and must be restored.
- **`syn_TAC2015_3rd`**: Jafari, Ioannou & Rudd (2015), IEEE Transactions on Automatic Control  
  证据: `dataset/export_all_synthetic.m`；使用方: demo3_Robust/JafariTAC_DiscreteAVC.m; demo3_Robust/demos；备注: Paper benchmark transcribed to a MAT artifact.
- **`syn_JVC2017_3rd`**: Jafari & Ioannou (2016), Journal of Vibration and Control  
  证据: `dataset/export_all_synthetic.m`；使用方: demo3_Robust/JafariJVC_Continuous.m; demo3_Robust/demos；备注: Continuous-time paper benchmark.
- **`syn_JVC2017_6th`**: Jafari & Ioannou (2016), Journal of Vibration and Control  
  证据: `dataset/export_all_synthetic.m`；使用方: demo3_Robust/JafariJVC_Continuous.m；备注: Higher-order F(s) flattening example.
- **`syn_Bai1997_4th`**: Bai & Lee (1997), IEEE Transactions on Speech and Audio Processing  
  证据: `dataset/export_all_synthetic.m`；使用方: demo3_Robust/Bai1997_Hinf.m; demo3_Robust/demos/Bai1997_MultiSine.m；备注: H-infinity active-headset benchmark.
- **`syn_Carmona2000_7th`**: Carmona & Alvarado (2000), ASME Journal of Vibration and Acoustics  
  证据: `dataset/export_all_synthetic.m`；使用方: demo1_RST/Carmona2000.m; demo1_RST/Landau2005.m；备注: Duct ANC benchmark; six-sample delay.
- **`syn_MassSpringDamper_2nd`**: Standard mass-spring-damper teaching model (m=1, k=1, b=0.5)  
  证据: `dataset/export_all_synthetic.m`；使用方: demo2_LQG/LQRdemo.m；备注: Not an acoustic plant and not suitable for ANC performance claims.
- **`syn_Ho2020_ALE`**: Ho, Shyu, Chang & Kuo (2020), IEEE/ACM TASLP  
  证据: `dataset/export_all_synthetic.m; sysid_models/20260708_Ho/build_syn_Ho2020_models.m`；使用方: planned ALE-FxLMS comparison；备注: Analytic synthetic paths; sysid_models contains equivalent FIR/ARMAX/SS forms.
- **`syn_RSTtoy_2nd`**: Repository-authored MOEA/RST benchmark  
  证据: `dataset/export_all_synthetic.m`；使用方: demo4_Robust/benchmark_MOEAs.m; demo4_Robust/run_eMOPSO.m；备注: Low-order numerical benchmark, not a measured acoustic model.
- **`syn_whitenoise`**: Repository-authored band-limited white-noise example  
  证据: `demo2_LQG/demo.m`；使用方: demo2_LQG/LQGdemo.m; demo2_LQG/LQGdemo2.m；备注: Legacy MAT layout containing disturbance-model matrices only.
- **`syn_bpf`**: Repository-authored band-pass-filter example  
  证据: `demo2_LQG/demo.m`；使用方: demo2_LQG/LQGdemo.m; demo2_LQG/LQGdemo2.m；备注: Legacy MAT layout used for early LQG experiments.
- **`lms_sysid`**: dSPACE primary/secondary, error/reference identification captures  
  证据: `sysid_models/20260120`；使用方: identification experiments；备注: Data bundle rather than a single controlled-object model.
- **`raw_dspace`**: dSPACE acoustic-duct primary-path capture  
  证据: `dataset/Environment.jpg`；使用方: dataset/armax_identification.m；备注: Upstream evidence for the legacy ARMAX benchmark.
- **`marino_tomei_2011_inline`**: Marino & Tomei adaptive output-regulation example (2011)  
  证据: `demo5_MarinoTomei/About_MarinoTomei.md`；使用方: demo5_MarinoTomei/MarinoTomei_2011_regulator_benchmark.m；备注: Parameters are embedded in two scripts and are not DataManager-loadable.
- **`marino_tomei_2016_inline`**: Marino & Tomei (2016), adaptive disturbance rejection for unknown stable linear systems  
  证据: `demo5_MarinoTomei/docs/MarinoTomei_2016_Derivation.md`；使用方: demo5_MarinoTomei；备注: P(s)=1/(s+1)^2 is defined inside the simulation.
- **`marino_tomei_2023_inline`**: Marino & Tomei 2023 unstable minimum-phase examples; citation details require verification  
  证据: `demo5_MarinoTomei/README.md`；使用方: demo5_MarinoTomei；备注: Inline switch selects one of four plants; exact paper-to-case mapping is not recorded.
- **`demo5_realistic_secondary_path_inline`**: Repository-authored realistic feedback-ANC comparison scenario  
  证据: `demo5_MarinoTomei/SIMULATION_REPORT.md`；使用方: demo5_MarinoTomei；备注: Useful scenario model, but not a measured path or literature benchmark.
- **`demo5_akhtar_paths_inline`**: Akhtar-style IMC-FxLMS benchmark; exact bibliographic source pending  
  证据: `demo5_MarinoTomei/results_AkhtarBenchmark_compare/AkhtarBenchmark_Marino_vs_IMCFxLMS.md`；使用方: demo5_MarinoTomei；备注: Path coefficients and scenario settings remain embedded in the comparison script.

## 附录 C：系统辨识模型动态清单

`sysid_models/` 默认不纳入 Git，因而不同机器的数量可能不同。下表由文件名与同名 Markdown 报告自动推导，不把“文件存在”等同于“模型已验证”。

| 资产 | 类型 | 日期 | 路径角色 | 结构/阶次 | 采样 | 来源判断 | 报告 | 状态 |
|:---|:---|:---|:---|:---|:---|:---|:---:|:---|
| `sysid_models/20251116/LMS_20251116_primpath_N1024_fs48k.mat` | FIR/LMS | 20251116 | primpath | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20251116/LMS_20251116_secpath_N1024_fs48k.mat` | FIR/LMS | 20251116 | secpath | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20251122/LMS_20251122_primpath_N1024_fs48k.mat` | FIR/LMS | 20251122 | primpath | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20251122/LMS_20251122_secpath_N1024_fs48k.mat` | FIR/LMS | 20251122 | secpath | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260120/LMS_20260120_pri-err_N1024_fs48k.mat` | FIR/LMS | 20260120 | pri-err | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260120/LMS_20260120_pri-ref_N1024_fs48k.mat` | FIR/LMS | 20260120 | pri-ref | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260120/LMS_20260120_sec-err_N1024_fs48k.mat` | FIR/LMS | 20260120 | sec-err | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260120/LMS_20260120_sec-ref_N1024_fs48k.mat` | FIR/LMS | 20260120 | sec-ref | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260120/armax_30303022_2026-01-20.mat` | ARMAX | unknown | unspecified | unknown | unknown | measured identification | 无 | undocumented |
| `sysid_models/20260120/lms_sysid_2026-01-20.mat` | FIR/LMS | unknown | unspecified | unknown | unknown | measured identification | 无 | undocumented |
| `sysid_models/20260323/LMS_20260323_primpath_N-unknown_fs48k_legacy.mat` | FIR/LMS | 20260323 | primpath | N=-unknown | 48 kHz | measured identification | 无 | legacy; metadata incomplete |
| `sysid_models/20260323/LMS_20260323_primpath_N1024_fs48k.mat` | FIR/LMS | 20260323 | primpath | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260323/LMS_20260323_secpath_N-unknown_fs48k_legacy.mat` | FIR/LMS | 20260323 | secpath | N=-unknown | 48 kHz | measured identification | 无 | legacy; metadata incomplete |
| `sysid_models/20260323/LMS_20260323_secpath_N1024_fs48k.mat` | FIR/LMS | 20260323 | secpath | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260329/LMS_20260329_primpath_N1024_fs48k.mat` | FIR/LMS | 20260329 | primpath | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260329/LMS_20260329_secpath_N1024_fs48k.mat` | FIR/LMS | 20260329 | secpath | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260411/LMS_20260411_pri-err_N1024_fs48k.mat` | FIR/LMS | 20260411 | pri-err | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260411/LMS_20260411_pri-ref_N1024_fs48k.mat` | FIR/LMS | 20260411 | pri-ref | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260411/LMS_20260411_sec-err_N1024_fs48k.mat` | FIR/LMS | 20260411 | sec-err | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260411/LMS_20260411_sec-ref_N1024_fs48k.mat` | FIR/LMS | 20260411 | sec-ref | N=1024 | 48 kHz | measured identification | 无 | undocumented |
| `sysid_models/20260706/ARMAX_20260706_secpath_[16-60-4-43]_fs12k_sim90.mat` | ARMAX | 20260706 | secpath | [16-60-4-43] | 12 kHz | measured identification | 同名 | documented candidate |
| `sysid_models/20260706/ARMAX_RETRY_20260706_secpath_fs12k.mat` | ARMAX | 20260706 | secpath | unknown | 12 kHz | measured identification | 同名 | provisional/retry |
| `sysid_models/20260706/LMS_20260706_secpath_N1024_fs48k.mat` | FIR/LMS | 20260706 | secpath | N=1024 | 48 kHz | measured identification | 同名 | documented candidate |
| `sysid_models/20260706/OE_20260706_secpath_[10-20-43]_fs12k_sim89.mat` | OE | 20260706 | secpath | [10-20-43] | 12 kHz | measured identification | 同名 | documented candidate |
| `sysid_models/20260706/OE_20260706_secpath_[10-40-43]_fs12k_sim90.mat` | OE | 20260706 | secpath | [10-40-43] | 12 kHz | measured identification | 同名 | documented candidate |
| `sysid_models/20260708_Ho/ARMAX_20260708_secpath_[0-3-0-2]_fs4k_syn.mat` | ARMAX | 20260708 | secpath | [0-3-0-2] | 4 kHz | literature synthetic | 批次 | reference equivalent |
| `sysid_models/20260708_Ho/FIR_20260708_primpath_L8_fs4k_syn.mat` | FIR/LMS | 20260708 | primpath | unknown | 4 kHz | literature synthetic | 批次 | reference equivalent |
| `sysid_models/20260708_Ho/FIR_20260708_secpath_L5_fs4k_syn.mat` | FIR/LMS | 20260708 | secpath | unknown | 4 kHz | literature synthetic | 批次 | reference equivalent |
| `sysid_models/20260708_Ho/SS_20260708_secpath_n4_fs4k_syn.mat` | state space | 20260708 | secpath | n=4 | 4 kHz | literature synthetic | 批次 | reference equivalent |
| `sysid_models/20260710_cylinder1dm/LMS_20260710_secpath_N1024_fs48k.mat` | FIR/LMS | 20260710 | secpath | N=1024 | 48 kHz | measured identification | 同名 | documented candidate |
| `sysid_models/20260710_cylinder1dm/LMS_20260710_secpath_N16_fs4k.mat` | FIR/LMS | 20260710 | secpath | N=16 | 4 kHz | measured identification | 同名 | documented candidate |
| `sysid_models/20260710_cylinder1dm/LMS_20260710_secpath_N32_fs4k.mat` | FIR/LMS | 20260710 | secpath | N=32 | 4 kHz | measured identification | 同名 | documented candidate |
| `sysid_models/20260710_cylinder1dm/LMS_20260710_secpath_N8_fs2k.mat` | FIR/LMS | 20260710 | secpath | N=8 | 2 kHz | measured identification | 同名 | documented candidate |
| `sysid_models/20260710_cylinder2dm/LMS_20260710_secpath_N16_fs1k.mat` | FIR/LMS | 20260710 | secpath | N=16 | 1 kHz | measured identification | 同名 | documented candidate |
| `sysid_models/20260710_cylinder2dm/LMS_20260710_secpath_N8_fs2k.mat` | FIR/LMS | 20260710 | secpath | N=8 | 2 kHz | measured identification | 同名 | documented candidate |
| `sysid_models/20260713_cylinder1dm/ARMAX_20260713_pripath_[4-8-2-2]_fs2k_sim87.mat` | ARMAX | 20260713 | pripath | [4-8-2-2] | 2 kHz | measured identification | 同名 | candidate; residual whiteness failed |
| `sysid_models/20260713_cylinder1dm/ARMAX_20260713_secpath_[4-8-2-1]_fs2k_sim90.mat` | ARMAX | 20260713 | secpath | [4-8-2-1] | 2 kHz | measured identification | 同名 | candidate; residual whiteness failed |
| `sysid_models/20260713_cylinder1dm/LMS_20260713_pripath_N16_fs2k.mat` | FIR/LMS | 20260713 | pripath | N=16 | 2 kHz | measured identification | 同名 | documented candidate |
| `sysid_models/20260713_cylinder1dm/LMS_20260713_secpath_N16_fs2k.mat` | FIR/LMS | 20260713 | secpath | N=16 | 2 kHz | measured identification | 同名 | documented candidate |

## 附录 D：完整性审计

### 注册资产缺失

注册表中的主资产在本机不存在；DataManager 对应入口会失败。

- 无。

### `dataset/` 未登记 MAT 文件

文件存在但没有进入统一注册表。

- 无。

### 注册模型内容问题

路径存在，但模型内容不足以重新计算动力学或复现实验。

- `armax_30303022: ARMAXmodel.model is empty; restore A/B/C polynomials`

### `sysid_models/` 缺少可关联报告

只能从文件名推断来源和结构，无法审查采集与验证条件。

- `sysid_models/20251116/LMS_20251116_primpath_N1024_fs48k.mat`
- `sysid_models/20251116/LMS_20251116_secpath_N1024_fs48k.mat`
- `sysid_models/20251122/LMS_20251122_primpath_N1024_fs48k.mat`
- `sysid_models/20251122/LMS_20251122_secpath_N1024_fs48k.mat`
- `sysid_models/20260120/LMS_20260120_pri-err_N1024_fs48k.mat`
- `sysid_models/20260120/LMS_20260120_pri-ref_N1024_fs48k.mat`
- `sysid_models/20260120/LMS_20260120_sec-err_N1024_fs48k.mat`
- `sysid_models/20260120/LMS_20260120_sec-ref_N1024_fs48k.mat`
- `sysid_models/20260120/armax_30303022_2026-01-20.mat`
- `sysid_models/20260120/lms_sysid_2026-01-20.mat`
- `sysid_models/20260323/LMS_20260323_primpath_N-unknown_fs48k_legacy.mat`
- `sysid_models/20260323/LMS_20260323_primpath_N1024_fs48k.mat`
- `sysid_models/20260323/LMS_20260323_secpath_N-unknown_fs48k_legacy.mat`
- `sysid_models/20260323/LMS_20260323_secpath_N1024_fs48k.mat`
- `sysid_models/20260329/LMS_20260329_primpath_N1024_fs48k.mat`
- `sysid_models/20260329/LMS_20260329_secpath_N1024_fs48k.mat`
- `sysid_models/20260411/LMS_20260411_pri-err_N1024_fs48k.mat`
- `sysid_models/20260411/LMS_20260411_pri-ref_N1024_fs48k.mat`
- `sysid_models/20260411/LMS_20260411_sec-err_N1024_fs48k.mat`
- `sysid_models/20260411/LMS_20260411_sec-ref_N1024_fs48k.mat`

### 其他脚本内嵌模型候选

脚本直接构造 tf/ss/zpk，且没有加载注册模型；需人工判断是对象、控制器还是辅助滤波器。

- `demo3_Robust/demos/HinfSynTest.m`
- `demo5_MarinoTomei/internal_model.m`

## 附录 E：治理规则

1. 新的长期使用模型先登记到 `dataset/model_registry.m`，再由 `DataManager` 或明确的加载函数暴露。
2. 实测辨识模型必须保留采集标识、路径角色、原始/辨识采样率、模型结构、验证指标和同名报告；缺一项即保持 candidate 状态。
3. 论文模型必须记录完整文献来源，并区分“论文直接给定”“按图表重建”“本工程为对比而构造”。
4. 同一物理路径的 FIR、ARMAX、OE、SS 是不同表示，不是独立实验来源；报告时必须共享同一采集 lineage。
5. 跨算法性能表只允许使用同一物理对象、采样率、测试信号、起止窗口和指标定义。
6. 每次新增、删除或移动模型后运行 `generate_model_report`，并审查“完整性审计”是否出现新条目。

## 附录 F：复现命令

```matlab
addpath('dataset', 'tools');
DataManager()               % 查看统一可加载入口
generate_model_report()     % 重新生成本报告
```
