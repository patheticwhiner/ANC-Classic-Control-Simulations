# dataset — 模型数据与导入工具

此目录存放主动噪声控制仿真所需的**模型数据文件**以及**数据导入工具**。所有`.mat`文件采用统一命名规范：`<类型>_<特征标识>_<日期>.mat`。

## 快速使用

```matlab
% 列出所有可用数据源
DataManager()

% 加载指定数据源
info = DataManager('armax_30303022');   % ARMAX辨识模型
info = DataManager('syn_whitenoise');   % 合成带限白噪声干扰模型
info = DataManager('syn_bpf');          % 合成带通滤波器被控对象
info = DataManager('lms_sysid');        % LMS系统辨识实验数据
info = DataManager('raw_dspace');       % dSPACE采集原始信号
```

所有数据源返回统一格式的struct，用`.type`字段区分数据类型。

## 模型清单

### 实测模型

| 数据源 | 文件名 | 类型 | 描述 | 使用方 |
|---|---|---|---|---|
| `armax_30303022` | `armax_30303022_2026-01-20.mat` | ARMAX辨识模型 | ARMAX(30,30,30,22), 实测声学管道 @48kHz | demo2_LQG, demo3_Robust |
| `lms_sysid` | `lms_sysid_2026-01-20.mat` | LMS辨识数据 | 4通道系统辨识数据 (pri/sec × err/ref) | 待编写辨识脚本 |
| `raw_dspace` | `raw_dspace_primpath.mat` | 原始信号 | dSPACE控制器导出的初级路径录音 | armax_identification.m |

### 合成理论模型 (syn_*.mat)

| 数据源 | 文件名 | 来源 | 阶数 | 域 | 使用方 |
|---|---|---|---|---|---|
| `syn_TAC2015_3rd` | `syn_TAC2015_3rd.mat` | Jafari et al. (2015) IEEE TAC | 3rd | 离散 Ts=1/480 | demo3_Robust ×3 |
| `syn_JVC2017_3rd` | `syn_JVC2017_3rd.mat` | Jafari & Ioannou (2016) JVC | 3rd | 连续 fs=10k | demo3_Robust ×3 |
| `syn_JVC2017_6th` | `syn_JVC2017_6th.mat` | Jafari & Ioannou (2016) JVC | 6th | 连续 | demo3_Robust |
| `syn_Bai1997_4th` | `syn_Bai1997_4th.mat` | Bai & Lee (1997) IEEE SAP | 4th | 连续 fs=4k | demo3_Robust ×2 |
| `syn_Carmona2000_7th` | `syn_Carmona2000_7th.mat` | Carmona & Alvarado (2000) ASME | 7th | 离散 Fs=2k | demo1_RST ×2 |
| `syn_MassSpringDamper_2nd` | `syn_MassSpringDamper_2nd.mat` | 教材标准模型 | 2nd | 连续 | demo2_LQG |
| `syn_Ho2020_ALE` | `syn_Ho2020_ALE.mat` | Ho et al. (2020) IEEE/ACM TASLP | P:7 S:4 | 离散 fs=4k | 待编写 |
| `syn_RSTtoy_2nd` | `syn_RSTtoy_2nd.mat` | MOEA benchmark | 2nd | 离散 Ts=1 | demo4_Robust ×2 |
| `syn_whitenoise` | `syn_whitenoise_ssmodel.mat` | 合成噪声模型 | — | SS | demo2_LQG |
| `syn_bpf` | `syn_bpf_ssmodel.mat` | 合成滤波器模型 | — | SS | demo2_LQG |

> **命名规范**: `syn_{来源}_{阶数}.mat`，来源含论文年份（如 TAC2015）作为文献标识，非数据采集日期。

## 工具脚本

| 文件 | 功能 |
|---|---|
| `DataManager.m` | 统一数据加载入口，case-switch选择模型 |
| `armax_identification.m` | 从原始信号数据辨识ARMAX模型并保存 |

## ARMAX 模型分析

### 性能上限理论分析

**[`About_Plant_Performance.md`](About_Plant_Performance.md)** — 被控对象性能上限的完整理论推导，涵盖 Bode 灵敏度积分、Poisson 积分约束、延迟约束、内模原理、前馈上限，以及四种 ANC 架构的定量对比。

### 分析工具

ARMAX 模型的结构分析工具位于 `dataset/analysis/analyze_armax_model.m`：

```matlab
% 分析 ARMAX 模型的能观标准型、互质性、零极点、灵敏度等
results = analyze_armax_model('armax_30303022');
```

该函数输出：
- 能观标准型状态空间矩阵 (Ao, Bo, Co, Do)
- 互质性检验结果 (is_coprime)
- 强可镇定判断 (is_strongly_stabilizable)
- 零极点位置与开环系统图形

## 目录结构

```
dataset/
├── README.md                           ← 本文件
├── DataManager.m                       ← 数据加载统一入口
├── analysis/analyze_armax_model.m      ← ARMAX 结构与控制性质分析
├── armax_identification.m              ← ARMAX辨识脚本
├── restructure_data.m                  ← 数据重构/合并脚本 (一次性)
├── armax_30303022_2026-01-20.mat       ← ARMAX(30,30,30,22)辨识模型
├── syn_whitenoise_ssmodel.mat          ← 合成带限白噪声干扰模型
├── syn_bpf_ssmodel.mat                 ← 合成带通滤波器模型
├── lms_sysid_2026-01-20.mat            ← LMS系统辨识数据 (4通道合并)
├── raw_dspace_primpath.mat             ← dSPACE原始录音信号
└── Environment.jpg                     ← 实验管道平台照片 (参考布局)
```

## 实验平台参考

管道ANC实验平台设备布局 (参见 `Environment.jpg`)：

| 设备 | 位置 (cm) |
|---|---|
| 主扬声器 (Primary Speaker) | 0 |
| 参考麦克风 (Reference Mic) | 21 |
| 次扬声器 (Secondary Speaker) | 50 |
| 误差麦克风 (Error Mic) | 73 |

扬声器-麦克风距离 (cm/采样点 @48kHz, 声速≈349m/s):
| | 参考麦克风 | 误差麦克风 |
|---|---|---|
| 主扬声器 | 21cm/29samp | 73cm/100samp |
| 次扬声器 | -29cm/-41samp | 23cm/32samp |
