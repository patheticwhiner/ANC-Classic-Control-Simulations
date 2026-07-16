# 基于经典控制律的主动噪声控制

在研究过程中，笔者发现单通道的主动噪声控制（Active Noise Cancellation, ANC）为经典自动控制理论提供了非常适于教学的实验设置。并且在实践中，也的确有许多研究探究控制理论与ANC信号控制领域的结合。本工程旨在探究经典控制理论在ANC中的适用性，同时完成研究中实验复现的整理。

**尽管基于RST的研究同样可以归于loop-shaping，而被间接归于鲁棒控制，这里仍然将其作为一种较为传统的控制理论单独分类，以符合一般的经典控制教材讲授的思路。*

## 项目初始化

在运行 demo、测试或模型分析前，从仓库根目录执行：

```matlab
run('project_init.m');
```

公共运行期函数位于 `shared/`，由 `project_init.m` 按职责统一加载。各 demo
自己的 `utils/` 仍由对应入口显式加入，`tools/` 中的维护命令不会自动加载。

## 模型统一入口

- [`MODEL_ANALYSIS_REPORT.md`](MODEL_ANALYSIS_REPORT.md)：模型族地图、动力学画像、逐模型解读、可比性边界、来源链路与完整性审计。
- [`dataset/model_registry.m`](dataset/model_registry.m)：模型身份与来源的可机读单一事实源，`DataManager` 从这里派生加载入口。
- [`dataset/ModelReference.md`](dataset/ModelReference.md)：论文/教学模型的公式、零极点和关键性质详解。

模型新增、删除或移动后，在 MATLAB 中运行：

```matlab
addpath('dataset', 'tools');
generate_model_report();
```

## 工程目录

```
ANC-Classic-Control-Simulations/
├── project_init.m 统一加载 shared/ 与 dataset/ 公共路径
├── shared/ 跨 demo 公共运行期函数（RST、DSP、控制分析）
├── tools/ 仓库维护与报告生成工具
├── assets/ 存放md文档所需的图片
│   
├── signal_excitation/ 此目录下存放与工程相关的噪声模型说明
│   └── AboutExcitation.md 此文档提供子目录下文件的说明
├── dataset/ 此目录下存放模型数据文件与数据导入工具 (详见 dataset/README.md)
│   ├── README.md                    目录说明与模型清单
│   ├── model_registry.m             模型身份、来源、用途和成熟度注册表
│   ├── DataManager.m                数据加载统一入口
│   ├── analysis/                    模型专用分析入口
│   ├── armax_identification.m       ARMAX辨识脚本
│   ├── armax_30303022_2026-01-20.mat ARMAX(30,30,30,22)辨识模型
│   ├── syn_whitenoise_ssmodel.mat   合成带限白噪声干扰模型
│   ├── syn_bpf_ssmodel.mat          合成带通滤波器模型
│   ├── lms_sysid_2026-01-20.mat     LMS系统辨识数据
│   └── raw_dspace_primpath.mat      dSPACE原始录音信号
│   
├── demo1_RST/ 此目录下存放一些用于调试与学习的临时文件
│   ├── 
│   ├── 
│   └── AboutRST.md 此文档提供子目录下文件的说明
|
├── demo2_LQG/ 此目录下存放一些用于调试与学习的临时文件
│   ├── LQRdemo.m   LQR控制学习入门
│   ├── LQGdemo.m   使用dlqr, dlqe完成控制器设计
│   ├── LQGdemo2.m  试图使用lqg函数完成控制器设计
│   ├── LQG_idmodel.m    将LQG设计拓展至辨识模型
│   └── AboutLQG.md 此文档对工程中设计的LQG控制方法作了相对完整的说明
|
├── demo3_Robust/ 鲁棒与自适应ANC (Jafari-Ioannou + H∞)
│   ├── Theory_Foundations.md         公共理论基础
│   ├── Derivation_Problem1~4.md      四道题的考试答题级推导
│   ├── ExperimentReport_TimeVarying.md  时变频率实验报告
│   ├── JafariANC_RealAcoustic.m      统一仿真入口
│   └── demos/                        开发变体与实验脚本
│
├── demo4_Robust/ ε-MOPSO RST控制器参数整定
│   ├── run_RST_eMOPSO.m              主入口脚本
│   ├── benchmark_MOEAs.m             参数分析
│   ├── utils/                         算法核心与工具函数
│   ├── output/                        运行输出
│   └── About_PSO_Foundations.md       PSO完整推导 (含定理证明)
│
└── README.md 此文档提供了工程的概述
```

## 实验记录

在上述仿真实验一直出现卡壳的状况下，建议引入一些辅助实验，帮助完成问题的解决。这些实验情况应该得到良好的记录，而不是随意进行。完整的项目进展记录在下节中，而每个项目各自进行的详细情况，记录在各自文档中。文档链接见下表。

| 目录   | 文档链接                                                     |
| ------ | ------------------------------------------------------------ |
| RST    | [AboutRST](.\demo1_RST\AboutRST.md)：极点配置法应用的实验记录文档 |
| LQG    | [AboutLQR](.\demo2_LQG\AboutLQR.md)：LQR仿真记录，用于快速了解LQR的基本思想和应用<br />[AboutLQGpt1](.\demo2_LQG\AboutLQGpt1.md)：针对带通滤波器模型 的 LQG初步仿真记录<br />[AboutLQGpt2](.\demo2_LQG\AboutLQGpt2.md)：针对声管道辨识模型 的 LQG仿真记录 |
| Robust | [Theory_Foundations](.\demo3_Robust\Theory_Foundations.md)：公共理论基础（自适应控制 + H∞鲁棒控制）<br />[题1-连续CC](.\demo3_Robust\Derivation_Problem1_Jafari_ContinuousCC.md) · [题2-离散AVC](.\demo3_Robust\Derivation_Problem2_Jafari_DiscreteAVC.md)<br />[题3-H∞综合](.\demo3_Robust\Derivation_Problem3_HinfSynthesis.md) · [题4-真实声学路径](.\demo3_Robust\Derivation_Problem4_RealAcousticPath.md)<br />[时变频率实验](.\demo3_Robust\ExperimentReport_TimeVarying.md) |
|        | 可供使用的[系统模型](.\dataset\README.md)                |
|        | 可供使用的[激励信号以及生成原则](.\signal_excitation\AboutExcitation.md) |

## 项目进展

对于不同类型的控制器，总是有从入门到进阶的两种版本。对于ANC中的控制器设计也是同理。

### （1）入门部分：固定控制器设计

|              | 基于RST的ANC                                                 | 基于LQG的ANC                                                 | 基于鲁棒控制的固定ANC1                                       | 基于鲁棒控制的ANC2                                           |
| ------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| **理论基础** | （1）能够按照需求为被控系统选取主导极点/辅助极点，设计预先指定环节；<br />（2）掌握灵敏度函数与Bezout多项式的基本概念，能够使用函数求解Diophantine方程；<br />（3）仿真验证所设计的控制器是否符合要求 | （1）能够对带通滤波器的输入-输出数据辨识为ARMAX模型，并进一步转化为状态空间模型表示；<br />（2）针对状态空间模型选取合适的Q,R参数，并使用工具完成Riccati方程的求解；<br />（3）仿真验证所设计的控制器是否符合要求 | （1）了解H∞鲁棒控制理论中的1984 Approach以及1987 DGKF Approach，明确鲁棒控制问题的标准定义；<br />（2）针对被控对象传递函数模型，选取合适的灵敏度约束权值，并使用工具求解控制器；<br />（3）仿真验证所设计的控制器是否符合要求 | （1）了解Jafari& Ioannou(2015)研究中提出的控制器设计方案；<br />（2）设计合适的H∞优化问题，并使用数学工具求解控制器；<br />（3）仿真验证所设计的控制器是否符合要求 |
| **简化验证** | ✅在Carmona 2000年的论文中完整介绍了基于管道声学的ANC实验。使用该论文提供的模型，以及论文提供的参数设计，能够通过仿真验证。 | ✅在钱梵梵等人2022年的论文中介绍了LQG应用于管道ANC的实验，将其中被控对象模型简化为带通滤波器，能够通过仿真验证。 | ❓使用白明宪等人1997年论文中介绍的耳机系统，以及论文提供的控制器设计，在存在白噪声干扰的情况下使用ANC算法滤除噪声。没有使用MATLAB hinfsyn工具完整控制器设计 | ✅针对简化模型，能够通过仿真验证复现                          |
| **仿真验证** | ❓针对实际辨识的被控对象所设计的控制器暂且不能完全满足需求，仿真结果不理想。 | ✅对于辨识的30阶次级通路模型，完成了LQG控制器设计以及仿真验证。对于参数选取等未作进一步优化。 | ❌暂未针对实际辨识的被控对象作仿真验证                        | ✅对ARMAX(30,30,30,22)真实声学模型，F(z)+FIR固定控制器实现23dB抑制，NLMS自适应时变跟踪22dB |
| **实验验证** | ❌                                                            | ❌暂未将仿真结果迁移至实物实验平台。                          | ❌                                                            | ❌                                                            |

### （2）进阶部分：自适应控制器设计

|              | 基于RST的ANC                                                 | 基于LQG的ANC                                                 | 基于鲁棒控制的ANC2                                           |
| ------------ | ------------------------------------------------------------ | ------------------------------------------------------------ | ------------------------------------------------------------ |
| **理论基础** | 需在Carmona等人研究的基础上，增加Landau等人研究中提供的自适应环节 | 需在已完成的LQG实验基础上，增加钱梵梵等人(2022)研究中的自适应环节 | 需在已完成的固定控制器设计基础上，增加Jafari& Ioannou(2015)研究中的自适应环节 |
| **简化验证** | ✅针对Carmona 2000论文中的模型与参数设计，能够通过仿真验证    |                                                              |                                                              |
| **仿真验证** | ❓暂未针对实际辨识的被控对象作仿真验证                        | ❌仿真发散，实验未通过                                        | ✅固定FIR 23dB抑制 + NLMS自适应时变跟踪22dB（详见 [题4](.\demo3_Robust\Derivation_Problem4_RealAcousticPath.md) 和 [实验报告](.\demo3_Robust\ExperimentReport_TimeVarying.md)） |
| **实验验证** | ❌                                                            | ❌                                                            | ❌                                                            |

---

## Git 工作流规范

**原则**：仓库管理者可直接在 `master` 提交并推送，无需 PR 流程。

```bash
# === 管理者工作流 (推荐) ===
# 1. 拉取最新
git checkout master
git pull --rebase

# 2. 选择性 stage 并提交 (conventional commits)
git add <修改的文件>          # 不要用 git add . 或 git add -A
git commit -m "fix: <简述>"

# 3. 直接推送 master
git push origin master


# === 协作者工作流 (需 PR) ===
# 1. 从最新 master 创建 feature 分支
git checkout -b feature/<简短描述>

# 2. 选择性 stage 并提交
git add <修改的文件>
git commit -m "fix: <简述>"

# 3. 推送 feature 分支 → 创建 PR → Review → Merge
git push -u origin feature/<简短描述>

# 4. 清理本地分支
git checkout master && git pull --rebase
git branch -D feature/<简短描述>
```

**Commit Message 格式**：`<type>: <简述>`，常用 type：`fix`、`feat`、`refactor`、`docs`。

**反模式**：❌ `git add -A` · ❌ force push 已共享的分支 · ❌ 一个 commit 包含多个不相关修改
