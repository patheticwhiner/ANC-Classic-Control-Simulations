# <font color="#1a1a2e">审稿意见</font>

## <font color="#1a1a2e">demo2_LQG —— 学术发表可行性评估</font>

> <font color="#555555"><b>审稿人视角 · 内部评审</b>　|　基于 LQG 最优控制的管道 ANC 研究（含 ARMAX 辨识模型 + Q 参数化自适应扩展）</font>

---

## §1 总评

<font color="#e67e22"><b>████████████████████████████████████</b></font>
<font color="#e67e22"><b>　Major Revision（大修）　</b></font>
<font color="#e67e22"><b>████████████████████████████████████</b></font>

本工作实现了基于 LQG（LQR + Kalman 滤波）的管道 ANC 控制器设计，从入门级仿真（LQR 弹簧-质量-阻尼 → LQG 固定控制器 → ARMAX 辨识模型 → Q 参数化自适应扩展）构成了一条完整的递进式技术路线。工作覆盖了三种被控对象（合成 BPF 模型 $_\text{syn}$、ARMAX 辨识声学模型 $_\text{armax}$、弹簧-质量-阻尼教学模型），并对 LQG 在 ANC 中的关键问题——分离原理、鲁棒性损失、Youla 参数化的时变算子非对易性——逐一进行了探讨。

| <font color="#27ae60"><b>✅ 主要贡献</b></font> | <font color="#c0392b"><b>❌ 主要不足</b></font> |
|:---|:---|
| 在 ANC 场景中系统地实现了 LQG 从<u>教学模型 → 合成模型 → 真实辨识模型</u>的三级递进验证，每一级的技术决策（Q 矩阵选择、R 参数调节、Kalman 增益特征）都有明确记录 | ❶ <b>LQG 的核心缺陷未被解决</b>：文档明确承认"LQG 丧失 LQR 的鲁棒性裕度"但未提出任何补偿措施（LTR、H∞ 混合、鲁棒 Kalman）——这使工作停留在"发现问题"而非"解决问题" |
| <u>时变算子非对易性</u>分析（demo.m）是一个理论深度远超工程论文平均水平的观察——$T_{12}Q_k \neq Q_kT_{12}$ 的仿真验证 + 与 $d\theta/dt$ 的相关性分析构成了一个有发表潜力的微贡献 | ❷ <b>Q 参数化自适应模块未完成</b>：LQG_idmodel.m 仿真（三）的 RLS 更新律存在明显问题（$\epsilon$ 计算式看起来是占位符，分母包含 $P\phi$ 项但分子缺少正确误差信号），且无收敛性验证 |
| 在三种不同模型上的参数调优经验（合成模型 R=10e-4 → 辨识模型 R=1e-8）为 <u>ANC 工程实践</u>提供了可操作的参考 | ❸ <b>无任何基线方法对比</b>：LQG 相对于 FxLMS（标准 ANC 方案）、RST（本项目的 demo1）或 H∞（demo3）的性能优势未被量化。审稿人会追问："在 ANC 中，LQG 比简单方法好在哪里？" |

---

## §2 逐项评估

### §2.1 新颖性（Novelty）

| 技术要素 | 来源 | <font color="#1a1a2e">审稿意见</font> |
|:---|:---|:---|
| LQR + Kalman 滤波器分离原理 | <font color="#7f8c8d">标准教材（Åström, Anderson-Moore）</font> | <font color="#c0392b"><b>✗ 非原创</b></font>。LQG 的基础理论自 1960 年代以来已完全确立。 |
| LQG 应用于 ANC 增广系统 | <font color="#7f8c8d">钱梵梵等 (2022) 博士论文</font> | <font color="#c0392b"><b>✗ 非原创</b></font>。钱梵梵的工作已明确提出了在 ANC 中使用 LQG + Youla 参数化的框架。 |
| 时变算子 $T_{12}Q_k \neq Q_kT_{12}$ 的非对易性验证 | <font color="#2980b9"><b>作者独立完成</b></font> | <font color="#27ae60"><b>★ 主要新颖性</b></font>。虽然非对易性本身是线性时变系统理论中的已知事实，但在 <u>ANC + Q 参数化 + 有限阶自适应滤波器</u>这个具体交叉点上做定量仿真验证（含 $d\theta/dt$ 相关性分析）是本工作的原创贡献。这个微贡献适合扩展为 Signal Processing Letters 短文的核心论点。 |
| 三级递进式模型验证策略（教学→合成→真实辨识） | <font color="#2980b9"><b>作者独立完成</b></font> | <font color="#27ae60"><b>◈ 有指导价值</b></font>。这不是算法贡献，但作为<u>方法论/教学框架</u>有发表价值——它系统展示了如何从简单到复杂逐步验证一个 ANC 控制器的有效性。适合投 Control Systems Magazine 的教育专栏或作为长文的 Methodology 章节。 |
| LQR 权重 $Q = C^TC$ 在 ANC 增广系统中的使用 | <font color="#7f8c8d">标准 LQR 技巧</font> | <font color="#e67e22"><b>△ 增量贡献</b></font>。Q = C'C 是"使输出为零"的经典技巧，但在 ANC 增广系统（plant + disturbance 联合状态）中的具体效果和参数灵敏度分析有一定的新信息量。 |

<font color="#c0392b"><b>▶ 审稿意见：</b></font>本工作有两个可发表的"原子"——(a) 非对易性分析 → SPL 短文；(b) 三级递进验证框架 → 教学/方法论文章。但如果试图把 LQG for ANC 本身作为新颖性卖点发表，审稿人会指出钱梵梵 (2022) 已有先例。作者需要明确：这篇论文是在已有 LQG-for-ANC 框架下做增量贡献，还是提出一个新的分析视角？

### §2.2 理论完备性（Theoretical Rigor）

> <font color="#c0392b"><b>❶ 严重</b></font> &nbsp; <b><u>LQG 鲁棒性损失的补偿机制缺失</u></b>
>
> AboutLQGpt1.md 中正确指出了 LQG 的核心理论缺陷——"分离原理设计的 LQG 控制器不再享有 LQR 的 -6 dB 增益裕度和 60° 相位裕度保证"。但随后未讨论任何补偿方案。在 ANC 场景中，次级通路建模误差是不可避免的，这意味着纯 LQG 在真实系统中很可能不稳定或性能严重退化。
>
> <font color="#2980b9"><b>→</b></font> 至少应讨论以下方案之一并实现：(a) LTR (Loop Transfer Recovery) — 通过设计 Kalman 滤波器使回路传递函数恢复到 LQR 的健壮特性；(b) 将建模不确定性纳入 Kalman 滤波器的过程噪声协方差（鲁棒 Kalman）；(c) 明确标注本工作为"理想模型下的性能上限研究"，在 Limitations 中限制适用范围。

> <font color="#e67e22"><b>❷ 中等</b></font> &nbsp; <b><u>Q 参数化自适应模块的数学推导缺失</u></b>
>
> LQG_idmodel.m 仿真（三）的 Youla 参数化尝试是工作的重要扩展方向，但当前代码中的 RLS 更新律：
>
> $\epsilon = (y + 0) / \text{denom}$
>
> 明显是占位符（分子中的 +0 暗示误差项尚未定义）。正确的 Youla-Q 自适应推导应包含：Youla 参数化控制器参数化形式、等效误差模型的构建、以及 $\epsilon$ 与残差信号 $e(k)$ 之间的正确映射关系。
>
> <font color="#2980b9"><b>→</b></font> 补全 Q-参数化自适应 RLS 的完整推导（3-5 页），并验证：固定参数时自适应应收敛到 LQG 最优解。如果短期内无法完成，建议将仿真（三）从主文中移除，作为"Future Work"章节的内容。

> <font color="#e67e22"><b>❸ 中等</b></font> &nbsp; <b><u>Kalman 滤波器的噪声协方差缺乏物理标定</u></b>
>
> 当前 $Q_n$ 和 $R_n$ 的取值（$Q_n=2000$, $R_n=10^{-4}$ 或 $Q_n=2$, $R_n=10^{-10}$）是纯数值调参的结果，未与物理量建立联系。在 ANC 中，$Q_n$ 应反映扰动源的实际功率谱密度，$R_n$ 应反映传声器的本底噪声。
>
> <font color="#2980b9"><b>→</b></font> 至少应在仿真脚本注释中说明：实际传声器本底噪声约为 [X] Pa，对应 $R_n$ 约为 [Y]；或添加一个参数灵敏度曲线（$Q_n$ / $R_n$ vs. 抑制比）。

> <font color="#7f8c8d"><b>❹ 轻微</b></font> &nbsp; <b><u>连续/离散混淆风险</u></b>
>
> AboutLQGpt1.md 的理论推导部分使用连续时间符号 $\int$，但实际实现全部是离散时间（`dlqr`, `dlqe`, $T_s=1/1000$）。这种混搭可能导致读者对推导与实现的对应关系产生困惑。
>
> <font color="#2980b9"><b>→</b></font> 统一使用离散时间符号推导，或在理论章节开头声明"以连续时间推导原理，以离散时间实现"。

### §2.3 实验验证（Experimental Validation）

<font color="#27ae60"><b>✅ 已有内容：</b></font>
- ✅ LQR 在弹簧-质量-阻尼教学模型上的完整验证（稳定性 + 稳态误差 + 控制能量）
- ✅ LQG 在合成 BPF 模型上的两种实现方式对比（`dlqr+dlqe` vs `lqg()` 一体函数）
- ✅ LQG 在 ARMAX 辨识模型上的验证（含分离 plant/disturbance 状态的两种仿真结构）
- ✅ 时变算子非对易性的定量验证（RMS 差 + 相关系数 + 频谱分析）
- ✅ 控制器设计工具函数的参数化（`design_observer_controller.m` 接受可调 Q/R 权值）

<font color="#c0392b"><b>❌ 缺失内容</b></font>（按重要性排序）：

> <font color="#c0392b"><b>❶ 必需</b></font> &nbsp; <b><u>与至少两种基线方法的定量对比</u></b>
>
> 管道 ANC 的事实标准是 FxLMS（自适应前馈）及其变种（NLMS、Filtered-X RLS）。本项目的 demo3_Robust 已有 NLMS 实现。建议对比方案：
>
> | 方法 | 来源 | 对比维度 |
> |:---|:---|:---|
> | FxLMS (N=64, μ=0.05) | demo3 Part 1 | 稳态抑制、收敛速度 |
> | RST(K2) 鲁棒控制器 | demo1 | 控制器阶数、对 ±10% 频偏的鲁棒性 |
> | LQG (本工作) | — | 基准线 |
>
> 至少完成一张对比表。

> <font color="#e67e22"><b>❷ 强烈建议</b></font> &nbsp; <b><u>参数灵敏度系统性分析</u></b>
>
> 当前代码中存在大量手工调参（$Q_{lqr}$, $R_{lqr}$, $Q_n$, $R_n$, $\lambda$）但未给出参数灵敏度曲线。建议对 ARMAX 模型场景生成：
> - $R_{lqr}$ （控制代价）vs. 抑制比 (dB) 的 Pareto 前沿
> - $Q_n$ / $R_n$ 比值 vs. 估计误差 RMS
>
> 这 2 张图本身就是有价值的发表材料——"LQG 在 ANC 中的参数调优指南"。

> <font color="#7f8c8d"><b>❸ 可选</b></font> &nbsp; <b><u>Simulink 模型与实际物理系统的对比</u></b>
>
> 项目根目录有 `simLQG.slx`，但当前所有仿真均在 `.m` 脚本中手工实现。如果 Simulink 模型代表的是连续时间/混合信号设计，应至少对比 Simulink 仿真与离散脚本仿真的差异（例如量化效应、连续/离散间采样保持效应）。

### §2.4 表述与结构（Presentation）

| <font color="#27ae60"><b>✅ 优势</b></font> | <font color="#c0392b"><b>⚠ 从学术论文视角看的问题</b></font> |
|:---|:---|
| 三级递进式结构（教学→合成→真实）清晰展示了验证方法论，每个级别都有独立可运行的脚本 | • <b>核心文档碎片化</b>：AboutLQR / AboutLQGpt1 / AboutLQGpt2 三份文档之间存在信息孤岛——LQR 文档讲弹簧-质量、pt1 讲 LQG 理论、pt2 只有两行占位。读者需要在三份文档间跳跃才能理解完整图景 |
| 代码注释包含大量操作层面的工程判断（"注意：R 取太大系统无抑制效果，R 取太小则控制信号饱和"） | • <b>缺少 Related Work 章节</b>：未定位钱梵梵 (2022) 的具体贡献边界——本工作相对于钱梵梵的增量是什么？复现？扩展时变部分？还是引入新的非对易性分析？ |
| 对 LQG 的局限性保持诚实——多处代码注释和文档明确提到鲁棒性损失、计算量问题、Matlab 版本兼容性问题 | • <b>Q 参数化部分的不完整破坏了叙事的完整性</b>：读者沿着 LQR → LQG → LQG+ARMAX → LQG+Adaptive 的自然递进路线阅读，在最后一级遇到一个未完成的仿真（三）——这在投稿论文中是不可接受的 |
|  | • <b>$Q = C^TC$ 的使用没有被形式化</b>：这是一个在 LQR 文献中广为人知的技巧，但第一次接触的读者会困惑"为什么状态权重恰好是 $C^TC$？"应给出推导：$\|y\|^2 = x^TC^TCx$ 等价于状态权重 $Q = C^TC$ |

---

## §3 按期刊层级的可行性评估

| 期刊 / 会议 | <font color="#1a1a2e">可行性</font> | 所需追加工作 | 预期篇幅 |
|:---|:---:|:---|:---:|
| <b>IEEE TAC / Automatica</b> | <font color="#c0392b"><b>✗ 不可行</b></font> | — | — |
| <b>IEEE/ACM TASLP</b> | <font color="#c0392b"><b>✗ 不可行</b></font> | — | — |
| <b>MSSP</b> | <font color="#e67e22"><b>△ 困难</b></font> | 完成 Q 参数化自适应 + 基线对比 + 参数灵敏度 + 鲁棒性补偿 (LTR) | 12–14 页 |
| <b><u>Applied Acoustics</u></b> | <font color="#e67e22"><b>△ 有可能</b></font> | 基线对比表 + 参数灵敏度 + 整合三份 About 文档为完整论文 | 8–10 页 |
| <b><u>IEEE Signal Processing Letters</u></b> | <font color="#2980b9"><b>◈ 最可行的期刊目标</b></font> | 聚焦非对易性分析，精简为一个清晰的论点 + 1 个仿真实验 | 4–5 页 |
| <b>InterNoise / ICSV</b> | <font color="#27ae60"><b>● 可直接投稿</b></font> | 补全 LQG vs FxLMS 对比实验 + 整合文档 | 4–6 页 |

> <font color="#2980b9"><b>▸ 推荐策略：</b></font>将工作切为两篇：(1) 非对易性分析 → SPL 短文（快速发表）；(2) 三级递进 LQG 验证方法论 + 参数灵敏度 → Applied Acoustics 长文。不要试图把所有内容塞进一篇——LQG 基础理论不是原创，非对易性分析才是。

---

## §4 修改建议（按优先级排列）

### <font color="#c0392b">优先级 1</font> &nbsp;·&nbsp; <font color="#c0392b"><b>必须完成</b></font> <small>（否则不建议投稿任何期刊）</small>

> <font color="#c0392b"><b>①</b></font> &nbsp; <b><u>完成或移除 Q 参数化自适应模块</u></b>
>
> 仿真（三）的 RLS 更新律需完整推导和正确实现。如果短期内无法完成，将其从主文中移除，放入 Future Work 章节——一篇投稿论文中不能有未完成的实验。

> <font color="#c0392b"><b>②</b></font> &nbsp; <b><u>明确相对于钱梵梵 (2022) 的增量贡献</u></b>
>
> 钱梵梵的博士论文已提出 LQG + Youla 参数化用于 ANC。本工作必须明确声明：哪些是复现、哪些是独立贡献。如果非对易性分析是唯一的独立贡献，则应聚焦于此。

> <font color="#c0392b"><b>③</b></font> &nbsp; <b><u>合并三份 About 文档为一份完整的技术报告</u></b>
>
> AboutLQR + AboutLQGpt1 + AboutLQGpt2 → 一份 `AboutLQG.md`，按标准论文结构（Introduction → Problem Formulation → Method → Experiments → Discussion）组织。

### <font color="#e67e22">优先级 2</font> &nbsp;·&nbsp; <font color="#e67e22"><b>强烈建议</b></font> <small>（Applied Acoustics / SPL 级别需要）</small>

> <font color="#e67e22"><b>④</b></font> &nbsp; <b><u>新增 LQG vs FxLMS 定量对比表</u></b>
>
> 同一 ARMAX 模型 + 同一扰动（固定频率 + 扫频），对比稳态抑制 (dB)、收敛时间 (s)、计算量/采样点。

> <font color="#e67e22"><b>⑤</b></font> &nbsp; <b><u>补充参数灵敏度分析</u></b>
>
> 生成 $R_{lqr}$ vs. 抑制比的 Pareto 前沿（展示控制代价与性能的权衡）、$Q_n/R_n$ vs. 估计误差 RMS。

> <font color="#e67e22"><b>⑥</b></font> &nbsp; <b><u>讨论 LTR 或等效鲁棒性恢复方案</u></b>
>
> 不一定要实现，但至少需要在 Discussion 中讨论——因为这是 LQG 在工程应用中最大的障碍。

### <font color="#7f8c8d">优先级 3</font> &nbsp;·&nbsp; <font color="#7f8c8d"><b>锦上添花</b></font> <small>（更高级别期刊需要）</small>

> <font color="#7f8c8d"><b>⑦</b></font> Simulink 连续/离散对比（量化效应、采样保持效应）
>
> <font color="#7f8c8d"><b>⑧</b></font> 控制器降阶（`design_observer_controller.m` 的 30 阶实现 → 平衡截断 → 性能损失量化）
>
> <font color="#7f8c8d"><b>⑨</b></font> 硬件实现可行性讨论（定点运算、实时性约束下的 Kalman 增益计算复杂度）

---

## §5 最终建议

<font color="#1a1a2e"><b>给作者的中肯总结：</b></font>

本工作最有价值的部分<font color="#c0392b"><b>不是</b></font> LQG 在 ANC 中的应用本身——钱梵梵 (2022) 已覆盖了这一框架——
也<font color="#c0392b"><b>不是</b></font> LQR/Kalman 的参数调优——那是工程实践而非可发表的研究贡献——
<font color="#27ae60"><b>而是</b></font> 时变算子 $T_{12}Q_k \neq Q_kT_{12}$ 的非对易性分析，以及三级递进式验证方法论。

> <i>大多数 ANC 论文将 Q 参数化自适应中的非对易性视为数值误差忽略不计。本工作用量化的相关性分析（$d\theta/dt$ vs. 误差 RMS）证明这是一个可测量、可分析的系统性效应——并且在参数快速变化时不能忽略。这是一个微小但真实的学术发现。</i>

但要从目前的<font color="#e67e22"><b>"方法复现 + 部分原创分析"</b></font>升级为<font color="#2980b9"><b>"有清晰学术贡献的论文"</b></font>，你需要做两个决定：

<font color="#c0392b"><b>　切分 + 聚焦。</b></font>

不要把 LQG 基础理论、三级递进验证、非对易性分析、Q 参数化自适应全部塞进一篇论文。切分为：(a) 非对易性分析 → SPL 短文（1 个清晰论点 + 1 个实验）；(b) 三级 LQG 验证方法论 → Applied Acoustics 长文（补基线对比 + 参数灵敏度）。

一旦你做出了这个切分，你就不只是在说：

> <font color="#7f8c8d">"我们用 LQG 方法在三个模型上做了 ANC，效果还不错，我们还发现了一个非对易性的问题"</font>

你是在说：

> <font color="#27ae60"><b>"在 Youla 参数化的 ANC 自适应框架中，当自适应速率与系统动态可比时，时变算子的非对易性是一个不可忽略的误差源——我们通过定量分析和仿真给出了这个效应的首次系统表征。"</b></font>

<font color="#c0392b"><b>这才是审稿人会认可的贡献。</b></font>

<font color="#e67e22"><b>████████████████████████████████████████████████████</b></font>
<font color="#e67e22"><b>　审稿决定　</b></font><b><u><font color="#d35400">Major Revision</font></u></b>
<font color="#e67e22"><b>████████████████████████████████████████████████████</b></font>

> <font color="#555555"><b>复审要求：</b>(1) 完成或移除 Q 参数化自适应模块；(2) 明确相对于钱梵梵 (2022) 的增量贡献；(3) 将非对易性分析扩展为可独立投稿的 SPL 短文草稿（如有意投稿 SPL）。InterNoise 可降低至仅需 LQG vs FxLMS 对比 + 整合文档。</font>

<p align="right"><font color="#95a5a6"><small>本审稿意见基于对 LQRdemo.m、LQGdemo.m、LQGdemo2.m、LQG_idmodel.m、demo.m、design_observer_controller.m、AboutLQR.md、AboutLQGpt1.md、AboutLQGpt2.md 及项目 README.md 的全文阅读。</small></font></p>
