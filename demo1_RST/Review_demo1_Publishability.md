# <font color="#1a1a2e">审稿意见</font>

## <font color="#1a1a2e">demo1_RST —— 学术发表可行性评估</font>

> <font color="#555555"><b>审稿人视角 · 内部评审</b>　|　基于 RST 极点配置与灵敏度塑形的管道 ANC 研究</font>

---

## §1 总评

<font color="#e67e22"><b>████████████████████████████████████</b></font>
<font color="#e67e22"><b>　Major Revision（大修）　</b></font>
<font color="#e67e22"><b>████████████████████████████████████</b></font>

本工作实现了基于 RST 多项式控制器结构的管道 ANC 系统，包含固定控制器（Carmona 2000 灵敏度函数塑形）与自适应控制器（Landau 2005 自适应极点配置）两条技术路线。固定控制器展示了从标称性能（K1, ~55 dB 抑制）到鲁棒性能（K2, ~18 dB 抑制 @±10% 频偏）的合理退化过程；自适应方案验证了 RLS/LMS 对时变频率的跟踪能力。

| <font color="#27ae60"><b>✅ 主要贡献</b></font> | <font color="#c0392b"><b>❌ 主要不足</b></font> |
|:---|:---|
| 将鲁棒性约束<u>转化为灵敏度频率模板</u>并嵌入极点配置预指定多项式——这是 Carmona 2000 的核心工程智慧被忠实地复现和记录 | ❶ <b>无基线对比</b>：未与任何其他控制方法（H∞、LQG、自适应 FIR 前馈）做定量比较 |
| 在自适应 RST 框架下同时实现了 <u>RLS 和 LMS</u> 两种参数适应律，并对比了它们在时变频率场景下的行为差异 | ❷ <b>理论贡献声明不清</b>：目前是"忠实复现"而非"提出新方法"。对本工作在 Carmona/Landau 原论文基础上的<u>增量贡献</u>没有明确界定 |
| <u>非对易性分析</u>（demo.m）切入了一个控制理论中真实存在的细粒度问题——时变算子的非交换性——这是大多数工程论文会忽略的理论观察 | ❸ <b>实验验证仅停留在仿真复现层面</b>：三种控制器（K1/K2/K3）均未在真实的 ARMAX 辨识模型上验证，自适应部分也未在实测数据上测试 |

---

## §2 逐项评估

### §2.1 新颖性（Novelty）

| 技术要素 | 来源 | <font color="#1a1a2e">审稿意见</font> |
|:---|:---|:---|
| RST 极点配置 + 灵敏度函数塑形 | <font color="#7f8c8d">Carmona & Alvarado (2000)</font> | <font color="#c0392b"><b>✗ 非原创</b></font>。完整的设计流程（Bezout 方程求解、三控制器设计 K1/K2/K3）均在原论文中给出。 |
| 灵敏度频率模板与鲁棒性"软化" | <font color="#7f8c8d">Carmona 2000 + Landau 2005</font> | <font color="#e67e22"><b>△ 增量贡献</b></font>。本工作用清晰的文档和代码对"鲁棒约束→灵敏度模板→预指定多项式"这一设计链路做了<u>可复现的拆解</u>，原论文中分散的工程考量被系统地串联起来。 |
| RLS vs LMS 自适应律对比 | <font color="#2980b9"><b>作者独立完成</b></font> | <font color="#27ae60"><b>★ 主要新颖性</b></font>。在自适应 RST 框架内对 RLS 和 LMS 做并排对比，特别是在时变频场景（扫频 + FM 调制）下的行为差异分析，这在已有文献中不常见。 |
| RST 控制器设计的数值稳定性分析 | <font color="#7f8c8d">Landau 2005 提及</font> | <font color="#e67e22"><b>△ 增量贡献</b></font>。对角化 Sylvester 矩阵病态化、RLS 协方差矩阵奇异的诊断具有工程参考价值，但未形成系统性分析。 |

<font color="#c0392b"><b>▶ 审稿意见：</b></font>本工作的核心新颖性不在于 RST 方法本身——这些都是已有经典——而在于 (a) 鲁棒约束嵌入极点配置的工程链路的可复现拆解，以及 (b) RLS vs LMS 在自适应 RST 框架下的并排对比。当前文档未明确聚焦这两个点，导致阅读时难以判断"什么是作者的原创贡献"。

### §2.2 理论完备性（Theoretical Rigor）

> <font color="#c0392b"><b>❶ 严重</b></font> &nbsp; <b><u>自适应 RST 的稳定性分析缺失</u></b>
>
> Landau 2005 自适应极点配置的理论保证依赖于"被控对象参数可辨识"和"持续激励条件"两个假设。本工作在时变频率场景下运行自适应算法，但未检查持续激励条件是否满足，也未分析参数收敛是否保证闭环稳定。
>
> <font color="#2980b9"><b>→</b></font> 至少需要：(a) 分析回归向量 $\Phi_k$ 的谱特性以验证持续激励；(b) 在时变频场景中画出闭环极点随时间的轨迹，验证其始终在单位圆内。

> <font color="#e67e22"><b>❷ 中等</b></font> &nbsp; <b><u>LMS 替代 RLS 的理论理由不充分</u></b>
>
> 文档中提到"RLS 协方差矩阵容易奇异"因此尝试 LMS，但 LMS 在有色噪声（ANC 场景的典型特征）下的收敛速度远慢于 RLS——这正是 RLS 被提出的原因。代码中 LMS 步长需缩小到 $10^{-7}$ 才不发散，暗示其在该场景下可能不适用。
>
> <font color="#2980b9"><b>→</b></font> 应给出 LMS 收敛步长的理论上界（与回归向量最大特征值的关系），或改用归一化 LMS (NLMS) 并解释理由。

> <font color="#7f8c8d"><b>❸ 轻微</b></font> &nbsp; <b><u>灵敏度函数与 H∞ 范数的等价关系未形式化</u></b>
>
> 文档中提到模量裕度 $\Delta_M = 1/\|S_{ye}\|_\infty$，以及延时裕度与灵敏度模板的关系，但引用的关键推导（Ref [29]）未被包含在文档中。这使得"软化鲁棒性"这一核心设计决策的理论依据不完整。
>
> <font color="#2980b9"><b>→</b></font> 在 AboutRST.md 或 Derivation 文档中补充该等价关系的推导（2–3 页即可），使其成为自包含的理论基础。

### §2.3 实验验证（Experimental Validation）

<font color="#27ae60"><b>✅ 已有内容：</b></font>
- ✅ 三控制器设计（K1 标称 / K2 鲁棒 / K3 宽带）均在简化模型上仿真成功
- ✅ RLS 与 LMS 自适应律在时变频率下的行为对比（扫频 + FM 调制）
- ✅ 灵敏度函数频率响应、Nyquist 图等频域分析完备
- ✅ 时变频场景的 spectrogram 可视化

<font color="#c0392b"><b>❌ 缺失内容</b></font>（按重要性排序）：

> <font color="#c0392b"><b>❶ 必需</b></font> &nbsp; <b><u>与至少一种基线方法的定量对比</u></b>
>
> 这是论文被接收的最低门槛。对于管道 ANC，公认的基线包括：自适应 FIR 前馈（FxLMS 或其变种）、H∞ 鲁棒控制（Bai 1997 方法，本项目 demo3 中已有实现）。没有对比表（抑制比 × 收敛时间 × 控制器阶数 × 计算量），审稿人无法判断 RST 方法相对于更简单或更现代方法是否有优势。
>
> 所需实验：同一被控对象 + 同一扰动场景（至少固定频率 + 扫频两种），对比 RST(K2) vs FxLMS vs H∞。指标：稳态抑制 (dB)、收敛时间 (s)、控制器阶数、浮点运算量/采样点。

> <font color="#e67e22"><b>❷ 强烈建议</b></font> &nbsp; <b><u>在真实 ARMAX 辨识模型上验证</u></b>
>
> 当前固定控制器（K1/K2/K3）使用的是 Carmona 原文的 7 阶 duct 模型（$f_s=2000$ Hz）。本项目已有真实声学管道的 30 阶 ARMAX 模型（$f_s=48000$ Hz）。在这两个模型上测试同一设计流程——差异有多大？鲁棒性"软化"在真实模型上是否仍然有效？这个问题直接关系到方法的工程可信度。
>
> 建议：在 `armax_30303022_2026-01-20.mat` 上重复 K2（鲁棒）控制器设计，报告两个模型上的抑制差异。

> <font color="#7f8c8d"><b>❸ 可选</b></font> &nbsp; <b><u>带建模不确定性的 Monte Carlo 鲁棒性测试</u></b>
>
> 对 plant 多项式的系数施加随机扰动（如 ±5% 均匀分布），运行 100 次 Monte Carlo，统计闭环稳定率和抑制比分布。这是验证鲁棒性"软化"方法有效性的直接证据。

### §2.4 表述与结构（Presentation）

| <font color="#27ae60"><b>✅ 优势</b></font> | <font color="#c0392b"><b>⚠ 从学术论文视角看的问题</b></font> |
|:---|:---|
| 代码内的 `%[text]` 注释富含工程洞察（如"控制器阶数 20-30 已是实时系统上限"、"Nyquist 图中模量裕度的几何含义"） | • <b>文档是模板空壳</b>：AboutRST.md 中大量章节（理论基础、仿真验证、实物实验验证、改进建议、参考资料）为空白占位符，无法作为独立技术文档阅读 |
| 三个控制器的设计链路（K1 → K2 → K3）形成自然的叙事弧：标称 → 鲁棒 → 宽带，教学意图清晰 | • <b>缺少 Problem Formulation 章节</b>：读者需要从代码逻辑反推"我们要解决什么问题"。应先声明"给定管道系统 $G(z)$ 和扰动频谱特征，设计控制器 $K(z)$ 使得..." |
| 数值诊断（Sylvester 矩阵条件数、RLS 协方差奇异性检测）穿插在代码中，显示作者对实现细节有扎实理解 | • <b>贡献声明分散</b>：RLS vs LMS 对比、非对易性分析、鲁棒性软化——这些有潜力的贡献点散落在不同文件的注释中，缺少一处集中声明 |
|  | • <b>三控制器各自独立成章但缺少横向比较</b>：K1/K2/K3 的数据（抑制比、控制器阶数、计算时间）未汇总为一张对比表 |

---

## §3 按期刊层级的可行性评估

| 期刊 / 会议 | <font color="#1a1a2e">可行性</font> | 所需追加工作 | 预期篇幅 |
|:---|:---:|:---|:---:|
| <b>IEEE TAC / Automatica</b> | <font color="#c0392b"><b>✗ 不可行</b></font> | — | — |
| <b>IEEE/ACM TASLP</b> | <font color="#c0392b"><b>✗ 不可行</b></font> | — | — |
| <b>MSSP (Mechanical Systems & Signal Processing)</b> | <font color="#e67e22"><b>△ 困难</b></font> | 需完整对比实验（≥3 基线）+ 真实模型验证 + Monte Carlo 鲁棒性 | 12–14 页 |
| <b><u>Applied Acoustics</u></b> | <font color="#2980b9"><b>◈ 最可行的期刊目标</b></font> | RLS vs LMS 对比正式化 + ARMAX 模型验证 + 1 个基线对比 | 8–10 页 |
| <b><u>InterNoise / ICSV</u></b> | <font color="#27ae60"><b>● 可直接投稿</b></font> | 精简为 1 个清晰对比（RLS vs LMS in RST for ANC），去除未完成章节 | 4–6 页 |

> <font color="#2980b9"><b>▸ 推荐策略：</b></font>先投 InterNoise 2026 短文（聚焦 RLS vs LMS 自适应 RST 对比），获得同行反馈后，补充基线实验扩展为 Applied Acoustics 长文。

---

## §4 修改建议（按优先级排列）

### <font color="#c0392b">优先级 1</font> &nbsp;·&nbsp; <font color="#c0392b"><b>必须完成</b></font> <small>（否则不建议投稿任何期刊）</small>

> <font color="#c0392b"><b>①</b></font> &nbsp; <b><u>填补 AboutRST.md 所有空白章节</u></b>
>
> 至少填充：理论基础（含关键公式）、仿真验证结果（定量数据）、改进建议。当前的空模板状态意味着没有独立的学术文档可用。

> <font color="#c0392b"><b>②</b></font> &nbsp; <b><u>明确本工作的增量贡献声明</u></b>
>
> 在 AboutRST.md 开头新增一段不超过 200 字的 Contribution 声明，回答："Carmona/Landau 做了 X，我们在此基础上做了 Y，Y 的独特之处在于 Z。" 如果 Y = 忠实复现，则投稿定位为"教学案例报告"而非"研究论文"——这在某些教学期刊（如 IEEE Control Systems Magazine 的教育专栏）是有价值的，但需要换投稿策略。

### <font color="#e67e22">优先级 2</font> &nbsp;·&nbsp; <font color="#e67e22"><b>强烈建议</b></font> <small>（Applied Acoustics 级别需要）</small>

> <font color="#e67e22"><b>③</b></font> &nbsp; <b><u>新增基线方法对比表</u></b>
>
> 至少对比 FxLMS（本项目的 demo3_Robust 中已有 NLMS 实现）。定量指标：稳态抑制 (dB)、收敛时间 (s)、控制器阶数。

> <font color="#e67e22"><b>④</b></font> &nbsp; <b><u>在真实 ARMAX 模型上验证 K2 控制器</u></b>
>
> 使用 `armax_30303022_2026-01-20.mat` 重跑 K2 设计并报告抑制差异。

> <font color="#e67e22"><b>⑤</b></font> &nbsp; <b><u>补充持续激励条件验证</u></b>
>
> 对自适应 RST 仿真，添加回归向量条件数或最小特征值的监控图。

### <font color="#7f8c8d">优先级 3</font> &nbsp;·&nbsp; <font color="#7f8c8d"><b>锦上添花</b></font> <small>（更高级别期刊需要）</small>

> <font color="#7f8c8d"><b>⑥</b></font> Monte Carlo 鲁棒性测试（100 次随机 plant 参数扰动）
>
> <font color="#7f8c8d"><b>⑦</b></font> 控制器降阶研究：30 阶 RST 控制器 → 低阶近似 → 性能损失量化
>
> <font color="#7f8c8d"><b>⑧</b></font> 硬件在环 (HIL) 验证：在 dSPACE 平台上部署并测试最简控制器 K1

---

## §5 最终建议

<font color="#1a1a2e"><b>给作者的中肯总结：</b></font>

本工作最有价值的部分<font color="#c0392b"><b>不是</b></font> RST 极点配置本身——那是 Carmona 和 Landau 的经典贡献——
也<font color="#c0392b"><b>不是</b></font> 灵敏度函数塑形——那是 1980 年代以来鲁棒控制的标配工具——
<font color="#27ae60"><b>而是</b></font> 用可复现的代码将"鲁棒约束 → 灵敏度模板 → 预指定多项式"这条设计链路做了完整的工程拆解，同时在自适应环节中并排对比了 RLS 与 LMS 的行为差异。

> <i>"Carmona 2000 的论文告诉你灵敏度模板长什么样，但没有告诉你怎么从零设计出这个模板。本工作补上了中间的工程步骤——而且是可运行的，不是只画了框图。"</i>

但要从目前的<font color="#e67e22"><b>"复现验证 + 教学案例"</b></font>升级为<font color="#2980b9"><b>"可发表的研究论文"</b></font>，你需要补上那一个缺失的环节：

<font color="#c0392b"><b>　定量对比。</b></font>

不与其他方法对比，审稿人永远可以问同一个问题："为什么我要用 RST 而不是 FxLMS（更简单）或 H∞（更现代）？" 你不能只在结论里说"RST 有工程上的优势"——你需要一张表，同一对象、同一扰动、不同方法，数字说话。

一旦你完成了对比实验，你就不只是在报告：

> <font color="#7f8c8d">"我们成功地复现了 Carmona 和 Landau 的方法并写成了可运行的代码"</font>

你是在说：

> <font color="#27ae60"><b>"在管道 ANC 的特定场景下，自适应 RST 在 [X] 指标上优于 FxLMS，在 [Y] 条件下接近 H∞ 但控制器阶数降低 [Z]%——并且我们通过实验验证了这些结论。"</b></font>

<font color="#c0392b"><b>这才是审稿人会认可的贡献。</b></font>

<font color="#e67e22"><b>████████████████████████████████████████████████████</b></font>
<font color="#e67e22"><b>　审稿决定　</b></font><b><u><font color="#d35400">Major Revision</font></u></b>
<font color="#e67e22"><b>████████████████████████████████████████████████████</b></font>

> <font color="#555555"><b>复审要求：</b>补充 AboutRST.md 全部空白章节 + 明确增量贡献声明 + 至少一项基线方法对比实验（FxLMS 或 H∞）。InterNoise 短文可降低要求至仅需 RLS vs LMS 对比 + 完整文档。</font>

<p align="right"><font color="#95a5a6"><small>本审稿意见基于对 Carmona2000.m (715行)、Landau2005.m (491行)、AboutRST.md、相关 functions/ 及项目 README.md 的全文阅读。</small></font></p>
