# <font color="#1a1a2e">审稿意见</font>

## <font color="#1a1a2e">题4「真实声学次级通路的鲁棒自适应 ANC」—— 学术发表可行性评估</font>

> <font color="#555555"><b>审稿人视角 · 内部评审</b>　|　本文档以期刊审稿人的语言风格，对题4的技术内容进行系统性评估，指出其在新颖性、理论完备性、实验验证三个维度上的<u>优势</u>与<u>不足</u>，并给出按期刊层级分级的修改建议。</font>

---

## <font color="#1a1a2e">1. 总评</font>

<font color="#e67e22"><b>████████████████████████████████████████████████████████████</b></font>

<font color="#e67e22">　<b>推荐意见</b>　</font><b><u><font color="#d35400">Major Revision</font></u></b><font color="#e67e22">　（大修后重投）</font>

<font color="#e67e22"><b>████████████████████████████████████████████████████████████</b></font>

<br>

本工作将 Jafari & Ioannou (2015/2017) 的连续时间鲁棒自适应 ANC 框架推广至离散时间域，针对一条实测声学管道的高阶 ARMAX 模型（<b>30 阶，9 个 NMP 零点，$f_s=48$ kHz</b>）完成了 $F(z)$ 频谱展平滤波器设计、固定 FIR 控制器的 $H_\infty$ 插值综合、以及自适应 RLS 控制器的实施与调试。作者在调试过程中发现了<font color="#c0392b"><b>高采样率下 RLS 回归向量的共线性问题</b></font>，提出了抽取式 RLS 的解决方案，并在时变频率实验中对比了 RLS 与 NLMS 的跟踪性能。

| <font color="#27ae60"><b>✅ 主要贡献</b></font> | <font color="#c0392b"><b>❌ 主要不足</b></font> |
|:---|:---|
| 集中在<u>工程诊断层面</u>——对 RLS 在 $f_s \gg f_0$ 条件下的数值失效给出了清晰的因果链分析，并提供了一个简单有效的补救措施（抽取 + 对角加载 + 参数重调）。这一<b>诊断逻辑的严密性</b>是本文最值得肯定的部分。 | ❶ 核心方法论（$F+K$ 分解、$H_\infty$ 插值）<b>非作者原创</b><br>❷ 抽取式 RLS <b>缺乏收敛性理论分析</b><br>❸ 实验局限于<b>单条声学路径</b>和<b>单频扰动</b>，且未与 ANC 领域的标准基线方法（FXLMS 族）进行系统对比 |

---

## <font color="#1a1a2e">2. 逐项评估</font>

### <font color="#2980b9">2.1 新颖性 &nbsp;·&nbsp; <small>NOVELTY</small></font>

| 技术要素 | 来源 | <font color="#1a1a2e">审稿意见</font> |
|:---|:---|:---|
| $F(z) + K(z,\theta)$ 控制架构 | <font color="#7f8c8d">Jafari & Ioannou (2015, 2017)</font> | <font color="#c0392b"><b>✗ 非原创</b></font>。应明确标注为「已有架构的应用」，不可声称原创。 |
| $F(z)$ 离散时间零点反射设计 | <font color="#7f8c8d">连续时间→离散时间的自然推广</font> | <font color="#e67e22"><b>△ 增量贡献</b></font>。$z_i \to 1/\bar{z}_i$ 是本科 DSP 标准内容。若作为独立贡献，不足以支撑一篇期刊论文。 |
| FIR + $H_\infty$ 插值固定控制器 | <font color="#7f8c8d">题2 已完整建立（SOCP + LS 最小范数）</font> | <font color="#c0392b"><b>✗ 非新贡献</b></font>。仅为将同一方法应用于不同的 $G_{\text{eff}} = G_0F$。 |
| <b>高采样率下 RLS 回归向量共线性诊断</b> | <font color="#2980b9"><b>作者独立完成</b></font> | <font color="#27ae60"><b>★ 主要新颖性所在</b></font>。对 $\phi$ 相邻元素相关度 $>0.999$ 导致信息矩阵近乎秩亏的识别，通过<u>逐一排除 $G_0F$ 近似误差、闭环反馈、投影过紧</u>三个替代假说的诊断流程，构成有价值的工程洞察。 |
| 抽取式 RLS 方案 + $M$ 选择准则 | <font color="#2980b9"><b>作者提出</b></font> | <font color="#e67e22"><b>⚠ 有保留</b></font>。实用性强，但 $M$ 依赖于 $f_0$——在自适应场景中 $f_0$ 正是未知量，构成<u>方法论的循环依赖</u>。 |
| NLMS vs RLS 扫频对比实验 | <font color="#2980b9"><b>作者独立完成</b></font> | <font color="#2980b9"><b>◈ 有指导价值</b></font>。RLS→稳态辨识，NLMS→非平稳跟踪，但 NLMS 是成熟算法，非贡献点。 |

> <font color="#c0392b"><b>▶ 审稿意见：</b></font>建议将论文核心贡献<u>重新聚焦为<font color="#c0392b"><b>一项</b></font></u>——<b>高采样率下自适应 FIR 控制的回归向量共线性问题及其解决方案</b>。$F(z)$ 设计、$H_\infty$ 插值等内容应压缩为背景（<i>Problem Formulation</i>）和已有方法回顾（<i>Related Work</i>），<font color="#c0392b"><b>不作为贡献声明</b></font>。

---

### <font color="#8e44ad">2.2 理论完备性 &nbsp;·&nbsp; <small>THEORETICAL RIGOR</small></font>

<font color="#c0392b"><b>缺失的分析</b></font>（按严重程度排序）：

> <font color="#c0392b"><b>❶ 严重</b></font> &nbsp; <b><u>信息矩阵条件数的定量分析</u></b>
>
> §5.2 定性地描述了"相邻元素相关度 $>0.999$"，但未给出 $\kappa(P^{-1})$ 作为系统参数 $(f_s, f_0, N, M)$ 的函数的定量表达式。
>
> 对于纯正弦扰动 $d(k) = \sin(\omega_0 k)$，回归向量 $\phi(k) = [d(k-d-1), \dots, d(k-d-N)]^T$ 的渐近 Gram 矩阵是 <b>Toeplitz 矩阵</b>，其元素为 $\Phi_{ij} = \frac{1}{2}\cos(\omega_0(i-j))$。该矩阵的<u>秩为 2</u>（由 $\sin(\omega_0 k)$ 和 $\cos(\omega_0 k)$ 两个基函数张成），$N-2$ 个特征值为<font color="#c0392b"><b>零</b></font>。在有限精度下，条件数 $\kappa(\Phi + \delta I) \sim O(1/\delta)$，随 $N$ 增大而恶化。
>
> <font color="#2980b9"><b>→</b></font> 这一分析<u>可以且应该</u>在论文中给出。它是从「工程经验」升级为「学术贡献」的<font color="#c0392b"><b>关键一步</b></font>。

> <font color="#c0392b"><b>❷ 严重</b></font> &nbsp; <b><u>抽取式 RLS 的收敛性分析</u></b>
>
> 引入抽取因子 $M$ 后，等效采样率为 $f_s/M$，$\phi_{\text{dec}}$ 的 Gram 矩阵条件数随 $M$ 增大而改善。但论文<u>未回答</u>：
> - $M$ 与收敛速率的关系？
> - $M$ 的增大是否导致参数估计的渐近方差增大（数据利用率下降为 $1/M$）？
> - 是否存在一个<font color="#2980b9"><b>最优 $M^*$</b></font> 平衡条件数与估计效率？

> <font color="#e67e22"><b>❸ 中等</b></font> &nbsp; <b><u>PE 条件在数值意义上的退化</u></b>
>
> 附录 A 指出专著中的 PE 条件在 $f_s=48$ kHz 下"数值意义上不成立"，但未给出数值 PE 的判据。建议定义 $\alpha_{\text{eff}} = \lambda_{\min}(\frac{1}{T_0}\sum \phi\phi^T)$，并分析 $\alpha_{\text{eff}}$ 与机器精度 $\epsilon_{\text{mach}}$ 的关系——当 $\alpha_{\text{eff}} < O(\epsilon_{\text{mach}})$ 时，<font color="#e67e22"><b>PE 在数值意义上失效</b></font>。

> <font color="#7f8c8d"><b>❹ 轻微</b></font> &nbsp; <b><u>$G_0F$ 展平不完全的定量刻画</u></b>
>
> §2.5 发现 $G_0F$ 在扰动频率处有 <b>−18.9 dB</b> 的残余衰减，归因于 9 个 NMP 零点的反射。建议给出 NMP 零点反射对通带增益影响的<u>上界估计</u>（如基于 Poisson-Jensen 公式），使这一观察从数值偶然上升为一般性结论。

---

### <font color="#16a085">2.3 实验验证 &nbsp;·&nbsp; <small>EXPERIMENTAL VALIDATION</small></font>

<font color="#27ae60"><b>✅ 已有内容：</b></font>

- ✅ 一条真实声学管道的 30 阶 ARMAX 辨识模型
- ✅ 固定控制器的 <b>23.1 dB</b> 抑制结果
- ✅ RLS 自适应从失败到成功的完整调试记录
- ✅ 时变扫频实验中 NLMS vs RLS 的系统对比（`ExperimentReport_TimeVarying.md`）

<br>

<font color="#c0392b"><b>❌ 缺失内容</b></font>（按重要性排序）：

> <font color="#c0392b"><b>❶ 必需</b></font> &nbsp; <b><u>与标准基线的对比</u></b>
>
> ANC 领域的事实标准是 <b>FXLMS</b>（filtered-X LMS）及其变体。论文<font color="#c0392b"><b>必须</b></font>包含在同一声学路径、同一扰动条件下与至少以下方法之一的系统对比：
> - <b>FXLMS</b>（最基础的基线）
> - <b>Filtered-U RLS</b>（RLS 在 ANC 中的标准应用形式）
> - 文献中已有的其他鲁棒 ANC 方法（如 Bai 1997 的 $H_\infty$ 方法，即题3）
>
> 对比维度应包括：<u>收敛速度</u>（多长时间达到 −20 dB）、<u>稳态抑制</u>（dB）、<u>计算量</u>（乘法/加法次数/采样周期）、<u>对频率变化的鲁棒性</u>。

> <font color="#e67e22"><b>❷ 强烈建议</b></font> &nbsp; <b><u>多声学路径验证</u></b>
>
> 当前仅测试了<u>一条</u>管道的 ARMAX 模型。建议至少包含 <b>3 条</b>不同的声学路径（不同长度/直径的管道、不同传声器-扬声器布置），以验证 $F(z)$ 设计方法和 $M$ 选择准则的泛化性。
>
> 作为<u>最低要求</u>，应通过对原始 ARMAX 系数的参数扰动（如对零极点施加 <b>±5% 的随机偏移</b>）生成一组变体模型进行消融测试。

> <font color="#2980b9"><b>❸ 建议</b></font> &nbsp; <b><u>多频扰动测试</u></b>
>
> 当前仅测试了单频（334.6 Hz）。实际 ANC 场景常涉及多个谐波分量（如发动机噪声基频 + 2 次/3 次谐波）。多频扰动下，$\phi$ 的 Gram 矩阵秩增大（$r = 2n_\omega$），共线性问题可能缓解——但 $N$ 的需求也相应增大。论文应包含至少一个<u>双频或三频</u>扰动的测试案例。

> <font color="#7f8c8d"><b>❹ 可选（终极目标）</b></font> &nbsp; <b><u>实时硬件验证</u></b>
>
> ANC 期刊（JASA、IEEE TASLP）的标准配置包含<u>实时 DSP 实现</u>。纯仿真的论文（即使使用实测辨识模型）在审稿中会被质疑：ARMAX 是<u>线性</u>的，而真实声学路径包含非线性（扬声器失真、传声器饱和、时变温度梯度）——这些因素在实时环境中可能导致仿真中<u>不存在</u>的失效模式。
>
> 如果现阶段无法完成硬件验证，建议在论文中增加一节 <b>"Limitations of Simulation-Based Validation"</b> 坦诚讨论这些因素。

---

### <font color="#d35400">2.4 表述与结构 &nbsp;·&nbsp; <small>PRESENTATION</small></font>

| <font color="#27ae60"><b>✅ 优势</b></font> | <font color="#c0392b"><b>⚠ 从学术论文视角看的问题</b></font> |
|:---|:---|
| 文档在<u>教学方法论</u>层面（"考试答题级题解"）是优秀的——推导清晰，步骤完整，代码映射精确。 | • 以<b>"题目-解答"</b>的习题格式组织，而非<br>&nbsp;&nbsp;<font color="#2980b9"><i>Motivation → Problem Formulation → Proposed Method →</i></font><br>&nbsp;&nbsp;<font color="#2980b9"><i>Theoretical Analysis → Experimental Validation → Conclusion</i></font><br>• 诊断叙述风格（"怀疑 X→排除→怀疑 Y→排除→定位 Z"）在工程报告中有效，但在论文中应重构为<u>"对已有方法在此场景中失效的系统性分析"</u><br>• <b>缺少正式的 <i>Related Work</i> 章节</b>，未将本工作置于 ANC 自适应控制文献的坐标系中 |

---

## <font color="#1a1a2e">3. 按期刊层级的可行性评估</font>

| 期刊 / 会议 | <font color="#1a1a2e">可行性</font> | 所需追加工作 | 预期篇幅 |
|:---|:---:|:---|:---:|
| <b>IEEE Trans. Automatic Control</b> | <font color="#c0392b"><b>✗ 不可行</b></font> | — | — |
| <b>Automatica</b> | <font color="#c0392b"><b>✗ 不可行</b></font> | — | — |
| <b>IEEE Trans. Audio, Speech, Language Processing</b> | <font color="#e67e22"><b>△ 困难</b></font> | 硬件实验 + FXLMS 对比 + 多路径 + 理论分析 | 10–12 页 |
| <b>Journal of Vibration and Control</b> <small>（Jafari 原论文发表期刊）</small> | <font color="#e67e22"><b>△ 有可能</b></font> | FXLMS 对比 + 多频扰动 + 条件数分析 | 10–14 页 |
| <b>Mechanical Systems and Signal Processing</b> | <font color="#e67e22"><b>△ 有可能</b></font> | 同上 + 可能需硬件验证 | 10–14 页 |
| <b><u>IEEE Signal Processing Letters</u></b> | <font color="#2980b9"><b>◈ 最可行的期刊目标</b></font> | 信息矩阵秩分析 + FXLMS 对比 + 多路径（2–3 条） | 4–5 页 |
| <b><u>InterNoise / ICSV</u></b> <small>（声学工程会议）</small> | <font color="#27ae60"><b>● 可直接投稿</b></font> | 几乎不需追加 | 4–6 页 |

<br>

> <font color="#2980b9"><b>▸ 推荐策略：</b></font>先投 <b>InterNoise</b> 或 <b>ICSV</b> 会议论文（快速获得出版物和反馈），然后在此基础上扩展为 <b>SPL</b> 短文。SPL 的核心卖点是<font color="#2980b9"><b>"一个清晰的观察 + 一个有效的解决方案"</b></font>——你的<u>共线性诊断 + 抽取式 RLS</u> 恰好符合这一格式。

---

## <font color="#1a1a2e">4. 修改建议（按优先级排列）</font>

### <font color="#c0392b">优先级 1</font> &nbsp;·&nbsp; <font color="#c0392b"><b>必须完成</b></font> <small>（否则不建议投稿任何期刊）</small>

> <font color="#c0392b"><b>①</b></font> &nbsp; <b><u>重新定义论文贡献</u></b>
>
> 将故事聚焦为一句话：
>
> <font color="#1a1a2e"><i>"我们发现在高采样率 ANC 应用中，标准 RLS 自适应 FIR 控制器因回归向量共线性而失效。我们分析了失效的数值根源（信息矩阵条件数与 $f_s/f_0$ 的关系），提出了抽取式 RLS 作为解决方案，并给出了抽取因子 $M$ 的选择准则。"</i></font>

> <font color="#c0392b"><b>②</b></font> &nbsp; <b><u>补充信息矩阵条件数的理论分析</u></b>
>
> 推导纯正弦输入下 $\phi$ 的 Gram 矩阵结构，给出条件数作为 $N$ 和 $\Delta\phi$ 的函数的表达式（至少是上界或渐近估计）。
>
> <font color="#c0392b">这是从"工程经验"升级为"学术贡献"的<u>关键一步</u>。</font>

> <font color="#c0392b"><b>③</b></font> &nbsp; <b><u>与至少一个基线方法对比</u></b>
>
> 最低要求：<b>FXLMS</b>。在同一声学路径和扰动条件下，画出收敛曲线（抑制 vs 时间）和计算量对比表。

---

### <font color="#e67e22">优先级 2</font> &nbsp;·&nbsp; <font color="#e67e22"><b>强烈建议</b></font> <small>（SPL 级别需要）</small>

> <font color="#e67e22"><b>④</b></font> &nbsp; <b><u>多声学路径测试</u></b>。至少 3 条不同的管道 ARMAX 模型（可通过参数扰动生成变体，但需在论文中说明生成方式）。

> <font color="#e67e22"><b>⑤</b></font> &nbsp; <b><u>解决 $M$ 选择准则的循环依赖</u></b>
>
> 当前的 $\Delta\phi = 2\pi f_0 M / f_s \approx 30^\circ \sim 60^\circ$ 需要已知 $f_0$。建议：
> - <b>(a)</b> 在固定控制器设计阶段，$f_0$ 已知（离线测量），$M$ 可据此设定
> - <b>(b)</b> 在自适应阶段，通过 FFT 粗估 $\hat{f}_0$ 后自适应调整 $M$
>
> 至少应在论文中<u>讨论这一限制</u>并给出实用建议。

> <font color="#e67e22"><b>⑥</b></font> &nbsp; <b><u>补充 $M$ 对收敛速率和稳态精度的定量影响分析</u></b>。通过仿真扫描 $M \in \{1, 5, 10, 15, 20, 30\}$，绘制收敛时间与稳态抑制随 $M$ 的变化曲线。

---

### <font color="#7f8c8d">优先级 3</font> &nbsp;·&nbsp; <font color="#7f8c8d"><b>锦上添花</b></font> <small>（JVC / MSSP 级别需要）</small>

> <font color="#7f8c8d"><b>⑦</b></font> &nbsp; <b>多频扰动测试</b>（至少双频）。
>
> <font color="#7f8c8d"><b>⑧</b></font> &nbsp; <b>抽取式 RLS 的收敛性证明</b>（或至少给出收敛条件的猜想并辅以数值证据）。
>
> <font color="#7f8c8d"><b>⑨</b></font> &nbsp; <b>实时硬件验证</b>。

---

## <font color="#1a1a2e">5. 最终建议</font>

<font color="#1a1a2e"><b>给作者的中肯总结：</b></font>

本工作最有价值的部分<font color="#c0392b"><b>不是</b></font> $F(z)$ 设计——那是 <b>Jafari & Ioannou</b> 的——也<font color="#c0392b"><b>不是</b></font> $H_\infty$ 插值——那是标准凸优化在 ANC 中的一次常规应用——<font color="#27ae60"><b>而是</b></font>你在调试 RLS 失败时展现的那种<u>诊断逻辑</u>：

<font color="#1a1a2e">
　<i>"从'自适应发散'的表象出发，逐一排除 $G_0F$ 近似误差、闭环反馈结构、投影参数选择<br>
　 三个替代假说，最终定位到采样率与信号带宽失配导致的数值秩亏。"</i>
</font>

这种<font color="#8e44ad"><b>"剥洋葱"式</b></font>的诊断过程，加上抽取式 RLS 的简洁补救方案，构成了一个很好的 <b>4 页 IEEE SPL 短文</b>的素材。

但要从目前的<font color="#e67e22"><b>"工程诊断报告"</b></font>升级为<font color="#2980b9"><b>"学术论文"</b></font>，你需要补上那一个缺失的环节：

<font color="#c0392b"><b>　理论。</b></font>

具体来说，就是 <b>§2.2 第 ❶ 条</b> 中提到的 <u>Gram 矩阵分析</u>。这不需要高深的数学——$\phi$ 在正弦输入下的 Gram 矩阵是一个 <b>Toeplitz 矩阵</b>，其元素是 $\frac{1}{2}\cos(\omega_0(i-j))$，特征值分布在数学上是<u>可以解析刻画的</u>。

一旦你写出了 $\kappa(\Phi)$ 作为系统参数的函数，你就不只是在报告：

> <font color="#7f8c8d">"我们发现 RLS 发散了，然后用抽取修好了。"</font>

你是在说：

> <font color="#27ae60"><b>"对于任何满足 $f_s/f_0 > C(N)$ 的系统，RLS 的信息矩阵条件数将超过 $10^{10}$，<br>
> 因此<u>必然</u>在有限精度下失效；抽取因子 $M$ 需要满足 $M > M^*(f_s, f_0, N)$<br>
> 才能将条件数降至可接受范围。"</b></font>

<font color="#c0392b"><b>这才是审稿人会认可的贡献。</b></font>



<font color="#e67e22"><b>████████████████████████████████████████████████████</b></font>

<font color="#e67e22"><b>　审稿决定　</b></font><b><u><font color="#d35400">Major Revision</font></u></b><font color="#e67e22"><b>　（大修）</b></font>

<font color="#e67e22"><b>████████████████████████████████████████████████████</b></font>

> <font color="#555555"><b>复审要求：</b>请作者在修订稿中明确标注新增的理论分析部分和补充实验，并对每一条审稿意见给出逐条回复（<i>point-by-point response</i>）。</font>

---

<p align="right"><font color="#95a5a6"><small>本审稿意见基于对题4文档及其附录、配套仿真代码、以及 ExperimentReport_TimeVarying.md 的全文阅读。</small></font></p>
