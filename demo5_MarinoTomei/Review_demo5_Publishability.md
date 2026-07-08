# <font color="#1a1a2e">审稿意见</font>

## <font color="#1a1a2e">demo5_MarinoTomei —— 学术发表可行性评估</font>

> <font color="#555555"><b>审稿人视角 · 内部评审</b>　|　基于 Marino-Tomei 自适应输出调节理论的非线性 ANC 控制研究（含调节器方程 + 参数投影 + 频率自适应估计）</font>

---

## §1 总评

<font color="#e67e22"><b>████████████████████████████████████</b></font>
<font color="#e67e22"><b>　Major Revision（大修）　</b></font>
<font color="#e67e22"><b>████████████████████████████████████</b></font>

本工作实现了 Marino & Tomei 系列论文（2011/2016/2023）中自适应输出调节理论的 MATLAB 仿真复现，覆盖了从基础内模原理到非线性自适应参数估计的完整技术链条。核心贡献在于：(a) 对 Marino-Tomei 理论进行了独立的数值验证，(b) 通过两种求解器（ode45 / Euler-RK4 混合）交叉验证了仿真的数值可信度，(c) 利用调节器方程求解给出了理论最优轨迹作为控制器性能的上界参考。

| <font color="#27ae60"><b>✅ 主要贡献</b></font> | <font color="#c0392b"><b>❌ 主要不足</b></font> |
|:---|:---|
| 完整复现了 Marino-Tomei 自适应输出调节框架——从状态滤波器 ($\xi_1, \xi_2$) → 参数投影算子 → 参考控制器 $R_c$ 的控制律全链路 | ❶ <b>理论复现为主，缺乏独立创新</b>：所有核心算法（状态滤波器、投影算子、调节器方程求解）均来自 Marino-Tomei 原始论文，仿真的贡献仅限于数值验证 |
| <u>调节器方程的理论轨迹对比</u>（实际轨迹 vs $\Pi w$ 理论最优）是一个有价值的分析视角——它把自适应控制器的性能与频域设计的最优解直接对标 | ❷ <b>无 ANC 声学场景验证</b>：所有仿真均在抽象的二阶非线性系统上运行，未使用项目 `dataset/` 中的真实声学传递函数或次级通路模型。审稿人会追问："这个控制器在真实 ANC 管道上能工作吗？" |
| 双求解器交叉验证（ode45 vs Euler/RK4）和扰动状态空间建模的正确性自检体现了良好的数值实验规范 | ❸ <b>参数投影算子的实现不一致</b>：`demos/demo_simple_adaptive.m` 使用简化版投影（硬边界 `[low, high]`），而主仿真使用完整的平滑投影算子。两种实现的等价性未在文档中证明，可能导致读者困惑 |
| 三区间分段扰动场景设计（含频率切换瞬态）针对性地测试了自适应机制的鲁棒性 | ❹ <b>2023 实验与主线的理论关联不清晰</b>：`MarinoTomei_2023_unstable.m` 使用了完全不同的系统模型（不稳定传函）和控制策略（极点配置风格），与 Marino-Tomei 输出调节框架的理论联系未建立 |

---

## §2 逐项评估

### §2.1 新颖性（Novelty）

| 技术要素 | 来源 | <font color="#1a1a2e">审稿意见</font> |
|:---|:---|:---|
| 自适应输出调节控制律 $u = -ke - \hat{\chi}_1 - \mu_1\hat{\theta}_1 - \mu_2\hat{\theta}_2$ | <font color="#7f8c8d">Marino & Tomei (2011)</font> | <font color="#c0392b"><b>✗ 非原创</b></font>。控制律结构与原始论文完全一致。 |
| 参数投影算子 (平滑版 + 简化硬边界版) | <font color="#7f8c8d">Marino & Tomei (2011), Pomet & Praly (1992)</font> | <font color="#c0392b"><b>✗ 非原创</b></font>。投影算子是自适应控制的标准工具。 |
| 调节器方程 $\Pi R = A\Pi + \frac{1}{\beta}b\Gamma + P$ 的数值求解 | <font color="#7f8c8d">Francis & Wonham (1976), Isidori (1995)</font> | <font color="#e67e22"><b>△ 增量贡献</b></font>。将调节器方程的解作为自适应控制器性能的理论上界是一个好的分析框架，但这个想法在输出调节文献中并非全新。 |
| 双求解器交叉验证 + 扰动状态空间正确性自检 | <font color="#2980b9"><b>作者独立完成</b></font> | <font color="#27ae60"><b>◈ 有方法论价值</b></font>。这不是算法贡献，但作为<u>数值实验规范</u>——在复杂非线性自适应仿真中使用两种独立求解器互相验证——是良好的工程实践，适合作为 SIAM 或 CSC 会议短文的方法论说明。 |
| 分段扰动场景（含频率切换）下的瞬态分析 | <font color="#7f8c8d">Marino & Tomei (2011) 原文有此建议</font> | <font color="#e67e22"><b>△ 增量贡献</b></font>。原文建议了分段场景，但具体的三区间设计（双频→单频→零参考）和切换点附近的瞬态放大分析是作者的独立工作。 |

<font color="#c0392b"><b>▶ 审稿意见：</b></font>本工作的新颖性不足——所有核心算法均来自已有文献。可作为教育性文章（如 IEEE Control Systems Magazine 的教育专栏）发表，强调"如何从零实现自适应输出调节控制器"的教学价值，而非声称方法论创新。或者，如果作者能补充在真实 ANC 声学系统上的实验结果（使用项目 `dataset/` 中的实测传递函数），则可以提升为应用验证型的短文。

### §2.2 理论完备性（Theoretical Rigor）

> <font color="#e67e22"><b>❶ 中等</b></font> &nbsp; <b><u>投影算子两种实现的理论等价性缺失</u></b>
>
> `About_MarinoTomei.md` 中详细推导了完整的平滑投影算子（含 $p_r(\hat{\theta})$ 和 $\text{grad}\,p_r$），但 `demos/demo_simple_adaptive.m` 使用了简化版：
>
> ```matlab
> if theta_hat <= low && update < 0, theta_dot = 0;
> elseif theta_hat >= high && update > 0, theta_dot = 0;
> else, theta_dot = update; end
> ```
>
> 这是硬边界投影，而非平滑投影。两种实现的关系（平滑投影在 $\epsilon_r \to 0$ 时趋近硬边界投影）未在文档中说明。在 $g=100$ 的大增益下，硬边界投影可能导致参数估计在边界处的 chattering。
>
> <font color="#2980b9"><b>→</b></font> 在 `About_MarinoTomei.md` 中补充一节讨论两种投影算子的关系、适用场景和 $\epsilon_r$ 的选择准则。

> <font color="#e67e22"><b>❷ 中等</b></font> &nbsp; <b><u>Lyapunov 稳定性分析的仿真验证缺失</u></b>
>
> `About_MarinoTomei.md` 给出了 Lyapunov 函数 $V = \frac{1}{2}e^2 + \frac{1}{2g}\sum(\hat{\theta}_i - \theta_i^*)^2$ 及其导数的理论分析，但仿真中未计算和绘制 $V(t)$ 的时间演化。在自适应控制论文中，展示 $V(t)$ 单调递减（或至少在瞬态后递减）是验证稳定性分析的标准做法。
>
> <font color="#2980b9"><b>→</b></font> 在主仿真的输出图中增加 Lyapunov 函数 $V(t)$ 及其导数 $\dot{V}(t)$ 的 subplot。

> <font color="#27ae60"><b>❸ 轻微</b></font> &nbsp; <b><u>参考控制器 $R_c$ 矩阵的构造仅适用于已知频率</u></b>
>
> 当前 $R_c$ 矩阵的构造显式依赖 $\theta_1^* = \omega_1^2 + \omega_2^2$ 和 $\theta_2^* = \omega_1^2 \cdot \omega_2^2$——这些是真实参数值。在 `MarinoTomei_ode45.m` 中，$R_c$ 在 `system_dynamics` 内使用 `theta1_true` 和 `theta2_true` 构造。这引入了"已知真实参数"的假设，与自适应控制器"参数未知"的前提相矛盾。实际实现应使用参数估计值 $\hat{\theta}_1, \hat{\theta}_2$ 来构造 $R_c$，并讨论由此引入的时变性和稳定性影响。
>
> <font color="#2980b9"><b>→</b></font> 讨论是否应该用 $\hat{\theta}$ 替代 $\theta^*$ 构造 $R_c$，并分析这对闭环稳定性的影响（这本身就是 Marino-Tomei 理论的一个开放问题）。

---

### §2.3 实验验证（Experimental Validation）

| 验证维度 | 现状 | 评分 |
|:---|:---|:---|
| 数值求解器交叉验证 | ode45 vs Euler/RK4，两种独立实现产生一致结果 | <font color="#27ae60"><b>★★★★★</b></font> |
| 扰动模型正确性自检 | `demo_ss_disturbance.m` 对比状态空间模型输出与解析表达式 | <font color="#27ae60"><b>★★★★★</b></font> |
| 参数收敛性验证 | 展示了 $\hat{\theta}_1, \hat{\theta}_2$ 的时间曲线及与真实值的对比 | <font color="#27ae60"><b>★★★★☆</b></font> |
| 调节器方程残差检验 | 求解后打印了残差范数 | <font color="#27ae60"><b>★★★★☆</b></font> |
| 与基线方法的对比 | <font color="#c0392b"><b>无</b></font> — 未与任何其他 ANC 方法（FxLMS, RST, H∞, LQG）对比 | <font color="#c0392b"><b>☆☆☆☆☆</b></font> |
| 真实声学数据验证 | <font color="#c0392b"><b>无</b></font> — 仅使用抽象二阶系统 | <font color="#c0392b"><b>☆☆☆☆☆</b></font> |
| 参数灵敏度分析 | 未系统扫描 $k, g, D$ 参数对性能的影响 | <font color="#e67e22"><b>★★☆☆☆</b></font> |
| Monte Carlo 统计 | 未进行多次随机初始化运行以评估成功率/方差 | <font color="#e67e22"><b>★★☆☆☆</b></font> |

<font color="#c0392b"><b>▶ 审稿意见：</b></font>数值实验的**内部一致性验证**（求解器交叉验证、扰动模型自检）做得好，但**外部有效性验证**（基线对比、真实数据）完全缺失。最低要求：(a) 至少与一种简单方法（如高增益 PID 或 FxLMS 的简化版）对比跟踪性能和控制代价；(b) 将系统参数替换为项目 `dataset/` 中的真实声学传递函数进行至少一组仿真。

---

### §2.4 表述结构（Presentation）

| 文档 | 质量评估 |
|:---|:---|
| `README.md` | <font color="#27ae60"><b>★★★★☆</b></font> 结构清晰，包含目录树、脚本说明表、阅读路径和快速开始。可改进：增加仿真输出的代表性截图 |
| `About_MarinoTomei.md` | <font color="#27ae60"><b>★★★★☆</b></font> 完整的理论推导，从系统模型到 Lyapunov 稳定性分析。可改进：与代码的交叉引用（例如标注"此公式对应 MarinoTomei_ode45.m 第 X 行"） |
| 代码注释 | <font color="#e67e22"><b>★★★☆☆</b></font> 中文注释覆盖了关键变量，但缺少：(a) 每个函数的输入/输出维度文档，(b) 关键公式与代码行的对应关系，(c) 已知 bug/限制的标注 |
| 图注与坐标轴标签 | <font color="#e67e22"><b>★★★☆☆</b></font> 大部分图有标题和标签，但存在硬编码中文。审稿人建议统一为英文标签或使用 LaTeX 解释器 |

---

## §3 发表建议

### 推荐路径（按可行性排序）

| 优先级 | 目标期刊/会议 | 所需补充工作 | 预期篇幅 |
|:---|:---|:---|:---|
| <font color="#27ae60"><b>A</b></font> | IEEE Control Systems Magazine (Education Column) | 补充教学框架（从内模原理 → 简化控制器 → 完整自适应控制器的递进教学设计），补全 Lyapunov 函数可视化 | 6–8 页 |
| <font color="#27ae60"><b>A</b></font> | CSC / CCC (中国控制会议) | 补充与至少一种基线方法的定量对比，参数灵敏度分析 | 4–6 页 |
| <font color="#e67e22"><b>B</b></font> | IEEE/ASME TMECH (短文) | 必须在真实 ANC 声学管道上进行实验验证，使用 `dataset/` 中的实测传递函数 | 6–8 页 |
| <font color="#c0392b"><b>C</b></font> | Automatica (短文) | 需要理论创新（如用 $\hat{\theta}$ 替代 $\theta^*$ 构造 $R_c$ 的时变稳定性分析），当前纯复现工作不适合 | — |

### 不推荐的投稿方向

- **纯算法期刊**（IEEE TAC, Automatica, SICON）：当前工作无独立的理论贡献
- **纯信号处理期刊**（IEEE TSP, Signal Processing）：未涉及声学信号处理的具体挑战（延迟、非线性失真、时变路径）

---

## §4 修订优先级排序

<font color="#c0392b"><b>█ 阻塞性</b></font> (必须在下一版解决)

1. 补充至少一种基线方法的定量对比
2. 在 `About_MarinoTomei.md` 中说明两种投影算子实现的关系

<font color="#e67e22"><b>█ 重要</b></font> (应在终稿前解决)

3. 在主仿真中增加 Lyapunov 函数 $V(t)$ 的绘制
4. 讨论 $R_c$ 矩阵应使用 $\hat{\theta}$ 还是 $\theta^*$ 的理论影响
5. 系统扫描 $k, g$ 参数的灵敏度

<font color="#2980b9"><b>█ 增强性</b></font> (如有余力)

6. 在真实 ANC 声学数据上运行至少一组仿真
7. Monte Carlo 多初值统计分析
8. 代码中增加公式-代码交叉引用注释
