# 多目标优化算法 (MOEAs) 数学与算法实现手册
包含：多目标问题定义、$\varepsilon$-MOPSO、MOPSO、NSGA-II、MODE 的核心数学推导与伪代码逻辑。

---

## 一、 多目标优化问题 (MOP) 标准数学表达

在编写算法框架前，须统一定义多目标优化问题（Multi-Objective Problem, MOP）：

$$
\begin{cases}
\underset{\mathbf{x}}{\text{Minimize}} \quad F(\mathbf{x}) = \big( f_1(\mathbf{x}), f_2(\mathbf{x}), \dots, f_M(\mathbf{x}) \big)^T \\
\text{subject to:} \\
g_j(\mathbf{x}) \leq 0, \quad j = 1, \dots, J \\
h_q(\mathbf{x}) = 0, \quad q = 1, \dots, Q \\
\mathbf{x} \in \mathcal{S} = [\mathbf{x}_{\min}, \mathbf{x}_{\max}]^D
\end{cases}
$$

- $\mathbf{x} = (x_1, \dots, x_D)$ 为 $D$ 维决策变量。
- $M$ 为目标数量。
- **帕累托支配定义 (Pareto Dominance $\prec$)**：解 $\mathbf{x}_1$ 支配 $\mathbf{x}_2$ ($\mathbf{x}_1 \prec \mathbf{x}_2$)，当且仅当：
  $$
  \forall m \in \{1,\dots,M\}, f_m(\mathbf{x}_1) \leq f_m(\mathbf{x}_2)
  $$
  且
  $$
  \exists k \in \{1,\dots,M\}, f_k(\mathbf{x}_1) < f_k(\mathbf{x}_2)
  $$

---

## 二、 四种算法核心机制横向对比

| 算法名称 | 算法类别 | 核心寻优机制 | 档案/非支配解维护机制 | 多样性维护策略 |
| :--- | :--- | :--- | :--- | :--- |
| **MOPSO** | 群智能 (Swarm) | 速度与位置公式更新 | 严格帕累托支配，维护外部存档 (Archive) | 自适应网格 (Adaptive Grid) / 拥挤度距离 |
| **$\varepsilon$-MOPSO** | 群智能 (Swarm) | 速度与位置公式更新 | $\varepsilon$-支配，维护外部存档 | $\varepsilon$-Box 网格划分 / Sigma 方法 |
| **NSGA-II** | 进化算法 (EA) | 交叉 (Crossover) 与变异 (Mutation) | 快速非支配排序 (Fast Non-dominated Sorting) | 拥挤度距离 (Crowding Distance) |
| **MODE** | 进化算法 (EA) | 差分变异与二项交叉 | 竞争选择 / 帕累托非支配排序 | 缩放因子 $F$ / 交叉概率 $CR$ |

---

## 三、 算法详细数学推导与编码公式

### 3.1 基于粒子群的算法 (MOPSO & $\varepsilon$-MOPSO)

两者共用底层的粒子运动学方程，但在**全局最优领导者(gbest)选择**和**外部存档(Archive)维护**上存在差异。

#### 3.1.1 基础粒子运动学方程 (Particle Kinematics)
对于第 $i$ 个粒子在迭代 $k$ 时的速度 $\mathbf{v}_k^i$ 和位置 $\mathbf{x}_k^i$：
$$
\mathbf{v}_{k+1}^i = w \mathbf{v}_k^i + c_1 r_1 \big(\mathbf{pbest}_k^i - \mathbf{x}_k^i\big) + c_2 r_2 \big(\mathbf{gbest}_k^i - \mathbf{x}_k^i\big)
$$
$$
\mathbf{x}_{k+1}^i = \mathbf{x}_k^i + \mathbf{v}_{k+1}^i
$$
- $w = w_{\max} - (w_{\max} - w_{\min})\frac{k}{k_{\max}}$ 为自适应惯性权重。
- $r_1, r_2 \sim U(0,1)$ 为随机矩阵。
- $\mathbf{pbest}_k^i$：个体历史最优（若新位置支配原 $\mathbf{pbest}$ 则更新；若互不支配则随机保留其一）。

#### 3.1.2 MOPSO (经典多目标粒子群)
- **Archive 更新**：使用严格帕累托支配逻辑。如果新粒子不被 Archive 中任何解支配，则加入 Archive，并剔除 Archive 中被新粒子支配的解。
- **$\mathbf{gbest}$ 选择 (网格轮盘赌)**：将目标空间划分为超立方体网格。统计每个网格内包含的粒子数 $N_j$。网格被选中的概率 $P_j \propto 1 / N_j$（或类似反比机制，为了提升多样性）。选中网格后，在网格内随机挑选一个粒子作为 $\mathbf{gbest}$。
- **截断机制**：若 Archive 超出最大容量，优先删除最拥挤网格内的粒子。

#### 3.1.3 $\varepsilon$-MOPSO (本文改进版算法)
- **$\varepsilon$-支配判断**：解 $\mathbf{x}_1$ $\varepsilon$-支配 $\mathbf{x}_2$ ($\mathbf{x}_1 \prec_\varepsilon \mathbf{x}_2$) 需满足：
  $$
  \frac{f_m(\mathbf{x}_1)}{1 + \varepsilon} \leq f_m(\mathbf{x}_2), \forall m \quad \text{且} \quad \exists k, \frac{f_k(\mathbf{x}_1)}{1 + \varepsilon} < f_k(\mathbf{x}_2)
  $$
- **Box 网格索引计算**：对于解的目标向量 $\mathbf{f}$，其所在的超立方体索引向量 $\mathcal{B} = (B_1, \dots, B_M)$ 计算如下：
  $$
  B_m = \left\lfloor \frac{\log(f_m)}{\log(1 + \varepsilon)} \right\rfloor
  $$
- **Archive 更新机制**：
  1. 若新解的 Box 索引与 Archive 中某解的 Box 索引完全相同，则保留目标值更小的解（或直接替换）。
  2. 若新解的 Box $\varepsilon$-支配现有 Archive 中某些解的 Box，新解入列，旧解删除。
  3. 天然限制容量上限，不需要额外的截断函数计算。

---

### 3.2 进化类算法 (NSGA-II & MODE)

进化类算法不使用速度概念，直接对候选解集合(Population)进行数学操作。

#### 3.2.1 NSGA-II (非支配排序遗传算法 II)
LLM 编码需实现三大核心模块：

**模块 A：快速非支配排序 (Fast Non-dominated Sorting)**
1. 对种群中每个解 $p$，计算两个参数：
   - $n_p$：支配解 $p$ 的其他解的数量。
   - $S_p$：解 $p$ 所支配的解的集合。
2. $n_p = 0$ 的所有解归为第一前沿面 $Front_1$ (Rank=1)。
3. 对于 $Front_1$ 中的每个解 $p$，遍历其 $S_p$ 中的解 $q$，执行 $n_q = n_q - 1$。如果 $n_q = 0$，则将 $q$ 加入 $Front_2$。依此类推划分所有 Rank。

**模块 B：拥挤度距离 (Crowding Distance, CD)**
为了保持多样性，对同一前沿面 $Front_k$ 内的解进行距离排序：
$$
CD_i = \sum_{m=1}^{M} \frac{f_m^{(i+1)} - f_m^{(i-1)}}{f_m^{\max} - f_m^{\min}}
$$
- 边界解（极值解）的距离设为无穷大 $CD = \infty$。
- **选择机制 (拥挤比较算子 $\prec_n$)**：$i \prec_n j$ 当且仅当 $Rank_i < Rank_j$，或者 ($Rank_i = Rank_j$ 且 $CD_i > CD_j$)。

**模块 C：SBX交叉与多项式变异 (数值优化标准算子)**
- 模拟二进制交叉 (SBX):
  $$
  x_1^{child} = 0.5 \big[ (1+\beta)x_1^{parent} + (1-\beta)x_2^{parent} \big]
  $$
  $$
  x_2^{child} = 0.5 \big[ (1-\beta)x_1^{parent} + (1+\beta)x_2^{parent} \big]
  $$
  $\beta$ 是根据特定概率分布由随机数 $u \in (0,1)$ 计算出的散布因子。
- 多项式变异:
  $$
  x^{new} = x + \delta_{q} (x_{\max} - x_{\min})
  $$

#### 3.2.2 MODE (多目标差分进化算法)
结合了差分进化（DE/rand/1/bin）与多目标选择机制。

**模块 A：差分变异 (Differential Mutation)**
对当前种群中的每个目标解 $\mathbf{x}_{i,G}$，随机选择三个不同的个体 $\mathbf{x}_{r1}, \mathbf{x}_{r2}, \mathbf{x}_{r3}$，生成变异向量 $\mathbf{v}_i$：
$$
\mathbf{v}_i = \mathbf{x}_{r1} + F \cdot (\mathbf{x}_{r2} - \mathbf{x}_{r3})
$$
- $F$：缩放因子，论文实验设为 $0.5$。

**模块 B：二项式交叉 (Binomial Crossover)**
将目标解 $\mathbf{x}_i$ 与变异向量 $\mathbf{v}_i$ 混合，生成试验向量 $\mathbf{u}_i$（针对第 $j$ 个维度）：
$$ 
u_{i,j} = 
\begin{cases} 
v_{i,j}, & \text{if } rand_j(0,1) \leq CR \text{ or } j = j_{rand} \\ 
x_{i,j}, & \text{otherwise} 
\end{cases} 
$$
- $CR$：交叉概率，论文实验设为 $0.2$。
- $j_{rand}$：随机选择的一个维度，保证至少有一个基因来自变异向量。

**模块 C：多目标选择 (Pareto-based Selection)**
计算试验向量 $\mathbf{u}_i$ 的多目标适应度：
1. 若 $\mathbf{u}_i \prec \mathbf{x}_i$ ($\mathbf{u}_i$ 严格支配 $\mathbf{x}_i$)：直接用 $\mathbf{u}_i$ 替换 $\mathbf{x}_i$。
2. 若 $\mathbf{x}_i \prec \mathbf{u}_i$：保留原解 $\mathbf{x}_i$。
3. 若互不支配：将 $\mathbf{u}_i$ 丢入外部 Archive（或并入混合种群）；在迭代末尾执行类似 NSGA-II 的非支配排序或基于拥挤度的截断机制，以维持 $N_{pop}$ 的种群规模。
