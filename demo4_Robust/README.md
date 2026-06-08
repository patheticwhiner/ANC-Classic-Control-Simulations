# demo4_Robust — 基于 ε-MOPSO 的 RST 控制器参数整定

## 概述

本目录实现了基于 **ε-MOPSO**（ε-多目标粒子群优化）算法的 RST 数字控制器自动整定框架。核心思想：将灵敏度函数整形问题转化为多目标优化，利用 Zames-Francis 积分等式构建目标函数，通过 ε-MOPSO 搜索满足频域约束（模值裕度、延迟裕度）的最优 X/Y 多项式零点，最终通过 Bezout 方程反解 R/S/T 控制器多项式。

## 目录结构

```
demo4_Robust/
├── run_RST_eMOPSO.m              ← 主入口脚本
├── benchmark_MOEAs.m             ← 参数敏感性分析
│
├── utils/                        ← 算法核心与工具函数
│   ├── eMOPSO_core.m               ε-MOPSO 多目标粒子群优化
│   ├── RST_objective.m             Zames-Francis 积分目标函数
│   ├── RST_constraints.m           频域模值/延迟裕度约束
│   ├── penalty_function.m          约束违例惩罚映射
│   ├── postprocess_RST.m           反解 R/S/T 控制器
│   ├── testFunctions.m             标准多目标测试函数
│   ├── export_for_validator.m      控制器导出工具
│   └── validate_controller.m       控制器频域/时域验证
│
├── output/                       ← 运行输出
│   ├── RST_eMOPSO_results.mat      Pareto 前沿 / Archive
│   ├── RST_model.mat               系统模型数据
│   ├── RST_controller.mat          控制器多项式
│   └── *.txt                       运行日志
│
├── About_PSO_Foundations.md       ← 📖 PSO 完整数学推导（含定理证明）
├── About_Solution.md              ← 🔬 RST 案例推导（柔性传动系统）
├── About_RST_eMOPSO_spec.md       ← 📐 数学建模与算法规范
├── About_MOEA_algorithms.md       ← 📋 四种 MOEA 算法横向对比
└── README.md                      ← 本文件
```

## 阅读路径

初学者建议按以下顺序阅读：

1. 📖 **[About_PSO_Foundations.md](About_PSO_Foundations.md)** — 从优化问题建模 → 标准 PSO → 收敛性定理证明 → MOPSO → ε-MOPSO 的完整数学推导
2. 📋 **[About_MOEA_algorithms.md](About_MOEA_algorithms.md)** — 四种算法（MOPSO / ε-MOPSO / NSGA-II / MODE）的横向对比与伪代码
3. 📐 **[About_RST_eMOPSO_spec.md](About_RST_eMOPSO_spec.md)** — RST 控制器整定的数学建模与 ε-MOPSO 算法规范
4. 🔬 **[About_Solution.md](About_Solution.md)** — 柔性传动系统案例的完整推导（物理模型 → 优化命题 → 控制器反解）
5. 💻 本文件 — 代码使用手册

## 快速开始

1. 在 MATLAB 中将工作目录设为 `demo4_Robust/`
2. 将 `functions/` 和 `utils/` 添加到路径：`addpath('../functions', 'utils')`
3. 运行主仿真：`run('run_RST_eMOPSO')`
4. 查看结果：`output/RST_eMOPSO_results.mat` 包含最优解

### 仿真参数

| 参数 | 值 | 说明 |
|------|-----|------|
| 被控对象 | `G(z) = z⁻¹·(0.2+0.15z⁻¹)/(1−1.2z⁻¹+0.45z⁻²)` | 离散二阶系统，极点 0.6±0.3i |
| 主导极点 | ζ=0.8, ωₙ=0.25π | 期望闭环动态 |
| H_R | `1+z⁻¹` | 控制器分子固定因子 |
| H_S | `1−z⁻¹` | 控制器分母固定因子（积分作用） |
| 种群大小 | 40 | ε-MOPSO 粒子数 |
| 最大迭代 | 80 | 优化终止代数 |
| ε 值 | 0.02 | ε-支配精度 |
| 目标数 | 3 | 零点匹配 + 延迟积分 + Bezout 残差 |

### 预期结果

- **闭环稳定**：最大极点模值 < 0.95
- **控制器可实现**：R、S 多项式根均在单位圆内
- **Bezout 残差**：≈ 0（闭环极点精确匹配设计值）
- **运行时间**：< 20 秒（Intel i7 / 40 粒子 × 80 迭代）

## 依赖

- MATLAB R2020b 或更高版本
- 项目 `functions/` 目录下的工具函数（`bezoutd.m`、`trimPolynomial.m` 等）
- Signal Processing Toolbox（`freqz` 用于频域计算）

## 参考

- Zames, G. and Francis, B.A., "Feedback, minimax sensitivity, and optimal robustness", IEEE TAC, 1983
- Laumanns, M. et al., "Combining Convergence and Diversity in Evolutionary Multiobjective Optimization", Evolutionary Computation, 2002
- Landau, I.D. et al., "Digital Control Systems: Design, Identification and Implementation", Springer, 2005
