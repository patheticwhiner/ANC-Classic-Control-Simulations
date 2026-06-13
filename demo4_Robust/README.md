# demo4_Robust — 基于 ε-MOPSO 的 RST 控制器参数整定

## 概述

本目录实现了基于 **ε-MOPSO**（ε-多目标粒子群优化）算法的 RST 数字控制器自动整定框架。核心思想：将灵敏度函数整形问题转化为多目标优化，利用 Zames-Francis 积分等式构建目标函数，通过 ε-MOPSO 搜索满足频域约束的最优 X/Y 多项式零点，最终通过 Bezout 方程反解 R/S/T 控制器多项式。

统一入口脚本 [`run_eMOPSO.m`](run_eMOPSO.m) 支持 **6 种运行模式**，通过 case-switch 分派：

| case_key | 描述 | 引擎 |
|----------|------|------|
| `'RST_toy'` | RST 控制器优化（二阶教学模型） | RST_engine |
| `'RST_armax'` | RST 控制器优化（ARMAX 实测声学管道, 完整30阶） | RST_engine |
| `'F1'` | 多目标测试函数 F1（2目标, Rastrigin） | TestFunction_engine |
| `'F2'` | 多目标测试函数 F2（2目标, ZDT型, 30维） | TestFunction_engine |
| `'F3'` | 多目标测试函数 F3（3目标, 球坐标） | TestFunction_engine |
| `'F4'` | 多目标测试函数 F4（3目标, 正弦振荡） | TestFunction_engine |

## 目录结构

```
demo4_Robust/
├── run_eMOPSO.m                   ← ★ 统一入口脚本（case-switch 分派）
├── run_RST_eMOPSO.m               ← 向后兼容层（thin wrapper → run_eMOPSO('RST_toy')）
├── benchmark_MOEAs.m              ← 参数敏感性分析
│
├── utils/                         ← 算法核心与工具函数
│   ├── eMOPSO_core.m                ε-MOPSO 多目标粒子群优化
│   ├── RST_objective.m              Zames-Francis 积分目标函数
│   ├── RST_constraints.m            频域模值/延迟裕度约束
│   ├── penalty_function.m           约束违例惩罚映射
│   ├── postprocess_RST.m            反解 R/S/T 控制器
│   ├── testFunctions.m              标准多目标测试函数（F1–F4）
│   ├── evaluateFrontier.m           前沿质量定量评估（GD/IGD/Spread/HV）
│   ├── export_for_validator.m       控制器导出工具
│   └── validate_controller.m        控制器频域/时域验证
│
├── output/                        ← 运行输出
│   ├── RST_eMOPSO_results.mat       RST_toy Pareto 前沿
│   ├── RST_eMOPSO_armax_results.mat RST_armax Pareto 前沿
│   ├── test_F1_results.mat ~ F4     测试函数优化结果
│   ├── RST_model.mat                系统模型数据
│   └── RST_controller.mat           控制器多项式
│
├── plans/
│   └── unified_eMOPSO_architecture.md ← 统一入口架构设计文档
│
├── About_PSO_Foundations.md        ← 📖 PSO 完整数学推导（含定理证明）
├── About_Solution.md               ← 🔬 RST 案例推导（柔性传动系统）
├── About_RST_eMOPSO_spec.md        ← 📐 数学建模与算法规范
├── About_MOEA_algorithms.md        ← 📋 四种 MOEA 算法横向对比
├── About_ARMAX_Debugging.md        ← 🐛 ARMAX 高阶系统调试全记录
├── Review_demo4_Publishability.md  ← 📝 学术发表可行性评估
└── README.md                       ← 本文件
```

## 阅读路径

初学者建议按以下顺序阅读：

1. 📖 **[About_PSO_Foundations.md](About_PSO_Foundations.md)** — 从优化问题建模 → 标准 PSO → 收敛性定理证明 → MOPSO → ε-MOPSO 的完整数学推导
2. 📋 **[About_MOEA_algorithms.md](About_MOEA_algorithms.md)** — 四种算法（MOPSO / ε-MOPSO / NSGA-II / MODE）的横向对比与伪代码
3. 📐 **[About_RST_eMOPSO_spec.md](About_RST_eMOPSO_spec.md)** — RST 控制器整定的数学建模与 ε-MOPSO 算法规范
4. 🔬 **[About_Solution.md](About_Solution.md)** — 柔性传动系统案例的完整推导
5. 🐛 **[About_ARMAX_Debugging.md](About_ARMAX_Debugging.md)** — ARMAX 高阶系统调试全记录（Bezout 零约束、ε-支配失效分析）
6. 📝 **[Review_demo4_Publishability.md](Review_demo4_Publishability.md)** — 学术发表可行性评估
7. � 本文件 — 代码使用手册

## 快速开始

### 基础用法

```matlab
% 在 demo4_Robust/ 目录下运行:

% RST 控制器优化（二阶教学模型）
run_eMOPSO('RST_toy')

% RST 控制器优化（ARMAX 实测声学管道模型, 完整30阶）
run_eMOPSO('RST_armax')

% 标准多目标测试函数
run_eMOPSO('F1')
run_eMOPSO('F2')
run_eMOPSO('F3')
run_eMOPSO('F4')

% 向后兼容（等价于 run_eMOPSO('RST_toy')）
run_RST_eMOPSO
```

### 仿真参数

| 参数 | RST_toy | RST_armax | F1–F4 |
|------|---------|-----------|-------|
| 被控对象 | 二阶教学模型 | ARMAX(30,30,30,22), 48kHz | 标准测试函数 |
| 决策变量 | nX+nY=4 | nX+nY=**62** (nX=50, nY=12) | 2–30 |
| 目标数 | 3 | 11 | 2–3 |
| 种群大小 | 40 | **200** | 100 |
| 最大迭代 | 80 | **400** | 200 |
| ε 值 | 0.02 | 0.02 | 0.01 |

> **⚠️ RST_armax 已弃用**：经完整调试，判定 ε-MOPSO + Zames-Francis 框架对 ARMAX 声学管道模型不具有可行性。详见 [`About_ARMAX_Debugging.md §7`](About_ARMAX_Debugging.md) 和 [`About_RST_eMOPSO_spec.md §B`](About_RST_eMOPSO_spec.md)。`run_eMOPSO('RST_armax')` 入口保留但不再维护。

### 预期结果

- **RST_toy**: 闭环稳定（max pole < 1），Bezout 方程可解，5/5 通过
- **RST_armax**: 闭环稳定但 GM 可能为负（条件稳定），Bezout 残差 < 1e-12，4/5 通过
- **F1–F4**: Pareto 前沿覆盖目标空间，收敛曲线上升至 Archive 容量饱和
- **运行时间**: RST_toy < 10s，F1–F4 < 30s，RST_armax ~10min（Intel i7，62维搜索空间）

## 依赖

- MATLAB R2020b 或更高版本
- Control System Toolbox（`balred`、`ss`、`tf`，仅 `benchmark_MOEAs.m` 使用）
- System Identification Toolbox（`idpoly`，仅 ARMAX 模型处理）
- 项目 `functions/` 目录下的工具函数（`bezoutd.m`、`trimPolynomial.m` 等）
- Signal Processing Toolbox（`freqz` 用于频域计算）

## 参考

- Zames, G. and Francis, B.A., "Feedback, minimax sensitivity, and optimal robustness", IEEE TAC, 1983
- Laumanns, M. et al., "Combining Convergence and Diversity in Evolutionary Multiobjective Optimization", Evolutionary Computation, 2002
- Landau, I.D. et al., "Digital Control Systems: Design, Identification and Implementation", Springer, 2005
- Sierra, M.R. and Coello Coello, C.A., "Improving PSO-Based Multi-objective Optimization Using Crowding, Mutation and ε-Dominance", EMO 2005
