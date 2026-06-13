---
name: demo4-armax-failure
description: demo4 ε-MOPSO+ZF 框架对 ARMAX 声学管道模型的三重结构性失败
metadata: 
  node_type: memory
  type: project
  originSessionId: e5895dbf-ac69-4e95-bb7a-4c6497c6c25b
---

# demo4 ε-MOPSO + Zames-Francis 框架对 ARMAX 模型的失败分析

## 三个互锁的结构性问题

| # | 问题 | 层 | 根因 |
|:---|:---|:---|:---|
| 1 | 伪 NMP 零点 | 被控对象 | ARMAX(30,30,30) 在 48kHz 过参数化，9 个 |z|>1 零点中 5 个 |z|<1.01 是辨识 artifact |
| 2 | Poisson 核退化 | 目标函数 | 对 |β|→1 的零点，Poisson 核退化为 δ-脉冲，500 点积分欠采样 |
| 3 | ε-支配组合爆炸 | 优化器 | 11 维目标空间 ε-支配概率 ~2⁻¹¹，Archive 在第 10 次迭代饱和后选择压力归零 |

## 因果链

```
ARMAX 过参数化 + 48kHz
  → 近单位圆伪 NMP 零点 (问题1)
  → ZF 为每个零点设 Poisson 积分目标 (问题2a)
  → 近单位圆零点使 Poisson 核退化 (问题2b)
  → 11 目标使 ε-支配概率 → 0 (问题3)
  → 优化器在第 10 次迭代后实质停止工作
```

## 框架适用性边界

此框架仅适用于：
- 存在真实物理 NMP 零点（|β| ≥ 1.2，数量 ≤ 3）
- 模型阶数合理（nah+nbh ≤ 20）
- 采样率适中（信号带宽/Nyquist ≥ 1:5）

## 相关分析文档

- `About_ARMAX_Debugging.md §7` — 结论与决策
- `About_RST_eMOPSO_spec.md §B` — 理论详细分析

## 相关记忆

- [[research-directions]] — 此后探索的三个替代方向
- [[demo3-architecture]] — Jafari 管线（替代方案）
