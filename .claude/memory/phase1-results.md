---
name: phase1-results
description: "Phase 1 被控对象性能上限分析结果 — ARMAX(30,30,30,22) 三层约束结构与变频率分析"
metadata: 
  node_type: memory
  type: project
  originSessionId: e5895dbf-ac69-4e95-bb7a-4c6497c6c25b
---

# Phase 1 结果：ARMAX 被控对象性能上限分析

## 三层约束

| 约束层 | 结果 | 含义 |
|:---|:---|:---|
| Bode 灵敏度积分 | = 0（无开环不稳定极点） | 不限制单频抑制深度，只限制抑制带宽 |
| NMP 零点 Poisson 约束 | 9 个零点但约束强度 0.005-0.212 | 近乎空洞 |
| 延迟 (458μs) | 反馈带宽上限 ~727 Hz | 仅约束纯反馈架构 |

## 核心纠正

- Bode 积分不限制单频抑制深度（零测度集合可任意抑制）→ **内模原理可突破**
- 前馈 ANC 不受 Bode 积分约束 → 受限于 S(z) 可逆性和因果性
- Youla Q 的局限性 = 参数化容量 (nQ 阶数)，非 Bode 积分

## 文档位置

- 正式报告: `dataset/About_Plant_Performance.md` (datasheet 格式)
- 分析脚本: `dataset/analyze_plant_performance.m`
- 全局计划: `plans/Global_Work_Plan.md`
