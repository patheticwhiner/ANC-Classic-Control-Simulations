---
name: demo3-architecture
description: Jafari-Ioannou 固定+自适应 ANC 管线架构（demo3_Robust 的实现细节）
metadata: 
  node_type: memory
  type: project
  originSessionId: e5895dbf-ac69-4e95-bb7a-4c6497c6c25b
---

# demo3 Jafari-Ioannou 固定+自适应 ANC 管线

## 管线架构

```
ARMAX(30,30,30,22) @ 48kHz
  ↓
① F(z) 频谱展平（全自动）
  - NMP 零点提取 → 单位圆内映射 (z→1/z̄)
  - B̃(z) = 镜像重建的最小相位分子
  - F(z) = A(z⁻¹)/B̃(z⁻¹) · L(z)  (L为低通)
  - 输出: G_eff = z^{-d}·B·F/A → 最小相位 + 平坦幅频
  ↓
② 固定 FIR K(z)（半手动）
  - 自动: findpeaks → 检测共振峰
  - 手动: 选第一共振峰 f₀, 设 FIR 长度 N=120
  - 约束: K(e^{jω₀}) = 1/G_eff(e^{jω₀})  (单频插值)
  - 输出: θ_fixed (120个 FIR 系数)
  ↓
③ 自适应 NLMS（全自动）
  - μ = 0.05, N_adapt = 64
  - 在线追踪时变干扰频谱
```

## 实测性能

| 测试场景 | 固定 FIR | 自适应 NLMS |
|:---|:---|:---|
| 固定频率 334.6Hz | 23 dB | 22 dB |
| 线性扫频 300→370Hz | 退化 | 维持 |
| 非线性调频 | 失效 | 维持 |

## 与 Youla (demo1) 的关键区别

| | Youla (demo1) | Jafari (demo3) |
|:---|:---|:---|
| 固定部分 | R₀/S₀（极点配置）| F(z) + K(z) FIR |
| 自适应部分 | Q 滤波器 (RLS) | NLMS 直接自适应 |
| 闭环极点 | 固定 | 移动（F 改变对象 + NLMS 在线调） |
| 稳定性保证 | Youla 引理（任何稳定 Q 都稳定）| 需要充分激励 + μ 选择 |
| 手动步骤 | P_D, H_R, H_S | f₀, N_fir |

## 代码入口

- `demo3_Robust/JafariANC_RealAcoustic.m` — 统一仿真脚本
- 固定 FIR 设计: lines 49–69
- 自适应 NLMS: lines 78–81

## 相关记忆

- [[research-directions]] — 方向 C：完善 Jafari 管线
- [[demo4-armax-failure]] — 被放弃的替代方案
