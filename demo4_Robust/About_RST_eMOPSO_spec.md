# RST 控制器整定与 ε-MOPSO 算法实现规范

本文档为基于 ε-MOPSO 的 RST 控制系统参数整定与优化算法的数学建模与程序实现指南。请严格依据以下数学定义与伪代码编写仿真程序。

---

## 一、 控制系统建模与决策变量定义 (System & Variables)

### 1.1 受控对象与固定参数
受控对象 (Plant) 离散传递函数：
$$
G(z^{-1}) = z^{-d} \frac{B(z^{-1})}{A(z^{-1})}
$$
- 必须已知 $A(z^{-1})$ 和 $B(z^{-1})$ 的系数，以及纯延迟步数 $d$。
- 系统开环不稳定极点记为 $\alpha_i$。
- 系统开环不稳定零点记为 $\beta_i = r_i e^{j\phi_i}$ （其中模长 $r_i > 1$）。

用户需预先指定的基准多项式：
- 闭环主导极点多项式 $P_D(z^{-1})$。
- RST控制器的固定部分：$H_R(z^{-1}) = 1 + z^{-1}$， $H_S(z^{-1}) = 1 - z^{-1}$。

### 1.2 决策变量与搜索空间
优化算法的搜索空间粒子（决策变量）为向量 $\boldsymbol{\theta}$，它代表待求多项式 $X$ 和 $Y$ 的实数零点：
$$
\boldsymbol{\theta} = (x_1, x_2, \dots, x_{n_X}, y_1, y_2, \dots, y_{n_Y})^T
$$
程序中对应多项式构建方式：
- $$
  X(z^{-1}, \boldsymbol{\theta}) = \prod_{i=1}^{n_X} (1 - x_i z^{-1})
  $$
- $$
  Y(z^{-1}, \boldsymbol{\theta}) = \prod_{i=1}^{n_Y} (1 - y_i z^{-1})
  $$

**待优化的期望灵敏度函数**：
$$
S_d(z^{-1}, \boldsymbol{\theta}) = \frac{A(z^{-1})}{P_D(z^{-1})} \frac{X(z^{-1}, \boldsymbol{\theta})}{Y(z^{-1}, \boldsymbol{\theta})}
$$

---

## 二、 多目标优化问题建模 (MOP Formulation)

### 2.1 目标函数 (Objective Functions)
依据 Zames-Francis 积分等式构建目标函数。每个开环不稳定零点 $\beta_i$ 对应一个目标 $f_i$，纯延迟 $d$ 对应目标 $f_d$。积分域为 $\omega \in [-\pi, \pi]$，对应 $z = e^{j\omega}$。

**目标 1：针对不稳定零点 $\beta_i = r_i e^{j\phi_i}$**
$$
f_i(\boldsymbol{\theta}) = \left| \int_{-\pi}^{\pi} \log|S_d(e^{j\omega}, \boldsymbol{\theta})| \frac{r_i^2 - 1}{1 - 2r_i \cos(\omega - \phi_i) + r_i^2} d\omega - \pi \log|B_a^{-1}(\beta_i)| \right|
$$
*(其中
$$
B_a(z) = \prod \frac{z - \alpha_k}{1 - \bar{\alpha}_k z}
$$
为系统不稳定极点构成的 Blaschke 乘积)*

**目标 2：针对纯时间延迟 (相当于无穷远零点)**
$$
f_d(\boldsymbol{\theta}) = \left| \int_{-\pi}^{\pi} \log|S_d(e^{j\omega}, \boldsymbol{\theta})| d\omega - \pi \sum_k \log|\alpha_k| \right|
$$

### 2.2 不等式约束 (Constraints)
系统实际的输出灵敏度函数为 $S_{yp}(z^{-1}, \boldsymbol{\theta})$（在本文寻优阶段，可直接令 $S_{yp} \approx S_d$ 进行频率模板验证）。对所有频率 $\omega \in [0, \pi]$ 必须满足以下绝对值约束（$\leq 0$）：

- **下界约束 (Modulus Margin):**
  $$
  g_1(e^{j\omega}, \boldsymbol{\theta}) = 1 - |1 - e^{-j\omega}|^{-1} - |S_d(e^{j\omega}, \boldsymbol{\theta})| \leq 0
  $$
- **上界约束 (Delay Margin):**
  $$
  g_2(e^{j\omega}, \boldsymbol{\theta}) = |S_d(e^{j\omega}, \boldsymbol{\theta})| - 1 - |1 - e^{-j\omega}|^{-1} \leq 0
  $$

---

## 三、 PSO 基本算法实现 (Standard PSO)

设种群大小为 $n_{pop}$，最大迭代次数为 $k_{max}$。第 $i$ 个粒子的位置为 $\mathbf{x}^i_k$（即参数 $\boldsymbol{\theta}$），速度为 $\mathbf{v}^i_k$。

**运动学更新方程：**
$$
\mathbf{x}^i_{k+1} = \mathbf{x}^i_k + \mathbf{v}^i_{k+1}
$$
$$
\mathbf{v}^i_{k+1} = w_k \mathbf{v}^i_k + c_1 r_{1,k}^i (\mathbf{pbest}^i_k - \mathbf{x}^i_k) + c_2 r_{2,k}^i (\mathbf{gbest}_k - \mathbf{x}^i_k)
$$

- $r_{1,k}^i, r_{2,k}^i \sim U(0,1)$ 为随机数。
- $c_1 = 1.5, c_2 = 1.7$。
- **自适应惯性权重 $w_k$：**
  $$
  w_{k+1} = w_{max} - (w_{max} - w_{min}) \frac{k}{k_{max}}
  $$
  (其中 $w_{max} = 0.9, w_{min} = 0.4$)

---

## 四、 约束处理机制 (Static Penalty Function)

使用外部静态惩罚技术处理 $g_1, g_2$ 约束。将违反约束的惩罚项 $\delta(\mathbf{x})$ 加到目标函数上：
$$
\varphi_m(\mathbf{x}) = f_m(\mathbf{x}) + \delta(\mathbf{x}) \quad (m \text{ 为目标维度索引})
$$

**惩罚项计算：**
$$
\delta(\mathbf{x}) = \sum_{j=1}^{q} \Lambda_j v_j^\kappa(\mathbf{x})
$$
其中：
- $q$ 为约束总数。
- $\Lambda_j$ 为大常数惩罚系数（如 $10^4$）。
- $\kappa = 2$。
- 约束惩罚函数
  $$
  v_j(\mathbf{x}) = \max(0, g_j(\mathbf{x}))
  $$
  *(如果在某些频点违反，求所有频点违反量之和或最大值)*

---

## 五、 ε-MOPSO 核心机制 (ε-Domination & Archive)

### 5.1 ε-支配定义 (ε-Domination Rule)
对于两个解（粒子） $\boldsymbol{\theta}_1$ 和 $\boldsymbol{\theta}_2$，存在 $m_f$ 个目标。我们称 **$\boldsymbol{\theta}_1$ $\varepsilon$-支配 $\boldsymbol{\theta}_2$** (记作 $\boldsymbol{\theta}_1 \prec_\varepsilon \boldsymbol{\theta}_2$)，当且仅当：
1. $$
   \forall m \in \{1, 2, \dots, m_f\}, \quad \frac{f_m(\boldsymbol{\theta}_1)}{1 + \varepsilon} \leq f_m(\boldsymbol{\theta}_2)
   $$
2. $$
   \exists m \in \{1, 2, \dots, m_f\}, \quad \frac{f_m(\boldsymbol{\theta}_1)}{1 + \varepsilon} < f_m(\boldsymbol{\theta}_2)
   $$

### 5.2 档案盒划分 (Algorithm 2: Box Indexing)
将目标空间按 $\varepsilon$ 对数划分网格。对于任意输入的目标向量 $\mathbf{f} = (f_1, f_2, \dots, f_{m_f})$，计算其所在的 Box 索引向量 $\mathcal{B}$：

```text
For all m = 1, 2, ..., m_f do:
    B_m = floor( log(f_m) / log(1 + ε) )
End For
Output B = (B_1, ..., B_{m_f})  // B 是该解的网格坐标
```

# 附录 A：多目标优化测试函数 (Test functions for multiobjective optimization)

## 1. 测试函数 F1 (Test Function F1)

$$
\begin{cases}
f_1(\boldsymbol{x}) = x_1 \\
f_2(\boldsymbol{x}) = g(\boldsymbol{x})h(\boldsymbol{x}) \\
n = 2
\end{cases} \quad (A.1)
$$

**变量约束与辅助函数：**
$$
\boldsymbol{x} = (x_1, x_2, \dots, x_n)
$$
$$
0 \leq x_1 \leq 1, \quad -30 \leq x_2 \leq 30
$$
$$
g(\boldsymbol{x}) = 11 + x_2^2 - 10\cos(2\pi x_2)
$$
$$
h(\boldsymbol{x}) = 
\begin{cases} 
1 - \sqrt{\frac{f_1(\boldsymbol{x})}{g(\boldsymbol{x})}}, & \text{if } f_1(\boldsymbol{x}) \leq g(\boldsymbol{x}) \\ 
0, & \text{otherwise} 
\end{cases}
$$

---

## 2. 测试函数 F2 (Test Function F2)

$$
\begin{cases}
f_1(\boldsymbol{x}) = x_1 \\
f_2(\boldsymbol{x}) = g(x_2, x_3, \dots, x_n)h(f_1, g) + 1 \\
n = 30
\end{cases} \quad (A.2)
$$

**变量约束与辅助函数：**
$$
0 \leq \boldsymbol{x} = (x_1, x_2, \dots, x_n) \leq 1
$$
$$
g(\boldsymbol{x}) = 1 + \frac{9}{n - 1} \sum_{i=2}^{n} x_i
$$
$$
h(f_1, g) = 1 - \sqrt{\frac{f_1}{g}} - \left(\frac{f_1}{g}\right)\sin(10\pi f_1)
$$

---

## 3. 测试函数 F3 (Test Function F3)

$$
\begin{cases}
f_1(\boldsymbol{x}) = (1 + x_3)\cos\left(\frac{\pi x_1}{2}\right)\cos\left(\frac{\pi x_2}{2}\right) \\
f_2(\boldsymbol{x}) = (1 + x_3)\cos\left(\frac{\pi x_1}{2}\right)\sin\left(\frac{\pi x_2}{2}\right) \\
f_3(\boldsymbol{x}) = (1 + x_3)\sin\left(\frac{\pi x_1}{2}\right) \\
n = 3
\end{cases} \quad (A.3)
$$

**变量约束：**
$$
0 \leq \boldsymbol{x} = (x_1, x_2, x_3) \leq 1
$$

---

## 4. 测试函数 F4 (Test Function F4)

$$
\begin{cases}
f_1(\boldsymbol{x}) = x_1 \\
f_2(\boldsymbol{x}) = x_2 \\
f_3(\boldsymbol{x}) = 3.5 - \sum_{i=1}^{n} x_i \sin(n\pi x_i) \\
n = 3
\end{cases} \quad (A.4)
$$

**变量约束：**
$$
0 \leq \boldsymbol{x} = (x_1, x_2, x_3) \leq 1
$$
*(注：根据上下文与标准多目标测试函数文献还原，对原 PDF OCR 乱码处的连乘/累加表示进行了规范化纠正整理。)*
