# About_MarinoTomei — 自适应输出调节：从系统模型到 Lyapunov 稳定性

> **本文档是 `demo5_MarinoTomei` 下所有仿真的公共理论地基。** 各仿真脚本直接引用本文档中的公式编号，不再重复推导。

## 适用仿真

| 仿真脚本 | 引用的理论章节 |
|:---|:---|
| `MarinoTomei_adaptive_regulation.m (integrator='ode45')` | §系统模型, §状态滤波器, §参数估计器, §控制律, §投影算子 |
| `MarinoTomei_adaptive_regulation.m (integrator='euler')` | 同上（Euler 定步长实现） |
| `MarinoTomei_regulator_benchmark.m (integrator='rk4')` | 同上 + §调节器方程 |
| `MarinoTomei_regulator_benchmark.m (integrator='ode45')` | 同上（ode45 实现） |
| `demos/demo_simple_adaptive.m` | §系统模型, §控制律（简化 D 矩阵 + 硬边界投影） |
| `MarinoTomei_2023_unstable.m` | 独立理论框架（不稳定传函 + 极点配置自适应） |
| `MarinoTomei_2016_adaptive_freq_est.m` | 独立理论框架（频率自适应估计），完整推导见 [docs/MarinoTomei_2016_Derivation.md](docs/MarinoTomei_2016_Derivation.md) |

## 系统模型

考虑以下二阶非线性系统：

$$
\begin{align}
\dot{x}_1 &= x_2 + a_1 x_1 + \frac{1}{\beta}u + w_1(t) \\
\dot{x}_2 &= \frac{1}{\beta}b_2 u \\
y &= x_1
\end{align}
$$

其中：
- $x_1$, $x_2$ 是系统状态
- $a_1$, $b_2$, $\beta$ 是未知系统参数
- $u$ 是控制输入
- $w_1(t)$ 是外部扰动
- $y$ 是系统输出

扰动 $w_1(t)$ 可以表示为：

$$
w_1(t) = 
\begin{cases}
\sin(\omega_1 t) + 0.2\sin(\omega_3 t), & t < 50 \\
\sin(\omega_1 t), & t \geq 50
\end{cases}
$$

## 控制目标

控制目标是使系统输出 $y$ 跟踪参考信号 $y_r = -w_2(t)$，其中：

$$
w_2(t) = 
\begin{cases}
\sin(\omega_2 t), & t < 100 \\
0, & t \geq 100
\end{cases}
$$

定义输出调节误差：

$$
e = y - y_r = x_1 - (-w_2)
$$

## 状态滤波器设计

为了处理未知参数，设计两个状态滤波器 $\xi_1$ 和 $\xi_2$：

$$
\begin{align}
\dot{\xi}_1 &= D\xi_1 + [0, u, 0, 0]^T \\
\dot{\xi}_2 &= D\xi_2 + [0, 0, 0, u]^T
\end{align}
$$

其中 $D$ 是一个赫尔维茨矩阵：

$$
D = \begin{bmatrix}
-d_2 & 1 & 0 & 0 \\
-d_2 & 0 & 1 & 0 \\
-d_3 & 0 & 0 & 1 \\
-d_4 & 0 & 0 & 0
\end{bmatrix}
$$

定义 $\mu_1 = \xi_{1,1}$ 和 $\mu_2 = \xi_{2,1}$ 作为控制器的输入。

## 参数估计器设计

为估计未知参数对应的量，设计以下估计器：

$$
\dot{\hat{\chi}} = D\hat{\chi} - [d_2, d_3, d_4, d_5]^T u
$$

参数更新律（采用投影算法）：

$$
\dot{\hat{\theta}}_i = g \cdot \text{Proj}(\mu_i \cdot e, \hat{\theta}_i), \quad i = 1,2
$$

其中投影算子 $\text{Proj}(\varphi, \hat{\theta})$ 定义为：

$$
\text{Proj}(\varphi, \hat{\theta}) = 
\begin{cases} 
\varphi, & \text{若 } p_r(\hat{\theta}) \leq 0 \\
\varphi, & \text{若 } p_r(\hat{\theta}) > 0 \text{ 且 } \langle \text{grad}\, p_r(\hat{\theta}), \varphi \rangle \leq 0 \\
\left[ I - \dfrac{p_r(\hat{\theta})\, \text{grad}\, p_r(\hat{\theta})\, \text{grad}\, p_r(\hat{\theta})^T}{\|\text{grad}\, p_r(\hat{\theta})\|^2} \right] \varphi, & \text{若 } p_r(\hat{\theta}) > 0 \text{ 且 } \langle \text{grad}\, p_r(\hat{\theta}), \varphi \rangle > 0
\end{cases}
$$

其中，$p_r(\hat{\theta}) = \dfrac{\|\hat{\theta}\|^2 - r_{\Omega_2}^2}{\epsilon_r^2 + 2\epsilon_r r_{\Omega_2}}$，$\text{grad}\, p_r(\hat{\theta}) = \dfrac{2\hat{\theta}}{\epsilon_r^2 + 2\epsilon_r r_{\Omega_2}}$，且 $\epsilon_r > 0
$ 为投影算子的平滑参数，$r_{\Omega_2}$ 为参数约束半径。

## 参考控制器和参考输入

为匹配扰动信号，我们设计了参考控制器，其状态向量满足：

$$
\dot{w}_c = R_c w_c
$$

其中矩阵 $R_c$ 依赖于扰动信号的频率参数 $\omega_1$ 和 $\omega_2$：

$$
R_c = 
\begin{bmatrix}
0 & 1 & 0 & 0 & 0 \\
-\theta_1^* & 0 & 1 & 0 & 0 \\
0 & 0 & 0 & 1 & 0 \\
-\theta_2^* & 0 & 0 & 0 & 1 \\
0 & 0 & 0 & 0 & 0 
\end{bmatrix}
$$

其中 $\theta_1^* = \omega_1^2 + \omega_2^2$ 和 $\theta_2^* = \omega_1^2 \cdot \omega_2^2$。

参考输入 $u_r(t)$ 通过以下方式计算：

$$
u_r(t) = \gamma_c w_c(t)
$$

其中 $\gamma_c = [1, 0, \ldots, 0]$ 是一个长度为 $2m+1$ 的向量，$m$ 是系统中频率成分的数量。需要注意的是，参考输入 $u_r(t)$ 的计算应当考虑扰动函数 $w_1(t)$ 和 $w_2(t)$ 的分段特性。关于此问题的详细讨论，请参阅 `reference_input_explanation.md` 文件。

## 控制律

最终的控制律设计为：

$$
u = -ke - \hat{\chi}_1 - \mu_1\hat{\theta}_1 - \mu_2\hat{\theta}_2
$$

其中 $k$ 是控制增益。理论上，参数 $\hat{\theta}_1$ 和 $\hat{\theta}_2$ 应当收敛到：

$$
\begin{align}
\hat{\theta}_1 &\rightarrow \omega_1^2 + \omega_2^2 \\
\hat{\theta}_2 &\rightarrow \omega_1^2 \cdot \omega_2^2
\end{align}
$$

## 稳定性分析

该控制系统的稳定性基于李雅普诺夫稳定性理论。考虑李雅普诺夫函数：

$$
V = \frac{1}{2}e^2 + \frac{1}{2g}\sum_{i=1}^{2}(\hat{\theta}_i - \theta_i^*)^2
$$

对李雅普诺夫函数求导可以证明在适当增益 $k$ 和 $g$ 的选择下，系统是稳定的，且输出调节误差 $e$ 最终会收敛到零。

## 仿真参数

在test2.m的仿真中，以下参数被使用：

- 系统参数：$a_1 = 0.5$, $b_2 = 0.5$, $\beta = 0.5$
- 扰动频率：$\omega_1 = 1$, $\omega_2 = 0.5$, $\omega_3 = 4$
- 控制器参数：$k = 7.1$, $g = 100$
- 赫尔维茨矩阵参数：$d = [2, 3, 2, 1]$
- 投影算子参数：$r_{\Omega_2} = 2$, $\epsilon_r = 0.1$
- 初始状态：$x_0 = [0.5; 0]$, $\hat{\theta}_0 = [0.1; 0.1]$

## 性能评估

控制性能通过以下指标评估：

1. 输出调节误差 $e_{output} = y - y_r$
2. 输入调节误差 $e_{input} = u - (-ke_{output})$
3. 参数估计收敛性 $\hat{\theta}_1 \rightarrow \omega_1^2 + \omega_2^2$ 和 $\hat{\theta}_2 \rightarrow \omega_1^2 \cdot \omega_2^2$

**MATLAB代码实现：**

在我们的测试文件中，投影算子的实现如下：

```matlab
% 投影算子的实现
function phi_proj = projection(phi, theta_hat, r_omega, epsilon_r)
    % 检查是否为标量情况
    if isscalar(theta_hat) && isscalar(phi)
        % 标量情况的处理
        % 计算p_r(theta_hat)
        p_r = (theta_hat^2 - r_omega^2)/(epsilon_r^2 + 2*epsilon_r*r_omega);
        
        % 计算梯度 grad p_r(theta_hat)
        grad_p_r = 2*theta_hat/(epsilon_r^2 + 2*epsilon_r*r_omega);
        
        % 根据条件判断使用哪种投影方式
        if p_r <= 0
            % 情况1：p_r(theta_hat) <= 0
            phi_proj = phi;
        else
            % 计算内积 <grad p_r(theta_hat), phi>
            inner_product = grad_p_r * phi;
            
            if inner_product <= 0
                % 情况2：p_r(theta_hat) > 0 且 <grad p_r(theta_hat), phi> <= 0
                phi_proj = phi;
            else
                % 情况3：p_r(theta_hat) > 0 且 <grad p_r(theta_hat), phi> > 0
                % 计算投影矩阵并应用到phi
                grad_norm_squared = grad_p_r^2;
                projection_factor = p_r * (grad_p_r^2) / grad_norm_squared;
                phi_proj = (1 - projection_factor) * phi;
            end
        end
    else
        % 向量情况的处理
        % 计算p_r(theta_hat)
        p_r = (norm(theta_hat)^2 - r_omega^2)/(epsilon_r^2 + 2*epsilon_r*r_omega);
        
        % 计算梯度 grad p_r(theta_hat)
        grad_p_r = 2*theta_hat/(epsilon_r^2 + 2*epsilon_r*r_omega);
        
        % 根据条件判断使用哪种投影方式
        if p_r <= 0
            % 情况1：p_r(theta_hat) <= 0
            phi_proj = phi;
        else
            % 计算内积 <grad p_r(theta_hat), phi>
            inner_product = dot(grad_p_r, phi);
            
            if inner_product <= 0
                % 情况2：p_r(theta_hat) > 0 且 <grad p_r(theta_hat), phi> <= 0
                phi_proj = phi;
            else
                % 情况3：p_r(theta_hat) > 0 且 <grad p_r(theta_hat), phi> > 0
                % 计算投影矩阵并应用到phi
                grad_norm_squared = norm(grad_p_r)^2;
                projection_matrix = eye(length(theta_hat)) - p_r * (grad_p_r * grad_p_r') / grad_norm_squared;
                phi_proj = projection_matrix * phi;
            end
        end
    end
end
```


