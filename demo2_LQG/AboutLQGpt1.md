# LQGbasedANC

Originally Designed by [24s153046@stu.hit.edu.cn](mailto:24s153046@stu.hit.edu.cn)

---

excitation/bandlimitedNoise.m：根据采样频率、目标带限噪声的频带范围生成外系统模型（带通滤波器→状态空间模型）

sysId/systemidentification.m：对所设计的次级通路（带通滤波器）选取ARMAX模型辨识，并转换为状态空间模型

sysId/identificationDelaytest.m：由于使用MATLAB函数辨识ARMAX模型需要提前知道模型纯时延个数，可以使用程序测试最优的模型时延个数

LQG_ANC.m： 一个更老版本的LQGdemo.m，此时还试图仿真钱梵梵等人的论文，因此并未完全选用带限噪声控制

LQGdemo.m：将上述两部分程序生成的状态空间模型合并为增广状态空间模型，并据此设计LQR与LQE（Kalman Filter），完成ANC仿真

LQGdemo2.m: 尝试使用MATLAB内置的LQG函数一步到位设计LQG控制器，并完成ANC仿真

---

## 1 理论基础

### 1.1 LQ问题

若系统是线性的，且性能泛函是状态变量/控制变量的二次型函数的积分，则这样的最优控制问题称为线性二次型最优控制问题（Linear Quadratic Optimal Problem）。由于二次型性能指标具有鲜明的物理意义，它代表了大量工程实际问题中提出的性能指标要求，并且在数学处理上比较简单，易于通过状态线性反馈实现闭环最优控制，便于工程师先，因而在实际工程问题中得到了广泛应用。

+ **二次型最优控制问题**：对于以下**线性系统**，确定**最优控制律**$u^*(t)$，使**二次型性能指标**最小。就实际工程应用而言，此性能指标较全面地体现了对复杂控制系统的性能要求。
+   $L_x = \frac{1}{2} \boldsymbol{x}^\top(t)\boldsymbol{Q}(t)\boldsymbol{x}(t)$ 体现对动态过程的要求。
+   $L_u = \frac{1}{2} \boldsymbol{u}^\top(t)\boldsymbol{R}(t)\boldsymbol{u}(t)$ 体现对控制能量的限制。

$$
\dot{\boldsymbol{x}}(t) = \boldsymbol{A}(t)\boldsymbol{x}(t) + \boldsymbol{B}(t)\boldsymbol{u}(t), \quad \boldsymbol{x}(t_0) = \boldsymbol{x}_0\\
\boldsymbol{y}(t) = \boldsymbol{C}(t)\boldsymbol{x}(t)
$$

式中，
$\boldsymbol{x} \in \mathbb{R}^n$, $\boldsymbol{u} \in \mathbb{R}^r$, $\boldsymbol{y} \in \mathbb{R}^m$，
$\boldsymbol{A}(t)\in \mathbb{R}^{n \times n}$、
$\boldsymbol{B}(t)\in \mathbb{R}^{n \times r}$和
$\boldsymbol{C}(t)\in \mathbb{R}^{m \times n}$。

$$
J=\frac{1}{2}\boldsymbol{x}^\top(t_f)\boldsymbol{Q}_f\boldsymbol{x}(t_f)+\frac{1}{2}\int_{t_0}^{t_f}[\boldsymbol{x}^\top(t)\boldsymbol{Q}(t)\boldsymbol{x}(t)+\boldsymbol{u}^\top(t)\boldsymbol{R}(t)\boldsymbol{u}(t)]\text{d}t
$$


式中，
$\boldsymbol{Q}_f,\ \boldsymbol{Q}(t)\in \mathbb{R}^{n\times n} \geq 0$，
$\boldsymbol{R}(t)\in\mathbb{R}^{r\times r} >0$
均为对称矩阵。

*注：LQ性能指标最小的物理意义是：在整个时间区间 $[t_0,t_f]$ 内，综合考虑过程中偏差、控制消耗的能量和终值误差3个方面总的结果要最小。

+ **有限时间的线性最优调节器**

*注：有限时间最优调节器问题只考察控制系统由**任意初态**恢复到**平衡状态**的行为；

---

*注：工程上所关心的另一类更广泛的问题时：除保证有限时间内系统的**非零初态响应最优性**之外，还要求系统具有**保持平衡状态**的能力；既有**最优性要求**又有**稳定性要求**。此时如果将调节器问题推广到无限时间的情况，就可以在无限时间内既考察**实际上有限时间内的响应**，又考察**系统的稳定性**。

$$
\dot{\boldsymbol{x}}(t) = \boldsymbol{A}(t)\boldsymbol{x}(t) + \boldsymbol{B}(t)\boldsymbol{u}(t), \quad \boldsymbol{x}(t_0) = \boldsymbol{x}_0
$$

$$
\boldsymbol{x} \in \mathbb{R}^n,\ \boldsymbol{u} \in \mathbb{R}^r,\ \boldsymbol{A}(t)\in \mathbb{R}^{n \times n},\ \boldsymbol{B}(t)\in \mathbb{R}^{n \times r}
$$

$$
J = \frac{1}{2} \int_{t_0}^{\infty} \left[ \boldsymbol{x}^\top(t)\boldsymbol{Q}(t)\boldsymbol{x}(t) + \boldsymbol{u}^\top(t)\boldsymbol{R}(t)\boldsymbol{u}(t) \right] \text{d}t
$$

$$
\boldsymbol{u}^*(t)=\arg \min_u J
$$

+ **无限时间的线性最优调节器**

  与有限时间调节器的主要差别在于 $t_f$ 改为 $\infty$ 。

  + **定常线性最优调节器**

    $\boldsymbol{Q}\in \mathbb{R}^{n\times n}$，
    $\boldsymbol{R}\in\mathbb{R}^{r\times r} >0$均为常值对称矩阵，此时存在唯一的最优控制：
    $\boldsymbol{u}^*(t)=-\boldsymbol{R}^{-1}(t)\boldsymbol{B}^\top \boldsymbol{P}\boldsymbol{x}(t)$

### 1.2 LQR设计

[lqr](https://ww2.mathworks.cn/help/control/ref/lti.lqr.html)，[dlqr](https://ww2.mathworks.cn/help/control/ref/lti.dlqr.html)

MATLAB控制系统工具箱中提供了求解线性二次型（LQ）最优控制问题的函数及算法，其中 `lqr` 与 `lqry` 函数可直接求解定常连续系统LQR问题相关的Riccati代数方程。

+ `lqr` 连续时间状态空间系统，无限时间定常状态调节器

  提供连续时间状态空间模型 $A,\ B$ 以及二次代价函数权重参数 $Q,\ R,\ N$，可通过求解连续时间Riccati方程来获取最优的状态反馈权重矩阵。
+ `lqry` 连续时间状态空间系统，无限时间定常输出调节器

  提供连续时间状态空间模型 $A,\ B$ 以及二次代价函数权重参数 $Q,\ R,\ N$，可通过求解连续时间Riccati方程来获取最优的输出反馈权重矩阵。
+ `dlqr` 离散时间状态空间系统，无限时间定常状态调节器

  提供离散时间状态空间模型 $A,\ B$ 以及二次代价函数权重参数 $Q,\ R,\ N$，可通过求解离散时间Riccati方程来获取最优的状态反馈权重矩阵。

### 1.3 LQE/Kalman Filter观测器设计

[卡尔曼滤波](https://ww2.mathworks.cn/help/control/ug/kalman-filtering.html)

<img src="..\assets\kalmdemo_02.png" width = 70% />
$$
\begin{align}
x(k+1)=&Ax(k)+Bu(k)+Gw(k)\\
y(k)=&C x(k)+v(k)
\end{align}
$$

其中 $w(k),\ v(k)$ 均为零均值高斯白噪声，协方差为： $\mathbb{E}[w(k)w(k)^\top]=\boldsymbol{Q}$， $\mathbb{E}[v(k)v(k)^\top]=\boldsymbol{R}$， $\mathbb{E}[w(k)v(k)^\top]=0$。

状态估计：

$$
\begin{align}
\hat{x}^+(k)=&\hat{x}^-(k)+L[y(k)-C\hat{x}^-(k)-Du(k)]\\
\hat{x}^-(k+1)=&A\hat{x}^+(k)+Bu(k)
\end{align}
$$

---

为离散时间系统构造Kalman状态观测器可选择以下不同的函数：

+ `dlqe`  Kalman estimator design for discrete-time systems.

  `[M,P,Z,E] = dlqe(A,G,C,Q,R)`  returns the gain matrix M such that the discrete, stationary Kalman filter with observation and time update equations produces an optimal state estimate x[n|n] of x[n] given y[n] and the past measurements.  The resulting Kalman estimator can be formed with DESTIM.

+ `kalman`  Kalman filter design, Kalman estimator

  `[kest,L,P] = kalman(sys,Qn,Rn,Nn)` This MATLAB function creates a state-space model kest of the Kalman estimator given the plant model sys and the noise covariance data Qn, Rn, Nn (matrices Q, R, N described in Description).

### 1.4 用于调节问题的LQG设计

[线性二次高斯 (LQG) 设计](https://ww2.mathworks.cn/help/control/getstart/linear-quadratic-gaussian-lqg-design.html)

线性二次高斯（Linear Quadratic Gaussian）控制是一种用于设计伺服控制器的现代状态空间方法，允许我们权衡调节/过程扰动和测量噪声。LQR与Kalman Filter的设计可使用分离性原理，LQG控制设计遵循以下步骤：

1. 构造LQ最优增益；
2. 构造Kalman Filter/状态估计器；
3. 通过连接LQ最优增益与Kalman Filter构建LQG设计。

<img src="..\assets\regulator.png" width = 60% />
$$
\begin{align}
    x(k+1)=&Ax(k)+Bu(k)+Gw(k)\\
    y(k)=&Cx(k)+Du(k)+Hw(k)+v(k)
\end{align}
$$

对于上述系统设计分为两步：

+ LQR：$u(k)=-K{x}(k)$
+ Kalman Filter：${x}(k+1)={A}{x}(k)+{B}u(k)+L[y(k)-{C}x(k)]+Lv(k)$

---

基于 `lqg` 函数的一步到位的快速设计方法，需满足以下条件：

<img src="..\assets\lqg1.png" width =50% />



- 需要最优的 LQG 控制器，并且 $\mathbb{E}(wv^\top)$ 或 $H$ 为非零值。
- 所有已知（确定性）输入均为控制输入，所有输出均为测得值。
- 积分器状态的加权独立于被控对象的状态和控制输入。

被控系统为：
$$
\dot{x}=Ax+Bu+w\\
y=Cx+Du+v
$$
还需预先指定指标函数的**权重矩阵 $Q_{xu}$** 和过程噪声与测量噪声的**协方差矩阵 $Q_{wv}$**，具体而言有：
$$
J=\mathbb{E}\{\lim_{\tau\rightarrow\infty}\cfrac{1}{\tau}\int_0^\tau [x^\top,u^\top]^\top Q_{xu}\begin{bmatrix}x\\u\end{bmatrix}\text{d}t\}\\
Q_{wv} = \mathbb{E}\left(\begin{bmatrix}w\\v\end{bmatrix}[w^\top, v^\top]^\top\right)
$$
为了使LQG指令可用于设计伺服控制器，还可预先指定误差向量的**权重矩阵 $Q_i$**，此时指标函数为：
$$
J=\mathbb{E}\{\lim_{\tau\rightarrow\infty}\cfrac{1}{\tau}\int_0^\tau \left([x^\top,u^\top]^\top Q_{xu}\begin{bmatrix}x\\u\end{bmatrix}+x_i^\top Q_i x_i\right)\text{d}t\}\\
x_i=\int(r-y)
$$

---

若需要引入更多的设计灵活性以设计LQG调节器应当使用 `lqr` ， `kalman` 和 `lqgreg` 指令，相应的条件如下：

+ 可任意指定过程噪声矩阵G和测量噪声矩阵H
+ 可存在非控制输入的确定性干扰或未测量输出
+ 支持对积分器状态、系统状态和控制量的差异化加权

## 2 带限噪声生成

### 2.1 实验内容

1. 设置实验采取的采样频率（Hz），带通滤波器的高低截止频率（Hz），带通滤波器类型以及滤波器阶次，使用MATLAB函数设计带通滤波器；程序选用了常规的butterworth滤波器，根据选取的归一化频率，使用 `butter` 函数生成了相应的数字滤波器，输出结果为冲激响应函数/传递函数系数（离散时间模型）。
2. 使用 `tf2ss` 函数将上面得到的传递函数转换为相应的状态空间模型（不区分离散/连续）；
3. 输入信号选取全频带白噪声，状态空间模型仿真得到相应的带限噪声输出。

### 2.2 实验结果

* 绘制所生成带限噪声的时域波形，使用fft函数分析其频谱（幅度谱），验证频率范围与所设范围一致；
* 绘制所涉及的带通滤波器的Bode幅频曲线，以及最终所生成噪声的pwelch功率谱密度函数曲线，验证频率范围与所设范围一致。

<figure align = center>
    <img src="..\assets\bndltd1.svg" width = 49% />
    <img src="..\assets\bndltd2.svg" width = 49% />
</figure>


+ 将所得到的外系统（exosystem）状态空间模型导出为bandlimitedNoise.mat，可导入到其它脚本文件中用于设计相应的控制器。

$$
\begin{align}
x_w(k+1)=&A_wx_w(k)+B_we(k)\\
w(k)=&C_w x_w(k)
\end{align}
$$

## 3 ARMAX次级通路辨识

### 3.1 实验内容

1. 设置实验采取的采样频率（Hz），带通滤波器的高低截止频率（Hz），带通滤波器类型以及滤波器阶次，使用MATLAB函数设计带通滤波器；参考MATLAB的ActiveNoiseCancellation示例程序，选用了Chebyshev II型滤波器，输出结果为传递函数。
2. 输入信号选取全频带白噪声，使用上述带通滤波器（次级通路）滤波，得到辨识所需的输入、输出信号序列。
3. 设定ARMAX模型的阶次以及延迟，使用 `iddata` 以及 `armax` 函数辨识，`polydata` 函数提取多项式系数。并绘制系列图像验证辨识效果。
4. 使用 `tf2ss` 函数，将传递函数模型转换为状态空间模型。并绘制系列图像验证辨识、状态空间模型的频率响应是否等价。

### 3.2 实验结果

* 绘制所生成带限噪声的时域波形，使用fft函数分析其频谱（幅度谱），验证频率范围与所设范围一致；
* 绘制所涉及的带通滤波器的Bode幅频曲线，以及最终所生成噪声的pwelch功率谱密度函数曲线，验证频率范围与所设范围一致。

<figure align = center>
    <img src="..\assets\sysId1.svg" width = 49% />
    <img src="..\assets\sysId2.svg" width = 49% />
</figure>


+ 将所得到的次级通路/对象系统（plant）状态空间模型导出为systemIdentification.mat，可导入到其它脚本文件中用于设计相应的控制器。

$$
\begin{align}
x(k+1)=&Ax(k)+Bu(k)+Ew(k)\\
y(k)=&C x(k)+Du(k)+Fw(k)
\end{align}
$$

### 3.3 基于系统辨识工具箱

[利用Simulink实现系统模型辨识-连续扫频 - 知乎](https://zhuanlan.zhihu.com/p/598926217)

[MATLAB系统辨识工具箱（ARMAX模型） - 知乎](https://zhuanlan.zhihu.com/p/568419600)

## 4 基于LQR/LQG的主动控制

### 4.0 理论说明

LQR具有很好的稳定裕度，但在和Kalman Filter组成LQG之后，闭环系统的稳定裕度不复存在。

LQG系统虽然是“最优”的，担现实中极易因为建模误差和扰动而导致不稳定。因此，纯粹的LQG并不是一种很实用的控制策略，只适用于简单的、不确定性小的系统，或者适用于为了先实现控制随便找个算法凑合一下的情况。

为了解决LQG的鲁棒稳定性差的问题，有人在LQG的基础上提出了LTR（回路传递恢复）作为补充，通过频域分析，设计闭环零极点位置，从而保证了系统的鲁棒稳定性。也有人干脆放弃LQG，发展出了其他鲁棒控制算法，如H∞。

### 4.1 ANC问题制定

反馈控制结构的主动噪声控制问题本质上是控制中的调节问题（regulation），可以使用LQG控制器实现。

<img src="..\assets\LQGdemo.png" width = 60% />

### 4.2 增广系统模型

$$
\begin{bmatrix}
x(k+1)\\ x_w(k+1)
\end{bmatrix}=
\begin{bmatrix}
A& EC_w\\
0& A_w
\end{bmatrix}
\begin{bmatrix}
x(k)\\ x_w(k)
\end{bmatrix}+
\begin{bmatrix}
B\\ 0
\end{bmatrix}u(k)+
\begin{bmatrix}
0\\ B_w
\end{bmatrix}e(k)
$$

$$
y(k) = 
\begin{bmatrix}
C& FC_w
\end{bmatrix}
\begin{bmatrix}
x\\ x_w
\end{bmatrix}
$$

记作：

$$
\begin{align}
\tilde{x}(k+1)=&\tilde{A}\tilde{x}(k)+\tilde{B}u(k)+\tilde{G}e(k)\\
y(k)=&\tilde{C}\tilde{x}+v(k)
\end{align}
$$

+ LQR： $u(k)=-K\tilde{x}(k)$ 
+ Kalman Filter： $\tilde{x}(k+1)=\tilde{A}\tilde{x}(k)+\tilde{B}u(k)+L(y(k)-\tilde{C}x(k))+Lv(k)$ 

### 4.3 LQR+LQE设计

#### （1）LQR设计

根据参考文献[1]选取权重矩阵 $Q$ 与 $R$，控制性能权重：$ J = \sum_{k=0}^{\infty} (x_k^T Q x_k + u_k^T R u_k)$。其中 $Q$ 为状态误差惩罚，$R$ 为控制能量消耗的惩罚。选择 $Q=C^\top C$ 是一个经典而有效的技巧，它意味着使无噪声的理想输出 $\hat{y}^\top \hat{y}=x^\top C^\top Cx$趋于0，对应于ANC的最终目的；而控制输入权重$R=10^{-4}$则意味着我们不太在意控制输入的能量消耗，允许控制器使用较大的能量，只要能够有效地抵消噪声即可。

```matlab
%% 设计LQR控制器
% 状态权重矩阵 Q（半正定）
Q = C' * C;
% 输入权重矩阵 R（正定）
R = 1e-4;
% 离散 LQR 求解（状态反馈增益）
[K, ~, ~] = dlqr(A, B, Q, R);
```

#### （2）Kalman滤波器设计

参考文献[1]中选取过程噪声与测量噪声协方差矩阵（为保证收敛，测量噪声设置较小）如下。其中$Q_n$意味着过程噪声的强度，而$R_n$则代表了测量噪声的强度。换言之，$Q_n$代表着我们所认为的系统建模不确定性（或外部随机干扰），在仿真实验中对应着主噪声的幅值**（而在实际实验中可以通过事先测量得到该数值）**，$R_n$代表着我们所认为的传感器测量误差。在我们所应对的情况下，系统模型可能存在较大的不确定性（或存在较大的外部随机干扰）而测量噪声一般认为较小（Kalman滤波器将更多地新人它自己的模型预测）。

```matlab
%% 设计卡尔曼滤波器（状态估计器）
% 假设过程噪声协方差矩阵和测量噪声协方差矩阵
Qn = 2;  % 过程噪声协方差
Rn = 1e-10 * eye(size(C, 1));  % 测量噪声协方差
% 使用离散卡尔曼滤波器求解最优增益矩阵
[L, ~, ~] = dlqe(A, G, C, Qn, Rn);
```

#### （3）ANC实验仿真

仿真步骤设置为 `N=2000` 步，逐步更新状态反馈控制律、增广动态系统（含外系统噪声模型）、Kalman状态估计器。

+ 控制信号$u(k)$：记录仿真中状态反馈控制律提供的控制信号；
+ 反噪声信号$Cx(k)$：记录仿真中动态系统的状态向量$\tilde{x}(k)$；
+ 干扰信号$w(k)=C_wx_w(k)$：记录仿真中动态系统的状态向量$\tilde{x}(k)$.

<img src="..\assets\LQGdemo1.svg" width = 60% />

### 4.4 直接LQG设计（一步到位）

`lqg` 函数一步到位设计仅适用于常规的状态反馈控制，不适用于输出调节问题。

<font color = red>调节问题可看作r=0的伺服问题，那么lqg能否间接用于输出调节？</font>

### 4.5 基于Simulink的仿真

<img src="..\assets\simLQG.png" width=80% />



## 5 不足与改进建议

1. 所生成的传递函数究竟是连续复频域传递函数还是离散复频域传函？所生成的状态空间矩阵属于连续时间模型还是离散时间模型？仍有待理论说明。
2. 未对多种条件作重复试验，实验条件特殊，且设置比较随意，不确定推广到其它数据是否适用。

   例如，在仿真实验中，所涉及的测量误差 v 需远小于过程误差，才能保证具有较好的调节特性。若两者量级相对接近，则起不到任何调节的xiao
3. 未讨论LQR主动控制的局限性，未说明为什么要引入LQG；未探究LQR与LQG反馈控制的带宽限制。
4. 对于LQR，LQG设计的参数（权重矩阵选取）未提供理论参考，未探究超出建议值会发生什么

### 5.1 带限噪声的作用

### 5.2 辨识方法的选取

### 5.3 控制方法的限制

## 参考资料

[1] 钱梵梵. 基于Youla参数化的自适应输出调节及应用研究[D/OL]. 上海大学, 2022[2024-12-18]. [https://link.cnki.net/doi/10.27300/d.cnki.gshau.2022.000228](https://link.cnki.net/doi/10.27300/d.cnki.gshau.2022.000228). DOI:[10.27300/d.cnki.gshau.2022.000228](https://doi.org/10.27300/d.cnki.gshau.2022.000228).

[2] [Steve Brunton](https://www.youtube.com/@Eigensteve). Control Bootcamp: Linear Quadratic Gaussian (LQG)[视频/OL].  (2017-02-07). [2025-04-04]. [https://www.youtube.com/watch?v=H4_hFazBGxU.](https://www.youtube.com/watch?v=H4_hFazBGxU.)

[3] MathWorks. State-space control design: LQG/LQR and pole-placement algorithms[EB/OL]. (n.d.) [2025-04-05]. https://www.mathworks.com/help/control/ref/lqr.html.

[4] 三脚猫Frank. LQR控制器— 线性二次型调节器 Linear Quadratic Regulator[EB/OL]. (2021-12-17). [2025-04-07]. https://zhuanlan.zhihu.com/p/139145957

[5] [尚戈继](https://www.zhihu.com/people/moon-half). LQG输出调节控制Matlab仿真实例[EB/OL]. (2020-12-12). [2025-04-07]. https://zhuanlan.zhihu.com/p/338340273
