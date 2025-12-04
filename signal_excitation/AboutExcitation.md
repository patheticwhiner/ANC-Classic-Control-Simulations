# 常用激励信号

## 文件目录

```
excitation/
├── myPRBS.m 自定义的PRBS生成函数
│   
├── simBandLtdWN.slx
├── generatePRBS.m
├── testPRBS.m
├── testBandLtdWN.m
├── testPNmodel.m
│   
├── bandlimitedNoise.m
├── bandlimitedNoise2.m
├── dtmcTIDisturbance.m
├── dtmcTVDisturbance.m
│   
└── AboutExcitation.md 对此目录下的文件作了说明
```



## （1）MATLAB Simulink的Band-Limited White Noise模块

```
testBandLtdWN.slx
testBandLtdWN.m
```

以下是为深入理解 Band-Limited White Noise 模块设计的实验方案，通过调整参数（采样时间、噪声功率、随机种子）对比时域特性、频谱特性及数据统计性质。实验步骤包含建模、参数设置、数据分析和可视化。Band-Limited White Noise模块生成的是高斯白噪声。

**实验目标**

1. 理解 `Sample Time` 对信号带宽的影响。
2. 分析 `Noise Power` 对信号幅值和功率谱密度（PSD）的影响。
3. 验证随机种子（`Seed`）对信号可重复性的作用。
4. 探索数据长度与统计特性的关系。

**实验步骤**

创建一个简单模型，包含以下模块：

• Band-Limited White Noise：生成带限白噪声。

• To Workspace：保存输出数据到MATLAB工作区（变量名设为`noise_data`，格式设为`Array`）。

• Scope（可选）：实时观察时域波形。

• 求解器：`Fixed-Step`，步长设为 `0.001秒`（需与噪声模块的 `Sample Time` 协调）。

---

#### 实验 1：采样时间（Sample Time）的影响

**目标**：观察 `Sample Time` 对信号带宽和时域特性的影响。

<font color = red>带限白噪声的采样频率/白噪声带宽由采样时间决定，并受到系统采样频率的约束（Simulink中设置的求解器仿真步长）</font>

**参数设置**

  • `Noise Power`：固定为 `0.1`。

  • `Seed`：固定为 `23341`（保证可重复性）。

  • `Sample Time`：分别设为 `0.001`, `0.002`, `0.005` , `0.01` 秒。

  • `sim_time`：仿真步数固定设置为1e5，仿真时间为 1e5×<u>Sample Time</u>

**预期结果**

• 时域：`Sample Time` 越小，信号变化越快速（高频分量更多）。

• 自相关函数：`Sample Time` 越小，自相关衰减越快（更接近理想白噪声）。

• PSD：带宽由 `1/(2*Sample Time)` 决定。

<figure>
    <img src="subassets\bndLtdWNfig1-1.svg" width = 49% />
    <img src="subassets\bndLtdWNfig1-2.svg" width = 49% />
    <img src="subassets\bndLtdWNfig1-3.svg" width = 49% />
    <img src="subassets\bndLtdWNfig1-4.svg" width = 49% />
</figure>



**实验结果**

• 时域：在相同种子的条件下，相同仿真步数的几个白噪声信号序列完全相同。

• 频域：频谱形状也完全相同，带宽由 `1/(2*Sample Time)` 决定，对应于采样频率更高的情况带宽更高。

• 自相关函数：`Sample Time` 越小，自相关衰减越快，更接近理想白噪声。


---

#### 实验 2：噪声功率（Noise Power）的影响

• 目标：分析 `Noise Power` 对信号幅值和功率谱密度的影响。

• 参数设置：

  • `Sample Time`：固定为 `0.001s`。

  • `Seed`：固定为 `23341`。

  • `Noise Power`：分别设为 `0.001`,`0.01`, `0.1`。

**预期结果**

• 时域幅值：`Noise Power` 越大，信号幅值标准差越大（直方图更宽）。

• PSD总功率：实际功率应与 `Noise Power` 匹配（误差由频谱估计方法引起）。

<figure>
	<img src="subassets\bndLtdWNfig2-1.svg" width = 32% />
    <img src="subassets\bndLtdWNfig2-2.svg" width = 32% />
    <img src="subassets\bndLtdWNfig2-3.svg" width = 32% />
</figure>
**实验结果**

• PSD平均数值为 PSD = 2 × Noise_Power = 2e-3，2e-2，2e-1 （V^2/Hz）

• 信号幅值分别为 Amplitude = sqrt( PSD × Fs ) = sqrt(2)，sqrt(20)，sqrt(200)；
有效值为 RMS = sqrt( Noise_Power × Fs )

• 自相关函数完全相同


---

#### 实验 3：随机种子（Seed）的影响

• 目标：验证种子对信号可重复性的作用。

**实验结果**

• 使用相同的种子进行两次实验产生的信号完全相同，固定种子可确保实验可重复性，适合调试；

• 不同种子生成的信号完全不同，但统计特性（均值、方差、PSD）一致。


---

#### 实验 4：数据长度的影响

• 目标：分析数据长度对统计估计（如PSD、方差）的影响。

**参数设置**：

  • `Sample Time=0.01s`, `Noise Power=0.1`, `Seed=23341`。

  • 仿真序列长度分别设为 `1e2`, `1e3`, `1e4`。

**预期结果**
• 方差估计：数据越长，估计值越接近理论值（`Noise Power / Sample Time`）。

• PSD平滑度：数据越长，PSD曲线越平滑（减小估计方差）。

• PSD平滑度：数据越长，自相关函数越接近白噪声。

<figure>
	<img src="subassets\bndLtdWNfig4-0.svg" width = 32% />
    <img src="subassets\bndLtdWNfig4-1.svg" width = 32% />
    <img src="subassets\bndLtdWNfig4-2.svg" width = 32% />
</figure>
#### 汇总

- **直接使用白噪声**：Simulink的`Band-Limited White Noise`模块生成的噪声通常可以直接用于系统辨识，因为它近似高斯分布。
- **是否需要高斯分布**：在大多数情况下，高斯分布假设是系统辨识方法的理论基础，但实际应用中对噪声分布的要求并不严格。
- **特殊情况处理**：如果噪声分布偏离高斯分布较远或是有色噪声，可以考虑对噪声进行预处理。

#### 附：模块的Noise Power概念辨析

Simulink的`Band-Limited White Noise`模块的`Cov`参数定义的是信号的**总功率（方差）**，而不是直接的单边PSD值。为了与`pwelch`计算的单边PSD结果一致，需要将总功率转换为单边PSD值。

**1. Simulink模块的`Cov`定义**

- `Cov`表示信号的**总功率（方差）**，即信号在整个频率范围（双边谱）上的总能量。
- 对于带限白噪声，信号的总功率均匀分布在频率范围`[-fs/2, fs/2]`（双边谱）。

因此，`Cov`实际上是双边谱的总功率，而不是单边PSD。

**2. 单边PSD与双边PSD的关系**

对于实信号，单边PSD与双边PSD的关系如下：

- 单边PSD只显示正频率部分（`[0, fs/2]`），但包含了负频率部分的功率。
- 因此，单边PSD的值是双边PSD的两倍（除了直流分量和奈奎斯特频率）。

如果`Cov`表示双边谱的总功率，那么单边PSD的理论值为：`theoretical_psd = 2 * noise_cov;`

**3. `pwelch`计算的PSD是单边谱**

MATLAB的`pwelch`函数默认计算的是单边PSD，频率范围为`[0, fs/2]`。为了与`pwelch`的结果一致，需要将`Cov`转换为单边PSD值。

**4. 总结**

- Simulink模块的`Cov`参数定义的是信号的总功率（方差），对应的是双边谱的总功率。
- 为了与`pwelch`计算的单边PSD结果一致，需要将`Cov`乘以2，得到单边PSD的理论值：`theoretical_psd = 2 * noise_cov;`
- 这种转换是因为单边PSD包含了负频率部分的功率。

## （2）MATLAB的prbs函数

[matlab——simulink从工作空间导入数据作为输入信号进行仿真-CSDN博客](https://blog.csdn.net/binheon/article/details/90295138)



## （3）BandLtdWN和PRBS如何选取？

#### 1. 信号特性差异

| 特性         | PRBS白噪声             | Simulink带限白噪声       |
| ------------ | ---------------------- | ------------------------ |
| **幅值分布** | 离散值（通常仅±1）     | 连续的高斯分布           |
| **周期性**   | 有固定周期(2^n-1)      | 无周期性（真随机）       |
| **确定性**   | 完全确定的序列         | 随机过程                 |
| **可重复性** | 同参数总是产生相同序列 | 需要固定随机种子才能重复 |

#### 2. 频谱特性

**PRBS白噪声**：

- 在特定频率处存在**谱线凹陷**(尤其在fs/p处，p为分频系数)
- 频谱呈现明显的**周期结构**
- 低频处表现较好，高频处会有衰减
- 频谱形状受PRBS级数和分频系数影响

**Simulink带限白噪声**：

- 在整个带宽内功率谱密度更加**均匀平坦**
- 不存在特定频率的谱线凹陷
- 在奈奎斯特频率(fs/2)附近有平滑的截止特性
- 频谱形状主要受采样频率影响

#### 3. 统计特性

**PRBS白噪声**：

- 自相关函数在零时滞处有峰值，其他位置近似为零，但会在周期处重复
- 非高斯分布（实际是二值分布）
- 固定的峰值因子（crest factor）

**Simulink带限白噪声**：

- 自相关函数更接近理想的冲击函数
- 严格的高斯分布
- 峰值因子有一定变化

#### 4. 系统辨识应用中的差异

**PRBS适用场景**：

- 对系统进行**线性辨识**时更有效
- 可以通过平均多个周期降低噪声影响
- 二值信号对于一些系统具有更好的激励能力
- 对非线性系统可能不够敏感

**带限白噪声适用场景**：

- 更接近实际工作环境中的噪声
- 适合需要连续幅值分布的场景
- 对于噪声统计特性要求严格的应用更合适
- 适用于高斯噪声假设的建模方法

#### 5. 生成效率

**PRBS白噪声**：

- 计算效率高，生成简单
- 存储需求小（仅需存储一个周期）

**Simulink带限白噪声**：

- 计算量稍大（需要生成随机数和滤波）
- 存储需求大（需要存储全部样本）

PRBS在系统辨识中被广泛使用，主要是因为其良好的激励特性和可重复性，而带限白噪声则更接近真实环境中的随机扰动。

### PRBS相较于后者的优势

#### 1. 信号激励能力优势

**PRBS的持续激励性(Persistent Excitation)更强**：

- PRBS信号具有最大峰值因子(Crest Factor)，提供更大幅度的输入激励
- 二值信号使得信号能量集中在±1，而非分散在多个幅度值上
- PRBS的能量在整个频带内分布更均匀(除特定凹陷频率外)，提供更有效的频谱覆盖

#### 2. 信号处理和分析优势

**PRBS具有确定性结构**：

- 自相关函数接近理想的冲击函数(δ函数)，使得线性系统的卷积计算简化
- 周期性使得可以通过平均多个周期来提高信噪比
- 确定性结构允许精确预测输入信号的每个样本值

#### 3. 实际应用优势

**实现和重复性**：

- 二值信号更易于在实际系统中生成(如通过继电器或开关)
- 对于模拟系统，不需要精确的多级电压控制，降低硬件复杂度
- 实验的精确重复性使得结果验证和比对更可靠

#### 4. 线性系统特性匹配

**线性系统的理论匹配**：

- 线性系统辨识主要关注系统的频率响应和脉冲响应
- PRBS的自相关特性使得系统输出与脉冲响应直接相关
- 线性系统对二值输入的响应包含了与连续幅值输入相同的信息

#### 5. 数据效率

**数据处理效率**：

- PRBS需要更少的样本量就能提供足够的系统信息
- 二值数据处理计算量小，特别适合实时或嵌入式系统辨识
- 参数计算收敛更快，需要较少的迭代

#### 总结对比

虽然带限高斯白噪声从纯理论角度可能更"理想"(特别是对于非线性系统或需要高斯噪声假设的方法)，但PRBS在实际的线性系统辨识中往往能提供更高效、更可靠的结果。

实际上，许多经典的系统辨识方法(如最小二乘法)在开发时就考虑了PRBS等伪随机信号的特性，因此形成了天然的优势互补。

## 随机性干扰的生成

### bandlimitedNoise.m

#### 功能

- 生成带限噪声信号，限制在指定频率范围内。
- 使用 Butterworth 带通滤波器设计，并转换为状态空间模型。
- 绘制时域波形、滤波器频率响应和功率谱密度（PSD）。

#### 使用方法

1. **参数设置**：
   - `fs`：采样频率（Hz）。
   - `f_low` 和 `f_high`：低、高截止频率（Hz）。
   - `order`：滤波器阶数。
2. **运行脚本**：
   - 在 MATLAB 中运行 [bandlimitedNoise.m](vscode-file://vscode-app/d:/Microsoft VS Code/resources/app/out/vs/code/electron-sandbox/workbench/workbench.html)。
3. **输出结果**：
   - 带限噪声信号 `w`。
   - 滤波器状态空间模型参数 `Aw`、`Bw`、`Cw`、`Dw`。
   - 图像

#### 输出文件

- `bandlimitedNoise.mat`：保存滤波器的状态空间模型参数。

#### 注意事项

- 确保截止频率满足 `0 < f_low, f_high < fs/2`。
- 滤波器阶数影响频率响应的陡峭程度。

#### 实验内容

1. 设置实验采取的采样频率（Hz），带通滤波器的高低截止频率（Hz），带通滤波器类型以及滤波器阶次，使用MATLAB函数设计带通滤波器；程序选用了常规的butterworth滤波器，根据选取的归一化频率，使用 `butter` 函数生成了相应的数字滤波器，输出结果为冲激响应函数/传递函数系数（离散时间模型）。
2. 使用 `tf2ss` 函数将上面得到的传递函数转换为相应的状态空间模型（不区分离散/连续）；
3. 输入信号选取全频带白噪声，状态空间模型仿真得到相应的带限噪声输出。

#### 实验结果

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



## 确定性干扰的生成

### dtmcTIDisturbance.m

### dtmcTVDisturbance.m

### MATLAB 仿真示例

下面给出一个 MATLAB 代码示例，演示如何使用差分形式的离散状态空间模型生成含有两个频率（例如 5Hz 和 10Hz）分量的简谐信号。代码清晰展示了状态更新和输出计算过程：

```matlab
% 采样参数
fs = 100;            % 采样频率 (Hz)
dt = 1/fs;           % 采样周期
% 两个简谐成分的频率
f1 = 5; f2 = 10;     
omega1 = 2*pi*f1; omega2 = 2*pi*f2;
% 构造每个频率的离散转移矩阵
A1 = [cos(omega1*dt), -sin(omega1*dt);
      sin(omega1*dt),  cos(omega1*dt)];
A2 = [cos(omega2*dt), -sin(omega2*dt);
      sin(omega2*dt),  cos(omega2*dt)];
% 合并成块对角矩阵
Ad = blkdiag(A1, A2);
% 输出矩阵，将两路信号的第一个分量相加
C = [1, 0, 1, 0];    
% 初始状态：两路信号均从相位0 (cos=1) 开始
x = [1; 0; 1; 0];    
N = 1000;            % 仿真步数
Y = zeros(1, N);
% 迭代仿真
for k = 1:N
    x = Ad * x;          % 状态更新 x[k+1] = A_d * x[k]
    Y(k) = C * x;        % 输出 y[k] = C*x[k] (两路信号之和)
end
% 绘图显示结果
t = (0:N-1)*dt;
plot(t, Y);
title('Two-tone discrete-time signal');
xlabel('Time (s)'); ylabel('Amplitude');
```

在这段代码中，首先设定采样周期 $dt$ 和两个频率 $\omega_1, \omega_2$。然后分别构造对应的 $2\times2$ 离散转移矩阵 $A_1,A_2$。通过 `blkdiag` 将它们组合成一个 4×4 的 $Ad$，对应两路简谐振荡。状态向量 $x=[x_1,x_2,x_3,x_4]^\top$ 包含了两路信号的 $(\cos,\sin)$ 分量。输出矩阵 $C=[1,0,1,0]$ 表示将两段信号的余弦分量相加得到最终输出。循环迭代计算状态更新 $x_{k+1}=A_d x_k$ 并计算输出 $y_k=Cx_k$，最后绘制 $y$ 随时间的变化。结果将显示两个正弦成分叠加的复合波形。