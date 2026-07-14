# `shared/` 公共运行期代码

本目录存放被多个 demo、测试或模型分析入口复用的 MATLAB 函数，按职责分域。
由根目录 `project_init.m` 统一加入 MATLAB 路径；入口脚本不应再硬编码旧的
`functions` 相对路径或依赖当前工作目录。

## 子目录职责

| 子目录 | 职责 | 代表函数 |
|---|---|---|
| `rst_synth/` | RST / Diophantine 综合、多项式代数和稳定性判断 | `bezoutd`, `addPolynomials`, `isschur` |
| `dsp/` | 通用离散信号处理原语 | `FIR`, `IIR_filter`, `realtimeIIR`, `sat` |
| `control_analysis/` | 零极点、灵敏度和分析结果展示 | `analyze_zeros_poles`, `plotSensitivity` |

模型专用的 ARMAX 分析入口位于 `dataset/analysis/`，因为它依赖
`DataManager` 和 `dataset` 中的模型注册表，不属于通用公共库。

新增代码按本质职责归类，而不是按首次使用它的 demo 归类。只被单个复现使用的
算法实现应留在该复现自己的目录（例如 `demo4_Robust/utils/`）。

部分 Live Script 和测试文件保留了局部同名函数，作为原始复现或独立测试实现。
MATLAB 会优先使用这些局部函数，因此修改 `shared/` 不会自动改变这些自包含实现。

## 入口约定

```matlab
projectRoot = fileparts(mfilename('fullpath'));
run(fullfile(projectRoot, 'project_init.m'));
```

若入口位于子目录，先用 `fileparts` 得到自身目录，再按实际深度定位根目录。
`project_init.m` 不加载 `tools/`，维护脚本应由调用者显式加入。
