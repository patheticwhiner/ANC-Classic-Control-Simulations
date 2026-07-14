# `tools/` 仓库维护工具

`tools/` 存放显式执行的仓库维护、报告生成和 Git 工作流工具，不属于运行期
公共函数库，也不由 `project_init.m` 自动加入路径。

- `clearMlxOutputs.m`：清理 MATLAB Live Script 的保存输出。
- `generate_model_report.m`：扫描 `dataset/` 与 `sysid_models/`，生成模型来源报告。
- `pre-commit`：Git 提交前清理 `.mlx` 和 `.m` 的实时脚本输出标记。

使用报告生成器时显式加载：

```matlab
addpath('dataset', 'tools');
generate_model_report();
```
