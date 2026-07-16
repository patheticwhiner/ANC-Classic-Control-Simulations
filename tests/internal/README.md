# Tests Shared Helpers

本目录只保存多个测试入口共享的函数，每个文件对应一个独立职责，
单一调用方的辅助逻辑都已内联为各文件的局部函数，不存在更深的包装层。

- `controller_demo1.m`：低阶 RST + Youla-Q RLS（含 RST 设计与限幅局部函数）。
- `controller_demo2.m`：增广 LQG + Q-RLS/LMS。
- `controller_demo3.m`：F(z)+FIR + Filtered-X LMS。
- `controller_demo4.m`：eMOPSO 灵敏度整形中央 RST + Youla-Q RLS（优化接线、
  Pareto 选解与自适应仿真为局部函数，eMOPSO 核心复用 `demo4_Robust/utils/`）。
- `controller_demo5.m`：Marino-Tomei Case A 纯反馈窄带频率估计器；显式记录
  离散有效反馈符号、Case A 适用性、频率估计历史和未限幅需求；正式候选使用
  `output_timing='updated'`；`previous` 仅作为 ARMAX 离散相位边界诊断。
- `controller_imc_fxlms.m`：Demo5 的纯反馈 IMC-FxLMS 基线。
- `load_cylinder1dm_signals.m`：加载或重新生成校准/评价信号套件（含信号生成局部函数）。
- `compute_metrics.m`：标准化抑制/控制指标。
- `generate_demo1234_reports.m`：Demo1-5 T1/T2 主实验分析图和 Markdown 报告生成器；
  显式传入 `partFilter='t3'` 时只负责宽带图像。
- `plot_demo_overview.m`：五个 Demo 共用的冻结评价综合图，统一绘制 T1
  fixed/adaptive、在线参数状态、T2 局部分箱和未限幅控制需求。
- `generate_broadband_report.m`：独立 T3 宽带失败边界报告生成器。

人工运行入口位于 `tests/` 根目录。请先运行 `startup.m`，然后运行五个
`demo*_cylinder1dm.m` 脚本或 `run_cylinder1dm_stage()`。宽带实验单独运行
`run_cylinder1dm_stage('broadband')`。
