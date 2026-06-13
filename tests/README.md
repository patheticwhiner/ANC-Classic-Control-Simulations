# Test Suite — 统一测试基准

Phase 2 标准化测试框架 + Phase 3 各 demo 控制器接口。

## 文件结构

```
tests/
├── generate_test_signals.m   # Phase 2: 信号生成 (T1/T2/T3 + y_open)
├── compute_metrics.m         # Phase 2: 统一指标计算
├── controller_interface.m    # Phase 2: 接口规范 + init_result() 模板
│
├── run_demo1_test.m          # Phase 3.1: Youla 极点配置 + Q-RLS
├── run_demo2_test.m          # Phase 3.2: LQG (LQR + Kalman)
├── run_demo3_test.m          # Phase 3.3: Jafari F(z)+FIR+NLMS
├── run_demo4_test.m          # Phase 3.4: RST 自动设计 (降阶+Pole Placement)
│
├── README.md                 # 本文件
└── output/                   # Phase 4: 汇总结果
```

## 快速开始

### 1. 一次性路径配置

```matlab
run('tests/startup.m');
```

或手动添加所有必需路径（若 startup 不可用）。

### 2. 生成测试信号

```matlab
signals = generate_test_signals('armax_30303022');
```

输出: `dataset/test_signals_armax_30303022.mat`

### 3. 运行各 demo 测试

```matlab
data = load('dataset/test_signals_armax_30303022.mat');
s = data.signals;

r1 = run_demo1_test(s, 'T1', 'fixed');   % Youla Q (当前占位)
r2 = run_demo2_test(s, 'T1', 'fixed');   % LQG
r3 = run_demo3_test(s, 'T1', 'fixed');   % Jafari FIR
r3a = run_demo3_test(s, 'T2', 'adaptive'); % Jafari NLMS
r4 = run_demo4_test(s, 'T1', 'fixed');   % RST (当前占位)
```

### 4. 批量运行

以下脚本一次性跑完所有 demo × 测试 × variant 组合：

```matlab
run('tests/startup.m');
s = load('dataset/test_signals_armax_30303022.mat').signals;

% 定义各 demo 的可用 variant
trials = {
    @run_demo1_test, {'fixed'}
    @run_demo2_test, {'fixed'}
    @run_demo3_test, {'fixed', 'adaptive'}
    @run_demo4_test, {'fixed'}
};

fprintf('========== 批量测试 ==========\n');
for d = 1:size(trials,1)
    fh = trials{d,1};
    name = functions(fh).function;
    for t = {'T1','T2','T3'}
        for v = trials{d,2}
            try
                r = fh(s, t{1}, v{1});
                fprintf('%-6s %-9s %s: %6.1f dB | umax=%6.2f | conv=%4.1fs\n', ...
                    r.demo, r.variant, r.test, r.supp_db, r.u_max, r.t_conv_s);
            catch e
                fprintf('%-6s %-9s %s: ERROR — %s\n', ...
                    name, v{1}, t{1}, e.message);
            end
        end
    end
end
fprintf('========== 完成 ==========\n');
```

## 三组测试信号

| 编号 | 类型 | 参数 | 验证目标 |
|:---|:---|:---|:---|
| T1 | 固定频率正弦 | 334 Hz | 稳态抑制能力 |
| T2 | 线性扫频 | 300→500 Hz, 10s | 频率鲁棒性 (20Hz 分箱) |
| T3 | 带限白噪声 | 100–800 Hz | 宽带综合性能 |

统一参数: `Tsim=10s`, `rng(42)`, **RMS(d)=0.8**。

## 当前实验结果 (T1, 固定 variant)

| Demo | 方法 | 抑制 | 状态 |
|:---|:---|:---|:---|
| demo1 | Youla R₀=0 (占位) | 0.0 dB | 框架就位，R₀/S₀ 待调参 |
| demo2 | LQG (Q=1, R=1) | 4.8 dB | ✅ 稳定 |
| demo3 | FIR (N=120) | 21.2 dB | ✅ 稳定 |
| demo4 | RST 降阶 (balred→8) | 发散 | 框架就位，待调参 |

## 接口契约

见 `controller_interface.m` — 运行无参数调用可打印完整规范:

```matlab
controller_interface();
```

## 结果 struct

| 字段 | 说明 |
|:---|:---|
| `demo` | `'demo1'`–`'demo4'` |
| `variant` | `'fixed'` / `'adaptive'` |
| `test` | `'T1'` / `'T2'` / `'T3'` |
| `supp_db` | 稳态抑制 (dB) |
| `supp_breakdown` | T2 分箱抑制 |
| `t_conv_s` | 收敛时间 (s) |
| `u_max`, `u_rms` | 控制信号统计 |
| `y_open_rms`, `y_closed_rms` | 输出 RMS |
| `extra` | demo 特定字段 |
