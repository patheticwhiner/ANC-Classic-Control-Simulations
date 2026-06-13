function result = controller_interface(varargin)
% CONTROLLER_INTERFACE  Phase 3 控制器接口规范与模板
%
%   本文件定义各 demo 控制器在 Phase 3 中应遵循的统一接口契约。
%   各 demo 需实现一个形如以下签名的函数：
%
%     function result = run_demoX_test(signals, test_name, variant, params)
%
%   其中 X ∈ {1,2,3,4}，参数说明见下方。
%
% -------------------------------------------------------------------------
% 接口契约
% -------------------------------------------------------------------------
%
% 输入:
%   signals    — Phase 2 生成的测试信号 struct（load 自
%                dataset/test_signals_<model>.mat）
%                字段: .model, .fs, .rng_seed, .norm_target, .norm_value,
%                      .T1, .T2, .T3
%                其中 .T? 含: .t, .d, .y_open (, .f_inst for T2)
%
%   test_name  — 测试编号: 'T1' | 'T2' | 'T3'
%
%   variant    — 控制器模式: 'fixed' | 'adaptive'
%                'fixed'     — 仅运行控制器的固定（非自适应）部分
%                'adaptive'  — 运行完整的自适应控制器
%                注意: 若 demo 不支持某 variant，应 error() 并给出明确消息
%
%   params     — demo 专用参数 struct。各 demo 自行定义所需字段。
%                建议包含:
%                  .modelName — 模型名，用于 DataManager 加载
%                示例 (demo3):
%                  params.N_fir = 120;
%                  params.mu_nlms = 0.05;
%                  params.N_adapt = 64;
%
% 输出:
%   result — 标准化 struct，由 compute_metrics() 生成，字段:
%     .demo            — 'demo1'|'demo2'|'demo3'|'demo4'
%     .variant         — 'fixed'|'adaptive'
%     .test            — 'T1'|'T2'|'T3'
%     .fs              — 采样频率 (Hz)
%     .Tsim            — 仿真时长 (s)
%     .supp_db         — 稳态抑制 (dB)
%     .supp_breakdown  — T2 专用：频率分箱抑制 (struct array)
%     .t_conv_s        — 收敛时间 (s)，固定控制器 = 0
%     .u_max           — max|u|
%     .u_rms           — RMS(u)
%     .y_open_rms      — 开环 RMS (稳态窗口)
%     .y_closed_rms    — 闭环 RMS (稳态窗口)
%     .d_rms           — 扰动 RMS (稳态窗口)
%     .extra           — demo 特定字段（由调用方在 compute_metrics 后追加）
%
% -------------------------------------------------------------------------
% 典型调用流程 (各 demo 的实现参考)
% -------------------------------------------------------------------------
%
%   function result = run_demo3_test(signals, test_name, variant, params)
%       % 1. 加载模型
%       info = DataManager(signals.model);
%       ...
%       % 2. 设计控制器（如果是第一次调用、或 params 要求重新设计）
%       theta_fixed = design_fixed_controller(info, params);
%       ...
%       % 3. 选择测试信号
%       test_sig = signals.(test_name);
%       ...
%       % 4. 运行仿真
%       if strcmp(variant, 'fixed')
%           y_closed = sim_fixed(test_sig.d, test_sig.Tsim);
%       else
%           y_closed = sim_adaptive(test_sig.d, test_sig.Tsim);
%       end
%       ...
%       % 5. 计算指标
%       meta = struct('demo', 'demo3', 'variant', variant, 'test', test_name);
%       result = compute_metrics(test_sig.y_open, y_closed, u, test_sig, meta);
%       ...
%       % 6. 追加 demo 专用字段
%       result.extra.gamma = gamma_fixed;
%       result.extra.norm_theta = norm(theta_fixed);
%   end
%
% -------------------------------------------------------------------------
% init_result — 返回标准化空 struct 模板
% -------------------------------------------------------------------------
%
%   用法:
%     result = init_result()
%     result = init_result('demo', 'demo3', 'variant', 'fixed', 'test', 'T1')
%
%   用于确保字段顺序和默认值一致，方便 Phase 4 汇总脚本处理。

    if nargin == 0
        % 无参数调用：打印接口文档
        help controller_interface;
        return;
    end

    % 解析键值对参数
    p = inputParser;
    p.addParameter('demo',    '');
    p.addParameter('variant', '');
    p.addParameter('test',    '');
    p.addParameter('fs',      []);
    p.addParameter('Tsim',    []);
    p.parse(varargin{:});
    opts = p.Results;

    result = struct();
    result.demo           = opts.demo;
    result.variant        = opts.variant;
    result.test           = opts.test;
    result.fs             = opts.fs;
    result.Tsim           = opts.Tsim;
    result.supp_db        = [];
    result.supp_breakdown = [];
    result.t_conv_s       = [];
    result.u_max          = [];
    result.u_rms          = [];
    result.y_open_rms     = [];
    result.y_closed_rms   = [];
    result.d_rms          = [];
    result.extra          = [];
end
