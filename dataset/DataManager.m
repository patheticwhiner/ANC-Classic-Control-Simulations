function info = DataManager(dataSource)
% DataManager.m - 数据集加载与模型导入统一入口
%
% 用法:
%   info = DataManager()              % 列出所有可用数据源
%   info = DataManager(dataSource)    % 加载指定数据源
%
% 输入:
%   dataSource - 字符串，指定要加载的数据源名称
%
% 输出:
%   info - struct，包含:
%     .type    - 数据类型: 'armax_model' | 'ss_model' | 'lms_data' | 'raw_signal'
%     .name    - 数据源名称
%     .fs      - 采样频率 (Hz)
%     以及该类型对应的数据字段 (见下方各 case)
%
% 可用数据源:
%   'armax_30303022'      - ARMAX(30,30,30,22) 辨识模型，实测声学管道 @48kHz
%   'syn_TAC2015_3rd'     - Jafari TAC 2015 离散AVC (3阶, Ts=1/480)
%   'syn_JVC2017_3rd'     - Jafari JVC 2017 连续CC  (3阶, fs=10k)
%   'syn_JVC2017_6th'     - Jafari JVC 2017 高阶示例 (6阶, 连续)
%   'syn_Bai1997_4th'     - Bai 1997 耳机H∞模型   (4阶, fs=4k)
%   'syn_Carmona2000_7th' - Carmona 2000 管道ANC   (7阶, Fs=2k)
%   'syn_MassSpringDamper_2nd' - 质量-弹簧-阻尼器  (2阶, 连续)
%   'syn_RSTtoy_2nd'      - RST教学模型           (2阶, Ts=1)
%   'syn_whitenoise'      - 合成带限白噪声干扰模型 (状态空间)
%   'syn_bpf'             - 合成带通滤波器被控对象模型 (状态空间)
%   'lms_sysid'           - LMS系统辨识实验数据 (4通道, 2026-01-20)
%   'raw_dspace'          - dSPACE采集原始信号数据

    % 数据集根目录
    dsDir = fileparts(mfilename('fullpath'));

    % 可用数据源注册表
    sources = struct(...
        'armax_30303022',       struct('file', 'armax_30303022_2026-01-20.mat', 'type', 'armax_model'), ...
        'syn_TAC2015_3rd',      struct('file', 'syn_TAC2015_3rd.mat',           'type', 'syn_tf'), ...
        'syn_JVC2017_3rd',      struct('file', 'syn_JVC2017_3rd.mat',           'type', 'syn_tf'), ...
        'syn_JVC2017_6th',      struct('file', 'syn_JVC2017_6th.mat',           'type', 'syn_tf'), ...
        'syn_Bai1997_4th',      struct('file', 'syn_Bai1997_4th.mat',           'type', 'syn_tf'), ...
        'syn_Carmona2000_7th',  struct('file', 'syn_Carmona2000_7th.mat',       'type', 'syn_tf'), ...
        'syn_MassSpringDamper_2nd', struct('file', 'syn_MassSpringDamper_2nd.mat', 'type', 'syn_tf'), ...
        'syn_RSTtoy_2nd',       struct('file', 'syn_RSTtoy_2nd.mat',            'type', 'syn_tf'), ...
        'syn_whitenoise',       struct('file', 'syn_whitenoise_ssmodel.mat',    'type', 'ss_model'), ...
        'syn_bpf',              struct('file', 'syn_bpf_ssmodel.mat',           'type', 'ss_model'), ...
        'lms_sysid',            struct('file', 'lms_sysid_2026-01-20.mat',      'type', 'lms_data'), ...
        'raw_dspace',           struct('file', 'raw_dspace_primpath.mat',       'type', 'raw_signal') ...
    );

    % 无参数调用：列出所有可用数据源
    if nargin == 0 || isempty(dataSource)
        fprintf('=== dataset/ 可用数据源 ===\n');
        names = fieldnames(sources);
        for i = 1:length(names)
            nm = names{i};
            s = sources.(nm);
            fprintf('  %-18s [%s]  %s\n', nm, s.type, s.file);
        end
        fprintf('\n用法: info = DataManager(''<name>'')\n');
        if nargout == 0, return; end
        info = sources;
        return;
    end

    % 查找数据源
    if ~isfield(sources, dataSource)
        error('未知数据源: ''%s''。使用 DataManager() 查看可用列表。', dataSource);
    end

    src = sources.(dataSource);
    fprintf('正在加载数据源: %s (%s)...\n', dataSource, src.file);

    % 统一加载
    fullPath = fullfile(dsDir, src.file);
    if ~exist(fullPath, 'file')
        error('数据文件不存在: %s', fullPath);
    end

    data = load(fullPath);
    info.name = dataSource;
    info.type = src.type;

    switch src.type
        case 'armax_model'
            % ARMAX辨识模型
            info.model  = data.ARMAXmodel.model;   % idpoly 对象
            info.orders = data.ARMAXmodel.orders;  % [na nb nc nk]
            info.fs     = data.ARMAXmodel.fs;

            fprintf('  ARMAX(%d,%d,%d,%d), fs=%d Hz\n', ...
                info.orders(1), info.orders(2), info.orders(3), info.orders(4), info.fs);

        case 'syn_tf'
            % 合成传递函数模型 (论文/教材中的理论被控对象)
            info.G0     = data.model.G0;       % tf/zpk 对象
            info.domain = data.model.domain;   % 'continuous' | 'discrete'
            info.orders = data.model.orders;
            info.source = data.model.source;
            info.desc   = data.model.desc;
            if isfield(data.model, 'Ts'),      info.Ts = data.model.Ts; end
            if isfield(data.model, 'fs'),      info.fs = data.model.fs; end
            if isfield(data.model, 'fs_nominal'), info.fs = data.model.fs_nominal; end
            if isfield(data.model, 'G0_tf'),   info.G0_tf = data.model.G0_tf; end
            if isfield(data.model, 'G0_zpk'),  info.G0_zpk = data.model.G0_zpk; end

            fprintf('  %s (%s, %s)\n', data.model.name, data.model.domain, data.model.source);

        case 'ss_model'
            % 合成状态空间模型: 文件直接导出 (Af,Bf,Cf)/(Aw,Bw,Cw) 矩阵
            % 构建 ss 对象以便统一操作
            if isfield(data, 'Aw')
                % 含干扰模型的系统 (如 syn_whitenoise)
                % 构建增广系统: 控制通道Bf + 干扰通道Bw, 输出[Cf, Cw]
                n_plant = size(data.Af, 1);
                n_dist  = size(data.Aw, 1);
                info.ss_plant = ss(data.Af, data.Bf, data.Cf, 0, -1);
                info.ss_dist  = ss(data.Aw, data.Bw, data.Cw, 0, -1);
                info.fs = 1000;  % 合成模型默认采样频率
                fprintf('  合成模型 (含干扰): plant=%d阶, dist=%d阶\n', n_plant, n_dist);
                % 同时保留原始矩阵以便向后兼容
                info.Af = data.Af; info.Bf = data.Bf; info.Cf = data.Cf;
                info.Aw = data.Aw; info.Bw = data.Bw; info.Cw = data.Cw;
            else
                % 仅被控对象的系统 (如 syn_bpf)
                info.ss = ss(data.Af, data.Bf, data.Cf, 0, -1);
                info.fs = 1000;
                fprintf('  合成模型: %d阶\n', size(data.Af, 1));
                info.Af = data.Af; info.Bf = data.Bf; info.Cf = data.Cf;
            end

        case 'lms_data'
            % LMS系统辨识数据 (4通道)
            info.pri_err = data.lms_pri_err;
            info.pri_ref = data.lms_pri_ref;
            info.sec_err = data.lms_sec_err;
            info.sec_ref = data.lms_sec_ref;
            fprintf('  LMS辨识数据: pri(误差/参考) + sec(误差/参考)\n');

        case 'raw_signal'
            % 原始信号数据
            info.signal.x  = data.x;
            info.signal.y  = data.y;
            info.signal.t  = data.t;
            info.fs        = data.fs;
            fprintf('  原始信号: %.0f Hz, %.2f秒, %d采样点\n', ...
                info.fs, length(data.x)/info.fs, length(data.x));

        otherwise
            error('未知数据类型: %s', src.type);
    end

    fprintf('  加载完成。\n');
end
