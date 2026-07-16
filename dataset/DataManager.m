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
%     .type    - 数据类型: 'armax_model' | 'fir_model' | 'ss_model' | 'lms_data' | 'raw_signal'
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
%   'syn_Ho2020_ALE'      - Ho 2020 窄带ALE前馈ANC (P:7阶, S:4阶, fs=4k)
%   'syn_RSTtoy_2nd'      - RST教学模型           (2阶, Ts=1)
%   'syn_whitenoise'      - 合成带限白噪声干扰模型 (状态空间)
%   'syn_bpf'             - 合成带通滤波器被控对象模型 (状态空间)
%   'lms_sysid'           - LMS系统辨识实验数据 (4通道, 2026-01-20)
%   'cylinder1dm_2k_primary_fir'   - Cylinder 1 dm 主路径 FIR(16) @2kHz
%   'cylinder1dm_2k_secondary_fir' - Cylinder 1 dm 次级路径 FIR(16) @2kHz
%   'raw_dspace'          - dSPACE采集原始信号数据

    persistent loadedData;
    if isempty(loadedData)
        loadedData = struct();
    end

    % 数据集根目录
    dsDir = fileparts(mfilename('fullpath'));

    % 从统一模型注册表派生 DataManager 可加载数据源。
    registered = model_registry('loadable');
    sources = struct();
    datasetPrefix = ['dataset' filesep];
    projectDir = fileparts(dsDir);
    for i = 1:numel(registered)
        artifact = strrep(registered(i).artifact, '/', filesep);
        if strncmp(artifact, datasetPrefix, length(datasetPrefix))
            artifact = artifact(length(datasetPrefix)+1:end);
            baseDir = dsDir;
        else
            baseDir = projectDir;
        end
        sources.(registered(i).loader_id) = struct( ...
            'file', artifact, 'type', registered(i).loader_type, ...
            'base_dir', baseDir, 'registry', registered(i));
    end

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
    fullPath = fullfile(src.base_dir, src.file);
    if ~exist(fullPath, 'file')
        error('数据文件不存在: %s', fullPath);
    end

    cacheKey = matlab.lang.makeValidName(dataSource);
    if isfield(loadedData, cacheKey)
        data = loadedData.(cacheKey);
    else
        data = load(fullPath);
        loadedData.(cacheKey) = data;
    end
    info.name = dataSource;
    info.type = src.type;

    switch src.type
        case 'armax_model'
            % ARMAX辨识模型
            if ~isfield(data, 'ARMAXmodel') || ~isfield(data.ARMAXmodel, 'model')
                error('DataManager:InvalidARMAXModel', ...
                    '数据文件 %s 缺少 ARMAXmodel.model。', src.file);
            end
            info.model  = data.ARMAXmodel.model;   % idpoly 对象
            info.orders = data.ARMAXmodel.orders;  % [na nb nc nk]
            info.fs     = data.ARMAXmodel.fs;
            info.registry = src.registry;
            info.artifact = fullPath;

            if isempty(info.model)
                warning('DataManager:EmptyARMAXModel', ...
                    ['ARMAXmodel.model 为空；当前文件只有阶次和采样率，' ...
                     '不能重新计算零极点或设计控制器。']);
            end

            fprintf('  ARMAX(%d,%d,%d,%d), fs=%d Hz\n', ...
                info.orders(1), info.orders(2), info.orders(3), info.orders(4), info.fs);

        case 'fir_model'
            % LMS FIR artifacts are stored as struct s.w/s.fs/s.mu.
            if ~isfield(data, 's') || ~isfield(data.s, 'w') || ...
                    ~isfield(data.s, 'fs')
                error('DataManager:InvalidFIRModel', ...
                    'FIR 文件 %s 缺少 s.w 或 s.fs。', src.file);
            end
            w = double(data.s.w(:).');
            info.model = struct('A', 1, 'B', w);
            info.orders = [0, numel(w) - 1, 0, 0];
            info.fs = double(data.s.fs);
            info.mu = iff_field(data.s, 'mu', NaN);
            info.registry = src.registry;
            info.artifact = fullPath;
            fprintf('  FIR(%d), fs=%d Hz, mu=%.3g\n', ...
                numel(w), info.fs, info.mu);

        case 'syn_tf'
            % 文献/教学 LTI 模型可能保存为 tf、zpk 或 ss。
            if ~isfield(data, 'model')
                error('DataManager:InvalidSyntheticModel', ...
                    '数据文件 %s 缺少 model struct。', src.file);
            end
            model = data.model;
            if isfield(model, 'G0')
                info.G0 = model.G0;
            elseif isfield(model, 'G0_zpk')
                info.G0 = model.G0_zpk;
            elseif isfield(model, 'G0_tf')
                info.G0 = model.G0_tf;
            elseif isfield(model, 'sys')
                info.G0 = model.sys;
                info.sys = model.sys;
            else
                error('DataManager:MissingLTIRepresentation', ...
                    '模型 %s 不含 G0、G0_zpk、G0_tf 或 sys。', dataSource);
            end
            info.domain = model.domain;
            info.orders = model.orders;
            info.source = model.source;
            info.desc   = model.desc;
            if isfield(model, 'Ts'),      info.Ts = model.Ts; end
            if isfield(model, 'fs'),      info.fs = model.fs; end
            if isfield(model, 'fs_nominal'), info.fs = model.fs_nominal; end
            if isfield(model, 'G0_tf'),   info.G0_tf = model.G0_tf; end
            if isfield(model, 'G0_zpk'),  info.G0_zpk = model.G0_zpk; end
            if isfield(model, 'A'), info.A = model.A; end
            if isfield(model, 'B'), info.B = model.B; end
            if isfield(model, 'C'), info.C = model.C; end
            if isfield(model, 'D'), info.D = model.D; end
            if isfield(model, 'A_coeffs'), info.A_coeffs = model.A_coeffs; end
            if isfield(model, 'B_coeffs'), info.B_coeffs = model.B_coeffs; end
            if isfield(model, 'A_poly'), info.A_poly = model.A_poly; end
            if isfield(model, 'B_poly'), info.B_poly = model.B_poly; end
            if isfield(model, 'd_delay'), info.delay_samples = model.d_delay; end

            fprintf('  %s (%s, %s)\n', model.name, model.domain, model.source);

        case 'syn_tf_feedforward'
            % 前馈 ANC 合成模型 (含 P(z) + S(z) 双路径)
            if isfield(data, 'model')
                model = data.model;
                info.P       = model.P;
                info.S       = model.S;
                info.P_desc  = model.P_desc;
                info.S_desc  = model.S_desc;
                info.domain  = model.domain;
                info.orders  = model.orders;
                info.source  = model.source;
                info.desc    = model.desc;
                if isfield(model, 'Ts'), info.Ts = model.Ts; end
                if isfield(model, 'fs'), info.fs = model.fs; end
                fprintf('  %s (P:%s, S:%s)\n', model.name, model.P_desc, model.S_desc);
            elseif all(isfield(data, {'P_num', 'P_den', 'S_num', 'S_den'}))
                % 兼容仅保存 FIR 系数的轻量 MAT 版本。
                info.fs = sscanf(src.registry.sample_rate, '%f', 1);
                info.Ts = 1 / info.fs;
                info.P_num = data.P_num; info.P_den = data.P_den;
                info.S_num = data.S_num; info.S_den = data.S_den;
                info.P = tf(data.P_num, data.P_den, info.Ts, 'Variable', 'z^-1');
                info.S = tf(data.S_num, data.S_den, info.Ts, 'Variable', 'z^-1');
                info.P_desc = 'Primary-path FIR coefficients';
                info.S_desc = 'Secondary-path FIR coefficients';
                info.domain = src.registry.domain;
                info.orders = [numel(data.P_num)-1, numel(data.S_num)-1];
                info.source = src.registry.source;
                info.desc = src.registry.notes;
                fprintf('  %s (coefficient-only P/S MAT)\n', dataSource);
            else
                error('DataManager:InvalidFeedforwardModel', ...
                    '模型 %s 不含 model struct 或 P/S 系数。', dataSource);
            end

        case 'ss_model'
            % 合成状态空间模型: 文件直接导出 (Af,Bf,Cf)/(Aw,Bw,Cw) 矩阵
            % 构建 ss 对象以便统一操作
            if isfield(data, 'Aw')
                % 干扰状态空间模型；部分历史文件不含被控对象矩阵。
                n_dist  = size(data.Aw, 1);
                info.ss_dist  = ss(data.Aw, data.Bw, data.Cw, 0, -1);
                info.fs = 1000;  % 合成模型默认采样频率
                info.Aw = data.Aw; info.Bw = data.Bw; info.Cw = data.Cw;
                if isfield(data, 'Af')
                    n_plant = size(data.Af, 1);
                    info.ss_plant = ss(data.Af, data.Bf, data.Cf, 0, -1);
                    info.Af = data.Af; info.Bf = data.Bf; info.Cf = data.Cf;
                    fprintf('  合成模型: plant=%d阶, dist=%d阶\n', n_plant, n_dist);
                else
                    fprintf('  合成干扰模型: dist=%d阶\n', n_dist);
                end
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

function value = iff_field(s, fieldName, fallback)
if isfield(s, fieldName)
    value = s.(fieldName);
else
    value = fallback;
end
end
