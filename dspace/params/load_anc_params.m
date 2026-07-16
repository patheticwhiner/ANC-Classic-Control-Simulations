function entry = load_anc_params(demo, test, variant)
%LOAD_ANC_PARAMS 读取一个 (demo, test, variant) 的冻结部署参数
%
%   entry = load_anc_params('demo1', 'T1', 'fixed')
%   entry = load_anc_params('plant', 'secondary')     % 被控对象次级路径
%   entry = load_anc_params('meta')                   % 溯源元数据
%
%   数据来自 dspace/params/anc_frozen_params.mat（由 export_frozen_params
%   生成，*.mat 不入 git —— 缺失时报错并给出再生命令）。

persistent cached cachedMtime;

matFile = fullfile(fileparts(mfilename('fullpath')), 'anc_frozen_params.mat');
if ~isfile(matFile)
    error('load_anc_params:MissingMat', ...
        ['未找到 %s\n请先运行: run(''dspace/startup.m''); ' ...
         'export_frozen_params'], matFile);
end

info = dir(matFile);
if isempty(cached) || isempty(cachedMtime) || cachedMtime ~= info.datenum
    payload = load(matFile, 'frozen');
    cached = payload.frozen;
    cachedMtime = info.datenum;
end

if nargin == 1
    key = demo;                                  % 'meta'
elseif nargin == 2
    key = sprintf('%s__%s', demo, test);         % 'plant__secondary'
else
    key = sprintf('%s__%s__%s', demo, test, variant);
end

if ~isfield(cached, key)
    error('load_anc_params:UnknownKey', ...
        '冻结参数中没有键 %s。可用键: %s', key, ...
        strjoin(fieldnames(cached), ', '));
end
entry = cached.(key);
end
