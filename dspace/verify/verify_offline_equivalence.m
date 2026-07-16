function results = verify_offline_equivalence(caseFilter)
%VERIFY_OFFLINE_EQUIVALENCE 脚本仿真 vs Simulink 离线模型逐样本对拍
%
%   results = verify_offline_equivalence()            % 全部冻结案例
%   results = verify_offline_equivalence('demo1')     % 按前缀过滤
%
%   对每个 (demo, test, variant):
%     1. 用冻结参数重跑 tests/internal/controller_demoN（真值）
%     2. 重新生成离线模型并 sim（同一 d 向量注入）
%     3. 断言 max|Δe|, max|Δu|:
%          fixed ≤ 1e-12（期望恰好 0）, adaptive ≤ 1e-9
%        并用 compute_metrics 复算 supp_db，要求 |Δ| < 0.01 dB
%     4. 写 dspace/output/equivalence_results.csv
%
%   门禁语义: 任何超差是移植 bug，不放宽容差。
%
%   demo3 特例: 脚本引擎从 k_start 才开始记录 y（之前全 0 是记录习惯，
%   不是物理行为），Simulink 中 e(k)=d(k)。对拍从 k_start 起算，
%   稳态窗口 (≥0.5s) 的 supp_db 不受影响。

if nargin < 1, caseFilter = ''; end

if ~license('test', 'Simulink')
    error('verify_offline_equivalence:NoSimulink', ...
        '本机无 Simulink 许可证。请在有许可证的机器上运行 Phase B 验证。');
end

dspaceRoot = fileparts(fileparts(mfilename('fullpath')));
outputDir = fullfile(dspaceRoot, 'output');
if ~isfolder(outputDir), mkdir(outputDir); end

meta = load_anc_params('meta');
matFile = fullfile(dspaceRoot, 'params', 'anc_frozen_params.mat');
payload = load(matFile, 'frozen');
keys = setdiff(fieldnames(payload.frozen), {'meta', 'plant__secondary'});
if ~isempty(caseFilter)
    keys = keys(startsWith(keys, caseFilter));
end
if isempty(keys)
    error('verify_offline_equivalence:NoCases', '过滤 "%s" 后无案例。', caseFilter);
end

signals = load_cylinder1dm_signals('evaluation');

fprintf('\n===== Simulink 离线等效性验证 (%d 案例, git %s) =====\n\n', ...
    numel(keys), meta.git_describe);

rows = cell(numel(keys), 9);
for index = 1:numel(keys)
    key = keys{index};
    parts = strsplit(key, '__');
    [demo, test, variant] = deal(parts{1}, parts{2}, parts{3});
    entry = payload.frozen.(key);
    test_sig = signals.(test);

    % ---- 1. 脚本真值 ----
    params = jsondecode(entry.params_json);
    runner = str2func(['controller_' demo]);
    result = runner(signals, test, variant, params);
    e_script = result.extra.y_closed(:).';
    u_script = result.extra.u(:).';

    % ---- 2. Simulink 仿真 ----
    build_offline_model(demo, test, variant);
    mdl = sprintf('offline_%s_%s_%s', demo, test, variant);
    t = test_sig.t(:);
    d = test_sig.d(:);
    assignin('base', 'anc_d_in', [t, d]);
    if controller_registry(demo, variant).has_ref
        assignin('base', 'anc_ref_in', [t, d]);   % 脚本以 d 本身为参考
    end
    out = sim(mdl, 'StopTime', num2str(t(end), '%.17g'));
    e_slx = out.e_log(:).';
    u_slx = out.u_log(:).';

    N = numel(e_script);
    if numel(e_slx) ~= N
        error('verify_offline_equivalence:LengthMismatch', ...
            '%s: 脚本 %d 样本 vs Simulink %d 样本。', key, N, numel(e_slx));
    end

    % ---- 3. 对拍 ----
    if strcmp(demo, 'demo3')
        k0 = entry.k_start;    % 见头注: 脚本前导零是记录习惯
    else
        k0 = 1;
    end
    max_de = max(abs(e_script(k0:end) - e_slx(k0:end)));
    max_du = max(abs(u_script(k0:end) - u_slx(k0:end)));

    if strcmp(variant, 'fixed')
        tol = 1e-12;
    else
        tol = 1e-9;
    end

    metricsMeta = struct('demo', demo, 'variant', variant, 'test', test);
    m_slx = compute_metrics(test_sig.y_open(:).', e_slx, u_slx, ...
        test_sig, metricsMeta);
    d_supp = abs(m_slx.supp_db - result.supp_db);

    pass = (max_de <= tol) && (max_du <= tol) && (d_supp < 0.01);

    rows(index, :) = {key, variant, max_de, max_du, tol, ...
        result.supp_db, m_slx.supp_db, d_supp, pass};
    fprintf('  %-28s max|Δe|=%.3e max|Δu|=%.3e (tol %.0e) Δsupp=%.5f dB %s\n', ...
        key, max_de, max_du, tol, d_supp, pick(pass, 'PASS', 'FAIL'));

    close_system(mdl, 0);
end

results = cell2table(rows, 'VariableNames', { ...
    'case', 'variant', 'max_abs_de', 'max_abs_du', 'tolerance', ...
    'supp_db_script', 'supp_db_slx', 'supp_db_diff', 'pass'});
csvFile = fullfile(outputDir, 'equivalence_results.csv');
writetable(results, csvFile);
fprintf('\n结果 CSV: %s\n', csvFile);

write_equivalence_report(dspaceRoot, results, meta);

nFail = sum(~cell2mat(rows(:, 9)));
if nFail > 0
    error('verify_offline_equivalence:Failures', ...
        '%d/%d 案例超差 —— 按移植 bug 处理，见 %s。', ...
        nFail, numel(keys), csvFile);
end
fprintf('全部 %d 案例通过等效性门禁。\n', numel(keys));
end

function write_equivalence_report(dspaceRoot, results, meta)
%WRITE_EQUIVALENCE_REPORT 生成 dspace/SIMULINK_EQUIVALENCE_REPORT.md
reportFile = fullfile(dspaceRoot, 'SIMULINK_EQUIVALENCE_REPORT.md');
fid = fopen(reportFile, 'w');
cleaner = onCleanup(@() fclose(fid));

fprintf(fid, '# Simulink 离线等效性报告\n\n');
fprintf(fid, ['**目的**: 证明 dspace/controllers/ 的 MATLAB Function 移植与 ' ...
    'tests/internal/ 脚本引擎在同一冻结参数、同一评价扰动 d（seed %d）下' ...
    '逐样本等效，作为 RT 部署前的硬门禁。\n\n'], meta.evaluation_seed);
fprintf(fid, ['**设置**: 离线闭环模型（被控对象 = cylinder1dm 实测 FIR 路径，' ...
    '与脚本引擎算术逐位一致的库块），From Workspace 注入与脚本完全同一的 ' ...
    'd 向量；FixedStepDiscrete Ts=%.6g s。\n\n'], meta.Ts);
fprintf(fid, ['**判定**: fixed 容差 max|Δ| ≤ 1e-12（期望恰好 0）；' ...
    'adaptive ≤ 1e-9；supp_db 复算差 < 0.01 dB。超差 = 移植 bug，不放宽。\n' ...
    'demo3 从 k_start 起对拍（脚本前导零是记录习惯，非物理行为，' ...
    '见 verify_offline_equivalence.m 头注）。\n\n']);
fprintf(fid, '**溯源**: git %s，参数导出 %s，报告生成 %s。\n\n', ...
    meta.git_describe, meta.exported_at, ...
    char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));

fprintf(fid, '| 案例 | 变体 | max\\|Δe\\| | max\\|Δu\\| | 容差 | supp 脚本 (dB) | supp Simulink (dB) | Δsupp (dB) | 结论 |\n');
fprintf(fid, '|---|---|---:|---:|---:|---:|---:|---:|---|\n');
for index = 1:height(results)
    fprintf(fid, '| %s | %s | %.3e | %.3e | %.0e | %.2f | %.2f | %.5f | %s |\n', ...
        results.case{index}, results.variant{index}, ...
        results.max_abs_de(index), results.max_abs_du(index), ...
        results.tolerance(index), results.supp_db_script(index), ...
        results.supp_db_slx(index), results.supp_db_diff(index), ...
        pick(results.pass(index), '通过', '**失败**'));
end

fprintf(fid, '\n**结论**: %d/%d 案例通过。%s\n', ...
    sum(results.pass), height(results), ...
    pick(all(results.pass), ...
    '全部移植控制器与脚本引擎等效，可进入 RT 模型阶段。', ...
    '存在失败案例 —— 修复移植 bug 后重跑本验证，不得带差异上板。'));
fprintf('等效性报告: %s\n', reportFile);
end

function out = pick(cond, a, b)
if cond, out = a; else, out = b; end
end
