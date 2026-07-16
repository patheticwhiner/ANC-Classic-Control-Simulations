function reportPath = generate_model_report(reportPath)
%GENERATE_MODEL_REPORT Build the repository-wide model provenance report.
%   GENERATE_MODEL_REPORT() writes MODEL_ANALYSIS_REPORT.md at repo root.
%   The scan is static: it does not require Control System Toolbox or load
%   model objects, so missing/unsynchronised artifacts remain reportable.

toolDir = fileparts(mfilename('fullpath'));
rootDir = fileparts(toolDir);
addpath(fullfile(rootDir, 'dataset'));

if nargin < 1 || isempty(reportPath)
    reportPath = fullfile(rootDir, 'MODEL_ANALYSIS_REPORT.md');
end

registry = model_registry();
loadable = registry(~cellfun(@isempty, {registry.loader_id}));
inlineRecords = registry(cellfun(@isempty, {registry.loader_id}));
dynamicError = '';
try
    evalc('dynamics = analyze_model_dynamics();');
catch exception
    dynamics = struct([]);
    dynamicError = exception.message;
end

[missingArtifacts, unregisteredDataset] = audit_dataset(rootDir, registry);
sysid = scan_sysid_models(rootDir);
embeddedCandidates = scan_embedded_models(rootDir, inlineRecords);
contentIssues = model_content_issues(dynamics);

fid = fopen(reportPath, 'w', 'n', 'UTF-8');
if fid < 0
    error('generate_model_report:CannotOpen', ...
        'Cannot write model report: %s', reportPath);
end
cleanup = onCleanup(@() fclose(fid));

writeln(fid, '# 模型来源与统一分析报告');
writeln(fid, '');
generatedDate = datestr(now, 'yyyy-mm-dd'); %#ok<DATST,TNOW1>
writeln(fid, sprintf('> 生成日期: %s  |  生成器: `tools/generate_model_report.m`  |  注册表: `dataset/model_registry.m`', ...
    generatedDate));
writeln(fid, '>');
writeln(fid, '> 本报告是模型资产、来源链路和适用范围的统一入口。公式级说明仍保留在 `dataset/ModelReference.md`；辨识过程细节以各 MAT 文件同名 Markdown 报告为准。');
writeln(fid, '');

write_model_map(fid);
write_dynamics_table(fid, dynamics, dynamicError);
write_model_portraits(fid);
write_comparability(fid);

writeln(fid, '## 五、实验选型');
writeln(fid, '');
writeln(fid, '- 跨 demo 横向对比暂时沿用 `armax_30303022`，因为当前测试基准已围绕它建立；但当前 MAT 的模型对象为空，修复前只能复用既有信号和分析快照，不能重新设计控制器。');
writeln(fid, '- 论文方法复现应使用对应的 `syn_*` 文献基准，只能证明算法复现一致性，不能直接外推到真实声学平台。');
writeln(fid, '- `sysid_models/` 中的新模型按采集批次动态盘点。它们适合候选模型筛选，但在残差白度、采集来源或验证报告不完整时，不应替换统一基准。');
writeln(fid, '- demo5 仍有脚本内嵌对象。它们已登记来源与风险，但在抽取成独立模型并接入测试基准前，不宜与 `dataset/` 模型混为一类。');
writeln(fid, '');

writeln(fid, '### 建议使用矩阵');
writeln(fid, '');
writeln(fid, '| 任务 | 首选模型 | 理由 | 限制 |');
writeln(fid, '|:---|:---|:---|:---|');
writeln(fid, '| 四类控制器统一比较 | `armax_30303022` | tests 已统一信号、指标和加载入口 | 48 kHz、高阶、Bezout 数值病态 |');
writeln(fid, '| RST/Youla 论文复现 | `syn_Carmona2000_7th` | 低阶管道 ANC 文献基准 | 不能代表当前实验平台 |');
writeln(fid, '| Jafari 自适应/鲁棒复现 | `syn_TAC2015_3rd`, `syn_JVC2017_*` | 与推导和论文场景一致 | 合成理论模型 |');
writeln(fid, '| H-infinity 复现 | `syn_Bai1997_4th` | 不稳定、非最小相位参考对象 | 耳机模型，不是管道模型 |');
writeln(fid, '| 前馈 ALE/FXLMS 结构验证 | `syn_Ho2020_ALE` | 明确的 P(z)+S(z) 路径对 | 合成等价模型，尚无实测验证 |');
writeln(fid, '| 最新硬件路径评估 | `sysid_models/<批次>/` | 保留采集批次和同名验证报告 | 需逐批审查，不能自动视为 canonical |');
writeln(fid, '');

writeln(fid, '## 附录 A：资产概览');
writeln(fid, '');
writeln(fid, sprintf('- 注册记录: **%d**（DataManager 可加载 %d，脚本内嵌/来源记录 %d）', ...
    numel(registry), numel(loadable), numel(inlineRecords)));
writeln(fid, sprintf('- 本机 `sysid_models/` MAT 文件: **%d**（有可关联报告 %d，无可关联报告 %d）', ...
    numel(sysid), sum([sysid.has_report]), sum(~[sysid.has_report])));
writeln(fid, sprintf('- 注册路径缺失: **%d**；`dataset/` 未登记 MAT 文件: **%d**', ...
    numel(missingArtifacts), numel(unregisteredDataset)));
writeln(fid, sprintf('- 注册资产内容问题: **%d**', numel(contentIssues)));
writeln(fid, sprintf('- 静态审计发现的其他脚本内嵌模型候选: **%d**', numel(embeddedCandidates)));
writeln(fid, '');

writeln(fid, '### 来源分层');
writeln(fid, '');
writeln(fid, '| 来源层级 | 含义 | 可支持的结论 |');
writeln(fid, '|:---|:---|:---|');
writeln(fid, '| `measured` / `measured-identification` | 来自硬件采集及其辨识结果 | 对指定装置、采样率和采集批次有效 |');
writeln(fid, '| `literature` / `literature-reproduction` | 论文给定或按论文重建 | 支持方法复现，不自动支持工程泛化 |');
writeln(fid, '| `synthetic` / `repository-synthetic` | 本工程构造的测试对象或扰动 | 支持功能与边界测试 |');
writeln(fid, '| `educational` | 教材或低阶教学对象 | 支持概念演示，不支持 ANC 性能声明 |');
writeln(fid, '');

writeln(fid, '## 附录 B：统一注册表');
writeln(fid, '');
writeln(fid, '| ID | 角色 | 来源 | 表示 | 域/采样 | 成熟度 | 主资产 |');
writeln(fid, '|:---|:---|:---|:---|:---|:---|:---|');
for i = 1:numel(registry)
    item = registry(i);
    domainRate = sprintf('%s / %s', item.domain, item.sample_rate);
    writeln(fid, sprintf('| `%s` | %s | `%s` | %s | %s | `%s` | `%s` |', ...
        md(item.id), md(item.role), md(item.origin), md(item.representation), ...
        md(domainRate), md(item.maturity), md(item.artifact)));
end
writeln(fid, '');

writeln(fid, '### 来源与证据');
writeln(fid, '');
for i = 1:numel(registry)
    item = registry(i);
    writeln(fid, sprintf('- **`%s`**: %s  ', md(item.id), md(item.source)));
    writeln(fid, sprintf('  证据: `%s`；使用方: %s；备注: %s', ...
        md(item.evidence), md(item.consumers), md(item.notes)));
end
writeln(fid, '');

writeln(fid, '## 附录 C：系统辨识模型动态清单');
writeln(fid, '');
writeln(fid, '`sysid_models/` 默认不纳入 Git，因而不同机器的数量可能不同。下表由文件名与同名 Markdown 报告自动推导，不把“文件存在”等同于“模型已验证”。');
writeln(fid, '');
if isempty(sysid)
    writeln(fid, '> 本机没有同步 `sysid_models/` MAT 文件。');
else
    writeln(fid, '| 资产 | 类型 | 日期 | 路径角色 | 结构/阶次 | 采样 | 来源判断 | 报告 | 状态 |');
    writeln(fid, '|:---|:---|:---|:---|:---|:---|:---|:---:|:---|');
    for i = 1:numel(sysid)
        item = sysid(i);
        writeln(fid, sprintf('| `%s` | %s | %s | %s | %s | %s | %s | %s | %s |', ...
            md(item.path), md(item.type), md(item.date), md(item.path_role), ...
            md(item.order), md(item.sample_rate), md(item.origin), ...
            md(item.report_kind), md(item.status)));
    end
end
writeln(fid, '');

writeln(fid, '## 附录 D：完整性审计');
writeln(fid, '');
write_audit_list(fid, '注册资产缺失', missingArtifacts, ...
    '注册表中的主资产在本机不存在；DataManager 对应入口会失败。');
write_audit_list(fid, '`dataset/` 未登记 MAT 文件', unregisteredDataset, ...
    '文件存在但没有进入统一注册表。');
write_audit_list(fid, '注册模型内容问题', contentIssues, ...
    '路径存在，但模型内容不足以重新计算动力学或复现实验。');

undocumented = {sysid(~[sysid.has_report]).path};
write_audit_list(fid, '`sysid_models/` 缺少可关联报告', undocumented, ...
    '只能从文件名推断来源和结构，无法审查采集与验证条件。');
write_audit_list(fid, '其他脚本内嵌模型候选', embeddedCandidates, ...
    '脚本直接构造 tf/ss/zpk，且没有加载注册模型；需人工判断是对象、控制器还是辅助滤波器。');

writeln(fid, '## 附录 E：治理规则');
writeln(fid, '');
writeln(fid, '1. 新的长期使用模型先登记到 `dataset/model_registry.m`，再由 `DataManager` 或明确的加载函数暴露。');
writeln(fid, '2. 实测辨识模型必须保留采集标识、路径角色、原始/辨识采样率、模型结构、验证指标和同名报告；缺一项即保持 candidate 状态。');
writeln(fid, '3. 论文模型必须记录完整文献来源，并区分“论文直接给定”“按图表重建”“本工程为对比而构造”。');
writeln(fid, '4. 同一物理路径的 FIR、ARMAX、OE、SS 是不同表示，不是独立实验来源；报告时必须共享同一采集 lineage。');
writeln(fid, '5. 跨算法性能表只允许使用同一物理对象、采样率、测试信号、起止窗口和指标定义。');
writeln(fid, '6. 每次新增、删除或移动模型后运行 `generate_model_report`，并审查“完整性审计”是否出现新条目。');
writeln(fid, '');

writeln(fid, '## 附录 F：复现命令');
writeln(fid, '');
writeln(fid, '```matlab');
writeln(fid, "addpath('dataset', 'tools');");
writeln(fid, 'DataManager()               % 查看统一可加载入口');
writeln(fid, 'generate_model_report()     % 重新生成本报告');
writeln(fid, '```');

fprintf('Model report written: %s\n', reportPath);
end

function write_model_map(fid)
writeln(fid, '## 一、先建立模型地图');
writeln(fid, '');
writeln(fid, '这些文件并不是一组可以直接按“抑制 dB”排序的同类对象。它们分别用来隔离不同控制难题：');
writeln(fid, '');
writeln(fid, '| 模型族 | 代表模型 | 它刻意暴露的问题 | 它不能证明什么 |');
writeln(fid, '|:---|:---|:---|:---|');
writeln(fid, '| 实测反馈声学路径 | `armax_30303022` | 高阶、长延迟、NMP 零点、辨识不确定性同时存在 | 不能单独归因某一种算法机制 |');
writeln(fid, '| 管道 RST 文献对象 | `syn_Carmona2000_7th` | 离散管道、明显延迟、非最小相位 | 不能代表当前 48 kHz 实验平台 |');
writeln(fid, '| Jafari 理论对象 | `syn_TAC2015_3rd`, `syn_JVC2017_*` | 低增益、边界/RHP 零点、频谱展平和自适应律 | 不能给出真实声学工程性能 |');
writeln(fid, '| 不稳定鲁棒对象 | `syn_Bai1997_4th` | 开环镇定与性能塑形必须同时完成 | 不能与稳定模型只比抑制量 |');
writeln(fid, '| 前馈双路径对象 | `syn_Ho2020_ALE` | P/S 因果关系、次级通路 NMP、filtered-x 结构 | 不能直接验证反馈 ANC 控制器 |');
writeln(fid, '| 教学对象 | `syn_RSTtoy_2nd`, `syn_MassSpringDamper_2nd` | 验证方程、接口和优化器是否基本正确 | 不能支撑 ANC 有效性结论 |');
writeln(fid, '| 扰动/滤波模型 | `syn_whitenoise`, `syn_bpf` | 生成测试信号或早期 LQG 场景 | 不是与声学路径并列的被控对象 |');
writeln(fid, '| 非线性输出调节 | demo5 脚本内嵌模型 | 未知频率、外系统、非线性参数自适应 | 不能用同一套 LTI 零极点指标概括 |');
writeln(fid, '');
writeln(fid, '**最重要的分类轴不是算法名称，而是：单路径/双路径、反馈/前馈、稳定/不稳定、最小/非最小相位、解析模型/辨识模型。** 只有这些轴一致，性能数值才有横向意义。');
writeln(fid, '');
end

function write_dynamics_table(fid, dynamics, dynamicError)
writeln(fid, '## 二、动力学画像');
writeln(fid, '');
writeln(fid, '下表直接从当前模型计算；标为“分析快照”的行来自既有分析文档。连续模型频率按 `f = omega/(2*pi)` 换算，仅表示模型内部时间尺度，不能当成声学 Hz 与离散模型横比。');
writeln(fid, '');
if isempty(dynamics)
    writeln(fid, sprintf('> 动力学提取失败：%s', md(dynamicError)));
    writeln(fid, '');
    return;
end
writeln(fid, '| 模型/路径 | 域 | 极点与稳定性 | 零点约束 | 延迟 | 扫描内最大增益 | 分析依据 |');
writeln(fid, '|:---|:---|:---|:---|:---|:---|:---|');
for i = 1:numel(dynamics)
    item = dynamics(i);
    writeln(fid, sprintf('| `%s` | %s | %s | %s | %s | %s | %s |', ...
        md(item.id), md(format_domain(item)), md(format_poles(item)), ...
        md(format_zeros(item)), md(format_delay(item)), ...
        md(format_peak(item)), md(format_basis(item.analysis_basis))));
end
writeln(fid, '');
writeln(fid, '读表时注意：');
writeln(fid, '');
writeln(fid, '- “稳定”只说明自由响应不发散，不代表容易控制。ARMAX 与 Carmona 都稳定，但延迟和 NMP 零点仍限制可达带宽与稳定逆。');
writeln(fid, '- `z=1` 或虚轴零点列为“边界零点”，与严格位于单位圆外/RHP 的 NMP 零点分开。');
writeln(fid, '- FIR 路径的原点极点主要表示有限记忆，不应和接近单位圆的轻阻尼动态等量看待；纯延迟已单独列出。');
writeln(fid, '- “扫描内最大增益”用于认识增益尺度，不自动等于共振峰；例如 Carmona 的最大值在 DC。');
writeln(fid, '');
end

function write_model_portraits(fid)
writeln(fid, '## 三、逐模型认识');
writeln(fid, '');
writeln(fid, '### 3.1 `armax_30303022`：真实但不可重新计算的当前基准');
writeln(fid, '');
writeln(fid, '- **怎么看它**：这是稳定的实测声学路径，不需要先解决开环镇定；主要动态集中在约 334.6 Hz，最大极点模约 0.9958，说明共振衰减慢、模型记忆长。');
writeln(fid, '- **真正障碍**：30 阶、9 个 NMP 零点、22 点延迟。22 点在 48 kHz 下是 0.458 ms，经验反馈带宽上限约 727 Hz；全阶 Bezout 方程还会形成约 `1e16` 的病态矩阵。');
writeln(fid, '- **研究含义**：它适合检验控制器能否同时承受延迟、NMP 和高阶辨识误差，不适合用来单独证明某个理论环节。');
writeln(fid, '- **数据警告**：当前 MAT 中 `ARMAXmodel.model` 是空数组。报告里的 30/29 阶、9 个 NMP 零点和 334.6 Hz 来自 `About_Plant_Performance.md` 的既有快照；必须恢复 A/B/C 多项式后才能重新设计和独立复核。');
writeln(fid, '');
writeln(fid, '### 3.2 `syn_Carmona2000_7th`：点数不大，但物理延迟更大');
writeln(fid, '');
writeln(fid, '- 7 个动态极点、6 点纯延迟，采样率 2 kHz。6 点只看数字似乎小于 ARMAX 的 22 点，但对应 **3 ms**，约是 0.458 ms 的 6.5 倍。');
writeln(fid, '- 当前系数计算得到 13 个单位圆外零点，说明稳定逆不是可用策略；RST/Youla 必须做闭环塑形，而不是直接消去对象。');
writeln(fid, '- 它适合研究低阶多项式设计和延迟影响；不应因为同为“管道模型”就替代 48 kHz 实测对象。');
writeln(fid, '');
writeln(fid, '### 3.3 Jafari 系列：把单一理论障碍逐步放大');
writeln(fid, '');
writeln(fid, '- `syn_TAC2015_3rd` 全部极点稳定，`z=1` 是边界零点而非严格 NMP；DC 增益为零，整个频带增益也很低。控制难点是需要很大的补偿增益，同时避免在边界零点附近作不可靠求逆。');
writeln(fid, '- `syn_JVC2017_3rd` 只有一对稳定极点和一个 RHP 零点 `s=0.2`，是理解“为什么要反射 NMP 零点、构造 F(s)”的最干净入口。');
writeln(fid, '- `syn_JVC2017_6th` 保留同一 RHP 零点问题，同时加入多重极点、原点零点和虚轴零点，用于检查频谱展平在高阶/边界零点下是否仍成立。');
writeln(fid, '- 三者是机制模型：先在这里理解算法为什么有效，再到 ARMAX 检查数值和工程边界。');
writeln(fid, '');
writeln(fid, '### 3.4 `syn_Bai1997_4th`：唯一以“先镇定”为首要任务的对象');
writeln(fid, '');
writeln(fid, '- 4 个极点中有 2 个 RHP 极点，4 个零点中有 2 个 RHP 零点。开环不稳定和非最小相位同时存在。');
writeln(fid, '- 对它的第一项评价是闭环内部稳定性，其次才是灵敏度塑形和抑制；一个抑制较低但可靠镇定的控制器，不能简单判为差于稳定对象上的高抑制控制器。');
writeln(fid, '- 它是 H-infinity 权重设计和水床效应的基准，不是管道 ANC 的缩小版。');
writeln(fid, '');
writeln(fid, '### 3.5 `syn_Ho2020_ALE`：理解前馈 ANC 的因果窗口');
writeln(fid, '');
writeln(fid, '- 主通路 P 延迟 5 点、次级通路 S 延迟 2 点，4 kHz 下分别为 1.25 ms 和 0.5 ms；参考信号因此提供约 **3 点/0.75 ms** 的因果余量。');
writeln(fid, '- S 含一个单位圆外零点（最大模 2），所以 `S^-1` 不是稳定因果滤波器；filtered-x/ALE 的意义正是避免把控制问题简化成精确求逆。');
writeln(fid, '- 这是双路径前馈对象。若拿它与反馈 ARMAX 模型比较，改变的不仅是模型，还包括信息结构，因此控制器抑制量不能直接横比。');
writeln(fid, '');
writeln(fid, '### 3.6 教学与信号模型：知道它们证明到哪里为止');
writeln(fid, '');
writeln(fid, '- `syn_RSTtoy_2nd` 低阶、稳定、最小相位，只能证明 RST 方程和优化器在良态问题上可运行。');
writeln(fid, '- 质量-弹簧-阻尼器用于理解状态反馈、观测器和阻尼，不包含声传播延迟、NMP 次级通路或辨识误差。');
writeln(fid, '- `syn_whitenoise` 是扰动生成模型，`syn_bpf` 是早期 LQG 场景对象；二者 MAT 中采样时间原为 `-1`（未指定），本报告按 DataManager 的名义 1 kHz 解释。其绝对频率结论因此低于明确采样率模型的可信度。');
writeln(fid, '');
writeln(fid, '### 3.7 demo5：不是缺一张零极点表，而是模型范式不同');
writeln(fid, '');
writeln(fid, 'demo5 的核心是非线性输出调节、外系统和未知频率估计。它的“模型”包含被控对象、扰动生成器和自适应参数化，不能压缩成一个 SISO 传递函数后与前四个 demo 公平比较。正确对齐方式是统一扰动场景、可用测量和控制能量，而不是强行统一零极点阶数。');
writeln(fid, '');
end

function write_comparability(fid)
writeln(fid, '## 四、哪些结论可以比较');
writeln(fid, '');
writeln(fid, '| 想回答的问题 | 应使用的模型 | 必须保持一致 | 不应混入 |');
writeln(fid, '|:---|:---|:---|:---|');
writeln(fid, '| 哪个控制器更适合当前实测反馈路径 | 恢复后的 `armax_30303022` 或同一新采集模型 | 对象、采样率、测试信号、窗口、控制限幅 | 论文低阶对象、前馈双路径对象 |');
writeln(fid, '| 某篇论文方法是否复现正确 | 该论文对应 `syn_*` | 原始参数、扰动、初值、指标 | 其他论文模型上的更高 dB |');
writeln(fid, '| NMP 零点处理是否有效 | JVC 3rd -> JVC 6th -> ARMAX | 相同 F/K 分解和稳定性判据 | 只在 RST toy 上成功 |');
writeln(fid, '| 延迟如何限制反馈带宽 | Carmona 与 ARMAX 可作定性对照 | 比较毫秒而不是采样点数 | 直接比较阶数或抑制 dB |');
writeln(fid, '| 前馈算法是否利用参考预览 | Ho P/S 或同一实测 P/S 对 | 主/次通路、参考噪声、因果余量 | 单路径反馈模型 |');
writeln(fid, '| SYSID 表示哪一种更好 | 同一次采集的 FIR/ARMAX/OE/SS | 同一训练/验证段和采样率 | 跨日期、跨装置的 fit 数值 |');
writeln(fid, '');
writeln(fid, '### 推荐认识顺序');
writeln(fid, '');
writeln(fid, '1. 用 RST toy 和质量-弹簧模型确认基本方程、极点和状态反馈概念。');
writeln(fid, '2. 用 JVC 3rd 单独理解一个 RHP 零点为何阻止稳定精确求逆。');
writeln(fid, '3. 用 JVC 6th 与 TAC 观察边界零点、低增益和高阶频谱展平。');
writeln(fid, '4. 用 Bai 区分“镇定问题”和“性能优化问题”。');
writeln(fid, '5. 用 Carmona 理解离散延迟和多项式设计，再进入 ARMAX 的高阶与数值病态。');
writeln(fid, '6. 最后用 Ho 与 demo5 分别扩展到前馈双路径和非线性外系统，不再强求单一 LTI 排名。');
writeln(fid, '');
end

function value = format_domain(item)
if strcmp(item.domain, 'discrete')
    if item.sample_rate_hz >= 1000
        value = sprintf('离散 / %.3g kHz', item.sample_rate_hz/1000);
    else
        value = sprintf('离散 / %.3g Hz', item.sample_rate_hz);
    end
else
    value = '连续';
end
end

function value = format_poles(item)
if strcmp(item.domain, 'discrete')
    edge = sprintf('max|p|=%.4g', item.stability_edge);
else
    edge = sprintf('max Re(p)=%.4g', item.stability_edge);
end
value = sprintf('%d 极点；%s；%s', item.pole_count, ...
    translate_stability(item.stability), edge);
end

function value = format_zeros(item)
if item.nmp_zeros > 0
    classText = sprintf('%d NMP', item.nmp_zeros);
elseif item.boundary_zeros > 0
    classText = sprintf('%d 边界', item.boundary_zeros);
else
    classText = '最小相位';
end
if item.boundary_zeros > 0 && item.nmp_zeros > 0
    classText = sprintf('%s + %d 边界', classText, item.boundary_zeros);
end
value = sprintf('%d 零点；%s', item.zero_count, classText);
end

function value = format_delay(item)
if item.delay_samples == 0
    value = '0';
elseif isfinite(item.delay_seconds)
    if item.delay_seconds >= 1
        value = sprintf('%d 点 / %.3g s', ...
            item.delay_samples, item.delay_seconds);
    else
        value = sprintf('%d 点 / %.3g ms', ...
            item.delay_samples, item.delay_seconds*1000);
    end
else
    value = sprintf('%d 点', item.delay_samples);
end
end

function value = translate_stability(value)
switch value
    case 'stable'
        value = '稳定';
    case 'unstable'
        value = '不稳定';
    case 'marginal'
        value = '临界稳定';
end
end

function value = format_peak(item)
if ~isfinite(item.response_peak_db)
    value = 'not available';
elseif item.response_peak_frequency_hz < 1e-10
    value = sprintf('%.2f dB at DC', item.response_peak_db);
else
    value = sprintf('%.2f dB @ %.4g Hz', ...
        item.response_peak_db, item.response_peak_frequency_hz);
end
end

function value = format_basis(basis)
if contains(basis, 'documented snapshot')
    value = '分析快照（当前对象为空）';
else
    value = '当前模型计算';
end
end

function issues = model_content_issues(dynamics)
issues = {};
for i = 1:numel(dynamics)
    if contains(dynamics(i).analysis_basis, 'current model object is empty')
        issues{end+1} = sprintf('%s: ARMAXmodel.model is empty; restore A/B/C polynomials', ...
            dynamics(i).id); %#ok<AGROW>
    end
end
end

function [missing, unregistered] = audit_dataset(rootDir, registry)
missing = {};
registeredMat = {};
for i = 1:numel(registry)
    relative = normalize_path(registry(i).artifact);
    fullPath = fullfile(rootDir, strrep(relative, '/', filesep));
    if ~exist(fullPath, 'file')
        missing{end+1} = relative; %#ok<AGROW>
    end
    if startsWith(relative, 'dataset/') && endsWith(relative, '.mat')
        registeredMat{end+1} = relative; %#ok<AGROW>
    end
end

files = dir(fullfile(rootDir, 'dataset', '*.mat'));
unregistered = {};
for i = 1:numel(files)
    relative = ['dataset/' files(i).name];
    if ~any(strcmp(relative, registeredMat))
        unregistered{end+1} = relative; %#ok<AGROW>
    end
end
end

function records = scan_sysid_models(rootDir)
template = struct('path', '', 'type', '', 'date', '', 'path_role', '', ...
    'order', '', 'sample_rate', '', 'origin', '', 'has_report', false, ...
    'report_kind', '无', 'status', '');
records = repmat(template, 0, 1);
sysidDir = fullfile(rootDir, 'sysid_models');
if ~isfolder(sysidDir)
    return;
end

files = dir(fullfile(sysidDir, '**', '*.mat'));
paths = cell(1, numel(files));
for i = 1:numel(files)
    paths{i} = normalize_path(fullfile(files(i).folder, files(i).name));
end
[~, orderIndex] = sort(paths);
files = files(orderIndex);

for i = 1:numel(files)
    fullPath = fullfile(files(i).folder, files(i).name);
    relative = relative_path(rootDir, fullPath);
    [~, baseName] = fileparts(files(i).name);

    item = template;
    item.path = relative;
    item.type = parse_type(baseName);
    item.date = first_token(baseName, '(20\d{6})', 'unknown');
    item.path_role = first_token(lower(baseName), ...
        '(primpath|pripath|secpath|pri-err|pri-ref|sec-err|sec-ref)', 'unspecified');
    item.order = parse_order(baseName);
    item.sample_rate = parse_sample_rate(baseName);
    if contains(lower(relative), '_syn') || contains(lower(relative), '20260708_ho')
        item.origin = 'literature synthetic';
    else
        item.origin = 'measured identification';
    end

    reportPath = regexprep(fullPath, '\.mat$', '.md', 'ignorecase');
    item.has_report = exist(reportPath, 'file') == 2;
    if item.has_report
        item.report_kind = '同名';
    elseif strcmp(item.origin, 'literature synthetic')
        batchReports = dir(fullfile(files(i).folder, 'SYN*.md'));
        if ~isempty(batchReports)
            reportPath = fullfile(batchReports(1).folder, batchReports(1).name);
            item.has_report = true;
            item.report_kind = '批次';
        end
    end
    item.status = infer_sysid_status(baseName, reportPath, item.has_report, item.origin);
    records(end+1) = item; %#ok<AGROW>
end
end

function candidates = scan_embedded_models(rootDir, inlineRecords)
files = dir(fullfile(rootDir, 'demo*', '**', '*.m'));
registered = normalize_cell({inlineRecords.artifact});
candidates = {};
for i = 1:numel(files)
    fullPath = fullfile(files(i).folder, files(i).name);
    relative = relative_path(rootDir, fullPath);
    if any(strcmp(relative, registered))
        continue;
    end
    try
        source = fileread(fullPath);
    catch
        continue;
    end
    constructsModel = ~isempty(regexp(source, ...
        '(?<![A-Za-z0-9_])(tf|ss|zpk|tf2ss)\s*\(', 'once'));
    loadsModel = ~isempty(regexp(source, ...
        '(?<![A-Za-z0-9_])(load|DataManager)\s*\(', 'once'));
    if constructsModel && ~loadsModel
        candidates{end+1} = relative; %#ok<AGROW>
    end
end
candidates = sort(unique(candidates));
end

function value = parse_type(baseName)
upperName = upper(baseName);
if startsWith(upperName, 'ARMAX')
    value = 'ARMAX';
elseif startsWith(upperName, 'LMS') || startsWith(upperName, 'FIR')
    value = 'FIR/LMS';
elseif startsWith(upperName, 'OE')
    value = 'OE';
elseif startsWith(upperName, 'SS')
    value = 'state space';
else
    value = first_token(baseName, '^([^_]+)', 'unknown');
end
end

function value = parse_order(baseName)
value = first_token(baseName, '\[([^\]]+)\]', '');
if ~isempty(value)
    value = ['[' value ']'];
    return;
end
value = first_token(baseName, '_N([^_]+)', '');
if ~isempty(value)
    value = ['N=' value];
    return;
end
value = first_token(baseName, '_n(\d+)', 'unknown');
if ~strcmp(value, 'unknown')
    value = ['n=' value];
end
end

function value = parse_sample_rate(baseName)
tokens = regexp(lower(baseName), 'fs(\d+)(k?)', 'tokens', 'once');
if isempty(tokens)
    value = 'unknown';
elseif strcmp(tokens{2}, 'k')
    value = [tokens{1} ' kHz'];
else
    value = [tokens{1} ' Hz'];
end
end

function status = infer_sysid_status(baseName, reportPath, hasReport, origin)
lowerName = lower(baseName);
if contains(lowerName, 'legacy')
    status = 'legacy; metadata incomplete';
elseif contains(lowerName, 'retry')
    status = 'provisional/retry';
elseif ~hasReport
    status = 'undocumented';
elseif strcmp(origin, 'literature synthetic')
    status = 'reference equivalent';
else
    status = 'documented candidate';
    try
        report = fileread(reportPath);
        if ~isempty(regexpi(report, 'Residual Whiteness\s*\|\s*Failed', 'once'))
            status = 'candidate; residual whiteness failed';
        end
    catch
        status = 'report unreadable';
    end
end
end

function value = first_token(textValue, pattern, fallback)
token = regexp(textValue, pattern, 'tokens', 'once');
if isempty(token)
    value = fallback;
else
    value = token{1};
end
end

function write_audit_list(fid, heading, items, explanation)
writeln(fid, sprintf('### %s', heading));
writeln(fid, '');
writeln(fid, explanation);
writeln(fid, '');
if isempty(items)
    writeln(fid, '- 无。');
else
    for i = 1:numel(items)
        writeln(fid, sprintf('- `%s`', md(items{i})));
    end
end
writeln(fid, '');
end

function writeln(fid, value)
fprintf(fid, '%s\n', char(value));
end

function value = md(value)
value = char(value);
value = strrep(value, '|', '\|');
value = regexprep(value, '[\r\n]+', ' ');
end

function value = normalize_path(value)
value = strrep(char(value), '\', '/');
end

function values = normalize_cell(values)
for i = 1:numel(values)
    values{i} = normalize_path(values{i});
end
end

function relative = relative_path(rootDir, fullPath)
root = [normalize_path(rootDir) '/'];
full = normalize_path(fullPath);
if startsWith(full, root)
    relative = full(length(root)+1:end);
else
    relative = full;
end
end
