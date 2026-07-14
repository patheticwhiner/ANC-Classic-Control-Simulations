function controller_validator()
% CONTROLLER_VALIDATOR  RST 控制器综合验证工具
%
% 功能：
%   加载 PPMaster 导出的模型文件（.mat，含 A,B,Ts）和控制器文件（.mat，含 R,S,T,Am,Bm），
%   通过 5 个标签页对闭环系统进行控制理论经典场景的全面验证：
%
%   Tab 1 - 📋 系统概览 ：左侧可拖动理论公式 (Theory.png) + 右侧闭环极零图
%                         + 多项式/稳定性面板
%   （系统框图 Schema.png 位于所有标签页上方的控制面板右侧，始终可见）
%   Tab 2 - 📊 频域分析 ：Syp/Sup/Syb 灵敏度函数、Nyquist 图、S+T=1 恒等式
%   Tab 3 - ⏱️ 时域分析 ：扰动阶跃/脉冲响应、控制量展示
%   Tab 4 - 🔇 噪声抑制 ：谐波噪声开环 vs 闭环对比、各频率衰减 dB 表
%   Tab 5 - 🛡️ 鲁棒性   ：模裕度、延迟裕度、增益/相位裕度、灵敏度峰值
%
% 用法：
%   controller_validator    % 启动 GUI，手动选择模型和控制器文件
%
% 依赖：
%   需要加载本项目的 shared/ 公共函数，以使用 isschur、addPolynomials 等函数。
%   脚本启动时会自动尝试添加路径。

%#ok<*NASGU>
%#ok<*AGROW>

% =========================================================================
% 路径初始化
% =========================================================================
scriptDir = fileparts(mfilename('fullpath'));
run(fullfile(scriptDir, 'project_init.m'));
% 同时添加当前目录（demo4_Robust 等子目录需要）
addpath(scriptDir);

% =========================================================================
% 全局状态
% =========================================================================
state = struct();
state.modelLoaded = false;
state.ctlLoaded   = false;
state.A     = [];  % 对象分母多项式
state.B     = [];  % 对象分子多项式
state.Ts    = [];  % 采样时间
state.R     = [];  % 控制器分子
state.S     = [];  % 控制器分母
state.Tval  = [];  % 跟踪增益
state.Am    = [];  % 跟踪模型分母
state.Bm    = [];  % 跟踪模型分子
state.modelFile = '';
state.ctlFile   = '';
state.statusHistory = {};

% =========================================================================
% 创建 GUI 主窗口
% =========================================================================
fig = figure( ...
    'Name', '🏗️ RST Controller Validation Tool', ...
    'NumberTitle', 'off', ...
    'Units', 'normalized', ...
    'Position', [0.05 0.05 0.85 0.82], ...
    'MenuBar', 'none', ...
    'ToolBar', 'figure', ...
    'Tag', 'mainFig');

% 顶层控制面板 — 继续增大高度
topPanel = uipanel(fig, 'Units', 'normalized', ...
    'Position', [0 0.72 1 0.26], ...
    'BorderType', 'none', 'BackgroundColor', get(fig, 'Color'));

% 加载模型按钮（缩小高度，靠上方排布）
uicontrol(topPanel, 'Style', 'pushbutton', ...
    'String', '📂 加载模型 (.mat)', ...
    'Units', 'normalized', 'Position', [0.01 0.76 0.14 0.18], ...
    'Callback', @(src,evt) loadModelCallback(), ...
    'FontWeight', 'bold', 'FontSize', 9, ...
    'BackgroundColor', [0.7 0.9 0.7]);

% 加载控制器按钮
uicontrol(topPanel, 'Style', 'pushbutton', ...
    'String', '📂 加载控制器 (.mat)', ...
    'Units', 'normalized', 'Position', [0.01 0.55 0.14 0.18], ...
    'Callback', @(src,evt) loadControllerCallback(), ...
    'FontWeight', 'bold', 'FontSize', 9, ...
    'BackgroundColor', [0.7 0.8 0.9]);

% 模型文件信息（移至按钮下方，左对齐一列）
uicontrol(topPanel, 'Style', 'text', ...
    'String', '模型: (未加载)', ...
    'Units', 'normalized', 'Position', [0.01 0.42 0.14 0.10], ...
    'HorizontalAlignment', 'left', 'FontSize', 9, 'Tag', 'lblModel', ...
    'BackgroundColor', 'none');

% 控制器文件信息
uicontrol(topPanel, 'Style', 'text', ...
    'String', '控制器: (未加载)', ...
    'Units', 'normalized', 'Position', [0.01 0.30 0.14 0.10], ...
    'HorizontalAlignment', 'left', 'FontSize', 9, 'Tag', 'lblController', ...
    'BackgroundColor', 'none');

% 采样周期显示（透明背景）
uicontrol(topPanel, 'Style', 'text', ...
    'String', '采样周期 Ts: -- s', ...
    'Units', 'normalized', 'Position', [0.01 0.18 0.14 0.10], ...
    'HorizontalAlignment', 'left', 'FontSize', 9, 'Tag', 'lblTs', ...
    'BackgroundColor', 'none');

% 状态指示灯（透明背景）
uicontrol(topPanel, 'Style', 'text', ...
    'String', '⏳ 等待加载...', ...
    'Units', 'normalized', 'Position', [0.01 0.04 0.14 0.12], ...
    'HorizontalAlignment', 'left', 'FontSize', 9, ...
    'ForegroundColor', [0.7 0.7 0.7], 'Tag', 'lblStatus', ...
    'BackgroundColor', 'none');

% 系统框图 (Schema.png) — 高度拉满，为控件面板腾出更多宽度
ax_schema = axes('Parent', topPanel, 'Units', 'normalized', ...
    'Position', [0.17 0.04 0.52 0.94], 'XTick', [], 'YTick', [], ...
    'Box', 'on', 'Tag', 'ax_schema');
% 启动时加载 Schema.png
imgFileSchema = fullfile(scriptDir, 'Images', 'Schema.png');
if exist(imgFileSchema, 'file') == 2
    imSchema = imread(imgFileSchema);
    axes(ax_schema);
    image(imSchema);
    daspect(ax_schema, [1 1 1]);
    set(ax_schema, 'XTick', [], 'YTick', []);
    title(ax_schema, '控制系统结构框图', 'FontSize', 8, 'FontWeight', 'bold');
end

% 右侧控件面板（加宽至 28%，字号全部升级）
% 各标签页控件共用同一区域但互斥可见，故每组均可从顶部开始松散排布
ctlPanel = uipanel(topPanel, 'Units', 'normalized', ...
    'Position', [0.71 0.02 0.28 0.96], ...
    'BorderType', 'none', 'BackgroundColor', get(fig, 'Color'));

% --- Tab 2 控件（频域分析） ---
btnFreq = uicontrol(ctlPanel, 'Style', 'pushbutton', ...
    'String', '🔄 刷新频域图', ...
    'Units', 'normalized', 'Position', [0.04 0.82 0.92 0.14], ...
    'Callback', @(src,evt) runFrequencyAnalysis(), ...
    'FontWeight', 'bold', 'FontSize', 10, ...
    'BackgroundColor', [0.7 0.8 0.9], ...
    'Visible', 'off');

% --- Tab 3 控件（时域分析） ---
lblTD = uicontrol(ctlPanel, 'Style', 'text', ...
    'String', '仿真点数:', ...
    'Units', 'normalized', 'Position', [0.04 0.82 0.50 0.07], ...
    'FontSize', 10, 'HorizontalAlignment', 'left', 'Visible', 'off');
editTD = uicontrol(ctlPanel, 'Style', 'edit', ...
    'String', '200', ...
    'Units', 'normalized', 'Position', [0.56 0.82 0.40 0.08], ...
    'FontSize', 10, 'Tag', 'editTimeDomainN', 'Visible', 'off');
btnTD = uicontrol(ctlPanel, 'Style', 'pushbutton', ...
    'String', '▶ 运行时域仿真', ...
    'Units', 'normalized', 'Position', [0.04 0.68 0.92 0.12], ...
    'Callback', @(src,evt) runTimeDomain(), ...
    'FontWeight', 'bold', 'FontSize', 10, ...
    'BackgroundColor', [0.7 0.9 0.5], ...
    'Visible', 'off');

% --- Tab 4 控件（噪声抑制，左右两栏布局） ---
% 左侧放基础参数，右侧放可选噪声参数，保持整体高度不变
lblNF = uicontrol(ctlPanel, 'Style', 'text', ...
    'String', '目标频率 (Hz):', ...
    'Units', 'normalized', 'Position', [0.04 0.87 0.40 0.08], 'FontSize', 10, 'HorizontalAlignment', 'left', 'Visible', 'off');
editNF = uicontrol(ctlPanel, 'Style', 'edit', ...
    'String', '100 200', ...
    'Units', 'normalized', 'Position', [0.04 0.78 0.40 0.09], ...
    'FontSize', 11, 'Tag', 'editFreqs', 'Visible', 'off');
lblND = uicontrol(ctlPanel, 'Style', 'text', ...
    'String', '仿真时长 (s):', ...
    'Units', 'normalized', 'Position', [0.04 0.68 0.40 0.08], ...
    'FontSize', 10, 'HorizontalAlignment', 'left', 'Visible', 'off');
editND = uicontrol(ctlPanel, 'Style', 'edit', ...
    'String', '1.0', ...
    'Units', 'normalized', 'Position', [0.04 0.58 0.40 0.09], ...
    'FontSize', 11, 'Tag', 'editDuration', 'Visible', 'off');

lblNT = uicontrol(ctlPanel, 'Style', 'text', ...
    'String', '信号类型:', ...
    'Units', 'normalized', 'Position', [0.04 0.48 0.40 0.08], ...
    'FontSize', 10, 'HorizontalAlignment', 'left', 'Visible', 'off');
popupNT = uicontrol(ctlPanel, 'Style', 'popupmenu', ...
    'String', {'仅谐波', '谐波+白噪', '谐波+有色噪', '纯白噪', '纯有色噪', '扫频正弦', '多音+白噪'}, ...
    'Units', 'normalized', 'Position', [0.04 0.38 0.40 0.10], ...
    'FontSize', 10, 'Tag', 'popupNoiseType', 'Visible', 'off', ...
    'Callback', @(src,evt) onNoiseTypeChanged());

% SNR 参数（含噪类型时可见）
lblSNR = uicontrol(ctlPanel, 'Style', 'text', ...
    'String', 'SNR (dB):', ...
    'Units', 'normalized', 'Position', [0.52 0.87 0.40 0.08], ...
    'FontSize', 10, 'HorizontalAlignment', 'left', 'Visible', 'off');
editSNR = uicontrol(ctlPanel, 'Style', 'edit', ...
    'String', '20', ...
    'Units', 'normalized', 'Position', [0.52 0.78 0.40 0.09], ...
    'FontSize', 11, 'Tag', 'editSNR', 'Visible', 'off');

% 噪声增益（白噪/有色噪/混合模式共用）
lblGain = uicontrol(ctlPanel, 'Style', 'text', ...
    'String', '噪声增益:', ...
    'Units', 'normalized', 'Position', [0.52 0.68 0.40 0.08], ...
    'FontSize', 10, 'HorizontalAlignment', 'left', 'Visible', 'off');
editGain = uicontrol(ctlPanel, 'Style', 'edit', ...
    'String', '1.0', ...
    'Units', 'normalized', 'Position', [0.52 0.59 0.40 0.09], ...
    'FontSize', 11, 'Tag', 'editNoiseGain', 'Visible', 'off');

% 随机种子（便于复现白噪/有色噪）
lblSeed = uicontrol(ctlPanel, 'Style', 'text', ...
    'String', '随机种子:', ...
    'Units', 'normalized', 'Position', [0.52 0.49 0.40 0.08], ...
    'FontSize', 10, 'HorizontalAlignment', 'left', 'Visible', 'off');
editSeed = uicontrol(ctlPanel, 'Style', 'edit', ...
    'String', '1', ...
    'Units', 'normalized', 'Position', [0.52 0.40 0.40 0.09], ...
    'FontSize', 11, 'Tag', 'editNoiseSeed', 'Visible', 'off');

% AR 系数（仅有色噪声时可见）
lblRho = uicontrol(ctlPanel, 'Style', 'text', ...
    'String', 'AR 系数 ρ:', ...
    'Units', 'normalized', 'Position', [0.52 0.30 0.40 0.08], ...
    'FontSize', 10, 'HorizontalAlignment', 'left', 'Visible', 'off', ...
    'TooltipString', '0 < ρ < 1，越接近1则相关性越强');
editRho = uicontrol(ctlPanel, 'Style', 'edit', ...
    'String', '0.9', ...
    'Units', 'normalized', 'Position', [0.52 0.21 0.40 0.09], ...
    'FontSize', 11, 'Tag', 'editRho', 'Visible', 'off');

% 运行按钮
btnNS = uicontrol(ctlPanel, 'Style', 'pushbutton', ...
    'String', '▶ 运行噪声仿真', ...
    'Units', 'normalized', 'Position', [0.04 0.06 0.88 0.11], ...
    'Callback', @(src,evt) runNoiseSuppression(), ...
    'FontWeight', 'bold', 'FontSize', 11, ...
    'BackgroundColor', [0.9 0.7 0.5], ...
    'Visible', 'off');

% =========================================================================
% 标签页组
% =========================================================================
tabGroup = uitabgroup(fig, 'Units', 'normalized', ...
    'Position', [0 0.02 1 0.69], 'Tag', 'tabGroup');

% 5 个标签页
tab1 = uitab(tabGroup, 'Title', '📋 系统概览');
tab2 = uitab(tabGroup, 'Title', '📊 频域分析');
tab3 = uitab(tabGroup, 'Title', '⏱️ 时域分析');
tab4 = uitab(tabGroup, 'Title', '🔇 噪声抑制');
tab5 = uitab(tabGroup, 'Title', '🛡️ 鲁棒性');

% 每个标签页内嵌一个 axes（占据主体）+ 右侧/底部控件面板
% 使用嵌套函数方式共用状态

% 为每个标签页创建 axes
% Tab 1 布局：左侧可拖动 RST 公式 (Theory.png) + 右侧极零图 + 多项式面板
scroll_theory = uipanel(tab1, 'Units', 'normalized', ...
    'Position', [0.03 0.03 0.47 0.94], ...
    'Title', 'RST 理论公式 (可拖动查看)', ...
    'Scrollable', 'on');
ax1_theory = axes('Parent', scroll_theory, ...
    'Units', 'pixels', 'Position', [0 0 941 1672], ...
    'XTick', [], 'YTick', [], 'Box', 'on', 'Tag', 'ax1_theory');
ax1_pz  = axes('Parent', tab1, 'Units', 'normalized', ...
    'Position', [0.54 0.35 0.44 0.60]);
panel1  = uipanel(tab1, 'Units', 'normalized', ...
    'Position', [0.54 0.03 0.44 0.30], 'Title', '模型 & 控制器多项式');
txt1    = uicontrol(panel1, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.02 0.96 0.96], 'HorizontalAlignment', 'left', ...
    'FontName', 'FixedWidth', 'FontSize', 9, 'Tag', 'txtOverview');
ax2a = axes('Parent', tab2, 'Units', 'normalized', 'Position', [0.08 0.55 0.42 0.40]);
ax2b = axes('Parent', tab2, 'Units', 'normalized', 'Position', [0.55 0.55 0.42 0.40]);
ax2c = axes('Parent', tab2, 'Units', 'normalized', 'Position', [0.08 0.08 0.42 0.40]);
ax2d = axes('Parent', tab2, 'Units', 'normalized', 'Position', [0.55 0.08 0.42 0.40]);
ax3a = axes('Parent', tab3, 'Units', 'normalized', 'Position', [0.08 0.55 0.42 0.40]);
ax3b = axes('Parent', tab3, 'Units', 'normalized', 'Position', [0.55 0.55 0.42 0.40]);
ax3c = axes('Parent', tab3, 'Units', 'normalized', 'Position', [0.08 0.08 0.90 0.40]);
ax4a = axes('Parent', tab4, 'Units', 'normalized', 'Position', [0.08 0.55 0.42 0.40]);
ax4b = axes('Parent', tab4, 'Units', 'normalized', 'Position', [0.55 0.55 0.42 0.40]);
ax4c = axes('Parent', tab4, 'Units', 'normalized', 'Position', [0.08 0.08 0.42 0.40]);
% Tab 4 右侧文本面板
panel4 = uipanel(tab4, 'Units', 'normalized', ...
    'Position', [0.55 0.02 0.42 0.50], 'Title', '频率衰减表');
txt4 = uicontrol(panel4, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.02 0.96 0.96], 'HorizontalAlignment', 'left', ...
    'FontName', 'FixedWidth', 'FontSize', 11, 'Tag', 'txtAttenuation');
ax5a = axes('Parent', tab5, 'Units', 'normalized', 'Position', [0.08 0.55 0.42 0.40]);
ax5b = axes('Parent', tab5, 'Units', 'normalized', 'Position', [0.55 0.55 0.42 0.40]);
% Tab 5 数值面板
panel5 = uipanel(tab5, 'Units', 'normalized', ...
    'Position', [0.08 0.02 0.90 0.50], 'Title', '鲁棒性指标');
txt5 = uicontrol(panel5, 'Style', 'text', 'Units', 'normalized', ...
    'Position', [0.02 0.02 0.96 0.96], 'HorizontalAlignment', 'left', ...
    'FontName', 'FixedWidth', 'FontSize', 11, 'Tag', 'txtRobustness');

% 存储所有句柄的标签
set(fig, 'UserData', state);

% 启动时立即加载理论公式图（无需等待模型/控制器加载）
runSystemOverview();

% =========================================================================
% 嵌套回调函数
% =========================================================================

    function loadModelCallback()
        [file, path] = uigetfile('*.mat', '选择模型文件 (.mat，含 A,B,Ts)');
        if isequal(file, 0), return; end
        fullPath = fullfile(path, file);
        data = load(fullPath);

        if ~isfield(data, 'A') || ~isfield(data, 'B') || ~isfield(data, 'Ts')
            errordlg('模型文件必须包含变量: A, B, Ts', '格式错误');
            return;
        end

        state.A  = data.A(:).';
        state.B  = data.B(:).';
        state.Ts = data.Ts;
        state.modelFile = fullPath;
        state.modelLoaded = true;

        % 确保多项式为行向量
        if size(state.A,1) > 1, state.A = state.A'; end
        if size(state.B,1) > 1, state.B = state.B'; end

        updateUI();
        logStatus(sprintf('✅ 模型已加载: %s (degA=%d, degB=%d)', ...
            file, length(state.A)-1, length(state.B)-1));

        if state.ctlLoaded
            refreshAllTabs();
        end
    end

    function loadControllerCallback()
        [file, path] = uigetfile('*.mat', '选择控制器文件 (.mat，含 R,S,T,Am,Bm)');
        if isequal(file, 0), return; end
        fullPath = fullfile(path, file);
        data = load(fullPath);

        if ~isfield(data, 'R') || ~isfield(data, 'S')
            errordlg('控制器文件必须包含变量: R, S', '格式错误');
            return;
        end

        state.R    = data.R(:).';
        state.S    = data.S(:).';
        state.ctlFile = fullPath;
        state.ctlLoaded = true;

        % 可选字段带默认值
        if isfield(data, 'T'),  state.Tval = data.T;  else state.Tval = 1; end
        if isfield(data, 'Am'), state.Am   = data.Am; else state.Am = 1; end
        if isfield(data, 'Bm'), state.Bm   = data.Bm; else state.Bm = 1; end

        % 确保行向量
        if size(state.R,1) > 1, state.R = state.R'; end
        if size(state.S,1) > 1, state.S = state.S'; end

        updateUI();
            logStatus(sprintf('✅ 控制器已加载: %s (degR=%d, degS=%d, T=%.3g)', ...
                file, length(state.R)-1, length(state.S)-1, state.Tval(1)));

        if state.modelLoaded
            refreshAllTabs();
        end
    end

    function updateUI()
        % 更新顶部信息栏
        h = findobj(fig, 'Tag', 'lblModel');
        if state.modelLoaded
            [~,f,ext] = fileparts(state.modelFile);
            set(h, 'String', sprintf('模型: %s%s', f, ext));
        else
            set(h, 'String', '模型: (未加载)');
        end

        h = findobj(fig, 'Tag', 'lblController');
        if state.ctlLoaded
            [~,f,ext] = fileparts(state.ctlFile);
            set(h, 'String', sprintf('控制器: %s%s', f, ext));
        else
            set(h, 'String', '控制器: (未加载)');
        end

        h = findobj(fig, 'Tag', 'lblTs');
        if ~isempty(state.Ts)
            set(h, 'String', sprintf('采样周期 Ts: %.2e s (%.1f Hz)', ...
                state.Ts, 1/state.Ts));
        end

        % 更新状态灯
        h = findobj(fig, 'Tag', 'lblStatus');
        if state.modelLoaded && state.ctlLoaded
            if isStable()
                set(h, 'String', '🟢 就绪 — 闭环系统稳定', ...
                    'ForegroundColor', [0 0.5 0]);
            else
                set(h, 'String', '🔴 就绪 — 闭环系统不稳定!', ...
                    'ForegroundColor', [0.8 0 0]);
            end
        elseif state.modelLoaded
            set(h, 'String', '🟡 等待加载控制器...', ...
                'ForegroundColor', [0.8 0.6 0]);
        elseif state.ctlLoaded
            set(h, 'String', '🟡 等待加载模型...', ...
                'ForegroundColor', [0.8 0.6 0]);
        end
    end

    function stable = isStable()
        % 检查闭环系统是否稳定
        if ~state.modelLoaded || ~state.ctlLoaded
            stable = false; return;
        end
        try
            AS = conv(state.A, state.S);
            BR = conv(state.B, state.R);

            % 对齐长度
            len = max(length(AS), length(BR));
            AS(len) = 0;
            BR(len) = 0;
            clPoly = AS + BR;

            r = roots(clPoly);
            stable = all(abs(r) < 1 - 1e-10);

            % 如果有 isschur 函数，使用更精确的判断
            if exist('isschur', 'file') == 2
                stable = isschur(clPoly);
            end
        catch
            stable = false;
        end
    end

    function logStatus(msg)
        state.statusHistory{end+1} = msg;
        fprintf('[controller_validator] %s\n', msg);
    end

    function refreshAllTabs()
        % 所有标签页全量刷新
        if ~state.modelLoaded || ~state.ctlLoaded, return; end

        % 根据当前可见标签页选择性刷新（首次全部刷新）
        runSystemOverview();
        runFrequencyAnalysis();
        runTimeDomain();
        % Tab 4 和 Tab 5 在用户点击时刷新
        updateUI();
    end

% =========================================================================
% Tab 1: 系统概览
% =========================================================================
    function runSystemOverview()
        cla(ax1_pz);
        cla(ax1_theory);
        set(txt1, 'String', '');

        % --- 左侧可滚动面板：加载 Theory.png，宽度自适应面板 ---
        imgFile = fullfile(scriptDir, 'Images', 'Theory.png');
        if exist(imgFile, 'file') == 2
            im = imread(imgFile);
            [imH, imW, ~] = size(im);
            % 获取可滚动面板的实际像素宽度
            panelPixels = getpixelposition(scroll_theory);
            availWidth = max(panelPixels(3) - 25, 100);  % 预留滚动条宽度
            % 按宽度缩放图像，保持比例
            scale = availWidth / imW;
            scaledH = round(imH * scale);
            imScaled = imresize(im, [scaledH, availWidth]);
            axes(ax1_theory);
            image(imScaled);
            set(ax1_theory, 'Units', 'pixels', ...
                'Position', [0 0 availWidth scaledH], ...
                'XTick', [], 'YTick', []);
            daspect(ax1_theory, [1 1 1]);
        else
            cla(ax1_theory);
            text(ax1_theory, 0.5, 0.5, '(Theory.png 未找到)', ...
                'Units', 'normalized', 'HorizontalAlignment', 'center', 'FontSize', 10);
            set(ax1_theory, 'XTick', [], 'YTick', []);
        end

        if ~state.modelLoaded || ~state.ctlLoaded
            title(ax1_pz, '请先加载模型和控制器文件');
            return;
        end

        % --- 右侧上部：闭环极零图 ---
        Ts_val = state.Ts;
        G = tf(state.B, state.A, Ts_val, 'Variable', 'z^-1');
        K = tf(state.R, state.S, Ts_val, 'Variable', 'z^-1');
        L = series(G, K);
        [Pcl, Zcl] = pzmap(feedback(L, 1));

        theta = linspace(0, 2*pi, 500);
        hold(ax1_pz, 'on');
        plot(ax1_pz, cos(theta), sin(theta), 'k--', 'LineWidth', 1);

        hP = [];
        hZ = [];
        if ~isempty(Pcl)
            hP = plot(ax1_pz, real(Pcl), imag(Pcl), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
        end
        if ~isempty(Zcl)
            hZ = plot(ax1_pz, real(Zcl), imag(Zcl), 'bo', 'MarkerSize', 8, 'LineWidth', 1.5);
        end

        if ~isempty(hP) && ~isempty(hZ)
            legend(ax1_pz, [hP(1), hZ(1)], {'闭环极点', '闭环零点'}, 'Location', 'best');
        elseif ~isempty(hP)
            legend(ax1_pz, hP(1), {'闭环极点'}, 'Location', 'best');
        end

        hold(ax1_pz, 'off');
        axis(ax1_pz, 'equal');
        xlim(ax1_pz, [-1.5 1.5]); ylim(ax1_pz, [-1.5 1.5]);
        xlabel(ax1_pz, '实轴'); ylabel(ax1_pz, '虚轴');
        grid(ax1_pz, 'on');

        nPoles = length(Pcl);
        nUnstable = sum(abs(Pcl) >= 1 - 1e-10);
        title(ax1_pz, sprintf( ...
            '闭环极零图  |  阶次: %d  |  稳定极点: %d/%d', ...
            nPoles, nPoles - nUnstable, nPoles), ...
            'FontSize', 10);

        % --- 底部右侧：多项式信息面板 ---
        if nUnstable == 0
            stabStr = 'Yes (stable)';
        else
            stabStr = 'No (unstable)';
        end
        textStr = sprintf( ...
            ['=== Model ===\n', ...
             'A (deg=%d) : %s\n', ...
             'B (deg=%d) : %s\n', ...
             'Ts = %.2e s  (Fs = %.1f Hz)\n\n', ...
             '=== Controller ===\n', ...
             'R (deg=%d) : %s\n', ...
             'S (deg=%d) : %s\n', ...
             'T = %.3g\n\n', ...
             '=== Closed Loop ===\n', ...
             'Order: %d\n', ...
             'Stable: %s  (%d/%d poles in unit circle)'], ...
            length(state.A)-1, poly2strLimit(state.A, 35), ...
            length(state.B)-1, poly2strLimit(state.B, 35), ...
            Ts_val, 1/Ts_val, ...
            length(state.R)-1, poly2strLimit(state.R, 35), ...
            length(state.S)-1, poly2strLimit(state.S, 35), ...
            state.Tval(1), ...
            nPoles, stabStr, nPoles - nUnstable, nPoles);
        set(txt1, 'String', textStr);
    end


% =========================================================================
% Tab 2: 频域分析
% =========================================================================
    function runFrequencyAnalysis()
        % 清除所有 axes
        cla(ax2a); cla(ax2b); cla(ax2c); cla(ax2d);

        if ~state.modelLoaded || ~state.ctlLoaded
            title(ax2a, '请先加载模型和控制器文件');
            return;
        end

        Ts_val = state.Ts;

        % 构造传递函数
        G = tf(state.B, state.A, Ts_val, 'Variable', 'z^-1');
        K = tf(state.R, state.S, Ts_val, 'Variable', 'z^-1');
        L = series(G, K);

        % 设计频率范围：避免 Nyquist 频率处奇异
        fNyq = 0.5 / Ts_val;
        fMin = max(1e-3, fNyq * 1e-4);
        fMax = fNyq * 0.99;
        freqHz = logspace(log10(fMin), log10(fMax), 2000);
        w = 2 * pi * freqHz;

        % 频率响应（避免 bode 在 Nyquist 附近 NaN）
        [magG, ~]   = freqRespSafe(G, w);
        [magK, ~]   = freqRespSafe(K, w);
        [magSyp, ~] = freqRespSafe(feedback(1, L), w);
        [magSup, ~] = freqRespSafe(feedback(-K, G), w);
        [magSyb, ~] = freqRespSafe(feedback(L, 1), w);

        % --- 图 2a: 输出灵敏度 Syp ---
        plotFreqDomain(ax2a, freqHz, magSyp);
        hold(ax2a, 'on');
        % 模板约束线
        z = exp(1i * 2 * pi * freqHz * Ts_val);
        Wsypu = min(6, 20*log10(1 + 1./abs(1 - 1./z)));
        Wsypd = real(20*log10(1 - 1./abs(1 - 1./z)));
        plot(ax2a, freqHz, Wsypu, 'r--', 'LineWidth', 1, 'DisplayName', '上界 6dB');
        plot(ax2a, freqHz, Wsypd, 'r:', 'LineWidth', 1, 'DisplayName', '下界');
        hold(ax2a, 'off');
        title(ax2a, '输出灵敏度 S_{yp} (含模板约束)');
        legend(ax2a, 'Location', 'best');

        % --- 图 2b: 输入灵敏度 Sup ---
        plotFreqDomain(ax2b, freqHz, magSup);
        hold(ax2b, 'on');
        plot(ax2b, freqHz([1 end]), [15 15], 'r--', 'LineWidth', 1, 'DisplayName', '上界 15dB');
        hold(ax2b, 'off');
        title(ax2b, '输入灵敏度 S_{up} (含模板约束)');
        legend(ax2b, 'Location', 'best');

        % --- 图 2c: 互补灵敏度 Syb ---
        plotFreqDomain(ax2c, freqHz, magSyb);
        hold(ax2c, 'on');
        Wsyb = 20*log10(1./abs(1 - 1./z));
        plot(ax2c, freqHz, Wsyb, 'r--', 'LineWidth', 1, 'DisplayName', '约束线');
        hold(ax2c, 'off');
        title(ax2c, '互补灵敏度 S_{yb}');
        legend(ax2c, 'Location', 'best');

        % --- 图 2d: S+T=1 恒等式验证 + Nyquist ---
        STsum = abs(magSyp(:) + magSyb(:));
        semilogx(ax2d, freqHz, 20*log10(STsum), 'b-', 'LineWidth', 1.5);
        hold(ax2d, 'on');
        plot(ax2d, freqHz([1 end]), [0 0], 'k--');
        hold(ax2d, 'off');
        xlabel(ax2d, '频率 (Hz)'); ylabel(ax2d, 'dB');
        set(ax2d, 'XScale', 'log');
        grid(ax2d, 'on');
        title(ax2d, 'S_{yp} + S_{yb} ≡ 1 验证');
        legend(ax2d, {'|S_{yp} + S_{yb}| (期望 0 dB)'}, 'Location', 'best');

        drawnow;
    end

    function [mag, ph] = freqRespSafe(sys, w)
        % 安全的频率响应计算，处理 z^-1 变量
        try
            [mag, ph] = bode(sys, w);
            mag = squeeze(mag);
            ph = squeeze(ph);
            % 替换 NaN / Inf
            mag(~isfinite(mag)) = 1e-10;
            ph(~isfinite(ph)) = 0;
        catch
            mag = ones(size(w));
            ph = zeros(size(w));
        end
    end

    function plotFreqDomain(ax, freqHz, magData)
        semilogx(ax, freqHz, 20*log10(abs(magData)), 'b-', 'LineWidth', 1.5);
        xlabel(ax, '频率 (Hz)'); ylabel(ax, '幅度 (dB)');
        set(ax, 'XScale', 'log');
        grid(ax, 'on');
    end

% =========================================================================
% Tab 3: 时域分析
% =========================================================================
    function runTimeDomain()
        cla(ax3a); cla(ax3b); cla(ax3c);

        if ~state.modelLoaded || ~state.ctlLoaded
            title(ax3a, '请先加载模型和控制器文件');
            return;
        end

        Ts_val = state.Ts;
        hEdit = findobj(fig, 'Tag', 'editTimeDomainN');
        N = str2double(get(hEdit, 'String'));
        if isnan(N) || N < 10, N = 200; end

        % 构造传递函数
        G = tf(state.B, state.A, Ts_val, 'Variable', 'z^-1');
        K = tf(state.R, state.S, Ts_val, 'Variable', 'z^-1');

        % 闭环：扰动到输出 Syp = 1/(1+GK)
        Syp_cl = feedback(1, series(G, K));

        % 从扰动到控制量 Sup = -K/(1+GK)
        Sup_cl = feedback(-K, G);

        % 仿真时间向量
        t = (0:N-1)' * Ts_val;

        % --- 图 3a: 扰动阶跃响应 ---
        [yStep, tStep] = step(Syp_cl, t(end));
        % 重采样到 t
        if length(tStep) > 1
            yStepInterp = interp1(tStep, yStep, t, 'linear', yStep(end));
        else
            yStepInterp = yStep * ones(size(t));
        end

        stairs(ax3a, t, yStepInterp, 'b-', 'LineWidth', 1.5);
        xlabel(ax3a, '时间 (s)'); ylabel(ax3a, 'y(k)');
        title(ax3a, '扰动阶跃响应 (d(k)=1(k))');
        grid(ax3a, 'on');

        % 计算超调和调节时间
        [peakVal, peakIdx] = max(abs(yStepInterp));
        finalVal = yStepInterp(end);  % 注意：Syp 的稳态应为 0（扰动抑制）
        overshoot = max(0, (peakVal - abs(finalVal)) / max(abs(finalVal), 1e-6) * 100);

        % --- 图 3b: 扰动脉冲响应 ---
        [yImp, tImp] = impulse(Syp_cl, t(end));
        if length(tImp) > 1
            yImpInterp = interp1(tImp, yImp, t, 'linear', 0);
        else
            yImpInterp = yImp * ones(size(t));
        end
        stem(ax3b, t(1:min(30,N)), yImpInterp(1:min(30,N)), 'b', 'LineWidth', 1);
        xlabel(ax3b, '时间 (s)'); ylabel(ax3b, 'y(k)');
        title(ax3b, '扰动脉冲响应 (d(k)=δ(k))');
        grid(ax3b, 'on');

        % --- 图 3c: 扰动→控制量 阶跃响应 ---
        [uStep, tUStep] = step(Sup_cl, t(end));
        if length(tUStep) > 1
            uStepInterp = interp1(tUStep, uStep, t, 'linear', uStep(end));
        else
            uStepInterp = uStep * ones(size(t));
        end
        stairs(ax3c, t, uStepInterp, 'r-', 'LineWidth', 1.5);
        xlabel(ax3c, '时间 (s)'); ylabel(ax3c, 'u(k)');
        title(ax3c, ...
            sprintf('控制量响应  |  峰值: %.3g  |  稳态: %.3g', ...
            max(abs(uStepInterp)), abs(uStepInterp(end))));
        grid(ax3c, 'on');

        drawnow;
    end

% =========================================================================
% Tab 4: 噪声抑制
% =========================================================================
    function runNoiseSuppression()
        cla(ax4a); cla(ax4b); cla(ax4c);

        if ~state.modelLoaded || ~state.ctlLoaded
            title(ax4a, '请先加载模型和控制器文件');
            return;
        end

        Ts_val = state.Ts;

        % 读取用户配置
        hEditFreqs = findobj(fig, 'Tag', 'editFreqs');
        freqStr = strtrim(get(hEditFreqs, 'String'));
        targetFreqs = sscanf(freqStr, '%f')';
        if isempty(targetFreqs)
            targetFreqs = [100 200];
        end

        hEditDur = findobj(fig, 'Tag', 'editDuration');
        simDur = str2double(get(hEditDur, 'String'));
        if isnan(simDur) || simDur <= 0, simDur = 1.0; end

        hPopup = findobj(fig, 'Tag', 'popupNoiseType');
        noiseType = get(hPopup, 'Value');

        % 时间向量
        Ns = max(2, round(simDur / Ts_val));  % lsim 要求至少 2 个采样点
        t = (0:Ns-1)' * Ts_val;

        % --- 读取 SNR 和 AR 系数 ---
        hEditSNR = findobj(fig, 'Tag', 'editSNR');
        snrDB = str2double(get(hEditSNR, 'String'));
        if isnan(snrDB) || snrDB < -20, snrDB = 20; end

        hEditGain = findobj(fig, 'Tag', 'editNoiseGain');
        noiseGain = str2double(get(hEditGain, 'String'));
        if isnan(noiseGain) || noiseGain <= 0, noiseGain = 1.0; end

        hEditSeed = findobj(fig, 'Tag', 'editNoiseSeed');
        noiseSeed = str2double(get(hEditSeed, 'String'));
        if isnan(noiseSeed) || noiseSeed < 0, noiseSeed = 1; end
        rng(round(noiseSeed), 'twister');

        hEditGain = findobj(fig, 'Tag', 'editNoiseGain');
        signalGain = str2double(get(hEditGain, 'String'));
        if isnan(signalGain) || signalGain <= 0, signalGain = 1.0; end

        hEditRho = findobj(fig, 'Tag', 'editRho');
        rhoVal = str2double(get(hEditRho, 'String'));
        if isnan(rhoVal) || rhoVal <= 0 || rhoVal >= 1, rhoVal = 0.9; end
        innovStd = sqrt(1 - rhoVal^2);  % 创新项标准差，保持 AR(1) 方差=1

        % 生成扰动信号
        disturbance = zeros(size(t));
        switch noiseType
            case 1  % 仅谐波
                phaseRad = deg2rad(str2double(get(hEditSeed, 'String')));
                for f = targetFreqs
                    disturbance = disturbance + signalGain * sin(2*pi*f*t + phaseRad);
                end
            case 2  % 谐波 + 白噪声
                d_harmonic = zeros(size(t));
                phaseRad = deg2rad(str2double(get(hEditSeed, 'String')));
                for f = targetFreqs
                    d_harmonic = d_harmonic + signalGain * sin(2*pi*f*t + phaseRad);
                end
                noiseLevel = signalGain * rms(d_harmonic) * 10^(-snrDB/20);
                disturbance = d_harmonic + noiseLevel * randn(size(t));
            case 3  % 谐波 + 有色噪声 (AR(1))
                d_harmonic = zeros(size(t));
                phaseRad = deg2rad(str2double(get(hEditSeed, 'String')));
                for f = targetFreqs
                    d_harmonic = d_harmonic + signalGain * sin(2*pi*f*t + phaseRad);
                end
                colored = zeros(size(t));
                colored(1) = noiseGain * randn;
                for i = 2:length(t)
                    colored(i) = rhoVal * colored(i-1) + noiseGain * innovStd * randn;
                end
                noiseLevel = signalGain * rms(d_harmonic) * 10^(-snrDB/20);
                disturbance = d_harmonic + noiseLevel * colored(:);
            case 4  % 纯白噪声
                disturbance = noiseGain * randn(size(t));
            case 5  % 纯有色噪声 (AR(1))
                disturbance = zeros(size(t));
                disturbance(1) = noiseGain * randn;
                for i = 2:length(t)
                    disturbance(i) = rhoVal * disturbance(i-1) + noiseGain * innovStd * randn;
                end
            case 6  % 扫频正弦 (chirp)
                f0 = max(1, min(targetFreqs)*0.5);
                f1 = max(targetFreqs) * 1.5;
                phaseRad = deg2rad(str2double(get(hEditSeed, 'String')));
                k = (f1 - f0) / max(t(end), eps);
                disturbance = signalGain * sin(2*pi*(f0*t + 0.5*k*t.^2) + phaseRad);
            case 7  % 多音 + 白噪声
                d_harmonic = zeros(size(t));
                phaseRad = deg2rad(str2double(get(hEditSeed, 'String')));
                for f = targetFreqs
                    d_harmonic = d_harmonic + signalGain * sin(2*pi*f*t + phaseRad);
                end
                noiseLevel = signalGain * rms(d_harmonic) * 10^(-snrDB/20);
                disturbance = d_harmonic + noiseLevel * randn(size(t));
        end

        % 构造传递函数
        G = tf(state.B, state.A, Ts_val, 'Variable', 'z^-1');
        K = tf(state.R, state.S, Ts_val, 'Variable', 'z^-1');
        Syp = feedback(1, series(G, K));

        % 开环输出 (仅经对象，无控制)
        y_open = lsim(G, disturbance, t);

        % 闭环输出 (经灵敏度函数)
        y_closed = lsim(Syp, disturbance, t);

        % --- 图 4a: 开环输出 ---
        plot(ax4a, t, disturbance, 'Color', [0.7 0.7 0.7], 'DisplayName', '扰动 d(k)');
        hold(ax4a, 'on');
        plot(ax4a, t, y_open, 'r-', 'LineWidth', 1, 'DisplayName', '开环输出 y_{open}');
        hold(ax4a, 'off');
        xlabel(ax4a, '时间 (s)'); ylabel(ax4a, '幅值');
        title(ax4a, '开环响应 (无控制)');
        legend(ax4a, 'Location', 'best'); grid(ax4a, 'on');

        % --- 图 4b: 闭环输出对比 ---
        plot(ax4b, t, disturbance, 'Color', [0.7 0.7 0.7], 'DisplayName', '扰动 d(k)');
        hold(ax4b, 'on');
        plot(ax4b, t, y_closed, 'b-', 'LineWidth', 1, 'DisplayName', '闭环输出 y_{closed}');
        hold(ax4b, 'off');
        xlabel(ax4b, '时间 (s)'); ylabel(ax4b, '幅值');
        title(ax4b, '闭环响应 (有控制)');
        legend(ax4b, 'Location', 'best'); grid(ax4b, 'on');

        % --- 图 4c: 开环 vs 闭环 PSD ---
        [Pxx_open, f_psd]  = computePSD(y_open, 1/Ts_val);
        [Pxx_closed, ~]    = computePSD(y_closed, 1/Ts_val);
        [Pxx_dist, ~]      = computePSD(disturbance, 1/Ts_val);

        semilogy(ax4c, f_psd, Pxx_dist, 'Color', [0.7 0.7 0.7], 'DisplayName', '扰动');
        hold(ax4c, 'on');
        semilogy(ax4c, f_psd, Pxx_open, 'r-', 'LineWidth', 1.2, 'DisplayName', '开环');
        semilogy(ax4c, f_psd, Pxx_closed, 'b-', 'LineWidth', 1.5, 'DisplayName', '闭环');
        % 在目标频率处绘制竖直线标记
        yl = ylim(ax4c);
        for f = targetFreqs
            plot(ax4c, [f f], yl, 'k--', 'LineWidth', 0.8, 'HandleVisibility', 'off');
        end
        hold(ax4c, 'off');
        xlabel(ax4c, '频率 (Hz)'); ylabel(ax4c, 'PSD (V^2/Hz)');
        title(ax4c, '功率谱密度对比');
        legend(ax4c, 'Location', 'best'); grid(ax4c, 'on');
        % 自动缩放：以最大目标频率的 1.5 倍为上限，但不裁剪 psd 数据
        xlim(ax4c, [0, min(max(f_psd), max(targetFreqs) * 1.5)]);

        % --- 计算各频率衰减 ---
        attenText = sprintf('═══ 频率衰减表 ═══\n\n');
        attenText = [attenText, sprintf('%-12s %12s %12s\n', '频率 (Hz)', '开环 (dB)', '闭环 (dB)')];
        attenText = [attenText, sprintf('%s\n', repmat('-', 1, 40))];

        totalAtten = 0;
        for f = targetFreqs
            [~, idx] = min(abs(f_psd - f));
            if idx > 0 && idx <= length(f_psd) && Pxx_open(idx) > 0 && Pxx_closed(idx) > 0
                openDB  = 10 * log10(Pxx_open(idx));
                closedDB = 10 * log10(Pxx_closed(idx));
                attenDB = openDB - closedDB;
                totalAtten = totalAtten + max(0, attenDB);
                attenText = [attenText, ...
                    sprintf('%8.1f Hz   %12.2f   %12.2f   (衰减: %.1f dB)\n', ...
                    f, openDB, closedDB, attenDB)]; %#ok<AGROW>
            end
        end
        attenText = [attenText, sprintf('%s\n', repmat('-', 1, 40))];
        attenText = [attenText, sprintf('总衰减 (和):  %.1f dB\n', totalAtten)];

        set(txt4, 'String', attenText);
        drawnow;
    end

    function [Pxx, f] = computePSD(x, Fs)
        % 使用 Welch 法计算功率谱密度
        n = length(x);
        nfft = min(2^nextpow2(n), 4096);
        window = hamming(min(n, 512));
        noverlap = round(length(window) * 0.5);
        [Pxx, f] = pwelch(x, window, noverlap, nfft, Fs);
    end

% =========================================================================
% Tab 5: 鲁棒性分析
% =========================================================================
    function runRobustness()
        cla(ax5a); cla(ax5b);
        set(txt5, 'String', '计算中...');

        if ~state.modelLoaded || ~state.ctlLoaded
            title(ax5a, '请先加载模型和控制器文件');
            set(txt5, 'String', '');
            return;
        end

        Ts_val = state.Ts;

        % 构造传递函数
        G = tf(state.B, state.A, Ts_val, 'Variable', 'z^-1');
        K = tf(state.R, state.S, Ts_val, 'Variable', 'z^-1');
        L = series(G, K);

        % 频率范围
        fNyq = 0.5 / Ts_val;
        freqHz = logspace(log10(fNyq*1e-4), log10(fNyq*0.99), 5000);
        w = 2 * pi * freqHz;

        % --- 图 5a: 开环 Nyquist + 模裕度 ---
        [magL, phL] = freqRespSafe(L, w);
        % 到 (-1,0) 点的距离
        Lreal = magL .* cosd(phL);
        Limag = magL .* sind(phL);
        dist = sqrt((Lreal + 1).^2 + Lreal.^2);
        dist(isnan(dist)) = inf;

        plot(ax5a, Lreal, Limag, 'b-', 'LineWidth', 1.2);
        hold(ax5a, 'on');
        plot(ax5a, -1, 0, 'rx', 'MarkerSize', 12, 'LineWidth', 2);
        % 单位圆
        th = linspace(0, 2*pi, 500);
        plot(ax5a, cos(th), sin(th), 'k:', 'LineWidth', 0.5);
        hold(ax5a, 'off');
        xlabel(ax5a, '实轴'); ylabel(ax5a, '虚轴');
        title(ax5a, 'Nyquist 图');
        grid(ax5a, 'on'); axis(ax5a, 'equal');

        % --- 图 5b: 灵敏度函数 + 模裕度标注 ---
        Syp = feedback(1, L);
        [magSyp, ~] = freqRespSafe(Syp, w);
        semilogx(ax5b, freqHz, 20*log10(magSyp), 'b-', 'LineWidth', 1.5);
        hold(ax5b, 'on');

        [SypPeak_dB, idxPeak] = max(20*log10(magSyp));
        plot(ax5b, freqHz(idxPeak), SypPeak_dB, 'ro', 'MarkerSize', 8, 'LineWidth', 2);

        % 6dB 约束线
        plot(ax5b, freqHz([1 end]), [6 6], 'r--', 'LineWidth', 1);
        hold(ax5b, 'off');
        xlabel(ax5b, '频率 (Hz)'); ylabel(ax5b, '幅度 (dB)');
        set(ax5b, 'XScale', 'log');
        title(ax5b, sprintf('输出灵敏度 |S_{yp}|  (峰值: %.2f dB)', SypPeak_dB));
        grid(ax5b, 'on'); legend(ax5b, 'Location', 'best');

        % --- 鲁棒性指标计算 ---
        % 模裕度 ΔM = min|1+L|
        modMargin = min(dist);
        modMargin_dB = 20 * log10(max(modMargin, 1e-10));

        % 灵敏度峰值 Ms
        Ms = max(magSyp);
        Ms_dB = 20 * log10(max(Ms, 1e-10));

        % 增益裕度和相位裕度
        try
            [Gm_dB, Pm_deg, Wcg, Wcp] = margin(L);
        catch
            Gm_dB = NaN; Pm_deg = NaN; Wcg = NaN; Wcp = NaN;
        end

        % 延迟裕度近似
        if ~isnan(Pm_deg) && Pm_deg > 0 && ~isnan(Wcp)
            delayMargin = deg2rad(Pm_deg) / Wcp;
        else
            delayMargin = NaN;
        end

        % --- 输出结果 ---
        robustText = sprintf('═══ 鲁棒性分析结果 ═══\n\n');
        robustText = [robustText, sprintf('采样周期 Ts:              %.2e s (%.1f Hz)\n', ...
            Ts_val, 1/Ts_val)];
        robustText = [robustText, sprintf('\n── 灵敏度与裕度 ──\n')];
        robustText = [robustText, sprintf('模裕度 ΔM:                %.4f  (%.2f dB)  [>0.5 (-6dB) 为优]\n', ...
            modMargin, modMargin_dB)];
        robustText = [robustText, sprintf('灵敏度峰值 Ms:            %.2f dB  [<6 dB 为优]\n', Ms_dB)];

        if ~isnan(Gm_dB)
            robustText = [robustText, sprintf('增益裕度 GM:              %.2f dB  [>6 dB 可接受]\n', Gm_dB)];
        else
            robustText = [robustText, '增益裕度 GM:              无法计算\n'];
        end

        if ~isnan(Pm_deg)
            robustText = [robustText, sprintf('相位裕度 PM:              %.1f°  [>30° 可接受]\n', Pm_deg)];
        else
            robustText = [robustText, '相位裕度 PM:              无法计算\n'];
        end

        if ~isnan(delayMargin)
            robustText = [robustText, sprintf('延迟裕度 Δτ:              %.3e s  (%.1f Ts)\n', ...
                delayMargin, delayMargin / Ts_val)];
        else
            robustText = [robustText, '延迟裕度 Δτ:              无法计算\n'];
        end

        robustText = [robustText, sprintf('\n── 灵敏度峰值位置 ──\n')];
        robustText = [robustText, sprintf('Syp 峰值频率:             %.2f Hz\n', freqHz(idxPeak))];

        robustText = [robustText, sprintf('\n── 模型与控制器信息 ──\n')];
        robustText = [robustText, sprintf('deg(A)=%d  deg(B)=%d  deg(R)=%d  deg(S)=%d\n', ...
            length(state.A)-1, length(state.B)-1, ...
            length(state.R)-1, length(state.S)-1)];

        % 状态评估
        robustText = [robustText, sprintf('\n── 综合评估 ──\n')];
        passed = 0;
        total = 0;

        if ~isnan(Gm_dB)
            total = total + 1;
            if Gm_dB > 6, passed = passed + 1;
                robustText = [robustText, '✅ 增益裕度 > 6 dB\n'];
            else
                robustText = [robustText, '⚠️  增益裕度 ≤ 6 dB\n'];
            end
        end

        if ~isnan(Pm_deg)
            total = total + 1;
            if Pm_deg > 30, passed = passed + 1;
                robustText = [robustText, '✅ 相位裕度 > 30°\n'];
            else
                robustText = [robustText, '⚠️  相位裕度 ≤ 30°\n'];
            end
        end

        total = total + 1;
        if Ms_dB < 6, passed = passed + 1;
            robustText = [robustText, '✅ 灵敏度峰值 < 6 dB\n'];
        else
            robustText = [robustText, '⚠️  灵敏度峰值 ≥ 6 dB\n'];
        end

        total = total + 1;
        if modMargin > 0.5, passed = passed + 1;
            robustText = [robustText, '✅ 模裕度 > 0.5\n'];
        else
            robustText = [robustText, '⚠️  模裕度 ≤ 0.5\n'];
        end

        robustText = [robustText, sprintf('\n通过率: %d / %d\n', passed, total)];

        set(txt5, 'String', robustText);
        drawnow;
    end

% =========================================================================
% 标签页切换回调
% =========================================================================
    set(tabGroup, 'SelectionChangedFcn', @(src, evt) onTabChanged(evt));

    function onTabChanged(evt)
        % 隐藏所有控件
        set(btnFreq, 'Visible', 'off');
        set([lblTD, editTD, btnTD], 'Visible', 'off');
        set([lblNF, editNF, lblND, editND, lblSNR, editSNR, ...
             lblGain, editGain, lblSeed, editSeed, lblRho, editRho, ...
             lblNT, popupNT, btnNS], 'Visible', 'off');
        if ~state.modelLoaded || ~state.ctlLoaded, return; end
        switch evt.NewValue.Title
            case '📋 系统概览'
                runSystemOverview();
            case '📊 频域分析'
                set(btnFreq, 'Visible', 'on');
                runFrequencyAnalysis();
            case '⏱️ 时域分析'
                set([lblTD, editTD, btnTD], 'Visible', 'on');
                runTimeDomain();
            case '🔇 噪声抑制'
                set([lblNF, editNF, lblND, editND, lblNT, popupNT, btnNS], 'Visible', 'on');
                onNoiseTypeChanged();  % 根据当前类型决定 SNR/Rho 显隐
                runNoiseSuppression();
            case '🛡️ 鲁棒性'
                runRobustness();
        end
    end

% =========================================================================
% 噪声类型切换回调 — 根据所选类型动态显示/隐藏参数
% =========================================================================
    function onNoiseTypeChanged()
        noiseType = get(popupNT, 'Value');
        isDeterministicTone = ismember(noiseType, [1, 2, 3, 6, 7]);
        isStochasticNoise   = ismember(noiseType, [2, 3, 4, 5, 7]);
        needsSNR            = ismember(noiseType, [2, 3, 7]);
        needsRho            = ismember(noiseType, [3, 5]);
        needsFreq           = ismember(noiseType, [1, 2, 3, 6, 7]);

        if needsFreq
            set([lblNF, editNF], 'Visible', 'on');
            if noiseType == 6
                set(lblNF, 'String', '起止频率 (Hz):');
            else
                set(lblNF, 'String', '频率列表 (Hz):');
            end
        else
            set([lblNF, editNF], 'Visible', 'off');
        end

        if isDeterministicTone
            set(lblGain, 'String', '谐波幅值:');
            set(lblSeed, 'String', '初相位 (deg):');
            set([lblGain, editGain], 'Visible', 'on');
            set([lblSeed, editSeed], 'Visible', 'on');
        elseif isStochasticNoise
            set(lblGain, 'String', '噪声增益:');
            set(lblSeed, 'String', '随机种子:');
            set([lblGain, editGain], 'Visible', 'on');
            set([lblSeed, editSeed], 'Visible', 'on');
        else
            set([lblGain, editGain, lblSeed, editSeed], 'Visible', 'off');
        end

        set([lblSNR, editSNR], 'Visible', ternaryVis(needsSNR));
        set([lblRho, editRho], 'Visible', ternaryVis(needsRho));

        % 纯随机噪声模式下，频率输入不参与生成，保留但弱化存在感
        if ismember(noiseType, [4, 5])
            set(lblNF, 'String', '频率参数 (不适用):');
        end
    end

    function vis = ternaryVis(flag)
        if flag
            vis = 'on';
        else
            vis = 'off';
        end
    end

% =========================================================================
% 工具函数
% =========================================================================
    function s = poly2strLimit(p, maxLen)
        % 将多项式转换为字符串，限制长度
        if isempty(p), s = '[]'; return; end
        p = p(:)';
        terms = {};
        for i = 1:min(length(p), 6)
            if i == 1
                terms{end+1} = sprintf('%.4g', p(i));
            elseif p(i) >= 0
                terms{end+1} = sprintf('+ %.4g z^{-%d}', p(i), i-1);
            else
                terms{end+1} = sprintf('- %.4g z^{-%d}', abs(p(i)), i-1);
            end
        end
        s = strjoin(terms, ' ');
        if length(p) > 6
            s = [s, sprintf(' ... (+%d项)', length(p)-6)];
        end
        if length(s) > maxLen
            s = [s(1:maxLen-3), '...'];
        end
    end

end
