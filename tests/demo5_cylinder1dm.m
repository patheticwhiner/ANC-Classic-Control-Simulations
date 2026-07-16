% DEMO5_CYLINDER1DM  Marino-Tomei adaptive frequency estimator vs IMC-FxLMS.
%
%   Feedback-only ANC comparison on the Cylinder 1 dm acoustic model.
%   Both controllers use only the error microphone — no reference sensor.
%
%   Marino-Tomei (2016): adaptive internal model estimates narrowband
%   frequencies online and generates the cancellation signal. Assumes the
%   disturbance consists of persistent sinusoidal components.
%
%   IMC-FxLMS: classical feedback ANC baseline. Reconstructs an internal
%   disturbance reference through the IMC structure and adapts a control
%   FIR via filtered-x NLMS.
%
%   T1 (357 Hz fixed tone) tests the narrowband Case A design on the default
%   measured 2 kHz FIR path. The benchmark is allowed to fail the 20 dB and
%   actuator gates.
%   T2 (300→420 Hz chirp) checks Case A applicability; the secondary-path
%   real part crosses zero in this band, which is a structural boundary.

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath'));
run(fullfile(scriptDir, 'startup.m'));

%% Part 0 — Experiment setup
signals = load_cylinder1dm_signals('evaluation');

fprintf('========== Demo5: Marino-Tomei vs IMC-FxLMS (feedback ANC) ==========\n');
fprintf('  Model: Cylinder 1 dm, fs=%d Hz\n', signals.fs);
fprintf('  Paths: %s / %s\n', signals.primary_model, signals.model);
fprintf('  Evaluation split: seed %d, T1/T2 held out\n', signals.rng_seed);
fprintf('  Information: error microphone only (no reference sensor)\n');
fprintf('  T1: %.0f Hz fixed tone\n', signals.T1.f_hz);
fprintf('  T2: %.0f -> %.0f Hz %s chirp\n', ...
    signals.T2.f_range(1), signals.T2.f_range(2), signals.T2.direction);
fprintf('  Actuator: hard limit 5, tuning demand limit 4\n\n');
fprintf('  Frozen parameters come from calibration; this script only replays\n');
fprintf('  them on the held-out evaluation split.\n\n');

%% Part 1 — T1 fixed: Marino-Tomei, epsilon=0 (frozen initial guess)
%
%  The frozen candidate uses the true design frequency but no frequency
%  adaptation. Its result is retained as the fixed reference, including
%  any control-constraint failure.

fprintf('--- Part 1: T1 fixed (Marino-Tomei, epsilon=0) ---\n');
fallbackMTfixed = struct('q', 1, 'k', 3e-3, 'epsilon', 0, ...
    'freq_init_hz', 357, 'method', 'exact', ...
    'ramp_seconds', 0.5, 'dc_cancel', false, 'output_timing', 'updated');
[paramsMTfixed, sourceMTfixed] = frozen_params_or_default(scriptDir, signals, ...
    'demo5', 'T1', 'fixed', fallbackMTfixed);
fprintf('  Parameters: %s\n', sourceMTfixed);
resultT1f_MT = controller_demo5(signals, 'T1', 'fixed', paramsMTfixed);

%% Part 2 — T1 adaptive: Marino-Tomei online frequency estimation
%
%  The adaptive candidate starts close to the target because the published
%  convergence result is local; the evaluation signal is never used to
%  alter this initial condition.

fprintf('\n--- Part 2: T1 adaptive (Marino-Tomei, epsilon>0) ---\n');
fallbackMTadapt = struct('q', 1, 'k', 3e-3, 'epsilon', 1e-4, ...
    'freq_init_hz', 350, 'method', 'exact', ...
    'ramp_seconds', 0.5, 'dc_cancel', false, 'output_timing', 'updated');
[paramsMTadapt, sourceMTadapt] = frozen_params_or_default(scriptDir, signals, ...
    'demo5', 'T1', 'adaptive', fallbackMTadapt);
fprintf('  Parameters: %s\n', sourceMTadapt);
resultT1a_MT = controller_demo5(signals, 'T1', 'adaptive', paramsMTadapt);

%% Part 3 — T1 IMC-FxLMS baseline
%
%  Feedback IMC-FxLMS on the same 357 Hz tone. Uses IMC internal reference
%  reconstruction — same information structure as Marino-Tomei (no ref mic).

fprintf('\n--- Part 3: T1 IMC-FxLMS baseline ---\n');
paramsIMC_T1 = struct('N_fir', 64, 'mu', 0.01, 'delta', 1e-4, ...
    'ramp_seconds', 0.5);
resultT1a_IMC = controller_imc_fxlms(signals, 'T1', 'adaptive', paramsIMC_T1);

%% Part 4 — T2 adaptive: Marino-Tomei chirp tracking
%
%  T2 is a 300→420 Hz linear chirp over 10 s. Marino-Tomei's frequency
%  estimator is designed for constant frequencies — the Lyapunov stability
%  argument assumes d(theta)/dt = 0 for the true parameters.
%
%  During a chirp, the adaptive law theta_dot = epsilon * sgn * w2 * e
%  produces an estimate that lags behind the instantaneous frequency.
%  The lag depends on epsilon: larger epsilon tracks faster but amplifies
%  noise; smaller epsilon smooths noise but increases lag.
%
%  This section reports both the chirp suppression and the frequency
%  tracking error as a function of time — a diagnostic that Demo1–4
%  cannot provide because they do not estimate explicit frequency parameters.

fprintf('\n--- Part 4: T2 adaptive (Marino-Tomei chirp tracking) ---\n');
fallbackMTchirp = struct('q', 1, 'k', 3e-4, 'epsilon', 1e-4, ...
    'freq_init_hz', 420, 'method', 'exact', ...
    'ramp_seconds', 0.5, 'dc_cancel', false, 'output_timing', 'updated');
[paramsMTchirp, sourceMTchirp] = frozen_params_or_default(scriptDir, signals, ...
    'demo5', 'T2', 'adaptive', fallbackMTchirp);
fprintf('  Parameters: %s\n', sourceMTchirp);
resultT2_MT = controller_demo5(signals, 'T2', 'adaptive', paramsMTchirp);

%% Part 5 — T2 IMC-FxLMS baseline
%
%  IMC-FxLMS on the same chirp. The FIR filter has no explicit frequency
%  model — it adapts tap weights directly. Comparison shows whether an
%  explicit frequency parameterisation helps or hurts on non-stationary tones.

fprintf('\n--- Part 5: T2 IMC-FxLMS baseline ---\n');
paramsIMC_T2 = struct('N_fir', 64, 'mu', 0.02, 'delta', 1e-4, ...
    'ramp_seconds', 0.5, 'norm_limit', 10);
resultT2_IMC = controller_imc_fxlms(signals, 'T2', 'adaptive', paramsIMC_T2);

%% Part 6 — Comparison summary
fprintf('\n========== Demo5 comparison summary ==========\n');
fprintf('  Controller           | Test | Variant  |   Supp dB | max|u| | Clips | t_conv\n');
fprintf('  ---------------------|------|----------|-----------|---------|-------|-------\n');
print_result('Marino-Tomei', resultT1f_MT);
print_result('Marino-Tomei', resultT1a_MT);
print_result('IMC-FxLMS    ', resultT1a_IMC);
print_result('Marino-Tomei', resultT2_MT);
print_result('IMC-FxLMS    ', resultT2_IMC);

%% Part 7 — Frozen overview and Case A diagnostics
figureRoot = fullfile(scriptDir, 'figures', 'cylinder1dm_2k');
figureDir = fullfile(figureRoot, 'demo5');
if ~isfolder(figureDir), mkdir(figureDir); end
overview = struct('T1_fixed', resultT1f_MT, ...
    'T1_adaptive', resultT1a_MT, 'T2', resultT2_MT);
plot_demo_overview(overview, fullfile(figureDir, 'demo5_analysis.png'));
% Reuse the report-grade T1/T2 producers so the readable Demo5 script emits
% the same evidence panels as the other Demo scripts. Design/candidate
% analysis remains tied to the full calibration stage and is rendered by
% render_cylinder1dm('demo5').
directStage = struct();
directStage.selections = [ ...
    struct('demo', 'demo5', 'test', 'T1', 'variant', 'fixed'), ...
    struct('demo', 'demo5', 'test', 'T1', 'variant', 'adaptive'), ...
    struct('demo', 'demo5', 'test', 'T2', 'variant', 'adaptive')];
directStage.evaluation_results = {resultT1f_MT, resultT1a_MT, resultT2_MT};
generate_demo1234_reports(directStage, figureRoot, scriptDir, 'demo5', 't1');
generate_demo1234_reports(directStage, figureRoot, scriptDir, 'demo5', 't2');
fprintf('\n  Case A T1 applicable: %d (min ratio %.3f)\n', ...
    resultT1a_MT.extra.case_a_applicable, resultT1a_MT.extra.case_a_min_ratio);
fprintf('  Case A T2 applicable: %d (min ratio %.3f, sign-consistent %d)\n', ...
    resultT2_MT.extra.case_a_applicable, resultT2_MT.extra.case_a_min_ratio, ...
    resultT2_MT.extra.case_a_sign_consistent);

%% Part 8 — T2 frequency tracking analysis
%
%  For the chirp case, Marino-Tomei exposes its internal frequency estimate
%  f_hat(t). Plotting it against the instantaneous chirp frequency f_inst(t)
%  quantifies the tracking lag directly.

if isfield(resultT2_MT.extra, 'f_est_log') && ...
   isfield(signals.T2, 'f_inst')
    t_log = resultT2_MT.extra.t_log;
    f_est = resultT2_MT.extra.f_est_log(1, :);
    t_chirp = signals.T2.t(:);
    f_inst  = signals.T2.f_inst(:);

    f_err_rms = resultT2_MT.extra.frequency_tracking_error_rms_hz;

    figure('Name', 'Demo5 T2 frequency tracking', ...
        'Position', [100 100 900 420], 'Color', 'w');
    hold on;
    h1 = plot(t_chirp, f_inst, 'Color', [0.6 0.6 0.6], 'LineWidth', 1.0, ...
        'DisplayName', 'Instantaneous chirp frequency');
    h2 = plot(t_log, f_est, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.5, ...
        'DisplayName', sprintf('Marino-Tomei estimate (RMS err %.1f Hz)', f_err_rms));
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    title('Demo5 T2: Marino-Tomei frequency tracking vs chirp');
    legend([h1 h2], 'Location', 'best');
    grid on;

    saveas(gcf, fullfile(figureDir, 'demo5_t2_tracking.png'));
    fprintf('\n  T2 tracking: RMS frequency error = %.1f Hz (steady window)\n', f_err_rms);
end

%% Part 8 — Discussion
fprintf('\n========== Discussion ==========\n');
fprintf(['  Marino-Tomei is a narrowband, feedback-only method. Its core\n' ...
         '  assumption — that the disturbance consists of a small number of\n' ...
         '  constant-frequency sinusoids — is satisfied by T1 but violated\n' ...
         '  by T2 (chirp) and T3 (broadband noise).\n\n' ...
         '  T1: The frozen evaluation is judged against 20 dB suppression,\n' ...
         '      demand <= 4, and zero clipping; partial suppression remains\n' ...
         '      a failed controller design.\n\n' ...
         '  T2: The secondary-path real-part sign changes inside 300-420 Hz.\n' ...
         '      The current Case A implementation is therefore marked\n' ...
         '      inapplicable before frequency tracking is called successful.\n\n' ...
         '  T3 (broadband): Not tested. The method has no mechanism for\n' ...
         '      broadband noise — it models discrete sinusoidal components.\n\n' ...
         '  Comparison with IMC-FxLMS: Both use only the error microphone.\n' ...
         '      IMC-FxLMS makes no explicit frequency model — it adapts FIR\n' ...
         '      tap weights directly. For narrowband tones (T1), Marino-Tomei\n' ...
         '      can be more parameter-efficient (2 states per tone vs 64+ taps).\n' ...
         '      For non-stationary tones (T2), the explicit frequency model\n' ...
         '      introduces tracking lag that FIR-based methods avoid.\n']);

%% Local functions

function print_result(label, r)
    fprintf('  %-20s | %-4s | %-8s | %8.2f | %7.3f | %5d | %6.2f\n', ...
        label, r.test, r.variant, r.supp_db, r.u_max, ...
        r.extra.saturation_count, r.t_conv_s);
end

function [params, source] = frozen_params_or_default( ...
        scriptDir, signals, demo, testName, variant, fallback)
params = fallback;
source = 'FIR diagnostic fallback (run the full stage to freeze parameters)';
stageFile = fullfile(scriptDir, 'output', 'cylinder1dm_2k', ...
    'demo1234_stage', 'demo1234_stage_results.mat');
if ~isfile(stageFile)
    return;
end

payload = load(stageFile, 'stage');
if ~isfield(payload, 'stage') || ~isfield(payload.stage, 'suite') ...
        || ~isfield(payload.stage.suite, 'meta') ...
        || ~isfield(payload.stage.suite.meta, 'secondary_model') ...
        || ~strcmp(payload.stage.suite.meta.secondary_model, signals.model)
    return;
end

stage = payload.stage;
index = find(strcmp({stage.selections.demo}, demo) ...
    & strcmp({stage.selections.test}, testName) ...
    & strcmp({stage.selections.variant}, variant), 1);
if isempty(index)
    return;
end

params = stage.selections(index).params;
source = sprintf('frozen calibration candidate %s', ...
    stage.selections(index).candidate);
end
