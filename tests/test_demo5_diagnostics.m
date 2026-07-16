function test_demo5_diagnostics()
%TEST_DEMO5_DIAGNOSTICS Lock down the Marino-Tomei sign and Case A probes.
%
% This test deliberately uses a one-sample FIR before the measured ARMAX
% path. The FIR is a small seam that distinguishes a correct paper feedback
% direction from its inverse without depending on the report renderer.

testDir = fileparts(mfilename('fullpath'));
run(fullfile(testDir, 'startup.m'));

% Minimal secondary path: S(z) = z^(-1), 4 seconds at 357 Hz.
fs = 2000;
Tsim = 4;
t = (0:round(Tsim * fs) - 1)' / fs;
d = sin(2 * pi * 357 * t);
firSignals = struct();
firSignals.model = 'unit_delay_fir';
firSignals.fs = fs;
firSignals.model_data = struct('A', 1, 'B', [0 1], ...
    'orders', [0 1 0 0], 'fs', fs);
firSignals.T1 = struct('f_hz', 357, 'fs', fs, 'Tsim', Tsim, ...
    't', t, 'd', d, 'y_open', d);

base = struct('q', 1, 'k', 1e-4, 'epsilon', 0, ...
    'freq_init_hz', 357, 'method', 'exact', 'ramp_seconds', 0, ...
    'actuator_limit', 5, 'output_timing', 'updated');

paper = base;
paper.feedback_sign = 'paper';
benchmark = base;
benchmark.feedback_sign = 'benchmark';
paperResult = controller_demo5(firSignals, 'T1', 'fixed', paper);
benchmarkResult = controller_demo5(firSignals, 'T1', 'fixed', benchmark);
assert(paperResult.supp_db > 0, ...
    'Paper sign must reduce the one-delay FIR residual.');
assert(benchmarkResult.supp_db < 0, ...
    'Benchmark sign must increase the one-delay FIR residual.');

% The default suite must now use the measured LMS FIR path.
signals = load_cylinder1dm_signals('evaluation');
firParams = base;
firParams.feedback_sign = 'paper';
firParams.ramp_seconds = 0.5;
firParams.output_timing = 'updated';
firResult = controller_demo5(signals, 'T1', 'fixed', firParams);
assert(strcmp(firResult.extra.path_family, 'measured_lms_fir'), ...
    'Default Cylinder suite must use the measured LMS FIR path.');
assert(firResult.extra.case_a_applicable, ...
    'Measured FIR T1 should satisfy the discrete Case A diagnostic.');

firPrevious = firParams;
firPrevious.output_timing = 'previous';
firPreviousResult = controller_demo5(signals, 'T1', 'fixed', firPrevious);
assert(firResult.supp_db > firPreviousResult.supp_db, ...
    ['Updated-state output must remain the formal FIR timing; reverting to ' ...
     'the previous-state diagnostic weakens the measured-path T1 result.']);

% Retain the old ARMAX path as a diagnostic comparison. At 357 Hz it is a
% discrete Case A boundary for the exact, updated-output implementation.
armaxSignals = signals;
armaxInfo = DataManager('cylinder1dm_2k_secondary');
armaxSignals.model = 'cylinder1dm_2k_secondary';
armaxSignals.model_data = struct('A', armaxInfo.model.A(:).', ...
    'B', armaxInfo.model.B(:).', 'orders', armaxInfo.orders, ...
    'fs', armaxInfo.fs, 'path_family', 'measured_armax');
armax = base;
armax.feedback_sign = 'paper';
armax.ramp_seconds = 0.5;
armaxResult = controller_demo5(armaxSignals, 'T1', 'fixed', armax);
assert(armaxResult.extra.case_a_min_ratio < 0.1, ...
    'ARMAX T1 must expose the near-zero discrete effective Case A ratio.');
assert(~armaxResult.extra.case_a_applicable, ...
    'ARMAX T1 must not be marked Case A applicable at 357 Hz.');

armaxPrevious = armax;
armaxPrevious.output_timing = 'previous';
previousResult = controller_demo5(armaxSignals, 'T1', 'fixed', armaxPrevious);
assert(previousResult.extra.case_a_min_ratio > armaxResult.extra.case_a_min_ratio, ...
    'Previous-state output must move the ARMAX T1 effective response away from zero.');

t2 = firParams;
t2.freq_init_hz = 420;
t2.epsilon = 1e-4;
t2.k = 3e-4;
t2Result = controller_demo5(signals, 'T2', 'adaptive', t2);
assert(~t2Result.extra.case_a_applicable, ...
    'T2 sweep must be rejected when the effective Case A sign changes.');

fprintf('test_demo5_diagnostics: sign, discrete Case A, and T2 boundary passed.\n');
end
