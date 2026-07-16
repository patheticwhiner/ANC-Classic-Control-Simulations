% DEMO3_CYLINDER1DM F(z)+FIR and IMC-FxNLMS adaptive experiment.
%
% T1 compares the fixed FIR design with the frozen source-reference
% IMC-FxNLMS candidate. T2 evaluates the same structure on a reverse chirp.

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath'));
run(fullfile(scriptDir, 'startup.m'));

%% Part 0 - Experiment setup
signals = load_cylinder1dm_signals('evaluation');

fprintf('========== Demo3: F(z)+FIR + IMC-FxNLMS ==========\n');
fprintf('  Evaluation split: seed %d, T1/T2 held out\n', signals.rng_seed);
fprintf('  Information: source reference + secondary-path model\n');
fprintf('  Actuator: hard limit 5, tuning demand limit 4\n\n');

%% Part 1 - T1 fixed: F(z)+FIR at 357 Hz
paramsT1f = struct('N_fir', 8, 'f_design', 357);
resultT1f = controller_demo3(signals, 'T1', 'fixed', paramsT1f);

%% Part 2 - T1 adaptive: frozen IMC-FxNLMS candidate
paramsT1a = struct('N_fir', 8, 'N_nlms', 32, 'mu_nlms', 0.003, ...
    'f_design', 357, 'adaptive_structure', 'imc_fxnlms', ...
    'theta_norm_limit', 10);
resultT1a = controller_demo3(signals, 'T1', 'adaptive', paramsT1a);

%% Part 3 - T2 adaptive: IMC-FxNLMS on a reverse chirp
paramsT2 = struct('N_fir', 8, 'N_nlms', 32, 'mu_nlms', 0.01, ...
    'f_design', 357, 'adaptive_structure', 'imc_fxnlms', ...
    'theta_norm_limit', 10);
resultT2 = controller_demo3(signals, 'T2', 'adaptive', paramsT2);

%% Part 4 - Frozen evaluation overview
figureDir = fullfile(scriptDir, 'figures', 'cylinder1dm_2k', 'demo3');
overview = struct('T1_fixed', resultT1f, ...
    'T1_adaptive', resultT1a, 'T2', resultT2);
plot_demo_overview(overview, fullfile(figureDir, 'demo3_analysis.png'));

%% Part 5 - Summary
fprintf('\n========== Demo3 summary ==========\n');
fprintf('  T1 fixed F(z)+FIR:   %6.2f dB  demand %.3f  clips %d\n', ...
    resultT1f.supp_db, resultT1f.extra.u_demand_max, ...
    resultT1f.extra.saturation_count);
fprintf('  T1 adaptive NLMS:    %6.2f dB  demand %.3f  clips %d  t_conv %.2f s\n', ...
    resultT1a.supp_db, resultT1a.extra.u_demand_max, ...
    resultT1a.extra.saturation_count, resultT1a.t_conv_s);
fprintf('  T2 adaptive NLMS:    %6.2f dB  demand %.3f  clips %d  t_conv %.2f s\n', ...
    resultT2.supp_db, resultT2.extra.u_demand_max, ...
    resultT2.extra.saturation_count, resultT2.t_conv_s);
