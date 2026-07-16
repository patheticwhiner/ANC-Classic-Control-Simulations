% DEMO4_CYLINDER1DM Epsilon-MOPSO sensitivity-shaping RST experiment.
%
% T1 compares the optimised central RST with the same controller plus the
% Demo1 Youla-Q RLS law. T2 evaluates that combination on a reverse chirp.

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath'));
run(fullfile(scriptDir, 'startup.m'));

%% Part 0 - Experiment setup
signals = load_cylinder1dm_signals('evaluation');

fprintf('========== Demo4: eMOPSO sensitivity-shaping RST ==========\n');
fprintf('  Evaluation split: seed %d, T1/T2 held out\n', signals.rng_seed);
fprintf('  Adaptive law: same Youla-Q RLS as Demo1\n');
fprintf('  Actuator: hard limit 5, tuning demand limit 4\n\n');

%% Part 1 - T1 fixed: eMOPSO notch RST at 357 Hz
paramsT1f = struct('nX', 10, 'nY', 5, 'f_notch', 357, ...
    'n_pop', 60, 'k_max', 120);
resultT1f = controller_demo4(signals, 'T1', 'fixed', paramsT1f);

%% Part 2 - T1 adaptive: same eMOPSO centre + Youla-Q RLS
paramsT1a = struct('nX', 10, 'nY', 5, 'f_notch', 357, ...
    'n_pop', 60, 'k_max', 120, ...
    'nQ', 4, 'lambda1', 0.995, 'F_diag', 1e-3);
resultT1a = controller_demo4(signals, 'T1', 'adaptive', paramsT1a);

%% Part 3 - T2 adaptive: eMOPSO centre + Youla-Q RLS
paramsT2 = struct('nX', 10, 'nY', 5, 'f_notch', 357, ...
    'n_pop', 60, 'k_max', 120, ...
    'nQ', 8, 'lambda1', 0.98, 'F_diag', 1e-3);
resultT2 = controller_demo4(signals, 'T2', 'adaptive', paramsT2);

%% Part 4 - Frozen evaluation overview
figureDir = fullfile(scriptDir, 'figures', 'cylinder1dm_2k', 'demo4');
overview = struct('T1_fixed', resultT1f, ...
    'T1_adaptive', resultT1a, 'T2', resultT2);
plot_demo_overview(overview, fullfile(figureDir, 'demo4_analysis.png'));

%% Part 5 - Summary
fprintf('\n========== Demo4 summary ==========\n');
fprintf('  T1 fixed eMOPSO RST: %6.2f dB  demand %.3f  clips %d  Ms %.2f dB\n', ...
    resultT1f.supp_db, resultT1f.extra.u_demand_max, ...
    resultT1f.extra.saturation_count, resultT1f.extra.design.Ms_db);
fprintf('  T1 adaptive Q-RLS:   %6.2f dB  demand %.3f  clips %d  t_conv %.2f s\n', ...
    resultT1a.supp_db, resultT1a.extra.u_demand_max, ...
    resultT1a.extra.saturation_count, resultT1a.t_conv_s);
fprintf('  T2 adaptive Q-RLS:   %6.2f dB  demand %.3f  clips %d  t_conv %.2f s\n', ...
    resultT2.supp_db, resultT2.extra.u_demand_max, ...
    resultT2.extra.saturation_count, resultT2.t_conv_s);
