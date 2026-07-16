% DEMO1_CYLINDER1DM Low-order RST and Youla-Q RLS experiment.
%
% T1 compares fixed RST with the same central controller plus Youla-Q RLS.
% T2 evaluates frequency tracking on the held-out reverse chirp. Broadband
% T3 is run independently by run_cylinder1dm_stage('broadband').

clear; close all; clc;
scriptDir = fileparts(mfilename('fullpath'));
run(fullfile(scriptDir, 'startup.m'));

%% Part 0 - Experiment setup
signals = load_cylinder1dm_signals('evaluation');

fprintf('========== Demo1: low-order RST + Youla-Q RLS ==========\n');
fprintf('  Evaluation split: seed %d, T1/T2 held out\n', signals.rng_seed);
fprintf('  Actuator: hard limit 5, tuning demand limit 4\n\n');

%% Part 1 - T1 fixed: low-order RST at 357 Hz
fprintf('--- T1 fixed RST ---\n');
paramsT1f = struct('f_design', 357);
resultT1f = controller_demo1(signals, 'T1', 'fixed', paramsT1f);
fprintf('  suppression = %.2f dB | demand = %.3f | clips = %d\n', ...
    resultT1f.supp_db, resultT1f.extra.u_demand_max, ...
    resultT1f.extra.saturation_count);

%% Part 2 - T1 adaptive: same central R0/S0 + Youla-Q RLS
fprintf('--- T1 adaptive Youla-Q RLS ---\n');
paramsT1a = struct('nQ', 2, 'lambda1', 0.995, 'F_diag', 1e-3);
resultT1a = controller_demo1(signals, 'T1', 'adaptive', paramsT1a);
fprintf('  suppression = %.2f dB | demand = %.3f | clips = %d | t_conv = %.2f s\n', ...
    resultT1a.supp_db, resultT1a.extra.u_demand_max, ...
    resultT1a.extra.saturation_count, resultT1a.t_conv_s);

%% Part 3 - T2 adaptive: Youla-Q RLS on a held-out reverse chirp
fprintf('--- T2 adaptive Youla-Q RLS ---\n');
paramsT2 = struct('nQ', 4, 'lambda1', 0.98, 'F_diag', 1e-3);
resultT2 = controller_demo1(signals, 'T2', 'adaptive', paramsT2);
fprintf('  suppression = %.2f dB | demand = %.3f | clips = %d | t_conv = %.2f s\n', ...
    resultT2.supp_db, resultT2.extra.u_demand_max, ...
    resultT2.extra.saturation_count, resultT2.t_conv_s);

%% Part 4 - Frozen evaluation overview
figureDir = fullfile(scriptDir, 'figures', 'cylinder1dm_2k', 'demo1');
overview = struct('T1_fixed', resultT1f, ...
    'T1_adaptive', resultT1a, 'T2', resultT2);
plot_demo_overview(overview, fullfile(figureDir, 'demo1_analysis.png'));

%% Part 5 - Summary
fprintf('\n========== Demo1 summary ==========\n');
fprintf('  T1 fixed RST:         %6.2f dB  demand %.3f  clips %d\n', ...
    resultT1f.supp_db, resultT1f.extra.u_demand_max, ...
    resultT1f.extra.saturation_count);
fprintf('  T1 adaptive Youla-Q:  %6.2f dB  demand %.3f  clips %d  t_conv %.2f s\n', ...
    resultT1a.supp_db, resultT1a.extra.u_demand_max, ...
    resultT1a.extra.saturation_count, resultT1a.t_conv_s);
fprintf('  T2 adaptive Youla-Q:  %6.2f dB  demand %.3f  clips %d  t_conv %.2f s\n', ...
    resultT2.supp_db, resultT2.extra.u_demand_max, ...
    resultT2.extra.saturation_count, resultT2.t_conv_s);
