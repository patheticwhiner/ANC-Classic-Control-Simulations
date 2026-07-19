function result = run_anc_migration_controller(controller_id, signals, params)
%RUN_ANC_MIGRATION_CONTROLLER Dispatch one controller on a prepared signal.

switch char(controller_id)
    case 'demo1_rst_fixed'
        result = controller_demo1(signals, 'T1', 'fixed', params);
    case 'demo1_qrls'
        result = controller_demo1(signals, 'T1', 'adaptive', params);
    case 'demo2_lqg_fixed'
        result = controller_demo2(signals, 'T1', 'fixed', params);
    case {'demo2_lms','demo2_qrls'}
        result = controller_demo2(signals, 'T1', 'adaptive', params);
    case 'demo3_fir_fixed'
        result = controller_demo3(signals, 'T1', 'fixed', params);
    case 'demo3_imc_fxnlms'
        result = controller_demo3(signals, 'T1', 'adaptive', params);
    case 'demo4_emopso_rst'
        result = controller_demo4(signals, 'T1', 'fixed', params);
    case 'demo4_emopso_qrls'
        result = controller_demo4(signals, 'T1', 'adaptive', params);
    case 'demo5_marino_fixed'
        result = controller_demo5(signals, 'T1', 'fixed', params);
    case 'demo5_marino_adaptive'
        result = controller_demo5(signals, 'T1', 'adaptive', params);
    case 'imc_fxlms'
        result = controller_imc_fxlms(signals, 'T1', 'adaptive', params);
    otherwise
        error('run_anc_migration_controller:unknownController', ...
            'Unknown controller id: %s', controller_id);
end
end
