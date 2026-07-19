function candidates = build_anc_migration_candidates(controller_id, fs, noise, design)
%BUILD_ANC_MIGRATION_CANDIDATES Sampling-rate-aware controller candidates.

controller_id = char(controller_id);
profile = lower(char(design.profile));
if ~ismember(profile, {'quick','standard'})
    error('build_anc_migration_candidates:unknownProfile', ...
        'Unknown tuning profile: %s', profile);
end

scale = fs/2000;
frequency = double(noise.frequency_hz);
center = mean([double(noise.low_hz), double(noise.high_hz)]);
if any(strcmpi(noise.type, {'linear_chirp','chirp'})), frequency = center; end
lambda = @(value) value.^(2000/fs);
taps = @(value) max(2, round(value*scale));
mu = @(value) value*min(1, 2000/fs);
actuator = design.actuator_limit;
candidates = repmat(struct('name', '', 'params', struct()), 0, 1);

switch controller_id
    case 'demo1_rst_fixed'
        candidates = add(candidates, 'rst', struct( ...
            'f_design', frequency, 'actuator_limit', actuator));

    case 'demo1_qrls'
        orders = unique(taps([2 4 8]));
        lambdas = lambda([0.98 0.995]);
        if strcmp(profile, 'quick'), orders = orders(1:min(2,end)); lambdas = lambdas(1); end
        for nQ = orders
            for forgetting = lambdas
                name = sprintf('q%d_l%.5f', nQ, forgetting);
                candidates = add(candidates, name, struct('f_design', frequency, ...
                    'nQ', nQ, 'lambda1', forgetting, 'lambda2', 1, ...
                    'F_diag', 1e-3, 'actuator_limit', actuator));
            end
        end

    case 'demo2_lqg_fixed'
        r_values = [1e-4 1e-3 1e-2];
        if strcmp(profile, 'quick'), r_values = [1e-3 1e-2]; end
        for R = r_values
            candidates = add(candidates, sprintf('lqg_r%g', R), ...
                base_lqg(R, actuator));
        end

    case {'demo2_lms','demo2_qrls'}
        orders = unique(taps([16 32]));
        gains = [0.05 0.1];
        if strcmp(profile, 'quick'), orders = orders(end); gains = gains(end); end
        for nQ = orders
            for gain = gains
                params = base_lqg(0.01, actuator);
                params.nQ = nQ;
                params.bp_band = valid_band(noise, fs);
                params.adaptive_warmup_seconds = min(0.5, 0.2*noise.duration);
                params.theta_norm_limit = 10;
                params.adaptation_gain = gain;
                if strcmp(controller_id, 'demo2_lms')
                    params.adaptation_method = 'lms';
                    params.lms_epsilon = 1e-6;
                    params.lms_leakage = 0;
                    name = sprintf('lms_n%d_mu%.3g', nQ, gain);
                else
                    params.adaptation_method = 'rls';
                    params.lambda = lambda(0.98);
                    params.F_init = 0.1;
                    name = sprintf('rls_n%d_g%.3g', nQ, gain);
                end
                candidates = add(candidates, name, params);
            end
        end

    case 'demo3_fir_fixed'
        lengths = unique(taps([8 16 24]));
        if strcmp(profile, 'quick'), lengths = lengths(1:min(2,end)); end
        for length_value = lengths
            candidates = add(candidates, sprintf('fir%d', length_value), struct( ...
                'N_fir', length_value, 'f_design', frequency, ...
                'actuator_limit', actuator, 'ramp_seconds', 0.1));
        end

    case 'demo3_imc_fxnlms'
        lengths = unique(taps([16 32]));
        gains = mu([0.003 0.01]);
        if strcmp(profile, 'quick'), lengths = lengths(end); gains = gains(end); end
        for length_value = lengths
            for gain = gains
                params = struct('N_fir', taps(8), 'N_nlms', length_value, ...
                    'mu_nlms', gain, 'f_design', frequency, ...
                    'adaptive_structure', 'imc_fxnlms', ...
                    'theta_norm_limit', 10, 'actuator_limit', actuator, ...
                    'ramp_seconds', 0.1);
                candidates = add(candidates, ...
                    sprintf('imc_n%d_mu%.4g', length_value, gain), params);
            end
        end

    case {'demo4_emopso_rst','demo4_emopso_qrls'}
        if strcmp(profile, 'quick')
            population = 16; iterations = 20;
        else
            population = 40; iterations = 80;
        end
        notch = frequency;
        if any(strcmpi(noise.type, {'bandlimited_noise','band_limited','white','white_noise','file'}))
            notch = [];
        end
        central = struct('nX', 6, 'nY', 3, 'f_notch', notch, ...
            'n_pop', population, 'k_max', iterations, 'seed', 42, ...
            'score_band', valid_band(noise, fs), ...
            'P_D_omega', min(0.45*pi, 0.25*pi*2000/fs), ...
            'actuator_limit', actuator);
        if strcmp(controller_id, 'demo4_emopso_rst')
            candidates = add(candidates, 'emopso_rst', central);
        else
            q_orders = unique(taps([2 4]));
            if strcmp(profile, 'quick'), q_orders = q_orders(end); end
            for nQ = q_orders
                params = central;
                params.nQ = nQ;
                params.lambda1 = lambda(0.98);
                params.lambda2 = 1;
                params.F_diag = 1e-3;
                candidates = add(candidates, sprintf('emopso_q%d', nQ), params);
            end
        end

    case {'demo5_marino_fixed','demo5_marino_adaptive'}
        frequencies = marino_frequencies(noise);
        gains = [3e-4 1e-3 3e-3];
        if strcmp(profile, 'quick'), gains = gains([1 end]); end
        for gain = gains
            params = struct('q', numel(frequencies), 'k', gain, ...
                'freq_init_hz', frequencies, 'method', 'exact', ...
                'ramp_seconds', min(0.5, 0.1*noise.duration), ...
                'dc_cancel', false, 'output_timing', 'updated', ...
                'actuator_limit', actuator, 'freq_min_hz', 20, ...
                'freq_max_hz', min(0.45*fs, max(500, 1.2*max(frequencies))));
            if strcmp(controller_id, 'demo5_marino_fixed')
                params.epsilon = 0;
            else
                params.epsilon = 1e-4*min(1, 2000/fs);
            end
            candidates = add(candidates, sprintf('mt_k%.3g', gain), params);
        end

    case 'imc_fxlms'
        lengths = unique(taps([32 64]));
        gains = mu([0.005 0.01 0.02]);
        if strcmp(profile, 'quick'), lengths = lengths(end); gains = gains([1 end]); end
        for length_value = lengths
            for gain = gains
                candidates = add(candidates, ...
                    sprintf('imc_n%d_mu%.4g', length_value, gain), struct( ...
                    'N_fir', length_value, 'mu', gain, 'delta', 1e-4, ...
                    'ramp_seconds', min(0.5, 0.1*noise.duration), ...
                    'norm_limit', 10, 'actuator_limit', actuator));
            end
        end

    otherwise
        error('build_anc_migration_candidates:unknownController', ...
            'Unknown controller id: %s', controller_id);
end

if ~design.auto_tune && ~isempty(candidates)
    candidates = candidates(ceil(numel(candidates)/2));
end

if isfield(design, 'manual_params') && ~isempty(fieldnames(design.manual_params))
    candidates = apply_manual_parameters(candidates, controller_id, ...
        design.manual_params);
end
end

function candidates = apply_manual_parameters(candidates, controller_id, manual)
definitions = anc_controller_parameter_definitions(controller_id);
allowed = {definitions.name};
names = fieldnames(manual);
unknown = setdiff(names, allowed);
if ~isempty(unknown)
    error('build_anc_migration_candidates:unknownManualParameter', ...
        'Unsupported manual parameter for %s: %s.', ...
        controller_id, strjoin(unknown, ', '));
end
for field_index = 1:numel(names)
    name = names{field_index};
    definition = definitions(strcmp(allowed, name));
    value = parse_anc_parameter_value(manual.(name), definition);
    for candidate_index = 1:numel(candidates)
        candidates(candidate_index).params.(name) = value;
        if strcmp(name, 'freq_init_hz')
            candidates(candidate_index).params.q = numel(value);
        end
    end
end
if any(strcmp(names, 'freq_min_hz')) || any(strcmp(names, 'freq_max_hz'))
    for candidate_index = 1:numel(candidates)
        p = candidates(candidate_index).params;
        if p.freq_max_hz <= p.freq_min_hz
            error('build_anc_migration_candidates:invalidFrequencyBounds', ...
                'freq_max_hz must be greater than freq_min_hz.');
        end
    end
end
end

function candidates = add(candidates, name, params)
entry = struct('name', name, 'params', params);
candidates(end+1,1) = entry;
end

function params = base_lqg(R, actuator)
params = struct('Q_lqr', 1, 'R_lqr', R, 'Qn_plant', 1e-4, ...
    'Qn_dist', 0.1, 'Rn', 0.1, 'control_scale', 1, ...
    'ramp_seconds', 0.1, 'actuator_limit', actuator);
end

function band = valid_band(noise, fs)
low = max(1, double(noise.low_hz));
high = min(0.45*fs, double(noise.high_hz));
if high <= low
    low = max(1, 0.05*fs);
    high = 0.4*fs;
end
band = [low high];
end

function frequencies = marino_frequencies(noise)
if strcmpi(noise.type, 'multisine')
    frequencies = unique(double(noise.frequencies_hz(:).'));
elseif any(strcmpi(noise.type, {'linear_chirp','chirp'}))
    frequencies = mean([double(noise.low_hz), double(noise.high_hz)]);
else
    frequencies = double(noise.frequency_hz);
end
end
