function signals = prepare_anc_migration_signals(models, noise, split_name, control_delay_samples)
%PREPARE_ANC_MIGRATION_SIGNALS Adapt imported models to tests/controllers.

if nargin < 3 || isempty(split_name), split_name = 'evaluation'; end
if nargin < 4 || isempty(control_delay_samples), control_delay_samples = 1; end
validateattributes(control_delay_samples, {'numeric'}, ...
    {'scalar','integer','nonnegative'});

fs = models.fs;
[source, noise_info] = generate_anc_migration_noise(noise, fs, split_name);
d = filter(models.primary.B, models.primary.A, source);
disturbance_rms = sqrt(mean(d.^2));
if disturbance_rms <= eps
    error('prepare_anc_migration_signals:zeroPrimaryOutput', ...
        'The selected primary model produced a zero-RMS disturbance.');
end
target_rms = double(noise.amplitude_rms);
validateattributes(target_rms, {'numeric'}, {'scalar','real','finite','positive'});
scale = target_rms/disturbance_rms;
source = source*scale;
d = d*scale;

if strcmp(models.secondary.kind, 'fir')
    B_secondary = [zeros(1, control_delay_samples), models.secondary.B];
    orders = [0, numel(B_secondary)-1, 0, control_delay_samples];
    path_family = 'measured_lms_fir';
else
    B_secondary = models.secondary.B;
    orders = models.secondary.orders;
    path_family = 'measured_armax';
end

test = make_test(noise_info, noise, source, d, fs);
signals = struct();
signals.model = models.id;
signals.primary_model = models.id;
signals.fs = fs;
signals.rng_seed = noise_info.seed;
signals.norm_target = 'rms_d';
signals.norm_value = target_rms;
signals.model_data = struct('A', models.secondary.A, ...
    'B', B_secondary, 'orders', orders, 'fs', fs, ...
    'path_family', path_family, ...
    'source_file', models.secondary.file, ...
    'control_delay_samples', control_delay_samples);
signals.T1 = test;
signals.meta = struct('split', split_name, 'source', source, ...
    'noise', noise_info, 'models', models, ...
    'primary_output_rms', sqrt(mean(d.^2)));
end

function test = make_test(info, noise, source, d, fs)
test = struct('type', info.type, 'fs', fs, 'Tsim', info.duration, ...
    't', info.time, 'source', source, 'd', d, 'y_open', d, ...
    'f_hz', double(noise.frequency_hz));

switch info.type
    case 'fixed_sine'
        test.f_hz = info.frequency_hz;
    case 'linear_chirp'
        test.f_range = [info.low_hz info.high_hz];
        test.f_inst = info.f_inst;
        test.direction = info.direction;
    case 'bandlimited_noise'
        test.f_range = [info.low_hz info.high_hz];
    case 'white_noise'
        test.f_range = [max(1, 0.01*fs), 0.45*fs];
    case 'multisine'
        test.frequencies_hz = info.frequencies_hz;
        test.f_range = [min(info.frequencies_hz), max(info.frequencies_hz)];
    case 'file'
        test.f_range = [max(1, 0.01*fs), 0.45*fs];
end
end
