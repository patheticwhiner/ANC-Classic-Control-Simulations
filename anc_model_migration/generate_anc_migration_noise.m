function [source, info] = generate_anc_migration_noise(noise, fs, split_name)
%GENERATE_ANC_MIGRATION_NOISE Create deterministic calibration/evaluation input.

if nargin < 3 || isempty(split_name), split_name = 'evaluation'; end
validateattributes(fs, {'numeric'}, {'scalar','real','finite','positive'});
required = {'type','duration','amplitude_rms'};
for index = 1:numel(required)
    if ~isfield(noise, required{index})
        error('generate_anc_migration_noise:missingSetting', ...
            'noise.%s is required.', required{index});
    end
end

split_name = validatestring(split_name, {'calibration','evaluation'});
if strcmp(split_name, 'calibration')
    seed = setting(noise, 'seed_calibration', 42);
else
    seed = setting(noise, 'seed_evaluation', 142);
end
previous_rng = rng;
restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
rng(seed, 'twister');

duration = double(noise.duration);
N = max(1, round(duration*fs));
t = (0:N-1)'/fs;
kind = lower(strrep(char(noise.type), '-', '_'));
frequency_hz = setting(noise, 'frequency_hz', 357);
low_hz = setting(noise, 'low_hz', 300);
high_hz = setting(noise, 'high_hz', 420);
frequencies_hz = double(setting(noise, 'frequencies_hz', frequency_hz));
frequencies_hz = frequencies_hz(:).';
direction = 'none';
f_inst = [];

switch kind
    case {'fixed_sine','sine','tone'}
        validate_frequency(frequency_hz, fs, 'frequency_hz');
        phase = 2*pi*rand();
        source = sin(2*pi*frequency_hz*t + phase) + 0.005*randn(N,1);
        kind = 'fixed_sine';

    case {'linear_chirp','chirp'}
        validate_band(low_hz, high_hz, fs);
        reverse = strcmp(split_name, 'evaluation') ...
            && logical(setting(noise, 'reverse_evaluation_chirp', true));
        if reverse
            start_hz = high_hz; stop_hz = low_hz; direction = 'down';
        else
            start_hz = low_hz; stop_hz = high_hz; direction = 'up';
        end
        rate = (stop_hz-start_hz)/duration;
        phase = 2*pi*(start_hz*t + 0.5*rate*t.^2) + 2*pi*rand();
        source = sin(phase);
        f_inst = start_hz + rate*t;
        kind = 'linear_chirp';

    case {'bandlimited_noise','band_limited','band'}
        validate_band(low_hz, high_hz, fs);
        filter_order = round(setting(noise, 'filter_order', 128));
        validateattributes(filter_order, {'numeric'}, {'scalar','integer','positive'});
        maximum_order = max(8, 2*floor((N-1)/6));
        filter_order = min(filter_order, maximum_order);
        if mod(filter_order, 2) ~= 0, filter_order = filter_order + 1; end
        b = fir1(filter_order, [low_hz high_hz]/(fs/2), 'bandpass');
        source = filter(b, 1, randn(N,1));
        kind = 'bandlimited_noise';

    case {'white','white_noise'}
        source = randn(N,1);
        low_hz = 0;
        high_hz = fs/2;
        kind = 'white_noise';

    case 'multisine'
        if isempty(frequencies_hz)
            error('generate_anc_migration_noise:emptyFrequencies', ...
                'At least one multisine frequency is required.');
        end
        for index = 1:numel(frequencies_hz)
            validate_frequency(frequencies_hz(index), fs, 'frequencies_hz');
        end
        source = zeros(N,1);
        phases = 2*pi*rand(size(frequencies_hz));
        for index = 1:numel(frequencies_hz)
            source = source + sin(2*pi*frequencies_hz(index)*t + phases(index));
        end
        low_hz = min(frequencies_hz);
        high_hz = max(frequencies_hz);

    case 'file'
        file_path = char(setting(noise, 'file', ''));
        if ~isfile(file_path)
            error('generate_anc_migration_noise:fileNotFound', ...
                'Noise file not found: %s', file_path);
        end
        [source, source_fs] = audioread(file_path);
        source = source(:,1);
        if source_fs ~= fs, source = resample(source, fs, source_fs); end
        if numel(source) < N
            error('generate_anc_migration_noise:fileTooShort', ...
                'Noise file contains %.3f s but %.3f s is required.', ...
                numel(source)/fs, duration);
        end
        source = source(1:N);
        low_hz = NaN; high_hz = NaN;

    otherwise
        error('generate_anc_migration_noise:unknownType', ...
            'Unsupported noise type: %s', noise.type);
end

source = double(source(:));
source = source-mean(source);
source_rms = sqrt(mean(source.^2));
if source_rms <= eps
    error('generate_anc_migration_noise:zeroSignal', ...
        'Generated source signal has zero RMS.');
end
source = source/source_rms;

info = struct('type', kind, 'split', split_name, 'fs', fs, ...
    'samples', N, 'duration', N/fs, 'seed', seed, ...
    'frequency_hz', frequency_hz, 'low_hz', low_hz, ...
    'high_hz', high_hz, 'frequencies_hz', frequencies_hz, ...
    'direction', direction, 'f_inst', f_inst, 'time', t);
end

function value = setting(settings, name, default_value)
if isfield(settings, name) && ~isempty(settings.(name))
    value = settings.(name);
else
    value = default_value;
end
end

function validate_frequency(value, fs, name)
if ~isscalar(value) || ~isfinite(value) || value <= 0 || value >= fs/2
    error('generate_anc_migration_noise:invalidFrequency', ...
        '%s must lie strictly between 0 and Nyquist.', name);
end
end

function validate_band(low_hz, high_hz, fs)
if ~isscalar(low_hz) || ~isscalar(high_hz) || ...
        ~isfinite(low_hz) || ~isfinite(high_hz) || ...
        low_hz <= 0 || high_hz <= low_hz || high_hz >= fs/2
    error('generate_anc_migration_noise:invalidBand', ...
        'The band must satisfy 0 < low_hz < high_hz < Nyquist.');
end
end
