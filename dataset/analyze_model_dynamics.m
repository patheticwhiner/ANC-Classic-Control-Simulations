function profiles = analyze_model_dynamics()
%ANALYZE_MODEL_DYNAMICS Compute comparable dynamics for registered LTI models.
%   The result separates a physical delay from numerator zeros and labels
%   marginal roots independently from unstable/non-minimum-phase roots.

ids = { ...
    'armax_30303022', ...
    'syn_TAC2015_3rd', ...
    'syn_JVC2017_3rd', ...
    'syn_JVC2017_6th', ...
    'syn_Bai1997_4th', ...
    'syn_Carmona2000_7th', ...
    'syn_MassSpringDamper_2nd', ...
    'syn_Ho2020_ALE', ...
    'syn_RSTtoy_2nd', ...
    'syn_whitenoise', ...
    'syn_bpf'};

profiles = repmat(empty_profile(), 0, 1);
for i = 1:numel(ids)
    info = DataManager(ids{i});
    switch ids{i}
        case 'armax_30303022'
            if isempty(info.model)
                profiles(end+1) = documented_armax_profile(info); %#ok<AGROW>
            else
                delay = info.orders(4);
                numerator = double(info.model.B);
                numerator = numerator(delay+1:end);
                system = tf(numerator, double(info.model.A), 1/info.fs, ...
                    'Variable', 'z^-1');
                profiles(end+1) = profile_lti(ids{i}, 'plant', system, delay); %#ok<AGROW>
            end

        case 'syn_Ho2020_ALE'
            delayP = leading_zeros(info.P_num);
            delayS = leading_zeros(info.S_num);
            systemP = tf(info.P_num(delayP+1:end), info.P_den, info.Ts, ...
                'Variable', 'z^-1');
            systemS = tf(info.S_num(delayS+1:end), info.S_den, info.Ts, ...
                'Variable', 'z^-1');
            profiles(end+1) = profile_lti([ids{i} ':P'], ...
                'primary path', systemP, delayP); %#ok<AGROW>
            profiles(end+1) = profile_lti([ids{i} ':S'], ...
                'secondary path', systemS, delayS); %#ok<AGROW>

        case 'syn_whitenoise'
            system = ss(info.Aw, info.Bw, info.Cw, 0, 1/info.fs);
            profiles(end+1) = profile_lti(ids{i}, ...
                'disturbance model', system, 0); %#ok<AGROW>

        case 'syn_bpf'
            system = ss(info.Af, info.Bf, info.Cf, 0, 1/info.fs);
            profiles(end+1) = profile_lti(ids{i}, ...
                'plant', system, 0); %#ok<AGROW>

        case 'syn_MassSpringDamper_2nd'
            profiles(end+1) = profile_lti(ids{i}, ...
                'plant (position output)', info.G0(1,1), 0); %#ok<AGROW>

        otherwise
            delay = 0;
            if isfield(info, 'delay_samples'), delay = info.delay_samples; end
            profiles(end+1) = profile_lti(ids{i}, 'plant', info.G0, delay); %#ok<AGROW>
    end
end
end

function result = empty_profile()
result = struct( ...
    'id', '', 'role', '', 'domain', '', 'sample_rate_hz', NaN, ...
    'analysis_basis', 'computed from current artifact', ...
    'poles', [], 'zeros', [], 'pole_count', 0, 'zero_count', 0, ...
    'unstable_poles', 0, 'marginal_poles', 0, ...
    'nmp_zeros', 0, 'boundary_zeros', 0, ...
    'stability', '', 'phase_class', '', 'delay_samples', 0, ...
    'delay_seconds', NaN, 'dominant_pole_value', NaN, ...
    'stability_edge', NaN, 'zero_edge', NaN, ...
    'dominant_pole_frequency_hz', NaN, 'dominant_pole_damping', NaN, ...
    'response_peak', NaN, 'response_peak_db', NaN, ...
    'response_peak_frequency_hz', NaN, 'response_scan', '', ...
    'dc_gain', NaN, 'dc_gain_db', NaN);
end

function result = profile_lti(id, role, system, delaySamples)
result = empty_profile();
result.id = id;
result.role = role;
result.delay_samples = delaySamples;

Ts = system.Ts;
isDiscrete = Ts > 0;
if isDiscrete
    result.domain = 'discrete';
    result.sample_rate_hz = 1/Ts;
    result.delay_seconds = delaySamples * Ts;
else
    result.domain = 'continuous';
end

result.poles = pole(system);
result.zeros = zero(system);
result.pole_count = numel(result.poles);
result.zero_count = numel(result.zeros);

tolerance = 1e-7;
if isDiscrete
    poleMetric = abs(result.poles);
    zeroMetric = abs(result.zeros);
    if ~isempty(poleMetric), result.stability_edge = max(poleMetric); end
    if ~isempty(zeroMetric), result.zero_edge = max(zeroMetric); end
    result.unstable_poles = sum(poleMetric > 1+tolerance);
    result.marginal_poles = sum(abs(poleMetric-1) <= tolerance);
    result.nmp_zeros = sum(zeroMetric > 1+tolerance);
    result.boundary_zeros = sum(abs(zeroMetric-1) <= tolerance);
else
    poleMetric = real(result.poles);
    zeroMetric = real(result.zeros);
    if ~isempty(poleMetric), result.stability_edge = max(poleMetric); end
    if ~isempty(zeroMetric), result.zero_edge = max(zeroMetric); end
    result.unstable_poles = sum(poleMetric > tolerance);
    result.marginal_poles = sum(abs(poleMetric) <= tolerance);
    result.nmp_zeros = sum(zeroMetric > tolerance);
    result.boundary_zeros = sum(abs(zeroMetric) <= tolerance);
end

if result.unstable_poles > 0
    result.stability = 'unstable';
elseif result.marginal_poles > 0
    result.stability = 'marginal';
else
    result.stability = 'stable';
end
if result.nmp_zeros > 0
    result.phase_class = 'non-minimum phase';
elseif result.boundary_zeros > 0
    result.phase_class = 'boundary zero';
else
    result.phase_class = 'minimum phase';
end

[result.dominant_pole_value, result.dominant_pole_frequency_hz, ...
    result.dominant_pole_damping] = dominant_pole(result.poles, Ts);

[result.response_peak, result.response_peak_frequency_hz, ...
    result.response_scan] = response_peak(system);
if isfinite(result.response_peak) && result.response_peak > 0
    result.response_peak_db = 20*log10(result.response_peak);
end

try
    result.dc_gain = max(abs(dcgain(system)), [], 'all');
    if isfinite(result.dc_gain) && result.dc_gain > 0
        result.dc_gain_db = 20*log10(result.dc_gain);
    elseif result.dc_gain == 0
        result.dc_gain_db = -Inf;
    end
catch
    result.dc_gain = NaN;
end
end

function result = documented_armax_profile(info)
% Current tracked MAT has an empty model object. Preserve the last audited
% dynamics from dataset/About_Plant_Performance.md and label the basis.
result = empty_profile();
result.id = 'armax_30303022';
result.role = 'plant';
result.domain = 'discrete';
result.sample_rate_hz = info.fs;
result.analysis_basis = 'documented snapshot; current model object is empty';
result.pole_count = 30;
result.zero_count = 29;
result.unstable_poles = 0;
result.marginal_poles = 0;
result.nmp_zeros = 9;
result.boundary_zeros = 0;
result.stability = 'stable';
result.phase_class = 'non-minimum phase';
result.delay_samples = info.orders(4);
result.delay_seconds = result.delay_samples / info.fs;
result.stability_edge = 0.9958;
result.zero_edge = 1.2402;
result.dominant_pole_frequency_hz = 334.6;
result.response_peak = 10^(2.5/20);
result.response_peak_db = 2.5;
result.response_peak_frequency_hz = 334.6;
result.response_scan = 'documented 0-5000 Hz analysis';
end

function [value, frequencyHz, damping] = dominant_pole(poles, Ts)
value = NaN;
frequencyHz = NaN;
damping = NaN;
if isempty(poles), return; end

if Ts > 0
    [~, index] = max(abs(poles));
    value = poles(index);
    if abs(value) < eps, return; end
    continuousPole = log(value) / Ts;
else
    [~, index] = max(real(poles));
    value = poles(index);
    continuousPole = value;
end

frequencyHz = abs(imag(continuousPole)) / (2*pi);
if abs(continuousPole) > eps
    damping = -real(continuousPole) / abs(continuousPole);
end
end

function [peak, frequencyHz, scanDescription] = response_peak(system)
peak = NaN;
frequencyHz = NaN;
Ts = system.Ts;
if Ts > 0
    nyquist = 1/(2*Ts);
    maximumHz = nyquist;
    if nyquist > 5000, maximumHz = 5000; end
    frequencyHzVector = linspace(0, maximumHz, 8193);
    omega = 2*pi*frequencyHzVector;
    scanDescription = sprintf('0-%g Hz', maximumHz);
else
    rootsScale = abs([pole(system); zero(system)]);
    rootsScale = rootsScale(isfinite(rootsScale) & rootsScale > 1e-8);
    if isempty(rootsScale)
        omega = logspace(-3, 3, 8192);
    else
        omega = logspace(log10(min(rootsScale)/100), ...
            log10(max(rootsScale)*100), 8192);
    end
    frequencyHzVector = omega/(2*pi);
    scanDescription = sprintf('%.3g-%.3g rad/s', omega(1), omega(end));
end

try
    response = freqresp(system, omega);
    magnitude = reshape(abs(response), [], numel(omega));
    magnitude = max(magnitude, [], 1).';
    [peak, index] = max(magnitude);
    frequencyHz = frequencyHzVector(index);
catch
    peak = NaN;
end
end

function count = leading_zeros(coefficients)
threshold = max(abs(coefficients)) * 1e-12;
first = find(abs(coefficients) > threshold, 1, 'first');
if isempty(first)
    count = numel(coefficients);
else
    count = first - 1;
end
end
