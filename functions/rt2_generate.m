function [rf_data, rf_data_full] = rt2_generate(P, mode)
    % [rf_data, rf_data_full] = rt2_generate(P, mode)
    % Generate aberrated RF data for a given set of parameters.
    %
    % Computes the time of flight through the aberrating layers,
    %   taking into account the speed of sound per layers, and
    %   imposing Snell's law on each of the tissue interfaces.
    %
    % The generated RF data amplitude per channel is based on:
    % - The incident angle at the transducer element, to take
    %       into account the element directivity.
    % - The transmission factor at each tissue interface.
    % - The attenuation in each layer of the medium, based
    %       on the frequency and the distance travelled in each layer.
    %
    % Optionally, the wavefront can be further aberrated by local speed
    % of sound variations, defined as a fraction of variation in P.medium_sos_variation.
    % To introduce 5% heterogeneity, set P.medium_sos_variation = 0.05.
    %
    if ~exist('mode', 'var'); mode = 'AC'; end

    %% compute rays and time of flight
    data = rt2_compute(P, mode);

    %% Element dictivity
    pitch = P.lambda;
    if isfield(P, 'pitch'); pitch = P.pitch; end
    W = pitch * 0.9; % small kerf
    directivity_fun = @(theta) cos(theta) .* sinc(pi * W / pitch .* sin(theta));

    %% Total attenuation along path
    att_fun = @(rays) sum(P.medium_attenuation .* [rays.length] * 1e2 * (P.Fc / 1e6));
    attenuation_tx = att_fun(data.rays_tx);
    attenuation_rx = cellfun(att_fun, data.rays_rx);
    attenuation_rt = attenuation_rx + attenuation_tx;
    attenuation_rt_mag = db2mag(-attenuation_rt);

    TF_tx = rt.raytheory.transmission_coeff(P, data.rays_tx, true);
    TF_rx = cellfun(@(rays) rt.raytheory.transmission_coeff(P, rays, false), data.rays_rx);
    TF_rt = TF_tx * TF_rx;

    %% add some wavespeed variation
    sos_vals = P.medium_soundspeeds .* (1 + P.medium_sos_variation .* randn(size(P.medium_soundspeeds)));
    sos_vals(:, 1) = P.medium_soundspeeds(1);

    ray_length_rx = cellfun(@(rays) [rays.length], data.rays_rx, 'UniformOutput', false);
    ray_length_rx = cat(1, ray_length_rx{:});
    sos_vals_rx = P.medium_soundspeeds .* (1 + P.medium_sos_variation .* randn(size(ray_length_rx)));
    sos_vals_rx(:, 1) = P.medium_soundspeeds(1);

    delay_in_lens = data.delay_in_lens;

    tof_tx = sum([data.rays_tx.length] ./ sos_vals) - delay_in_lens;
    tof_rx = sum(ray_length_rx ./ sos_vals_rx, 2);
    tof_round_trip_aberrated = tof_tx + tof_rx;
    tof_round_trip_aberrated = tof_round_trip_aberrated(:).';
    local_aberrration_time = tof_round_trip_aberrated - data.tof_round_trip;

    %% Generate RF
    P.Fs = 20 * P.Fc; % Sampling frequency
    P.num_cycles = 3;

    % Set the impulse response and excitation of the emit aperture
    image_pulse = sin(2 * pi * P.Fc * (0:1 / P.Fs:P.num_cycles / P.Fc)); % BP66 %,  sin(2*pi*fc*(0:1/fs:1/fc)); % BP100 %
    image_pulse = image_pulse .* hamming(length(image_pulse))';
    image_pulse_env = envelope(image_pulse);

    [peak_pulse, imax] = max(image_pulse_env);
    ttp = imax / P.Fs;

    % time vector for delaying
    t_min = 0;
    t_max = 20 / P.Fc;
    t_vec = t_min:1 / P.Fs:t_max;
    N = numel(t_vec);

    DT = 0; % min(local_aberrration_time);
    RF_pulse = repmat(image_pulse, P.num_elements, 1).';
    RF_pulse(N, 1) = 0; % pad with zeros
    if P.medium_sos_variation > 0
        RF_pulse = rt.util.pulse_delaying_RF(RF_pulse, local_aberrration_time + DT, P.Fs);
    end

    % apply 'real' attenuation etc.
    RF_pulse = RF_pulse .* attenuation_rt_mag; % attenuation
    RF_pulse = RF_pulse .* TF_rt; % transmission factors
    RF_pulse = RF_pulse .* directivity_fun(data.theta_rx); % element directivity

    rf_data.RF = RF_pulse; % RF data
    rf_data.ttp = ttp; % time to peak
    rf_data.t_vec = t_vec; % time vector
    rf_data.t_min = t_min; % min time
    rf_data.t_max = t_max; % max time

    %% full RF based on time of flight mode
    if P.medium_sos_variation > 0
        tof_round_trip = tof_round_trip_aberrated(:).';
    else
        tof_round_trip = data.tof_round_trip;
    end

    t_max = max(tof_round_trip) * 1.2;
    t_vec = t_min:1 / P.Fs:t_max;
    N = numel(t_vec);
    RF_pulse = repmat(image_pulse, P.num_elements, 1).';
    RF_pulse(N, 1) = 0; % pad with zeros

    % apply 'real' attenuation etc.
    RF_pulse = RF_pulse .* attenuation_rt_mag; % attenuation
    RF_pulse = RF_pulse .* TF_rt; % transmission factors
    RF_pulse = RF_pulse .* directivity_fun(data.theta_rx); % element directivity

    % delay with NC
    RF_delayed_full = rt.util.pulse_delaying_RF(RF_pulse, tof_round_trip - ttp, P.Fs);

    rf_data_full.RF = RF_delayed_full; % RF data
    rf_data_full.ttp = ttp; % time to peak
    rf_data_full.t_vec = t_vec; % time vector
    rf_data_full.t_min = t_min; % min time
    rf_data_full.t_max = t_max; % max time

end
