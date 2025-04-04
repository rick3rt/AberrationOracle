function [rf_data, rf_data_full] = rt2_generate(P, mode)

    if ~exist('mode', 'var'); mode = 'AC'; end

    %% % compute rays and time of flight
    data = rt2_compute(P, mode);

    %% Element dictivity
    lambda = P.lambda;
    W = lambda * 0.8; % small kerf
    directivity_fun = @(theta) cos(theta) .* sinc(pi * W / lambda .* sin(theta));

    %% Total attenuation along path
    att_fun = @(rays) sum(P.medium_attenuation .* [rays.length] * 1e2 * (P.Fc / 1e6));
    attenuation_tx = att_fun(data.rays_tx);
    attenuation_rx = cellfun(att_fun, data.rays_rx);
    attenuation_rt = attenuation_rx + attenuation_tx;
    attenuation_rt_mag = db2mag(-attenuation_rt);

    TF_tx = rt.raytheory.transmission_coeff(P, data.rays_tx, true);
    TF_rx = cellfun(@(rays) rt.raytheory.transmission_coeff(P, rays, false), data.rays_rx);
    TF_rt = TF_tx * TF_rx;

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

    RF_pulse = repmat(image_pulse, P.num_elements, 1).';
    RF_pulse(N, 1) = 0; % pad with zeros

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

    t_max = max(data.tof_round_trip) * 1.2;
    t_vec = t_min:1 / P.Fs:t_max;
    N = numel(t_vec);
    RF_pulse = repmat(image_pulse, P.num_elements, 1).';
    RF_pulse(N, 1) = 0; % pad with zeros

    % apply 'real' attenuation etc.
    RF_pulse = RF_pulse .* attenuation_rt_mag; % attenuation
    RF_pulse = RF_pulse .* TF_rt; % transmission factors
    RF_pulse = RF_pulse .* directivity_fun(data.theta_rx); % element directivity

    % delay with NC
    RF_delayed_full = rt.util.pulse_delaying_RF(RF_pulse, data.tof_round_trip - ttp, P.Fs);

    rf_data_full.RF = RF_delayed_full; % RF data
    rf_data_full.ttp = ttp; % time to peak
    rf_data_full.t_vec = t_vec; % time vector
    rf_data_full.t_min = t_min; % min time
    rf_data_full.t_max = t_max; % max time

end
