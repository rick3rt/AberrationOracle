function data = rt2_compute(P, mode)

    switch mode
        case 'NC' % no correction, homogeneous medium
            data = rt2_compute_homogeneous(P);
        case 'LC' % lens corrected
            data = rt2_compute_tof_trace(P, 2);
        case 'AC' % aberration corrected
            to_layer = numel(P.medium_soundspeeds);
            data = rt2_compute_tof_trace(P, to_layer);
        otherwise
            error('Unknown mode: %s', mode)
    end

end

function data = rt2_compute_tof_trace(P, to_layer)

    % point source compensation
    v_source = [P.x_source; P.z_source];
    v_piezo = [P.x_piezo; P.z_piezo];
    dist_se = vecnorm(v_piezo - v_source);
    tx_delay_lens = dist_se / P.medium_soundspeeds(1); % lens speed of sound
    delay_in_lens = min(tx_delay_lens);
    
    % % determine depth per interface??
    % z_ints = cellfun(@(p) rt.util.segeval(p, P.x_pixel), P.medium_interfaces);
    % z_dist = diff([0 z_ints P.z_pixel]);
    % if any(z_dist < 0)
    %     error('Pixel not under last interface!')
    % end
    % tof_tx_rt = sum(z_dist ./ P.medium_soundspeeds); % need to add lens term?
    % z_c0 = tof_tx_rt * P.c0;
    % v_pixel = [P.x_pixel; P.z_pixel];
    
    % DONT 
    % if to_layer < numel(P.medium_soundspeeds)
    %     tmp = rt2_compute_tof_trace(P, numel(P.medium_soundspeeds)); % final 
    %     tau = tmp.tof_round_trip;
    %     taum = min(tau);
    %     tau2 = taum - 2*P.lens_thickness / P.medium_soundspeeds(1);
    %     z2 = tau2/2*P.medium_soundspeeds(2);
    %     P.z_pixel = z2 + P.lens_thickness; 
    % end

    % trace transmit
    [rays_tx, tof_tx, theta_tx] = rt.ray_bending(P, P.x_source, P.z_source, P.x_pixel, P.z_pixel, to_layer);
    tof_tx = tof_tx - delay_in_lens; % correct propagation time in lens

    % trace receive
    tof_rx_all = zeros(1, P.num_elements);
    theta_rx_all = zeros(1, P.num_elements);
    rays_rx_all = cell(1, P.num_elements);
    for ke = 1:P.num_elements
        x_end = P.x_piezo(ke);
        z_end = P.z_piezo(ke);
        [rays_rx, tof_rx, theta_rx] = rt.ray_bending(P, x_end, z_end, P.x_pixel, P.z_pixel, to_layer);
        tof_rx_all(ke) = tof_rx;
        theta_rx_all(ke) = theta_rx;
        rays_rx_all{ke} = rays_rx;
    end

    % valid rays, tof = NaN if final refraction not valid (i.e. critical angle)
    msk_valid = find(~isnan(tof_rx_all)); % valid rays in reception.
    tof_round_trip = tof_tx + tof_rx_all; % aberration corrected time of flight

    % f number mask
    half_opening_angle_rad = atan(1/2 / P.f_number);
    f_number_msk = abs(theta_rx_all) < half_opening_angle_rad;

    % collect results
    data.rays_tx = rays_tx;
    data.rays_rx = rays_rx_all;
    data.tof_tx = tof_tx;
    data.delay_in_lens = delay_in_lens;
    data.tof_rx = tof_rx_all;
    data.tof_round_trip = tof_round_trip;
    data.theta_tx = theta_tx; %
    data.theta_rx = theta_rx_all;
    data.f_number_msk = f_number_msk;
    data.valid_refraction_msk = msk_valid;

end

function data = rt2_compute_homogeneous(P)

    % determine depth per interface
    z_ints = cellfun(@(p) rt.util.segeval(p, P.x_pixel), P.medium_interfaces);
    z_dist = diff([0 z_ints P.z_pixel]);
    if any(z_dist < 0)
        error('Pixel not under last interface!')
    end

    % tof_tx_rt = sum(z_dist ./ P.medium_soundspeeds); % need to add lens term?
    % z_c0 = tof_tx_rt * P.c0;

    v_pixel = [P.x_pixel; P.z_pixel]; % z_c0];
    v_source = [P.x_source; P.z_source];
    v_piezo = [P.x_piezo; P.z_piezo];

    dist_se = vecnorm(v_piezo - v_source);
    tx_delay_c0 = dist_se / P.c0;
    delay_in_lens = min(tx_delay_c0);

    tof_tx = vecnorm(v_pixel - v_source) / P.c0 - delay_in_lens; % correct for point source propagation in lens
    tof_rx = vecnorm(v_piezo - v_pixel) / P.c0;
    tof_round_trip = tof_tx + tof_rx;

    v_rx = v_pixel - v_piezo; % return vector
    theta_rx_all = atan2(v_rx(2, :), v_rx(1, :)) - pi / 2;

    half_opening_angle_rad = atan(1/2 / P.f_number);
    f_number_msk = abs(theta_rx_all) < half_opening_angle_rad;

    % collect results
    data.tof_tx = tof_tx;
    data.tof_rx = tof_rx;
    data.tof_round_trip = tof_round_trip;
    data.theta_rx = theta_rx_all;
    data.f_number_msk = f_number_msk;
    data.valid_refraction_msk = 1:numel(P.x_piezo); % all rays valid in homogeneous medium

end
