function data = rt2_compute_tof_grid(P, mode, Z, X)

    switch mode
        case 'NC' % no correction, homogeneous medium
            data = rt_compute_grid_homogeneous(P, Z, X);
        case 'LC' % lens corrected
            data = rt_compute_grid_tof_trace(P, 2, Z, X);
        case 'AC' % aberration corrected
            data = rt_compute_grid_tof_trace(P, 4, Z, X);
        otherwise
            error('Unknown mode: %s', mode)
    end

end

function data = rt_compute_grid_tof_trace(P, to_layer, Z, X)

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
    % v_pixel = [P.x_pixel; z_c0];

    Nx = size(X, 2);
    Nz = size(Z, 1);
    tof_round_trip_out = zeros(Nz, Nx, P.num_elements);
    half_opening_angle_rad = atan(1/2 / P.f_number);

    % use 3L recon
    % [data, Pars] = threeLayerReconWrapper(Pars, RF_I, RF_Q)

    % time_remaining_progbar_ui(0, Nx * Nz);
    % for kx = 1:Nx
    %     for kz = 1:Nz

    %         xp = X(kz, kx);
    %         zp = Z(kz, kx);

    %         % trace transmit
    %         [~, tof_tx, ~] = rt.ray_bending(P, P.x_source, P.z_source, xp, zp, to_layer);
    %         tof_tx = tof_tx - delay_in_lens; % correct propagation time in lens

    %         % trace receive
    %         tof_rx_all = zeros(1, P.num_elements);
    %         theta_rx_all = zeros(1, P.num_elements);
    %         rays_rx_all = cell(1, P.num_elements);
    %         for ke = 1:P.num_elements
    %             xe = P.x_piezo(ke);
    %             ze = P.z_piezo(ke);
    %             [rays_rx, tof_rx, theta_rx] = rt.ray_bending(P, xe, ze, xp, zp, to_layer);
    %             tof_rx_all(ke) = tof_rx;
    %             theta_rx_all(ke) = theta_rx;
    %             rays_rx_all{ke} = rays_rx;
    %         end

    %         % valid rays, tof = NaN if final refraction not valid (i.e. critical angle)
    %         % msk_valid = find(~isnan(tof_rx_all)); % valid rays in reception.
    %         tof_round_trip = tof_tx + tof_rx_all; % aberration corrected time of flight

    %         % f number mask
    %         f_number_msk = abs(theta_rx_all) < half_opening_angle_rad;
    %         tof_round_trip(~f_number_msk) = NaN; % mask fnumber for reconstruction
    %         tof_round_trip_out(kz, kx, :) = tof_round_trip; % collect results

    %         time_remaining_progbar_ui(kz + (kx - 1) * Nz, Nx * Nz);
    %     end
    % end

    data.Nz = Nz;
    data.Nx = Nx;
    data.tof_round_trip = tof_round_trip_out;
    data.X = X;
    data.Z = Z;

end

function data = rt_compute_grid_homogeneous(P, Z, X)

    % % determine depth per interface
    % z_ints = cellfun(@(p) rt.util.segeval(p, P.x_pixel), P.medium_interfaces);
    % z_dist = diff([0 z_ints P.z_pixel]);
    % if any(z_dist < 0)
    %     error('Pixel not under last interface!')
    % end

    % tof_tx_rt = sum(z_dist ./ P.medium_soundspeeds); % need to add lens term?
    % z_c0 = tof_tx_rt * P.c0;

    % v_pixel = [P.x_pixel; z_c0];
    v_source = [P.x_source; P.z_source];
    v_piezo = [P.x_piezo; P.z_piezo];

    dist_se = vecnorm(v_piezo - v_source);
    tx_delay_c0 = dist_se / P.c0;
    delay_in_lens = min(tx_delay_c0);
    t_lens_cor = 2 * (1 / P.medium_soundspeeds(1) - 1 / P.c0) * P.lens_thickness;
    Z = Z - P.lens_thickness;
    % P.medium_soundspeeds(1) = lens wavespeed

    % tof_tx = vecnorm(v_pixel - v_source) / P.c0 - delay_in_lens; % correct for point source propagation in lens
    % tof_rx = vecnorm(v_piezo - v_pixel) / P.c0;
    % tof_round_trip = tof_tx + tof_rx;

    % v_rx = v_pixel - v_piezo; % return vector
    % theta_rx_all = atan2(v_rx(2, :), v_rx(1, :)) - pi / 2;

    % half_opening_angle_rad = atan(1/2 / P.f_number);
    % f_number_msk = abs(theta_rx_all) < half_opening_angle_rad;

    % % collect results
    % data.tof_tx = tof_tx;
    % data.tof_rx = tof_rx;
    % data.tof_round_trip = tof_round_trip;
    % data.theta_rx = theta_rx_all;
    % data.f_number_msk = f_number_msk;
    % data.valid_refraction_msk = 1:numel(P.x_piezo); % all rays valid in homogeneous medium

    Nx = size(X, 2); Nz = size(Z, 1);
    tau_tx = hypot(X - P.x_source, Z - P.z_source) / P.c0 - delay_in_lens + t_lens_cor;
    tau_rx = hypot(X - reshape(P.x_piezo, 1, 1, []), Z - reshape(P.z_piezo, 1, 1, [])) / P.c0;
    theta_rx = atan2(Z - reshape(P.z_piezo, 1, 1, []), X - reshape(P.x_piezo, 1, 1, [])) - pi / 2;
    half_opening_angle_rad = atan(1/2 / P.f_number);

    f_number_msk = abs(theta_rx) < half_opening_angle_rad;

    tof_round_trip = tau_tx + tau_rx; % aberration corrected time of flight
    tof_round_trip(~f_number_msk) = NaN; % mask fnumber for reconstruction

    data.Nz = Nz;
    data.Nx = Nx;
    data.tof_round_trip = tof_round_trip;
    data.X = X;
    data.Z = Z;

end
