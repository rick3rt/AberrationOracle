function [tof_ac, tof_nc] = rt_trace_tof(P, BFC, xp, zp)


    % ray tracing from source to pixel, and pixel to all elements
    % ===========================================================
    to_layer = 4; % trace to brain layer
    [rays_tx, tof_tx, theta_tx] = rt.ray_bending(P.x_source, P.z_source, xp,zp, BFC, to_layer);

    tof_rx_all = zeros(1, P.num_elements);
    theta_rx_all = zeros(1, P.num_elements);
    rays_rx_all = cell(1, P.num_elements);
    for ke = 1:P.num_elements
        x_end = P.x_piezo(ke);
        z_end = P.z_piezo(ke);
        [rays_rx, tof_rx, theta_rx] = rt.ray_bending(x_end, z_end, xp, zp, BFC, to_layer);
        tof_rx_all(ke) = tof_rx;
        theta_rx_all(ke) = theta_rx;
        rays_rx_all{ke} = rays_rx;
    end

    ind_valid = find(~isnan(tof_rx_all)); % valid rays in reception.
    tof_round_trip = tof_tx + tof_rx_all;

    % calculate time of flight for homogenous medium
    % ===========================================================
    c0 = 1540;
    tof_round_trip_c0 = vecnorm([xp;zp] - [P.x_source; P.z_source]) / c0 + ...
        vecnorm([P.x_piezo; P.z_piezo] - [xp;zp]) / c0;

    rx_vec = [xp;zp] - [P.x_piezo; P.z_piezo];
    theta_rx_all_c0 = atan2(rx_vec(2, :), rx_vec(1, :)) - pi / 2;


    % mask outside f-number
    % half opening angle corresponding to f-number
    half_opening_angle_rad = atan(1/2 / P.f_number);
    f_number_idx = find(abs(theta_rx_all) > half_opening_angle_rad);
    f_number_idx_nc = find(abs(theta_rx_all_c0) > half_opening_angle_rad);

    
    % determine outlier in center; and fix
    
    % error_rt_tof = tof_round_trip-tof_round_trip_c0;
    % msk = error_rt_tof(error_rt_tof==0);
    % tof_round_trip(msk) = NaN;
    % tof_round_trip = fillmissing(tof_round_trip, 'spline');

    msk = isoutlier(diff(tof_round_trip),'mean');
    tof_round_trip(logical([msk 0])) = NaN;
    tof_round_trip = fillmissing(tof_round_trip, 'spline');
    

    tof_round_trip(f_number_idx) = NaN;
    tof_round_trip_c0(f_number_idx_nc) = NaN;
    



    tof_ac = tof_round_trip;
    tof_nc = tof_round_trip_c0;

end