function tof_ac = rt_trace_tof(P, BFC, xp, zp)

    % ray tracing from source to pixel, and pixel to all elements
    % ===========================================================
    to_layer = numel(P.medium_soundspeeds); % trace to final layer
    [rays_tx, tof_tx, theta_tx] = rt.ray_bending(P.x_source, P.z_source, xp, zp, BFC, to_layer);
    tof_tx = tof_tx - P.tx_add_to_ac;

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

    % mask outside f-number
    % half opening angle corresponding to f-number
    half_opening_angle_rad = atan(1/2 / P.f_number);
    f_number_idx = find(abs(theta_rx_all) > half_opening_angle_rad);

    % mask fnumber for reconstruction
    tof_round_trip(f_number_idx) = NaN;

    tof_ac = tof_round_trip;

end

% % calculate time of flight for homogenous medium
% % ===========================================================
% c0 = 1540;
% delta_tof_nc = (1 / P.bone_wavespeed -1/1540) * P.bone_thickness;
% dz = delta_tof_nc * c0;
%
% vec_pixel = [xp; zp + dz];
%
% tof_round_trip_c0 = vecnorm(vec_pixel - [P.x_source; P.z_source]) / c0 + ...
%     vecnorm([P.x_piezo; P.z_piezo] - vec_pixel) / c0;
%
% rx_vec = vec_pixel - [P.x_piezo; P.z_piezo];
% theta_rx_all_c0 = atan2(rx_vec(2, :), rx_vec(1, :)) - pi / 2;
% f_number_idx_nc = find(abs(theta_rx_all_c0) > half_opening_angle_rad);
% tof_round_trip_c0(f_number_idx_nc) = NaN;
% tof_nc = tof_round_trip_c0;
