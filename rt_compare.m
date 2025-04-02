function [out, BFC] = rt_compare(P)
    % compare multi layer ray tracing with homogenous assumption
    % P: parameters struct

    % prepare ray tracer
    % ===========================================================

    % define sound speed per layer
    BFC.medium_soundspeeds = [P.lens_wavespeed,
                              P.skin_wavespeed,
                              P.bone_wavespeed,
                              P.brain_wavespeed];

    BFC.medium_density = [1.0, 1.0, 1.0, 1.0]; % kg/m^3

    % define  interfaces between tissue layers
    BFC.medium_interfaces = {P.lens_thickness,
                             [P.bone_curvature, 0, P.lens_thickness + P.distance_trans_bone],
                             [P.bone_curvature, 0, P.lens_thickness + P.distance_trans_bone + P.bone_thickness]};

    BFC.medium_soundspeeds = BFC.medium_soundspeeds(:).'; % guarantee row vector
    BFC.medium_density = BFC.medium_density(:).'; % guarantee row vector
    BFC.medium_impendace = BFC.medium_soundspeeds .* BFC.medium_density; % kg/(m^2 s)

    BFC.medium_interfaces = BFC.medium_interfaces(:).'; % guarantee row vector

    % deremine interface derivatives for normal computation
    BFC.medium_interface_derivatives = cellfun(@rt.util.interface_derivative, BFC.medium_interfaces, 'UniformOutput', false);

    % put other variables in BFC
    BFC.XPiezo = P.x_piezo;
    BFC.ZPiezo = P.z_piezo;
    BFC.XRecon = P.x_piezo;
    BFC.ZRecon = 0:P.lambda / 2:100 * P.lambda;
    BFC.LensThickness = P.lens_thickness;

    % plot the medium
    % ===========================================================
    % figure(1); clf;
    % hold on
    % rt.plot.medium(BFC);
    % scatter(P.x_pixel * 1e3, P.z_pixel * 1e3, 'ko', 'filled');
    % set(gca, 'YDir', 'reverse')
    % daspect([1 1 1])
    % ===========================================================

    % ray tracing from source to pixel, and pixel to all elements
    % ===========================================================
    to_layer = 4; % trace to brain layer
    [rays_tx, tof_tx, theta_tx] = rt.ray_bending(P.x_source, P.z_source, P.x_pixel, P.z_pixel, BFC, to_layer);

    tof_rx_all = zeros(1, P.num_elements);
    theta_rx_all = zeros(1, P.num_elements);
    rays_rx_all = cell(1, P.num_elements);
    for ke = 1:P.num_elements
        x_end = P.x_piezo(ke);
        z_end = P.z_piezo(ke);
        [rays_rx, tof_rx, theta_rx] = rt.ray_bending(x_end, z_end, P.x_pixel, P.z_pixel, BFC, to_layer);
        tof_rx_all(ke) = tof_rx;
        theta_rx_all(ke) = theta_rx;
        rays_rx_all{ke} = rays_rx;
    end

    ind_valid = find(~isnan(tof_rx_all)); % valid rays in reception.
    tof_round_trip = tof_tx + tof_rx_all;

    % calculate time of flight for homogenous medium
    % ===========================================================
    c0 = 1540;
    delta_tof_nc = (1/P.bone_wavespeed - 1/1540)*P.bone_thickness;
    dz = delta_tof_nc*c0;
    
    vec_pixel = [P.x_pixel; P.z_pixel+dz];

    tof_round_trip_c0 = vecnorm(vec_pixel  - [P.x_source; P.z_source]) / c0 + ...
        vecnorm([P.x_piezo; P.z_piezo] - vec_pixel) / c0;

    rx_vec = vec_pixel - [P.x_piezo; P.z_piezo];
    theta_rx_all_c0 = atan2(rx_vec(2, :), rx_vec(1, :)) - pi / 2;

    % get rid of axial shift:
    tof_round_trip0 = tof_round_trip - min(tof_round_trip);
    tof_round_trip_c0_0 = tof_round_trip_c0 - min(tof_round_trip_c0);

    wave_period = 1 / P.Fc;

    % half opening angle corresponding to f-number
    half_opening_angle_rad = atan(1/2 / P.f_number);
    f_number_idx = find(abs(theta_rx_all) < half_opening_angle_rad);
    f_number_idx_nc = find(abs(theta_rx_all_c0) < half_opening_angle_rad);

    % collect results
    % ===========================================================
    out.ac_rays_tx = rays_tx;
    out.ac_tof_tx = tof_tx;
    % out.ac_theta_tx = theta_tx; %
    out.ac_rays_rx = rays_rx_all;
    out.ac_tof_rx = tof_rx_all;
    out.ac_theta_rx = theta_rx_all;
    out.ac_f_number_idx = f_number_idx;
    out.ac_valid_ray_idx = ind_valid;
    out.ac_tof_round_trip = tof_round_trip;

    out.nc_theta_rx = theta_rx_all_c0;
    out.nc_f_number_idx = f_number_idx_nc;
    out.nc_tof_round_trip = tof_round_trip_c0;

    out.error_tof_round_trip = tof_round_trip_c0_0 - tof_round_trip0;
    out.error_tof_round_trip_relative = out.error_tof_round_trip ./ wave_period;
    out.half_opening_angle_rad = half_opening_angle_rad;
end
