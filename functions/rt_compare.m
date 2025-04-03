function [P, BFC, out] = rt_compare(P, BFC)
    % compare multi layer ray tracing with homogenous assumption
    % P: parameters struct

    if ~exist('BFC','var'); BFC = struct(); end

    % prepare ray tracer
    % ===========================================================

    % TODO IMPLEMENT ATTENUATION LENS? (FIX POINT SOURCE STUFF)
    if ~isfield(BFC, 'medium_attenuation')
        BFC.medium_attenuation = [0.0, 0.54, 6.9, 0.6]; % dB/(MHz*cm)
    end

    % define sound speed per layer
    if ~isfield(BFC, 'medium_soundspeeds')
        BFC.medium_soundspeeds = [P.lens_wavespeed,
                                  P.skin_wavespeed,
                                  P.bone_wavespeed,
                                  P.brain_wavespeed];
    end

    if ~isfield(BFC, 'medium_density')
        BFC.medium_density = [1.2, 1.02, 2.0, 1.001]; % kg/m^3
    end

    % define  interfaces between tissue layers
    if ~isfield(BFC, 'medium_interfaces')
        BFC.medium_interfaces = {P.lens_thickness,
                                 [P.bone_curvature, 0, P.lens_thickness + P.distance_trans_bone],
                                 [P.bone_curvature, 0, P.lens_thickness + P.distance_trans_bone + P.bone_thickness]};
    end

    %% compute derived properties and pack other settings in BFC
    BFC.medium_soundspeeds = BFC.medium_soundspeeds(:).'; % guarantee row vector
    BFC.medium_density = BFC.medium_density(:).'; % guarantee row vector
    BFC.medium_interfaces = BFC.medium_interfaces(:).'; % guarantee row vector

    % derived properties
    BFC.medium_impendace = BFC.medium_soundspeeds .* BFC.medium_density; % kg/(m^2 s)
    % deremine interface derivatives for normal computation
    BFC.medium_interface_derivatives = cellfun(@rt.util.interface_derivative, BFC.medium_interfaces, 'UniformOutput', false);

    % put other variables in BFC
    BFC.XPiezo = P.x_piezo;
    BFC.ZPiezo = P.z_piezo;
    BFC.XRecon = P.x_piezo;
    BFC.ZRecon = 0:P.lambda / 2:100 * P.lambda;
    BFC.LensThickness = P.lens_thickness;

    %% ===========================================================

    % compute source location - homogeneous and lens speed of sound
    c0 = P.c0;
    xs = P.z_source * tan(P.theta_source);
    vs = [xs; P.z_source];
    ve = [P.x_piezo; P.z_piezo];
    ds = vecnorm(ve - vs);
    tx_delay_lens = ds / P.lens_wavespeed;
    tx_delay_c0 = ds / c0;
    P.x_source = xs;
    P.tx_add_to_ac = min(tx_delay_lens);
    P.tx_add_to_nc = min(tx_delay_c0);

    % plot the medium
    % ===========================================================
    % figure(1); clf;
    % hold on
    % rt.plot.medium(BFC);
    % scatter(P.x_pixel * 1e3, P.z_pixel * 1e3, 'ko', 'filled');
    % set(gca, 'YDir', 'reverse')
    % daspect([1 1 1])
    % ===========================================================

    %% ray tracing from source to pixel, and pixel to all elements
    % ===========================================================
    to_layer = numel(BFC.medium_soundspeeds); % trace to brain layer
    [rays_tx, tof_tx, theta_tx] = rt.ray_bending(P.x_source, P.z_source, P.x_pixel, P.z_pixel, BFC, to_layer);
    tof_tx = tof_tx - P.tx_add_to_ac; % correct propagation time in lens

    tof_rx_all = zeros(1, P.num_elements);
    theta_rx_all_ac = zeros(1, P.num_elements);
    rays_rx_all = cell(1, P.num_elements);
    for ke = 1:P.num_elements
        x_end = P.x_piezo(ke);
        z_end = P.z_piezo(ke);
        [rays_rx, tof_rx, theta_rx] = rt.ray_bending(x_end, z_end, P.x_pixel, P.z_pixel, BFC, to_layer);
        tof_rx_all(ke) = tof_rx;
        theta_rx_all_ac(ke) = theta_rx;
        rays_rx_all{ke} = rays_rx;
    end

    ind_valid = find(~isnan(tof_rx_all)); % valid rays in reception.
    tof_round_trip_ac = tof_tx + tof_rx_all; % aberration corrected time of flight

    %% calculate time of flight for homogenous medium
    % ===========================================================
    delta_tof_nc = (1 / P.bone_wavespeed -1/1540) * P.bone_thickness;
    dz = delta_tof_nc * c0;

    v_pixel = [P.x_pixel; P.z_pixel + dz];
    v_source = [P.x_source; P.z_source];
    v_piezo = [P.x_piezo; P.z_piezo];

    tof_tx_nc = vecnorm(v_pixel - v_source) / c0 - P.tx_add_to_nc;
    tof_rx_nc = vecnorm(v_piezo - v_pixel) / c0;
    tof_round_trip_nc = tof_tx_nc + tof_rx_nc;

    v_rx = v_pixel - v_piezo; % return vector
    theta_rx_all_nc = atan2(v_rx(2, :), v_rx(1, :)) - pi / 2;

    % get rid of axial shift:
    tof_round_trip_ac0 = tof_round_trip_ac - min(tof_round_trip_ac); % remove offset
    tof_round_trip_nc0 = tof_round_trip_nc - min(tof_round_trip_nc); % remove offset

    wave_period = 1 / P.Fc;

    % half opening angle corresponding to f-number
    half_opening_angle_rad = atan(1/2 / P.f_number);
    f_number_idx_ac = find(abs(theta_rx_all_ac) < half_opening_angle_rad);
    f_number_idx_nc = find(abs(theta_rx_all_nc) < half_opening_angle_rad);

    %% collect results
    % ===========================================================
    out.ac_rays_tx = rays_tx;
    out.ac_tof_tx = tof_tx;
    % out.ac_theta_tx = theta_tx; %
    out.ac_rays_rx = rays_rx_all;
    out.ac_tof_rx = tof_rx_all;
    out.ac_theta_rx = theta_rx_all_ac;
    out.ac_f_number_idx = f_number_idx_ac;
    out.ac_valid_ray_idx = ind_valid;
    out.ac_tof_round_trip = tof_round_trip_ac;

    out.nc_theta_rx = theta_rx_all_nc;
    out.nc_f_number_idx = f_number_idx_nc;
    out.nc_tof_round_trip = tof_round_trip_nc;

    out.error_tof_round_trip = tof_round_trip_nc0 - tof_round_trip_ac0;
    out.error_tof_round_trip_relative = out.error_tof_round_trip ./ wave_period;
    out.half_opening_angle_rad = half_opening_angle_rad;

end
