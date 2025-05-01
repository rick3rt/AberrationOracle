function P = rt2_derived_parameters(P)
    % utility function to compute derived parameters

    P.medium_soundspeeds = P.medium_soundspeeds(:).'; % guarantee row vector
    P.medium_density = P.medium_density(:).'; % guarantee row vector
    P.medium_interfaces = P.medium_interfaces(:).'; % guarantee row vector

    % derived properties
    P.medium_impendace = P.medium_soundspeeds .* P.medium_density; % kg/(m^2 s)

    % deremine interface derivatives for normal computation
    P.medium_interface_derivatives = cellfun(@rt.util.interface_derivative, P.medium_interfaces, 'UniformOutput', false);

    % reconstriuction grid
    P.max_depth = max(P.max_depth, 1.2 * P.z_pixel);

    P.x_recon = P.x_piezo;
    P.z_recon = 0:P.lambda / 2:P.max_depth;

    % source from angle
    P.x_source = P.z_source .* tan(P.theta_source);

    if ~isfield(P, 'medium_sos_variation')
        P.medium_sos_variation = 0;
    end

end
