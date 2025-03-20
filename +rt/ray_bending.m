function [ray, tof, theta_min] = ray_bending(x_start, z_start, x_target, z_target, BFC, to_layer)

    % cost fun
    cost_fun = @(theta) rt.ray_bending_tof(theta, x_start, z_start, x_target, z_target, BFC, to_layer);

    % determine angle range
    z_lens = BFC.LensThickness;
    x_extend = minmax(BFC.XRecon);
    r_left = [x_extend(1) - x_start; z_lens - z_start];
    r_right = [x_extend(2) - x_start; z_lens - z_start];
    angle_left = atan2(r_left(2), r_left(1)) - pi/2;
    angle_right = atan2(r_right(2), r_right(1)) - pi/2;

    % figure(80); clf; scatter(x_start, z_start);  hold on; scatter(x_extend(1), z_lens); % daspect([1 1 1])

    % coarse grid search
    NAngles = 30;
    theta_vals = linspace(angle_right, angle_left, NAngles);
    tof_grid = arrayfun(cost_fun, theta_vals);
    [~, imin] = min(tof_grid);
    theta0 = theta_vals(imin);

    % minimize:
    theta_min = fminsearch(cost_fun, theta0);
    [tof, ray] = cost_fun(theta_min);

end
