function [ray, tof, theta_min] = ray_bending(x_start, z_start, x_target, z_target, BFC, to_layer)

    % cost fun
    cost_fun = @(theta) rt.ray_bending_tof(theta, x_start, z_start, x_target, z_target, BFC, to_layer);

    % determine angle range
    z_lens = BFC.LensThickness;
    x_extend = rt.util.minmax(BFC.XRecon);
    r_left = [x_extend(1) - x_start; z_lens - z_start];
    r_right = [x_extend(2) - x_start; z_lens - z_start];
    angle_left = atan2(r_left(2), r_left(1)) - pi/2;
    angle_right = atan2(r_right(2), r_right(1)) - pi/2;

    % figure(80); clf; scatter(x_start, z_start);  hold on; scatter(x_extend(1), z_lens); % daspect([1 1 1])

    % coarse grid search
    NAngles = 50;
    theta_vals = linspace(angle_right, angle_left, NAngles);
    tof_grid = arrayfun(cost_fun, theta_vals, 'UniformOutput', true);
    [~, imin] = min(tof_grid);
    theta0 = theta_vals(imin);

    % minimize with fminsearch
    options = optimset('TolX',1e-8); % refine tolerance, especially important in TX
    % [theta_min,fval,exit_flag,out] = fminsearch(cost_fun, theta0,options);
    theta_min = fminsearch(cost_fun, theta0, options);
    [tof, ray] = cost_fun(theta_min);


    % final check, if angle of ray(k).dir and ray(k-1).dir_refracted
    % withing range
    theta = acos(dot(ray(end).dir, ray(end-1).dir_refracted));
    if theta > 0.1
        tof = NaN; theta_min = NaN; 
    end


    % minimize with fminbnd
    %theta_min = fminbnd(cost_fun, theta0-0.1, theta0+0.1);
    %[tof, ray] = cost_fun(theta_min);



end
