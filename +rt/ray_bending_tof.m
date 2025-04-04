function [tof, ray] = ray_bending_tof(P, theta, x_start, z_start, x_target, z_target, to_layer)

    % log_info = @(fmt, varargin) fprintf(['[%s] ' fmt], mfilename, varargin{:});
    log_info = @(varargin) [];

    log_info('[ENTER] theta %.3f to_layer %i\n', theta, to_layer);

    % Cast first ray
    ray.critical_angle = 0;
    ray.medium_idx = 1;
    ray.start = [x_start; z_start];
    ray.dir = rt.util.vec_rotate(theta);
    ray = rt.raytheory.ray_traverse(ray, P);

    if ray.length > 1e6; tof = 1e10 + theta; return; end
    if any(isnan(ray.dir_refracted)); tof = 1e10 + theta; return; end

    % cast subsequent rays traversing the scene if any object hit
    k = 1; % num hits
    while ray(k).hit_idx > 0 
        log_info('ray tracing k = %i\n', k);
        k = k + 1;
        ray(k).critical_angle = 0;
        ray(k).medium_idx = ray(k - 1).medium_idx + 1;
        ray(k).start = ray(k - 1).intersection;
        ray(k).dir = ray(k - 1).dir_refracted;

        if k == to_layer
            log_info('target layer reached k = %i\n', k);
            break; % if we reach the target layer
        end
        % if ray points upwards, break
        if (ray(k).dir(2) < 0); tof = 1e10 + theta; return; end

        ray(k) = rt.raytheory.ray_traverse(ray(k), P);
    end

    if ray(k).medium_idx ~= to_layer % k < to_layer +1
        log_info('Target layer not reached, no intersection or critical angles?\n');
        tof = 1e10 + theta; return;
    end

    % connect path
    ray(k).intersection = [x_target; z_target];
    tmp = ray(k).intersection - ray(k).start;
    ray(k).length = vecnorm(tmp);
    ray(k).normal = [NaN; NaN];
    ray(k).dir_refracted = [NaN; NaN];
    ray(k).dir = tmp ./ ray(k).length;

    % compute tof
    tof = sum([ray(1:to_layer).length] ./ P.medium_soundspeeds(1:to_layer));

end

% Rays struct fields:
%   critical_angle
%   medium_idx
%   start
%   dir
%   hit_idx
%   intersection
%   length
%   normal
%   dir_refracted
