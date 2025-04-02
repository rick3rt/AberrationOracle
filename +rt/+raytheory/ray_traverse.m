function ray = ray_traverse_spline(ray, BFC)
    
    import rt.raytheory.*

    % find closest idx
    dist = inf;
    ray.hit_idx = -1;

    if any(isnan(ray.dir)); ray.length = inf; return; end
    % test self intersection: 
    %   max(1,ray.medium_idx-1):numel(BFC.medium_interfaces)
    % otherwise:
    %   ray.medium_idx:numel(BFC.medium_interfaces)
    
    %for k = ray.medium_idx:numel(BFC.medium_interfaces(1:to_layer))
    k = ray.medium_idx;
        sp_test = BFC.medium_interfaces{k};
        if isstruct(sp_test)
            [~, x_intersect] = raylength_spline(sp_test, ray.start, ray.dir, 0);
            % x_intersect = intersect_spline_line(sp_test, ray.start, ray.dir);
            z_intersect = ppval(sp_test, x_intersect);
        else
            slope = ray.dir(2) ./ ray.dir(1);
            if isnearinf(slope)
                x_intersect = ray.start(1);
                p_ray = sp_test; % 
            else
                intercept = ray.start(2) - ray.start(1) * slope;
                p_ray = [slope, intercept];
                if any(isnan(p_ray)); ray.length = inf; return; end
                p_test = poly_comb(p_ray, sp_test);
                x_intersect = poly_roots_valid(p_test, ray);
                if numel(x_intersect) == 0; ray.length = inf; return; end
            end
            z_intersect = polyval(p_ray, x_intersect);
        end
        v_intersect = [x_intersect; z_intersect];
        dist_object = vecnorm(v_intersect - ray.start);

        if (dist_object < dist)
            dist = dist_object;
            ray.intersection = v_intersect;
            ray.length = dist;
            ray.hit_idx = k;
        end
    % end

    if isinf(dist)
        ray.length = inf;
        return;
    end
    

    % determine normal;
    N = sign(ray.dir(2)) * -1*[rt.util.segeval(BFC.medium_interface_derivatives{ray.hit_idx}, ray.intersection(1)); -1];
    % N =- [polyval(BFC.medium_interface_derivatives{ray.hit_idx}, ray.intersection(1)); -1];
    N = N ./ norm(N);
    ray.normal = N;
    R = ray.dir;
    c1 = BFC.medium_soundspeeds(ray.medium_idx);
    c2 = BFC.medium_soundspeeds(ray.medium_idx + 1);

    % ray.theta_in = acos(dot(ray.normal,ray.dir));

    mu = c2/c1;
    R2 = sqrt(1 - mu .^ 2 * (1 - dot(N, R) .^ 2)) * N + mu * (R - dot(N, R) * N);
    
    if imag(R2) == 0
        ray.dir_refracted = R2; 
    else
        ray.dir_refracted = [NaN; NaN];
        ray.critical_angle = 1;
    end

end

function val = isnearinf(slope)
    val = abs(slope) > 1e16;
end