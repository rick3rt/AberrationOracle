function [TF, TF_all] = transmission_coeff(BFC, rays, go_forward)

    TF_all = zeros(numel(rays)-1, 1);

    transmission_factor = @(Z1, Z2, ti, to) 2 * Z2 * cos(ti) / (Z2 * cos(ti) + Z1 * cos(to));
    % transmission_factor = @(Z1, Z2, ti, to) 2 * Z2 / (Z2 + Z1);
    
    if go_forward
        % TX
        for k = 1:numel(rays) - 1
            ray = rays(k);
            theta_i = acos(dot(ray.dir, ray.normal));
            theta_o = acos(dot(ray.dir_refracted, ray.normal));
            Z1 = BFC.medium_impendace(k);
            Z2 = BFC.medium_impendace(k + 1);
            TF_all(k) = transmission_factor(Z1, Z2, theta_i, theta_o);
        end
    else
        % RX
        for k = 1:numel(rays) - 1
            ray = rays(k);
            theta_o = acos(dot(ray.dir, ray.normal));
            theta_i = acos(dot(ray.dir_refracted, ray.normal));
            Z1 = BFC.medium_impendace(k + 1);
            Z2 = BFC.medium_impendace(k);
            TF_all(k) = transmission_factor(Z1, Z2, theta_i, theta_o);
        end
    end

    TF = prod(TF_all);

end

