% Syntax:
%   rays(rays)
%
% Description:
%   This function visualizes a set of rays by plotting their starting points,
%   directions, intersections, and normals. It uses various plotting utilities
%   from the `rt.plot` package to render the rays and their components.
%
% Input:
%   rays - An array of ray structures, where each structure contains the
%          following fields:
%          - start: A 3-element vector representing the starting point of the ray.
%          - dir: A 3-element unit vector representing the direction of the ray.
%          - length: A scalar representing the length of the ray.
%          - intersection: A 3-element vector representing the intersection point.
%          - normal: A 3-element vector representing the normal at the intersection.
%
% R. Waasdorp, 19-03-2025

function rays(rays, plot_normals, plot_true_refraction)

    if ~exist('plot_normals','var'); plot_normals = 0; end
    if ~exist('plot_true_refraction','var'); plot_true_refraction = 0; end

    ray_opts = {'color','k','linewidth',1.2};
    normal_opts = {'b', 'linewidth',3};
    scat_opts = {20,'ko','filled'};

    rt.plot.plot_point(rays(1).start.' * 1e3, scat_opts{:});
    for k = 1:numel(rays)
        ray = rays(k);
        if ray.hit_idx == -1; break; end
        rt.plot.plot_ray(ray.start * 1e3, ray.length * 1e3 * ray.dir, ray_opts{:});
        if k < numel(rays) && plot_true_refraction
            rt.plot.plot_ray(rays(k + 1).start * 1e3, ...
                rays(k + 1).length * 1e3 * ray.dir_refracted, normal_opts{:});
        end
        rt.plot.plot_point(ray.intersection.' * 1e3, scat_opts{:});

        if plot_normals
            rt.plot.plot_vec(ray.intersection * 1e3 -ray.normal, ...
                ray.intersection * 1e3 +ray.normal, normal_opts{:});
        end
    end
end
