clear
clc
addpath('functions\')

% This case is based on the paper:
%   Zhou, J., Guo, Y., Sun, Q., Lin, F., Jiang, C., Xu, K., & Ta, D. (2024).
%   Transcranial ultrafast ultrasound Doppler imaging: A phantom study.
%   Ultrasonics, 144, 107430. https://doi.org/10.1016/j.ultras.2024.107430

load('segmentations/segmentation_zhou2024.mat')
CASE_NAME = 'Zhou2024';

% set transducer frequency
P.Fc = 6.25e6;
P.Fs = 20 * P.Fc; % Sampling frequency
P.num_cycles = 3;
P.c0 = 1540; % reference sound speed, homogenous assumption
P.lambda = P.c0 / P.Fc;

% define transducer
P.lens_thickness = 5 * P.lambda; % m
P.num_elements = 128;
P.pitch = 0.3e-3;
P.x_piezo = (0:P.num_elements - 1) .* P.pitch - (P.num_elements - 1) / 2 * P.pitch;
P.z_piezo = 0 * P.x_piezo;

% define source and target point
P.z_source = -10;
P.theta_source = deg2rad(0);
% P.theta_source = deg2rad([-2 2]);
P.z_source = repelem(P.z_source, numel(P.theta_source));
% choose pixel
P.max_depth = 50 * P.lambda; % m
P.x_pixel = 5 * P.lambda; % center of image
% P.z_pixel = 40e-3; % center of image

% f-number in reconstruction
P.f_number = 2;

P.name = 'default';

% define medium
P.medium_soundspeeds = [1000 % lens
                        1540 % skin
                        2460 % bone
                        1540]; % brain
P.medium_density = [1.2, 1.0, 1.2, 1.00]; % kg/m^3
P.medium_attenuation = [0.0, 0.54, 3.0, 0.54]; % dB/(MHz*cm)

% relative distances
dist_trans_bone = 4e-3; % distance from transducer to bone
dist_bone_pix = 10e-3 +15e-3;
bone_thickness = 1.3e-3; % m
bone_curvature = 1; % 1/m


P.medium_interfaces = ...
    {P.lens_thickness,
     sp1,
     sp2};

P.z_pixel = P.lens_thickness + dist_trans_bone + bone_thickness + dist_bone_pix;

P = rt2_derived_parameters(P);
P_all = P;

% =============================================================================
% Plot medium
% =============================================================================

figure(1); clf;
rt.plot.medium(P);
yline(P.z_pixel * 1e3)
xline(P.x_pixel * 1e3)
% yline(z_c0*1e3)

%% ===================================================================
result = rt2_simulate(P);
% result = rt2_simulate_multi_source(P);

figs = rt2_plot_results(P, result);
rt2_plot_result_comparison(P, P_all, result)

%% MUTTI SOURCE

Ps = P;
% Ps.theta_source = deg2rad([-2 0 2]);
Ps.theta_source = deg2rad(linspace(-6, 6, 7));
Ps.z_source = repelem(Ps.z_source, numel(Ps.theta_source));
result = rt2_simulate_multi_source(Ps);
figs = rt2_plot_results_psf(Ps, result);

%% save figures
fig_output_path = fullfile('figs2', CASE_NAME, P.name);
[~, ~] = mkdir(fig_output_path);
save_fig_fun = @(f) exportgraphics(f, fullfile(fig_output_path, [strrep(f.Name, ' ', '_') '.png']), 'Resolution', 300);
cellfun(save_fig_fun, figs); % save figures

%% make variations
P_all = P;

% relative distances
Pn = P;
Pn.name = 'Distance Bone-Pixel +30L';
dist_trans_bone_v = dist_trans_bone; % 10 * P.lambda; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix +30 * P.lambda; % 20 * P.lambda;
bone_thickness_v = bone_thickness; % 4e-3; % m
bone_curvature_v = bone_curvature; % 2; % 1/m
Pn.medium_interfaces = ...
    {Pn.lens_thickness,
 [bone_curvature_v, 0.02, Pn.lens_thickness + dist_trans_bone_v],
 [bone_curvature_v, 0.03, Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v]};
Pn.z_pixel = Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v + dist_bone_pix_v;
Pn = rt2_derived_parameters(Pn);
P_all(end + 1) = Pn;

% More distance
Pn = P;
Pn.name = 'Distance Bone-Pixel +90L';
dist_trans_bone_v = dist_trans_bone; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix +90 * P.lambda;
bone_thickness_v = bone_thickness; % m
bone_curvature_v = bone_curvature; % 1/m
Pn.medium_interfaces = ...
    {Pn.lens_thickness,
 [bone_curvature_v, 0.02, Pn.lens_thickness + dist_trans_bone_v],
 [bone_curvature_v, 0.03, Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v]};
Pn.z_pixel = Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v + dist_bone_pix_v;
Pn = rt2_derived_parameters(Pn);
P_all(end + 1) = Pn;

% thicker skin
Pn = P;
Pn.name = 'Distance Trans-Bone +2mm';
dist_trans_bone_v = dist_trans_bone +2.0e-3; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix;
bone_thickness_v = bone_thickness; % m
bone_curvature_v = bone_curvature; % 1/m
Pn.medium_interfaces = ...
    {Pn.lens_thickness,
 [bone_curvature_v, 0.02, Pn.lens_thickness + dist_trans_bone_v],
 [bone_curvature_v, 0.03, Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v]};
Pn.z_pixel = Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v + dist_bone_pix_v;
Pn = rt2_derived_parameters(Pn);
P_all(end + 1) = Pn;

% curvy bone skin
Pn = P;
Pn.name = 'Curvature Bone +3 1/m';
dist_trans_bone_v = dist_trans_bone; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix;
bone_thickness_v = bone_thickness; % m
bone_curvature_v = bone_curvature + 3; % 1/m
Pn.medium_interfaces = ...
    {Pn.lens_thickness,
 [bone_curvature_v, 0.02, Pn.lens_thickness + dist_trans_bone_v],
 [bone_curvature_v, 0.03, Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v]};
Pn.z_pixel = Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v + dist_bone_pix_v;
Pn = rt2_derived_parameters(Pn);
P_all(end + 1) = Pn;

figure(2); clf;
rt.plot.medium(Pn);
yline(Pn.z_pixel * 1e3)
xline(Pn.x_pixel * 1e3)

%% Simulate scenarios
for kp = 1:numel(P_all)
    Pn = P_all(kp);
    result(kp) = rt2_simulate(Pn);
    figs = rt2_plot_results(Pn, result(kp));

    % save figures
    fig_output_path = fullfile('figs2', CASE_NAME, Pn.name);
    [~, ~] = mkdir(fig_output_path);
    save_fig_fun = @(f) exportgraphics(f, fullfile(fig_output_path, [strrep(f.Name, ' ', '_') '.png']), 'Resolution', 300);
    cellfun(save_fig_fun, figs); % save figures
end

%% Plot overview resolution
modes = {'AC', 'NC'};

psf_metrics = [result.data_psf];
getter_fun = @(field) arrayfun(@(pm) cellfun(@(m) pm.(m).(field), modes), psf_metrics, 'UniformOutput', 0);

res_x = getter_fun('res_x');
res_x = cat(1, res_x{:});

figs = {};

f = figure(1); clf;
f.Name = 'Resolution_lateral';
f.Position = [200 100 600 600];
figs{end + 1} = f;
subplot(211)
barh(res_x * 1e3)
set(gca, 'YDir', 'reverse')
yticklabels({P_all.name})

subplot(212)

barh(abs((res_x(:, 1) - res_x(:, 2)) ./ res_x(:, 2)) * 100)
set(gca, 'YDir', 'reverse')
yticklabels({P_all.name})

peak_val = getter_fun('peak');
peak_val = cat(1, peak_val{:});

f = figure(2); clf;
f.Name = 'Peak_pixel';
f.Position = [800 100 600 600];
figs{end + 1} = f;
subplot(211)
barh(peak_val)
set(gca, 'YDir', 'reverse')
yticklabels({P_all.name})

subplot(212)

barh(abs((peak_val(:, 1) - peak_val(:, 2)) ./ peak_val(:, 2)) * 100)
set(gca, 'YDir', 'reverse')
yticklabels({P_all.name})

% save figures
fig_output_path = fullfile('figs2', CASE_NAME);
[~, ~] = mkdir(fig_output_path);
save_fig_fun = @(f) exportgraphics(f, fullfile(fig_output_path, [strrep(f.Name, ' ', '_') '.png']), 'Resolution', 300);
cellfun(save_fig_fun, figs); % save figures

%%
return
%%
img = imread('zhou.png');
img = rgb2gray(img);
img = imbinarize(img, 0.8);

z1 = find2d(img, 2, 'first');
z2 = find2d(img, 2, 'last');

z1 = fillmissing(z1, 'linear');
z2 = fillmissing(z2, 'linear');

x_img = linspace(P.x_piezo(1), P.x_piezo(end), size(img, 2));
dx = x_img(2) - x_img(1);
z_img = (0:size(img, 1) - 1) * dx +1e-3;

figure(1); clf;
imagesc(x_img, z_img, img)
hold on
plot(z1)
plot(z2)

daspect([1 1 1])

z1p = z_img(z1);
z2p = z_img(z2) -0.5e-3;
sp1 = rt.util.spline_fit_knots(x_img, z1p, 6);
sp2 = rt.util.spline_fit_knots(x_img, z2p, 7);

figure(2); clf;
imagesc(x_img, z_img, img)
hold on
plot(x_img, ppval(sp1, x_img))
hold on
plot(x_img, ppval(sp2, x_img))
daspect([1 1 1])

save('seg_zhou.mat', 'sp1', 'sp2');

%%
runcor = [1.6 1.4 1.1]
rcorr = [1.1 1.0 0.8];

mean((runcor - rcorr) ./ runcor)
