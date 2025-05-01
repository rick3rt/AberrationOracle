clear
clc

addpath functions \

CASE_NAME = 'fat_layer';

% set transducer frequency
P.Fc = 8e6;
P.Fs = 20 * P.Fc; % Sampling frequency
P.num_cycles = 3;
P.c0 = 1540; % reference sound speed, homogenous assumption
P.lambda = P.c0 / P.Fc;

% define transducer
P.lens_thickness = 5 * P.lambda; % m
P.pitch = 0.15e-3;
P.num_elements = 128;
P.x_piezo = (0:P.num_elements - 1) .* P.pitch - (P.num_elements - 1) / 2 * P.pitch;
P.z_piezo = 0 * P.x_piezo;

% define source and target point
P.z_source = -10;
P.theta_source = deg2rad(0);

% choose pixel
P.max_depth = 100 * P.lambda; % m
P.x_pixel = 0 * P.lambda; % center of image
P.z_pixel = 50 * P.lambda; % center of image

% f-number in reconstruction
P.f_number = 2;

P.name = 'default';

% define medium
P.medium_soundspeeds = [1480 % lens
                        1480 % fat
                        1600]; % fat
% 1600]; % brain
% P.medium_soundspeeds = [1000% lens
% 1540 % fat
% 1540]; % fat
% 1600]; % brain
P.medium_density = [1.2, 1.02, 1.02]; %, 1.001]; % kg/m^3
P.medium_attenuation = [0.0, 0.54, 0.54]; % dB/(MHz*cm)

% define medium
P.medium_soundspeeds = [1000 % lens
                        1600 % skin
                        1400 % fat
                        1550]; % liver
% P.medium_soundspeeds = [1000% lens
% 1540 % fat
% 1540]; % fat
% 1600]; % brain
P.medium_density = [1.2, 1.02, 1.02, 1.02]; %, 1.001]; % kg/m^3
P.medium_attenuation = [0.0, 0.54, 0.54, 0.54]; % dB/(MHz*cm)

% introduce local aberrations:
P.medium_sos_variation = 0.002; % 0.1 %

% relative distances
skin_thickness = 2e-3; %  * P.lambda; % distance from transducer to bone
dist_bone_pix = 120 * P.lambda;
fat_thickness = 10e-3; % m
bone_curvature = 0; % 1/m

P.medium_interfaces = ...
    {P.lens_thickness,
 [bone_curvature, 0.0, P.lens_thickness + skin_thickness]
 [bone_curvature, 0.0, P.lens_thickness + skin_thickness + fat_thickness]};

P.z_pixel = P.lens_thickness + skin_thickness + fat_thickness + dist_bone_pix;

P = rt2_derived_parameters(P);

% =============================================================================
% Plot medium
% =============================================================================

figure(1); clf;
rt.plot.medium(P);
yline(P.z_pixel * 1e3)
xline(P.x_pixel * 1e3)
% yline(z_c0*1e3)
drawnow

% ===================================================================
result = rt2_simulate(P);
figs = rt2_plot_results(P, result);

%%

%%
return

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
dist_trans_bone_v = skin_thickness; % 10 * P.lambda; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix +30 * P.lambda; % 20 * P.lambda;
bone_thickness_v = fat_thickness; % 4e-3; % m
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
dist_trans_bone_v = skin_thickness; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix +90 * P.lambda;
bone_thickness_v = fat_thickness; % m
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
Pn.name = 'Distance Trans-Bone +0.5mm';
dist_trans_bone_v = skin_thickness +0.5e-3; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix;
bone_thickness_v = fat_thickness; % m
bone_curvature_v = bone_curvature; % 1/m
Pn.medium_interfaces = ...
    {Pn.lens_thickness,
 [bone_curvature_v, 0.02, Pn.lens_thickness + dist_trans_bone_v],
 [bone_curvature_v, 0.03, Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v]};
Pn.z_pixel = Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v + dist_bone_pix_v;
Pn = rt2_derived_parameters(Pn);
P_all(end + 1) = Pn;

% curvy bone
Pn = P;
Pn.name = 'Curvature Bone +10 1/m';
dist_trans_bone_v = skin_thickness; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix;
bone_thickness_v = fat_thickness; % m
bone_curvature_v = bone_curvature + 10; % 1/m
Pn.medium_interfaces = ...
    {Pn.lens_thickness,
 [bone_curvature_v, 0.02, Pn.lens_thickness + dist_trans_bone_v],
 [bone_curvature_v, 0.03, Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v]};
Pn.z_pixel = Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v + dist_bone_pix_v;
Pn = rt2_derived_parameters(Pn);
P_all(end + 1) = Pn;

% thick bone
Pn = P;
Pn.name = 'Bone thickness x2';
dist_trans_bone_v = skin_thickness; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix;
bone_thickness_v = fat_thickness * 2; % m
bone_curvature_v = bone_curvature; % 1/m
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
