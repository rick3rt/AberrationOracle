clear
clc

CASE_NAME = 'rat_spline';
load('seg_spline.mat')

% set transducer frequency
P.Fc = 15e6;
P.Fs = 20 * P.Fc; % Sampling frequency
P.num_cycles = 3;
P.c0 = 1540; % reference sound speed, homogenous assumption
P.lambda = P.c0 / P.Fc;

% define transducer
P.lens_thickness = 5 * P.lambda; % m
P.num_elements = 128;
P.x_piezo = (0:P.num_elements - 1) .* P.lambda - (P.num_elements - 1) / 2 * P.lambda;
P.z_piezo = 0 * P.x_piezo;

% crop aperture (optional)
% P.x_piezo(P.x_piezo < min(sp1.breaks)) = [];
% P.x_piezo(P.x_piezo > max(sp1.breaks)) = [];
% P.z_piezo = 0 * P.x_piezo;
% P.num_elements = numel(P.x_piezo);

% define source and target point
P.z_source = -10;
P.theta_source = deg2rad(0);

% choose pixel
P.max_depth = 80 * P.lambda; % m
P.x_pixel = -15 * P.lambda; % center of image
P.z_pixel = 50 * P.lambda; % center of image

% f-number in reconstruction
P.f_number = 2;

P.name = 'default';

% define medium
P.medium_soundspeeds = [1000 % lens
                        1600 % skin
                        3300 % bone
                        1570]; % brain
P.medium_density = [1.2, 1.02, 2.0, 1.001]; % kg/m^3
P.medium_attenuation = [0.0, 0.54, 6.9, 0.6]; % dB/(MHz*cm)

% relative distances
dist_trans_bone = 10 * P.lambda; % distance from transducer to bone
dist_bone_pix = 10 * P.lambda;
bone_thickness = 0.38e-3; % m
bone_curvature = 20; % 1/m

sp1_offset = rt2_offset_spline(P, sp1, P.lens_thickness + dist_trans_bone);
sp2_offset = rt2_offset_spline(P, sp2, P.lens_thickness + dist_trans_bone + bone_thickness);

P.medium_interfaces = ...
    {P.lens_thickness,
 sp1_offset,
 sp2_offset};

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
figs = rt2_plot_results(P, result);
rt2_plot_result_comparison(P, P_all, result);

disp('Test run over!')

%% Multiple sources

Ps = P;
% Ps.theta_source = deg2rad([-2 0 2]);
Ps.theta_source = deg2rad(linspace(-6, 6, 7));
Ps.z_source = repelem(Ps.z_source, numel(Ps.theta_source));
result = rt2_simulate_multi_source(Ps);
figs = rt2_plot_results_psf(Ps, result);
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
dist_trans_bone_v = dist_trans_bone; % 10 * P.lambda; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix +30 * P.lambda; % 20 * P.lambda;
bone_thickness_v = bone_thickness; % 4e-3; % m
bone_curvature_v = bone_curvature; % 2; % 1/m

Pn.medium_interfaces = ...
    {Pn.lens_thickness,
 rt2_offset_spline(P, sp1, Pn.lens_thickness + dist_trans_bone_v);
 rt2_offset_spline(P, sp2, Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v)};

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
 rt2_offset_spline(P, sp1, Pn.lens_thickness + dist_trans_bone_v);
 rt2_offset_spline(P, sp2, Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v)};
Pn.z_pixel = Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v + dist_bone_pix_v;
Pn = rt2_derived_parameters(Pn);
P_all(end + 1) = Pn;

% thicker skin
Pn = P;
Pn.name = 'Distance Trans-Bone +0.5mm';
dist_trans_bone_v = dist_trans_bone +0.5e-3; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix;
bone_thickness_v = bone_thickness; % m
bone_curvature_v = bone_curvature; % 1/m
Pn.medium_interfaces = ...
    {Pn.lens_thickness,
 rt2_offset_spline(P, sp1, Pn.lens_thickness + dist_trans_bone_v);
 rt2_offset_spline(P, sp2, Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v)};

Pn.z_pixel = Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v + dist_bone_pix_v;
Pn = rt2_derived_parameters(Pn);
P_all(end + 1) = Pn;

% thick bone
Pn = P;
Pn.name = 'Bone thickness x2';
dist_trans_bone_v = dist_trans_bone; % distance from transducer to bone
dist_bone_pix_v = dist_bone_pix;
bone_thickness_v = bone_thickness * 2; % m
bone_curvature_v = bone_curvature; % 1/m
Pn.medium_interfaces = ...
    {Pn.lens_thickness,
 rt2_offset_spline(P, sp1, Pn.lens_thickness + dist_trans_bone_v);
 rt2_offset_spline(P, sp2, Pn.lens_thickness + dist_trans_bone_v + bone_thickness_v)};
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

%% save results

% result_clean = result;
% result_clean = rmfield(result_clean, 'rf_data');
% result_clean = rmfield(result_clean, 'rf_data_full');

save(fullfile(fig_output_path, 'results.mat'), 'result', 'P', 'P_all')
