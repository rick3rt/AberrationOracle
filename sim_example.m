clear
clc

addpath('functions')

% set transducer frequency
P.Fc = 7.5e6;
P.Fs = 20 * P.Fc; % Sampling frequency
P.num_cycles = 3;
P.c0 = 1540; % reference sound speed, homogenous assumption
P.lambda = P.c0 / P.Fc;

% define transducer
P.lens_thickness = 5 * P.lambda; % m
P.num_elements = 128;
P.x_piezo = (0:P.num_elements - 1) .* P.lambda - (P.num_elements - 1) / 2 * P.lambda;
P.z_piezo = 0 * P.x_piezo;

% define source and target point
P.z_source = -10;
P.theta_source = deg2rad(0);

% choose pixel
P.max_depth = 80 * P.lambda; % m
P.x_pixel = -5 * P.lambda; % center of image
P.z_pixel = 50 * P.lambda; % center of image

% f-number in reconstruction
P.f_number = 2;

% P.name = 'default';

% define medium
P.medium_soundspeeds = [1000 % lens
                        1600 % skin
                        3000 % bone
                        1570]; % brain
% P.medium_soundspeeds(:) = 1540;
P.medium_density = [1.2, 1.02, 2.0, 1.001]; % kg/m^3
P.medium_attenuation = [0.0, 0.54, 6.9, 0.6]; % dB/(MHz*cm)

dist_trans_bone = 10 * P.lambda; % distance from transducer to bone
bone_thickness = 4e-3; % m
bone_curvature = 5; % 1/m

P.medium_interfaces = ...
    {P.lens_thickness,
 [bone_curvature, 0.01, P.lens_thickness + dist_trans_bone],
 [bone_curvature, 0.035, P.lens_thickness + dist_trans_bone + bone_thickness]};

P = rt2_derived_parameters(P);

% =============================================================================
% Plot medium
% =============================================================================

figure(1); clf;
rt.plot.medium(P);
yline(P.z_pixel * 1e3)
xline(P.x_pixel * 1e3)

%% Compute

data = struct();
fprintf('Computing Max Intensity...\n')

mode = 'NC'; % no correction, homogeneous medium
data.NC = rt2_compute(P, mode);

mode = 'LC'; % lens corrected
data.LC = rt2_compute(P, mode);

mode = 'AC'; % full corrected
data.AC = rt2_compute(P, mode);

[rf_data, rf_data_full] = rt2_generate(P);
metrics = rt2_compute_max_intensity(P, data, rf_data);

modes = fieldnames(metrics);
getter_fun = @(field) cellfun(@(m)metrics.(m).(field), fieldnames(metrics));

figure(11); clf;
bar(getter_fun('peak'))
xticklabels(modes)

rf_plot = rf_data_full;
% rf_plot = rf_data;

LOPTS = {'linewidth', 1.5};

figure(12); clf;
imagesc(P.x_piezo * 1e3, rf_plot.t_vec * 1e6, rf_plot.RF)
colormap bone
hold on
plot(P.x_piezo * 1e3, data.AC.tof_round_trip * 1e6, LOPTS{:});
plot(P.x_piezo * 1e3, data.NC.tof_round_trip * 1e6, LOPTS{:});
plot(P.x_piezo * 1e3, data.LC.tof_round_trip * 1e6, LOPTS{:});

%%

%% reconstruct PSF
fprintf('Computing PSFs...\n')
clear data_psf
m = 'AC';
data_psf.(m) = rt2_compute_psf(P, rf_data_full, m);
m = 'LC';
data_psf.(m) = rt2_compute_psf(P, rf_data_full, m);
m = 'NC';
data_psf.(m) = rt2_compute_psf(P, rf_data_full, m);

%%

modes = fieldnames(data_psf);
NModes = numel(modes);

figure(21); clf; hold on
for km = 1:NModes
    m = modes{km};
    plot(data_psf.(m).xv * 1e3, abs(data_psf.(m).IQ_line), 'DisplayName', m)
end
legend

%%

figure(22); clf;
for km = 1:NModes
    subplot(1, NModes, km)
    m = modes{km};
    imagesc(data_psf.(m).xv * 1e3, data_psf.(m).zv * 1e3, abs(data_psf.(m).IQ_grid))
    title(m)
    colormap bone
    daspect([1 1 1])
end
