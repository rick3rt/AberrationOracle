clc
clear all
% close all

% set transducer frequency
P.Fc = 7.5e6;
P.lambda = 1540 / P.Fc;
P.lens_wavespeed = 1000; % m/s
P.lens_thickness = 5 * P.lambda; % m
P.skin_wavespeed = 1600; % m/s
P.brain_wavespeed = 1570; % m/s

% bone properties
P.bone_wavespeed = 3000; % m/s
P.bone_thickness = 2e-3; % m
P.bone_curvature = 10; % m^-1

% distance to bone and pixel
P.distance_trans_bone = 10 * P.lambda;
P.distance_bone_pixel = 15 * P.lambda;

% define transducer
P.num_elements = 128;
P.x_piezo = (0:P.num_elements - 1) .* P.lambda - (P.num_elements - 1) / 2 * P.lambda;
P.z_piezo = 0 * P.x_piezo;

% define source and target point
P.x_source = P.x_piezo(64);
P.z_source = P.z_piezo(64);
P.x_pixel = 0; % center of image
P.z_pixel = P.distance_bone_pixel + P.distance_trans_bone + P.bone_thickness + P.lens_thickness;

% f-number in reconstruction
P.f_number = 1.5;

% from bone parameters, make realization of interfaces
% xax = linspace(-P.num_elements * P.lambda, P.num_elements * P.lambda, 100); %
% P.bone_outer_surface = (P.lens_thickness + P.distance_trans_bone) + P.bone_curvature * xax .^ 2;
% P.bone_inner_surface = (P.lens_thickness + P.distance_trans_bone + P.bone_thickness) + P.bone_curvature * xax .^ 2;

P.name = 'default';

zpix_fun = @(P) P.distance_bone_pixel + P.distance_trans_bone + P.bone_thickness + P.lens_thickness;

P2 = P;
P2.name = 'dist4_curv';
P2.distance_trans_bone = 10 * P2.lambda;
P2.distance_bone_pixel = 15 * P2.lambda;
P2.z_pixel = zpix_fun(P2);

%%
Ptest = P2;
[out, BFC] = rt_compare(Ptest);
figs = rt_compare_plot(Ptest, out, BFC);

% save figes
output_folder = 'figs';
[~, ~] = mkdir(output_folder);
save_fig_fun = @(f) exportgraphics(f, fullfile(output_folder, [f.Name '.png']), 'Resolution', 300);
save_figs_fun = @(fc) cellfun(save_fig_fun, fc);
save_figs_fun(figs)

%% Generate RF
P.Fs = 40 * P.Fc; % Sampling frequency
P.num_cycles = 3;

% Set the impulse response and excitation of the emit aperture
image_pulse = sin(2 * pi * P.Fc * (0:1 / P.Fs:P.num_cycles / P.Fc)); % BP66 %,  sin(2*pi*fc*(0:1/fs:1/fc)); % BP100 %
image_pulse = image_pulse .* hamming(length(image_pulse))';
image_pulse_env = envelope(image_pulse);

[~,imax] =max(image_pulse_env);
ttp = imax/P.Fs;

figure(99); 
plot(image_pulse)
hold on
plot(image_pulse_env )



t_min = 0;
t_max = max(out.ac_tof_round_trip(:)) * 1.5;
t_max = max(out.error_tof_round_trip) * 1.5;
t_vec = t_min:1 / P.Fs:t_max;
N = numel(t_vec);

RF_pulse = repmat(image_pulse, P.num_elements, 1).';
RF_pulse(N, 1) = 0; % pad with zeros

% apply apodization to mimick attenuation
RF_pulse = RF_pulse.*hamming(P.num_elements).';

% RF_delayed = pulse_delaying_RF(RF_pulse, out.ac_tof_round_trip-ttp, P.Fs);
% 
% figure(101); clf;
% % subplot(121)
% % imagesc(P.x_piezo * 1e3, t_vec * 1e6, RF_pulse)
% % subplot(122)
% imagesc(P.x_piezo * 1e3, t_vec * 1e6, RF_delayed)
% hold on
% plot(P.x_piezo * 1e3, out.ac_tof_round_trip * 1e6, 'k--', 'LineWidth', 2)
% plot(P.x_piezo * 1e3, out.nc_tof_round_trip * 1e6, 'r--', 'LineWidth', 2)



RF_delayed_nc = pulse_delaying_RF(RF_pulse, out.error_tof_round_trip, P.Fs);
figure(102);clf;
imagesc(P.x_piezo * 1e3, t_vec * 1e6, RF_delayed_nc)


%% 

RF_pulse_sum_ac = sum(RF_pulse, 2);
RF_pulse_sum_nc = sum(RF_delayed_nc, 2);

figure(1012);clf;
plot(t_vec*1e6, RF_pulse_sum_ac)
hold on 
plot(t_vec*1e6, RF_pulse_sum_nc)

RF_pulse_sum_ac = sum(RF_pulse(:,out.ac_f_number_idx), 2);
RF_pulse_sum_nc = sum(RF_delayed_nc(:,out.nc_f_number_idx), 2);

figure(1013);clf;
plot(t_vec*1e6, RF_pulse_sum_ac)
hold on 
plot(t_vec*1e6, RF_pulse_sum_nc)