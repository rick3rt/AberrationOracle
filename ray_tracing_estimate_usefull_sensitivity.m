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
output_folder= 'figs';
[~,~]=mkdir(output_folder);
save_fig_fun = @(f) exportgraphics(f, fullfile(output_folder, [f.Name '.png']),'Resolution',300);
save_figs_fun = @(fc) cellfun(save_fig_fun, fc);
save_figs_fun(figs)

%% Plot rays
