% clc
clear all
% close all

% set transducer frequency
P.Fc = 15e6;
P.lambda = 1540 / P.Fc;
P.lens_wavespeed = 1000; % m/s
P.lens_thickness = 5 * P.lambda; % m
P.skin_wavespeed = 1600; % m/s
P.brain_wavespeed = 1570; % m/s

% bone properties
P.bone_wavespeed = 3000; % m/s
P.bone_thickness = 0.4e-3; % m
P.bone_curvature = 15; % 1/m

% distance to bone and pixel
P.distance_trans_bone = 5 * P.lambda;
P.distance_bone_pixel = 10 * P.lambda;

% define transducer
P.num_elements = 128;
P.x_piezo = (0:P.num_elements - 1) .* P.lambda - (P.num_elements - 1) / 2 * P.lambda;
P.z_piezo = 0 * P.x_piezo;

% define source and target point
P.x_source = P.x_piezo(64);
P.z_source = P.z_piezo(64);
P.x_pixel = 0; % center of image
P.z_pixel = P.lens_thickness + P.distance_trans_bone + P.bone_thickness + P.distance_bone_pixel;

% f-number in reconstruction
P.f_number = 1.5;

% from bone parameters, make realization of interfaces
% xax = linspace(-P.num_elements * P.lambda, P.num_elements * P.lambda, 100); %
% P.bone_outer_surface = (P.lens_thickness + P.distance_trans_bone) + P.bone_curvature * xax .^ 2;
% P.bone_inner_surface = (P.lens_thickness + P.distance_trans_bone + P.bone_thickness) + P.bone_curvature * xax .^ 2;
 
P.name = 'default';
zpix_fun = @(P) P.lens_thickness + P.distance_trans_bone + P.bone_thickness + P.distance_bone_pixel;

% Run and show results 
[out, BFC] = rt_compare(P);
figs = rt_compare_plot(P, out, BFC);
metrics = rt_test_improvement(P, out, BFC);

% return 

%% Determine new parameter sets

n = 1; 
P_all = P;


Pn = P; n = n + 1;
Pn.name = 'dist tranducer-bone +10L';
Pn.distance_trans_bone = P.distance_trans_bone + 10 * P.lambda;
P_all(n) = Pn;

Pn = P; n = n + 1;
Pn.name = 'dist bone-pixel +10L';
Pn.distance_bone_pixel = P.distance_bone_pixel + 10 * P.lambda;
P_all(n) = Pn;

Pn = P; n = n + 1;
Pn.name = 'curved bone';
Pn.bone_curvature = 30; % 1/m
P_all(n) = Pn;

Pn = P; n = n + 1;
Pn.name = 'thicker bone (1.5x)';
Pn.bone_thickness = 1.5*P.bone_thickness;
P_all(n) = Pn;

Pn = P; n = n + 1;
Pn.name = 'f-number 2';
Pn.f_number = 2;
P_all(n) = Pn;


for k = 1:numel(P_all)
    P_all(k).z_pixel = zpix_fun(P_all(k));
end


%% Run for all; 

output_folder = 'figs';
[~, ~] = mkdir(output_folder);
save_fig_fun = @(f) exportgraphics(f, fullfile(output_folder, [strrep(f.Name,' ','_') '.png']), 'Resolution', 300);
save_figs_fun = @(fc) cellfun(save_fig_fun, fc);


for k = 1:numel(P_all)
    Ptest = P_all(k);
    [out, BFC] = rt_compare(Ptest);
    figs = rt_compare_plot(Ptest, out, BFC);
    metrics(k) = rt_test_improvement(P, out, BFC, false);
   
    % save_figs_fun(figs)
end

% Plot bar graphs of metrics
ac_peak = [metrics.ac_peak];
nc_peak = [metrics.nc_peak];
imp = [metrics.improvement];

f = figure(100);clf;
f.Position = [400 150 600 800];
subplot(211)
barh([nc_peak;ac_peak].')
yticklabels({P_all.name})
set(gca,'YDir','reverse')
legend('No Correction', 'Aberration Corrected', 'Orientation','horizontal',...
 'Location','northoutside')
subplot(212)
barh(imp.')
yticklabels({P_all.name})
set(gca,'YDir','reverse')
title('Relative improvement (peak AC / NC)')

return 
