% clc
clear all
% close all

% set transducer frequency
P.Fc = 15e6;
P.c0 = 1540; % reference sound speed, homogenous assumption
P.lambda = P.c0 /  P.Fc;
P.lens_wavespeed = 1000; % m/s
P.lens_thickness = 5 * P.lambda; % m
P.skin_wavespeed = 1600; % m/s
P.brain_wavespeed = 1570; % m/s

% bone properties
P.bone_wavespeed = 3000; % m/s
P.bone_thickness = 0.5e-3; % m
P.bone_curvature = 15; % 1/m

% distance to bone and pixel
P.distance_trans_bone = 5 * P.lambda;
P.distance_bone_pixel = 20 * P.lambda;

% define transducer
P.num_elements = 128;
P.x_piezo = (0:P.num_elements - 1) .* P.lambda - (P.num_elements - 1) / 2 * P.lambda;
P.z_piezo = 0 * P.x_piezo;

% define source and target point
ke = 64;
% P.x_source = P.x_piezo(ke);
% P.z_source = P.z_piezo(ke);
% P.x_source = 0;
P.z_source = -10;
P.theta_source = deg2rad(0);

P.x_pixel = -10*P.lambda; % center of image
P.z_pixel = P.lens_thickness + P.distance_trans_bone + P.bone_thickness + P.distance_bone_pixel;

% f-number in reconstruction
P.f_number = 2;

P.name = 'default';
zpix_fun = @(P) P.lens_thickness + P.distance_trans_bone + P.bone_thickness + P.distance_bone_pixel;

% fig saving
output_folder = 'figs_rat';
[~, ~] = mkdir(output_folder);
save_fig_fun = @(f) exportgraphics(f, fullfile(output_folder, [strrep(f.Name, ' ', '_') '.png']), 'Resolution', 300);
save_figs_fun = @(fc) cellfun(save_fig_fun, fc);

% Run and show results
f_outer = @(x) 20*x.^2 + sin(x *1e3)*5e-5 + sin(x *1e3/2-1)*3e-5;
f_inner = @(x) 20*x.^2 + sin(x *1e3)*3e-5 + sin(x *1e3*1.1+2.546)*3e-5 + 0.5e-3;

TEST_RESOLUTION = false;
BFC = rt_make_interfaces(P, f_outer, f_inner); % make custom spline interfaces
[P, BFC,out] = rt_compare(P, BFC);
figs = rt_compare_plot(P, BFC, out);
[metrics, data, figs2] = rt_test_improvement(P, BFC,out, TEST_RESOLUTION);
figs = [figs figs2];

% Compare plots
f = figure(13); clf;
f.Name = ['metrics_' P.name];
figs{end + 1} = f;
if TEST_RESOLUTION; subplot(131); end
b = bar([metrics.nc_peak; metrics.ac_peak], 'FaceColor', 'flat');
b.CData = lines(2);
xticklabels({'NC', 'AC'})
title(sprintf('Intensity: %.2fx', metrics.ac_peak / metrics.nc_peak))

if TEST_RESOLUTION
    subplot(132)
    b = bar([metrics.nc_res_x; metrics.ac_res_x], 'FaceColor', 'flat');
    b.CData = lines(2);
    xticklabels({'NC', 'AC'})
    title(sprintf('Resolution: %.2fx', metrics.ac_res_x / metrics.nc_res_x))
    
    subplot(133)
    b = bar([metrics.nc_IQ_x_peak; metrics.ac_IQ_x_peak], 'FaceColor', 'flat');
    b.CData = lines(2);
    xticklabels({'NC', 'AC'})
    title(sprintf('Peak Pixel Value: %.2fx', metrics.ac_IQ_x_peak / metrics.nc_IQ_x_peak))
end
drawnow

% save_figs_fun(figs);
% return

%% Determine new parameter sets

n = 1;
P_all = P;

Pn = P; n = n + 1;
Pn.name = 'dist tranducer-bone +20L';
Pn.distance_trans_bone = P.distance_trans_bone + 20 * P.lambda;
P_all(n) = Pn;

Pn = P; n = n + 1;
Pn.name = 'dist bone-pixel 10L';
Pn.distance_bone_pixel = 10 * P.lambda;
P_all(n) = Pn;

Pn = P; n = n + 1;
Pn.name = 'dist bone-pixel 100L';
Pn.distance_bone_pixel = 100 * P.lambda;
P_all(n) = Pn;

Pn = P; n = n + 1;
Pn.name = 'dist tranducer-bone +20L bone-pixel 80L';
Pn.distance_trans_bone = P.distance_trans_bone + 20 * P.lambda;
Pn.distance_bone_pixel = 80 * P.lambda;
P_all(n) = Pn;

Pn = P; n = n + 1;
Pn.name = 'curved bone';
Pn.bone_curvature = 30; % 1/m
P_all(n) = Pn;

Pn = P; n = n + 1;
Pn.name = 'thicker bone (2x)';
Pn.bone_thickness = 2 * P.bone_thickness;
P_all(n) = Pn;

for k = 1:numel(P_all)
    P_all(k).z_pixel = zpix_fun(P_all(k));
end

%% Run for all;

TEST_RESOLUTION = 1;

clear metrics
for k = 1:numel(P_all)

    Ptest = P_all(k);

    BFC = struct(); 
    BFC = rt_make_interfaces(Ptest, f_outer, f_inner); % make custom spline interfaces

    [Ptest, BFC, out]  = rt_compare(Ptest, BFC);
    figs = rt_compare_plot(Ptest, BFC, out);
    [metrics(k), ~, figs2] = rt_test_improvement(Ptest, BFC, out, TEST_RESOLUTION);
    figs = [figs figs2];

    % saving of figures
    output_folder = fullfile('figs_rat_spline', Ptest.name);
    [~, ~] = mkdir(output_folder);
    save_fig_fun = @(f) exportgraphics(f, fullfile(output_folder, [strrep(f.Name, ' ', '_') '.png']), 'Resolution', 300);
    save_figs_fun = @(fc) cellfun(save_fig_fun, fc);

    save_figs_fun(figs)
end

%% Plot bar graphs of metrics
ac_peak = [metrics.ac_peak];
nc_peak = [metrics.nc_peak];
imp = [metrics.improvement];

f = figure(100); clf;
f.Position = [200 150 600 800];
subplot(211)
barh([nc_peak; ac_peak].')
yticklabels({P_all.name})
set(gca, 'YDir', 'reverse')
xlabel('Absolute Intensity Difference')
legend('No Correction', 'Aberration Corrected', 'Orientation', 'horizontal', ...
    'Location', 'northoutside')
subplot(212)
barh(imp.')
yticklabels({P_all.name})
set(gca, 'YDir', 'reverse')
title('Relative improvement (peak AC / NC)')
xlabel('Relative Intenisity improvement')


ac_res = [metrics.ac_res_x];
nc_res = [metrics.nc_res_x];
imp = (nc_res - ac_res)./nc_res;

f = figure(101); clf;
f.Position = [800 150 600 800];
subplot(211)
barh([nc_res; ac_res].')
xlabel('Resolution (mm)')
yticklabels({P_all.name})
set(gca, 'YDir', 'reverse')
legend('No Correction', 'Aberration Corrected', 'Orientation', 'horizontal', ...
    'Location', 'northoutside')
subplot(212)
barh(100*imp.')
yticklabels({P_all.name})
set(gca, 'YDir', 'reverse')
title('Relative improvement Resolution (peak AC / NC)')
xlabel('Resolution Improvement (%)')

% ac_res = [metrics.ac_IQ_x_peak];
% nc_res = [metrics.nc_IQ_x_peak];
% imp = [nc_res ./ ac_res];
% f = figure(102); clf;
% f.Position = [400 150 600 800];
% subplot(211)
% barh([nc_res; ac_res].')
% yticklabels({P_all.name})
% set(gca, 'YDir', 'reverse')
% legend('No Correction', 'Aberration Corrected', 'Orientation', 'horizontal', ...
%     'Location', 'northoutside')
% subplot(212)
% barh(imp.')
% yticklabels({P_all.name})
% set(gca, 'YDir', 'reverse')
% title('Relative improvement Peak internsity (BFd)')

return

%% Make complexer shape


xv = P.x_piezo; 


fun1 = @(x) 20*x.^2 + sin(xv *1e3)*5e-5 + sin(xv *1e3/2-1)*3e-5;
fun2 = @(x) 20*x.^2 + sin(xv *1e3)*3e-5 + sin(xv *1e3*1.1+2.546)*3e-5 + 0.5e-3;

z1 = fun1(xv);
z2 = fun2(xv);

sp1 = spline_fit_knots(xv, z1, 8);
sp2 = spline_fit_knots(xv, z2, 8);


figure(1);clf;
plot(xv, fun1(xv))
hold on 
plot(xv, fun2(xv))
plot(xv, ppval(sp1, xv))
plot(xv, ppval(sp2, xv))


set(gca,'ydir','reverse')
axis equal