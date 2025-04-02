% clc
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
P.bone_curvature = 00; % 1/m

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
[metrics,data] = rt_test_improvement(P, out, BFC);

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
    metrics(k) = rt_test_improvement(P, out, BFC, true);
   
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

%% new parameter set 
P2 = P;
P2.name = 'dist4_curv';
P2.distance_trans_bone = 10 * P2.lambda;
P2.distance_bone_pixel = 15 * P2.lambda;
P2.z_pixel = zpix_fun(P2);

%%
Ptest = P2;
[out, BFC] = rt_compare(Ptest);
figs = rt_compare_plot(Ptest, out, BFC);



%% save figes
% output_folder = 'figs';
% [~, ~] = mkdir(output_folder);
% save_fig_fun = @(f) exportgraphics(f, fullfile(output_folder, [f.Name '.png']), 'Resolution', 300);
% save_figs_fun = @(fc) cellfun(save_fig_fun, fc);
% save_figs_fun(figs)


%% Total attenuation along path
BFC.medium_attenuation = [1.0, 0.54, 6.9, 0.6]; % dB/(MHz*cm)
att_fun = @(rays) sum(BFC.medium_attenuation .* [rays.length] * 1e2 * (P.Fc / 1e6));

attenuation_tx = att_fun(out.ac_rays_tx);
attenuation_rx = cellfun(att_fun, out.ac_rays_rx);
attenuation_rt = attenuation_rx + attenuation_tx;
attenuation_rt_mag = db2mag(-attenuation_rt);

figure(11);clf;
plot(P.x_piezo * 1e3, attenuation_rx, '-', 'LineWidth', 2)
hold on 
plot(P.x_piezo * 1e3, attenuation_rt, '-', 'LineWidth', 2)

%% Transmission coefficients

[TF_tx, TF_all] = rt.raytheory.transmission_coeff(BFC, out.ac_rays_tx, true);
TF_rx = cellfun(@(rays) rt.raytheory.transmission_coeff(BFC, rays, false), out.ac_rays_rx);
TF_rt = TF_tx * TF_rx;


%% Generate RF
P.Fs = 40 * P.Fc; % Sampling frequency
P.num_cycles = 3;

% Set the impulse response and excitation of the emit aperture
image_pulse = sin(2 * pi * P.Fc * (0:1 / P.Fs:P.num_cycles / P.Fc)); % BP66 %,  sin(2*pi*fc*(0:1/fs:1/fc)); % BP100 %
image_pulse = image_pulse .* hamming(length(image_pulse))';
image_pulse_env = envelope(image_pulse);

[~, imax] = max(image_pulse_env);
ttp = imax / P.Fs;

% figure(99);
% plot(image_pulse)
% hold on
% plot(image_pulse_env)


error_rt_tof = out.error_tof_round_trip;
error_rt_tof(error_rt_tof==0) = NaN
error_rt_tof = fillmissing(error_rt_tof, 'spline');
error_rt_tof = error_rt_tof-min(error_rt_tof);

t_min = 0;
% t_max = max(out.ac_tof_round_trip(:)) * 1.5; % full RF
t_max = max(error_rt_tof) * 5; % delay with error
t_vec = t_min:1 / P.Fs:t_max;
N = numel(t_vec);

RF_pulse = repmat(image_pulse, P.num_elements, 1).';
RF_pulse(N, 1) = 0; % pad with zeros

% apply apodization to mimick attenuation
% RF_pulse = RF_pulse .* hamming(P.num_elements).';
RF_pulse = RF_pulse .* attenuation_rt_mag;
RF_pulse = RF_pulse .* TF_rt;

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

RF_delayed_nc = pulse_delaying_RF(RF_pulse, error_rt_tof, P.Fs);
figure(102); clf;
imagesc(P.x_piezo * 1e3, t_vec * 1e6, RF_delayed_nc)
colormap bone

% RF_pulse_sum_ac = sum(RF_pulse, 2);
% RF_pulse_sum_nc = sum(RF_delayed_nc, 2);
% env_sum_ac = envelope(RF_pulse_sum_ac);
% env_sum_nc = envelope(RF_pulse_sum_nc);
% peak_ac = max(env_sum_ac);
% peak_nc = max(env_sum_nc);
% fprintf('Intensity improvement all elements: %.2fx\n', peak_ac/peak_nc);

figure(1012); clf;
% subplot(211)
% plot(t_vec * 1e6, RF_pulse_sum_ac)
% hold on
% plot(t_vec * 1e6, RF_pulse_sum_nc)
% title('Sum all elements')

RF_pulse_sum_ac = sum(RF_pulse(:, out.ac_f_number_idx), 2);
RF_pulse_sum_ac_nof = sum(RF_pulse(:, out.nc_f_number_idx), 2);
RF_pulse_sum_nc = sum(RF_delayed_nc(:, out.nc_f_number_idx), 2);

env_sum_ac = envelope(RF_pulse_sum_ac);
env_sum_nc = envelope(RF_pulse_sum_nc);
env_sum_ac_nof = envelope(RF_pulse_sum_ac_nof);
peak_ac = max(env_sum_ac);
peak_nc = max(env_sum_nc);
peak_ac_nof = max(env_sum_ac_nof);
fprintf('Intensity improvement f-masked:            %.2fx\n', peak_ac/peak_nc);
fprintf('Intensity improvement f-masked (DASmask):  %.2fx\n', peak_ac_nof/peak_nc);

% subplot(212)
plot(t_vec * 1e6, RF_pulse_sum_ac)
hold on
plot(t_vec * 1e6, RF_pulse_sum_nc)
title('Taking f-number into account')
legend('Aberration Corrected', 'No Correction')

return; 

%% Determine
% - theta_i: incident angle
% - theta_o: transmitted / refracted angle

ray = out.ac_rays_rx{end}(2)
% ray = out.ac_rays_tx(1);

rays = out.ac_rays_tx;
transmission_factor = @(Z1, Z2, ti, to) 2 * Z2 * cos(ti) / (Z2 * cos(ti) + Z1 * cos(to));
transmission_factor = @(Z1, Z2, ti, to) 2 * Z2 / (Z2 + Z1);

TF_factors = 0;

% ONLY FORWARD!
for k = 1:numel(rays) - 1
    ray = rays(k);
    theta_i = acos(dot(ray.dir, ray.normal));
    theta_o = acos(dot(ray.dir_refracted, ray.normal));
    Z1 = BFC.medium_impendace(k);
    Z2 = BFC.medium_impendace(k + 1);
    TF_factors(k) = transmission_factor(Z1, Z2, theta_i, theta_o);
end

% TODO: RETURN, need swapping of TF.


attenuation_dB = sum(BFC.medium_attenuation .* [rays.length] * 1e2 * (P.Fc / 1e6));
attenuation_mag = db2mag(-attenuation_dB);


%% determine PSF shape;
[out, BFC] = rt_compare(P);
[metrics, data] = rt_test_improvement(P, out, BFC);

span_lambda = 10;

xvec = P.x_pixel + (-span_lambda*P.lambda:P.lambda/2:span_lambda*P.lambda);
zvec = P.z_pixel + (-span_lambda*P.lambda:P.lambda/2:span_lambda*P.lambda);
Nx = numel(xvec); Nz = numel(zvec);

IQ_full = hilbert(data.RF);
IQ_xax_ac = zeros(1, Nx);
IQ_xax_nc = zeros(1, Nx);
IQ_zax = zeros(Nz, 1);

kx = find(xvec==0,1)+3;
% kx = 1;
[tof_ac, tof_nc] = rt_trace_tof(P, BFC, xvec(kx), P.z_pixel);
% delta_tof_nc = (2*P.bone_thickness/P.bone_wavespeed - 2*P.bone_thickness/1540);
delta_tof_nc = 2*(1/P.bone_wavespeed - 1/1540)*P.bone_thickness;

% delta_tof_nc = min(tof_ac) - min(tof_nc);

figure(99);clf;
imagesc(P.x_piezo*1e3, data.tvec*1e6, data.RF);
hold on 
plot(P.x_piezo*1e3, (tof_ac)*1e6)
plot(P.x_piezo*1e3, (tof_nc+delta_tof_nc)*1e6)

%%

% time_remaining_progbar_ui(0, Nx)
for kx = 1:Nx
    % kx=  find(xvec==0,1);
    [tof_ac, tof_nc] = rt_trace_tof(P, BFC, xvec(kx), P.z_pixel);
    IQ_interp = interp1_per_channel(data.tvec, IQ_full, tof_ac);
    IQ_xax_ac(kx) = sum(IQ_interp,'omitmissing');

    % delta_tof_nc = min(tof_ac) - min(tof_nc);

    IQ_interp = interp1_per_channel(data.tvec, IQ_full, tof_nc+delta_tof_nc);
    IQ_xax_nc(kx) = sum(IQ_interp,'omitmissing');
    % time_remaining_progbar_ui(kx, Nx)
end

% dont care about Z? 
% % time_remaining_progbar_ui(kz, Nz)
% for kz = 1:Nz
%     [tof_ac, tof_nc] = rt_trace_tof(P, BFC, P.x_pixel, zvec(kz));
%     IQ_interp = interp1_per_channel(data.tvec, IQ_full, tof_ac);
%     IQ_zax(kz) = sum(IQ_interp,'omitmissing');
%     % time_remaining_progbar_ui(kz, Nz)
% end


figure(14); clf
plot(xvec*1e3, abs(IQ_xax_nc))
hold on 
plot(xvec*1e3, abs(IQ_xax_ac))
legend('NC', 'AC')
% figure(15); plot(zvec, abs(IQ_zax))


% end



