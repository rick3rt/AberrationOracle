clc
clear all
close all

%% PARAMETERS

% set transducer frequency
P.Fc = 7.5e6;
P.lambda = 1540 / P.Fc;
P.lens_wavespeed = 1000; % m/s
P.lens_thickness = 5 * P.lambda; % m
P.skin_wavespeed = 1600; % m/s
P.brain_wavespeed = 1570; % m/s

% bone properties
P.bone_wavespeed = 2780; % m/s
P.bone_thickness = 2e-3; % m
P.bone_curvature = 0; % 1/m

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
P.z_source = P.z_piezo(ke);
% P.x_source = 0;
P.z_source = -20 * P.lambda;
P.theta_source = deg2rad(5);

P.x_pixel = 0; % center of image
P.z_pixel = P.lens_thickness + P.distance_trans_bone + P.bone_thickness + P.distance_bone_pixel;

% f-number in reconstruction
P.f_number = 2;

%% ==================

% prepare ray tracer
% ===========================================================
c0 = 1540;

% define sound speed per layer
BFC.medium_soundspeeds = [P.lens_wavespeed,
                          P.skin_wavespeed,
                          P.bone_wavespeed,
                          P.brain_wavespeed];

BFC.medium_density = [1.0, 1.0, 1.0, 1.0]; % kg/m^3

% define  interfaces between tissue layers
BFC.medium_interfaces = {P.lens_thickness,
                         [P.bone_curvature, 0, P.lens_thickness + P.distance_trans_bone],
                         [P.bone_curvature, 0, P.lens_thickness + P.distance_trans_bone + P.bone_thickness]};

BFC.medium_soundspeeds = BFC.medium_soundspeeds(:).'; % guarantee row vector
BFC.medium_density = BFC.medium_density(:).'; % guarantee row vector
BFC.medium_impendace = BFC.medium_soundspeeds .* BFC.medium_density; % kg/(m^2 s)

BFC.medium_interfaces = BFC.medium_interfaces(:).'; % guarantee row vector

% deremine interface derivatives for normal computation
BFC.medium_interface_derivatives = cellfun(@rt.util.interface_derivative, BFC.medium_interfaces, 'UniformOutput', false);

% put other variables in BFC
BFC.XPiezo = P.x_piezo;
BFC.ZPiezo = P.z_piezo;
BFC.XRecon = P.x_piezo;
BFC.ZRecon = 0:P.lambda / 2:100 * P.lambda;
BFC.LensThickness = P.lens_thickness;

% compute source location - homogeneous and lens speed of sound
xs = P.z_source * tan(P.theta_source);
vs = [xs; P.z_source];
ve = [P.x_piezo; P.z_piezo];
ds = vecnorm(ve - vs);
tx_delay_lens = ds / P.lens_wavespeed;
tx_delay_c0 = ds / c0;
P.x_source = xs;
P.tx_add_to_ac = min(tx_delay_lens);
P.tx_add_to_nc = min(tx_delay_c0);

% plot the medium
% ===========================================================
% figure(1); clf;
% hold on
% rt.plot.medium(BFC);
% scatter(P.x_pixel * 1e3, P.z_pixel * 1e3, 'ko', 'filled');
% set(gca, 'YDir', 'reverse')
% daspect([1 1 1])
% ===========================================================

% ray tracing from source to pixel, and pixel to all elements
% ===========================================================
to_layer = 4; % trace to brain layer
[rays_tx, tof_tx, theta_tx] = rt.ray_bending(P.x_source, P.z_source, P.x_pixel, P.z_pixel, BFC, to_layer);

tof_tx = tof_tx - P.tx_add_to_ac;

tof_rx_all = zeros(1, P.num_elements);
theta_rx_all = zeros(1, P.num_elements);
rays_rx_all = cell(1, P.num_elements);
for ke = 1:P.num_elements
    x_end = P.x_piezo(ke);
    z_end = P.z_piezo(ke);
    [rays_rx, tof_rx, theta_rx] = rt.ray_bending(x_end, z_end, P.x_pixel, P.z_pixel, BFC, to_layer);
    tof_rx_all(ke) = tof_rx;
    theta_rx_all(ke) = theta_rx;
    rays_rx_all{ke} = rays_rx;
end

ind_valid = find(~isnan(tof_rx_all)); % valid rays in reception.
tof_round_trip = tof_tx + tof_rx_all;

% calculate time of flight for homogenous medium
% ===========================================================

delta_tof_nc = (1 / P.bone_wavespeed -1/1540) * P.bone_thickness;
dz = delta_tof_nc * c0;

v_pixel = [P.x_pixel; P.z_pixel + dz];
v_source = [P.x_source; P.z_source];
v_piezo = [P.x_piezo; P.z_piezo];

tof_tx_nc = vecnorm(v_pixel - v_source) / c0 - P.tx_add_to_nc;
tof_rx_nc = vecnorm(v_piezo - v_pixel) / c0;
tof_round_trip_nc = tof_tx_nc + tof_rx_nc;

v_rx = v_pixel - v_piezo; % return vector
theta_rx_all_nc = atan2(v_rx(2, :), v_rx(1, :)) - pi / 2;

% get rid of axial shift:
tof_round_trip0 = tof_round_trip - min(tof_round_trip);
tof_round_trip_c0_0 = tof_round_trip_nc - min(tof_round_trip_nc);

wave_period = 1 / P.Fc;

% half opening angle corresponding to f-number
half_opening_angle_rad = atan(1/2 / P.f_number);
f_number_idx = find(abs(theta_rx_all) < half_opening_angle_rad);
f_number_idx_nc = find(abs(theta_rx_all_nc) < half_opening_angle_rad);

% collect results
% ===========================================================
out.ac_rays_tx = rays_tx;
out.ac_tof_tx = tof_tx;
% out.ac_theta_tx = theta_tx; %
out.ac_rays_rx = rays_rx_all;
out.ac_tof_rx = tof_rx_all;
out.ac_theta_rx = theta_rx_all;
out.ac_f_number_idx = f_number_idx;
out.ac_valid_ray_idx = ind_valid;
out.ac_tof_round_trip = tof_round_trip;

out.nc_theta_rx = theta_rx_all_nc;
out.nc_f_number_idx = f_number_idx_nc;
out.nc_tof_round_trip = tof_round_trip_nc;

out.error_tof_round_trip = tof_round_trip_c0_0 - tof_round_trip0;
out.error_tof_round_trip_relative = out.error_tof_round_trip ./ wave_period;
out.half_opening_angle_rad = half_opening_angle_rad;

%% Plot this output
P.name = 'test';
rt_compare_plot(P, out, BFC)

%%
apply_attenuation = 1;

% Elemenet directivity
lambda = P.lambda;
W = lambda * 0.9; % small kerf
directivity_fun = @(theta) cos(theta) .* sinc(pi * W / lambda .* sin(theta));

% theta_test = linspace(-pi/2, pi/2, 128);
% figure(10);clf;
% plot(rad2deg(theta_test), directivity_fun(theta_test))

%% Total attenuation along path

att_fun = @(rays) sum(BFC.medium_attenuation .* [rays.length] * 1e2 * (P.Fc / 1e6));

attenuation_tx = att_fun(out.ac_rays_tx);
attenuation_rx = cellfun(att_fun, out.ac_rays_rx);
attenuation_rt = attenuation_rx + attenuation_tx;
attenuation_rt_mag = db2mag(-attenuation_rt);

% figure(11);clf;
% plot(P.x_piezo * 1e3, attenuation_rx, '-', 'LineWidth', 2)
% hold on
% plot(P.x_piezo * 1e3, attenuation_rt, '-', 'LineWidth', 2)

%% Transmission coefficients
TF_tx = rt.raytheory.transmission_coeff(BFC, out.ac_rays_tx, true);
TF_rx = cellfun(@(rays) rt.raytheory.transmission_coeff(BFC, rays, false), out.ac_rays_rx);
TF_rt = TF_tx * TF_rx;

% figure(1);clf;
% plot(TF_rt)

%% Generate RF
P.Fs = 20 * P.Fc; % Sampling frequency
P.num_cycles = 3;

% Set the impulse response and excitation of the emit aperture
image_pulse = sin(2 * pi * P.Fc * (0:1 / P.Fs:P.num_cycles / P.Fc)); % BP66 %,  sin(2*pi*fc*(0:1/fs:1/fc)); % BP100 %
image_pulse = image_pulse .* hamming(length(image_pulse))';
image_pulse_env = envelope(image_pulse);

[peak_pulse, imax] = max(image_pulse_env);
ttp = imax / P.Fs;

figure(99);
plot(image_pulse)
hold on
plot(image_pulse_env)

error_rt_tof = out.error_tof_round_trip;
% BUG FXED NOT NEEDED ANYMORE?
% error_rt_tof(error_rt_tof==0) = NaN;
% error_rt_tof = fillmissing(error_rt_tof, 'spline');
% error_rt_tof = error_rt_tof-min(error_rt_tof);
%

t_min = 0;
%t_max = max(error_rt_tof) * 20; % delay with error
t_max = 20 / P.Fc;

t_vec = t_min:1 / P.Fs:t_max;
N = numel(t_vec);

RF_pulse = repmat(image_pulse, P.num_elements, 1).';
RF_pulse(N, 1) = 0; % pad with zeros

% apply apodization to mimick attenuation
% RF_pulse = RF_pulse .* hamming(P.num_elements).';

% apply 'real' attenuation etc.
if apply_attenuation
    RF_pulse = RF_pulse .* attenuation_rt_mag; % attenuation
    RF_pulse = RF_pulse .* TF_rt; % transmission factors
    RF_pulse = RF_pulse .* directivity_fun(out.ac_theta_rx); % element directivity
end
RF_delayed_nc = pulse_delaying_RF(RF_pulse, error_rt_tof, P.Fs);

cmap_lines = lines(2);
figure(102); clf;
subplot(121)
imagesc(P.x_piezo * 1e3, t_vec * 1e6, RF_delayed_nc)
colormap bone
xline(P.x_piezo(out.ac_f_number_idx([1 end])) * 1e3, 'color', cmap_lines(1, :));
xline(P.x_piezo(out.nc_f_number_idx([1 end])) * 1e3, 'color', cmap_lines(2, :));

RF_pulse_sum_ac = sum(RF_pulse(:, out.ac_f_number_idx), 2);
RF_pulse_sum_nc = sum(RF_delayed_nc(:, out.nc_f_number_idx), 2);
% RF_pulse_sum_nc = sum(RF_delayed_nc(:, out.ac_f_number_idx), 2);

env_sum_ac = envelope(RF_pulse_sum_ac);
env_sum_nc = envelope(RF_pulse_sum_nc);
peak_ac = max(env_sum_ac);
peak_nc = max(env_sum_nc);
fprintf('Intensity improvement f-masked:     %.2fx\n', peak_ac / peak_nc);

if numel(t_vec) ~= numel(RF_pulse_sum_ac)
    error('t_max is too short')
end

subplot(122)
plot(t_vec * 1e6, RF_pulse_sum_ac)
hold on
plot(t_vec * 1e6, RF_pulse_sum_nc)
title('Taking f-number into account')
legend('Aberration Corrected', 'No Correction')

%% collect  metrics
metrics.peak_pulse = peak_pulse;
metrics.ac_peak = peak_ac;
metrics.nc_peak = peak_nc;
metrics.improvement = peak_ac / peak_nc;
metrics.ac_RF_sum = RF_pulse_sum_ac;
metrics.nc_RF_sum = RF_pulse_sum_nc;

%% OPTIONAL SHOW FULL RF
test_resolution = 0;

% if nargout ~= 2 && ~test_resolution
%     return
% end

RF_full = RF_pulse;
t_max_full = max(out.ac_tof_round_trip(:)) * 1.5; % full RF
t_vec_full = t_min:1 / P.Fs:t_max_full;
RF_full(numel(t_vec_full), 1) = 0; % pad with zeros

if apply_attenuation
    RF_full = RF_full .* attenuation_rt_mag; % attenuation
    RF_full = RF_full .* TF_rt; % transmission factors
    RF_full = RF_full .* directivity_fun(out.ac_theta_rx); % element directivity
end
RF_delayed_full = pulse_delaying_RF(RF_full, out.ac_tof_round_trip - ttp, P.Fs);
cmap_lines = lines(2);
% cmap_lines = flip(cmap_lines);

figure(101); clf;
imagesc(P.x_piezo * 1e3, t_vec_full * 1e6, RF_delayed_full)
colormap bone
hold on

plot(P.x_piezo * 1e3, out.ac_tof_round_trip * 1e6, '--', 'color', cmap_lines(1, :), 'LineWidth', 1);
plot(P.x_piezo(out.ac_f_number_idx) * 1e3, out.ac_tof_round_trip(out.ac_f_number_idx) * 1e6, '-', 'color', cmap_lines(1, :), 'LineWidth', 2);
plot(P.x_piezo * 1e3, out.nc_tof_round_trip * 1e6, '--', 'color', cmap_lines(2, :), 'LineWidth', 1);
plot(P.x_piezo(out.nc_f_number_idx) * 1e3, out.nc_tof_round_trip(out.nc_f_number_idx) * 1e6, '-', 'color', cmap_lines(2, :), 'LineWidth', 2);

data.RF = RF_delayed_full;
data.tvec = t_vec_full;
data.ttp = ttp;

if ~test_resolution
    return
end

% keyboard

%% Test resolution
fprintf('Testing Resolution...\n');

span_lambda = 5;
xvec = P.x_pixel + (-span_lambda * P.lambda:P.lambda / 4:span_lambda * P.lambda);
zvec = P.z_pixel + (-span_lambda * P.lambda:P.lambda / 4:span_lambda * P.lambda);
Nx = numel(xvec);
Nz = numel(zvec);

IQ_full = hilbert(data.RF);
IQ_xax_ac = zeros(1, Nx);
IQ_xax_nc = zeros(1, Nx);

% Abberration corrected PSF
time_remaining_progbar_ui(0, Nx)
for kx = 1:Nx
    % if 1; kx =  ceil(Nx/2);
    tof_ac = rt_trace_tof(P, BFC, xvec(kx), P.z_pixel);
    IQ_interp = interp1_per_channel(data.tvec, IQ_full, tof_ac, 'cubic');
    IQ_xax_ac(kx) = sum(IQ_interp, 'omitmissing');
    %IQ_interp = interp1_per_channel(data.tvec, IQ_full, tof_nc + dtof,'cubic');
    %IQ_xax_nc(kx) = sum(IQ_interp, 'omitmissing');

    time_remaining_progbar_ui(kx, Nx)
end

%%  Homogeneous PSF
span_lambda_x = 10;
span_lambda_z = 20;
xvec_nc = P.x_pixel + (-span_lambda_x * P.lambda:P.lambda / 4:span_lambda_x * P.lambda);
zvec = P.z_pixel + (-span_lambda_z * P.lambda:P.lambda / 4:span_lambda_z * P.lambda);
Nx = numel(xvec_nc); Nz = numel(zvec);
[X, Z] = meshgrid(xvec_nc, zvec);

tau_tx = hypot(X - P.x_source, Z - P.z_source) / c0 - P.tx_add_to_nc;
tau_rx = hypot(X - reshape(P.x_piezo, 1, 1, []), Z - reshape(P.z_piezo, 1, 1, [])) / c0;
theta_rx = atan2(Z - reshape(P.z_piezo, 1, 1, []), X - reshape(P.x_piezo, 1, 1, [])) - pi / 2;
half_opening_angle_rad = atan(1/2 / P.f_number);

IQ_nc = zeros(Nz, Nx);
for kx = 1:Nx
    tau = squeeze(tau_tx(:, kx, :) + tau_rx(:, kx, :));
    f_mask = squeeze(abs(theta_rx(:, kx, :)) < half_opening_angle_rad);
    IQ_interp = interp1_per_channel(data.tvec, IQ_full, tau, 'cubic');
    IQ_line = sum(IQ_interp .* f_mask, 2);
    IQ_nc(:, kx) = IQ_line;
end

[v, kz, kx] = maxij(abs(IQ_nc));
% debug plot
figure(123); clf;
imagesc(xvec_nc * 1e3, zvec * 1e3, flogc(IQ_nc), [-60 0])
daspect([1 1 1])
hold on
scatter(xvec_nc(kx) * 1e3, zvec(kz) * 1e3)
colormap bone

% interpolate to same grid as AC
IQ_xax_nc = IQ_nc(kz, :);
IQ_xax_nc = interp1(xvec_nc, IQ_xax_nc, xvec);

%%

% dont care about Z?
% % time_remaining_progbar_ui(kz, Nz)
% for kz = 1:Nz
%     [tof_ac, tof_nc] = rt_trace_tof(P, BFC, P.x_pixel, zvec(kz));
%     IQ_interp = interp1_per_channel(data.tvec, IQ_full, tof_ac);
%     IQ_zax(kz) = sum(IQ_interp,'omitmissing');
%     % time_remaining_progbar_ui(kz, Nz)
% end

% determine resolution.
[w_ac, y50_ac, x1_ac, x2_ac] = fwhm_center(xvec * 1e3, abs(IQ_xax_ac));
[w_nc, y50_nc, x1_nc, x2_nc] = fwhm_center(xvec * 1e3, abs(IQ_xax_nc));
% [fwhm_val_ac, fwhm_pos_ac] = fwhm(abs(IQ_xax_ac), xvec * 1e3);
% [fwhm_val_nc, fwhm_pos_nc] = fwhm(abs(IQ_xax_nc), xvec * 1e3);

cmap = lines(2);

f = figure(69); clf
f.Name = ['resolution_' P.name];
% figs{end+1} = f;
plot(xvec * 1e3, abs(IQ_xax_nc))
hold on
plot(xvec * 1e3, abs(IQ_xax_ac))
legend('NC', 'AC')

plot([x1_nc x2_nc], [y50_nc y50_nc], 'Color', cmap(1, :), 'HandleVisibility', 'off')
plot([x1_ac x2_ac], [y50_ac y50_ac], 'Color', cmap(2, :), 'HandleVisibility', 'off')

title(sprintf('Lateral Resolution - NC: %.3f mm  - AC: %.3f mm', w_nc, w_ac))

% resolution metrics
metrics.ac_res_x = w_ac; % in mm
metrics.nc_res_x = w_nc; % in mm
metrics.ac_IQ_x = IQ_xax_ac;
metrics.nc_IQ_x = IQ_xax_nc;
metrics.ac_IQ_x_peak = max(abs(IQ_xax_ac));
metrics.nc_IQ_x_peak = max(abs(IQ_xax_nc));



%% DEBUG RAY TRACER 



    % compute source location - homogeneous and lens speed of sound
    c0 = P.c0;
    xs = P.z_source * tan(P.theta_source);
    vs = [xs; P.z_source];
    ve = [P.x_piezo; P.z_piezo];
    ds = vecnorm(ve - vs);
    tx_delay_lens = ds / P.lens_wavespeed;
    tx_delay_c0 = ds / c0;
    P.x_source = xs;
    P.tx_add_to_ac = min(tx_delay_lens);
    P.tx_add_to_nc = min(tx_delay_c0);



    x_start = P.x_source;
    z_start = P.z_source;

    x_target = P.x_pixel    ;
    z_target = P.z_pixel    ;
    to_layer = 4;

    % cost fun
    cost_fun = @(theta) rt.ray_bending_tof(theta, x_start, z_start, x_target, z_target, BFC, to_layer);

    % determine angle range
    z_lens = BFC.LensThickness;
    x_extend = rt.util.minmax(BFC.XRecon);
    r_left = [x_extend(1) - x_start; z_lens - z_start];
    r_right = [x_extend(2) - x_start; z_lens - z_start];
    angle_left = atan2(r_left(2), r_left(1)) - pi/2;
    angle_right = atan2(r_right(2), r_right(1)) - pi/2;

    figure(80); clf; scatter(x_start, z_start);  hold on; scatter(x_extend(1), z_lens); % daspect([1 1 1])

    % coarse grid search
    NAngles = 50;
    theta_vals = linspace(angle_right, angle_left, NAngles);
    [tof_grid, rays_theta] = arrayfun(cost_fun, theta_vals, 'UniformOutput', 0);
    tof_grid = [tof_grid{:}];

    [~, imin] = min(tof_grid);
    theta0 = theta_vals(imin);

    % minimize with fminsearch
    options = optimset('TolX',1e-6);
    [theta_min,fval,exit_flag,out] = fminsearch(cost_fun, theta0,options);
    [tof, rays] = cost_fun(theta_min);


    % final check, if angle of ray(k).dir and ray(k-1).dir_refracted
    % withing range
    theta = acos(dot(rays(end).dir, rays(end-1).dir_refracted));
    if theta > 0.1
        tof = NaN; theta_min = NaN; 
    end



figure(99); clf; hold on;
plot(theta_vals, tof_grid, '.-')
scatter(theta_min, tof);
% ylim([0 10e-6])
xlim(rt.util.minmax(theta_vals))

x_end = 0; z_end= 0;

% plot the rays
figure(2); clf;
hold on
rt.plot.medium(BFC, x_start, z_start, x_target, z_target, x_end, z_end);
% rt.plot.rays(rays);
rt.plot.rays(rays_theta{imin-1});
% rt.plot.rays(rays_theta{imin});
rt.plot.rays(rays_theta{imin+1});
rt.plot.rays(rays);

% rt.plot.rays(rays_rx);
xlim(rt.util.minmax(BFC.XRecon) * 1e3); ylim(rt.util.minmax(BFC.ZRecon) * 1e3)
daspect([1 1 1])
set(gca, 'ydir', 'reverse')

%% debug call
cost_fun = @(theta) rt.ray_bending_tof(theta, x_start, z_start, x_target, z_target, BFC, to_layer);

cost_fun(theta0);