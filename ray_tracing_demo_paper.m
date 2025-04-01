clc
clear
% close all

import rt.*
import rt.plot.*

% profile on
%% define transfucer and Reconstruction grid

wvl = 1540/15e6;
grid_size = wvl;

% the transducer
BFC.NElems = 96;
BFC.XPiezo = (1:BFC.NElems) .* wvl - BFC.NElems / 2 * wvl;
BFC.ZPiezo = 0 * BFC.XPiezo;
BFC.LensThickness = 10 * wvl;

% reconstruction grid
% BFC.XRecon = (1:BFC.NElems) .* wvl - BFC.NElems * grid_size/2;
BFC.XRecon = BFC.XPiezo(1):grid_size:BFC.XPiezo(end);
BFC.ZRecon = 0:grid_size:7e-3;
BFC.Nx = numel(BFC.XRecon);
BFC.Nz = numel(BFC.ZRecon);

% Some points sources
% BFC.NSources = 9;
% BFC.ZS = -10 * ones(1, BFC.NSources) * 1e-3;
% BFC.XS = linspace(-2, 2, BFC.NSources) * 1e-3;
% BFC.addToDelay = zeros(1, BFC.NSources);
BFC.NSources = 1;
BFC.ZS = -1; % -5e-3; BFC.ZPiezo(32 + 16);
BFC.XS = -20e-3; BFC.XPiezo(32 + 16);
BFC.addToDelay = zeros(1, BFC.NSources);

% setup the medium
% USED FOR FIRST FIGURE
BFC.medium_soundspeeds = [1000 1600 3200 1570];
BFC.medium_interfaces = {BFC.LensThickness, [-2500 25 0.05 2.4e-3], [-1600 25 0.01 4.5e-3]};
% ALSO NICE:
% BFC.medium_interfaces = {BFC.LensThickness, [1200 20 0.0 3e-3], [1000 20 -0.017 5.e-3]};
% BFC.medium_interfaces = {BFC.LensThickness, [-1200 20 0.0 3e-3], [-1000 20 0.017 5.e-3]};

% determine interface derivatives for normal computation
BFC.medium_interface_derivatives = cellfun(@rt.util.interface_derivative, BFC.medium_interfaces, 'UniformOutput', false);

%% Plot Transducer and Tissue Interfaces

% PLOT transducer and img grid

ks = 1;
x_start = BFC.XS(ks);
z_start = BFC.ZS(ks);
x_target = 2e-3;
z_target = 1.1e-3;

ke = 80;
x_end = BFC.XPiezo(ke);
z_end = BFC.ZPiezo(ke);

figure(1); clf;
hold on
rt.plot.medium(BFC, x_start, z_start, x_target, z_target, x_end, z_end);
xlim(rt.util.minmax(BFC.XRecon) * 1e3); ylim(rt.util.minmax(BFC.ZRecon) * 1e3)
daspect([1 1 1])
set(gca, 'ydir', 'reverse')

%% Find Path From Start to Pixel Position
to_layer = 2;
tic
[rays_tx, tof_tx, theta_tx] = rt.ray_bending(x_start, z_start, x_target, z_target, BFC, to_layer);
toc
tic
[rays_rx, tof_rx, theta_rx] = rt.ray_bending(x_end, z_end, x_target, z_target, BFC, to_layer);
toc
% cost fun
cost_fun = @(theta) rt.ray_bending_tof(theta, x_start, z_start, x_target, z_target, BFC, to_layer);
% cost_fun = @(theta) rt.ray_bending_tof(theta, x_end, z_end, x_target, z_target, BFC, to_layer);

% refine angle range
z_lens = BFC.LensThickness;
x_extend = rt.util.minmax(BFC.XRecon);
r_left = [x_extend(1) - x_start; z_lens - z_start];
r_right = [x_extend(2) - x_start; z_lens - z_start];
angle_left = atan2(r_left(2), r_left(1)) - pi / 2;
angle_right = atan2(r_right(2), r_right(1)) - pi / 2;
theta_vals = linspace(angle_right, angle_left, 1000);

% coarse grid search
% theta_vals = linspace(-pi / 2, pi / 2, 50);

tof_grid = arrayfun(cost_fun, theta_vals);

[~, imin] = min(tof_grid);
theta0 = theta_vals(imin);

% theta0 = -0.0002;
% [tof, rays_tx] = cost_fun(theta0);

figure(99); clf; hold on;
plot(theta_vals, tof_grid, '.-')
scatter(theta_tx, tof_tx);
% ylim([0 10e-6])
xlim(rt.util.minmax(theta_vals))

% plot the rays
figure(2); clf;
hold on
rt.plot.medium(BFC, x_start, z_start, x_target, z_target, x_end, z_end);
rt.plot.rays(rays_tx);
rt.plot.rays(rays_rx);
xlim(rt.util.minmax(BFC.XRecon) * 1e3); ylim(rt.util.minmax(BFC.ZRecon) * 1e3)
daspect([1 1 1])
set(gca, 'ydir', 'reverse')

return

%% Compute Refraction Corrected Time Delays for all Pixels

% make layer mask
[X, Z] = meshgrid(BFC.XRecon, BFC.ZRecon);
layer_mask = ones(BFC.Nz, BFC.Nx);
for ki = 1:numel(BFC.medium_interfaces)
    z_int = rt.util.segeval(BFC.medium_interfaces{ki}, BFC.XRecon);
    msk = Z > z_int;
    layer_mask(msk) = layer_mask(msk) + 1;
end

figure(3); clf;
imagesc(BFC.XRecon * 1e3, BFC.ZRecon * 1e3, layer_mask)
hold on
rt.plot.medium(BFC);
daspect([1 1 1])

%% Compute Refraction Corrected Time Delays for all Pixels

skip = 1; % if 1, compute every pixel, if higher, number of elements to skip
fprintf('Computing LUT for grid size %i x %i, compute points %i x %i\n\n\n', ...
    BFC.Nx, BFC.Nz, round(BFC.Nx / skip), round(BFC.Nz / skip));

xr = BFC.XRecon(1:skip:end);
zr = BFC.ZRecon(1:skip:end);
Nxr = numel(xr);
Nzr = numel(zr);
[Xr, Zr] = meshgrid(xr, zr);
layer_mask_r = layer_mask(1:skip:end, 1:skip:end);

LUT_T = NaN(Nzr, Nxr, BFC.NSources);
LUT_R = NaN(Nzr, Nxr, BFC.NElems);
LUT_T_interp = NaN(BFC.Nz, BFC.Nx, BFC.NSources);
LUT_R_interp = NaN(BFC.Nz, BFC.Nx, BFC.NElems);

rt.util.progressbar_ui(0, Nxr);
for kx = 1:Nxr
    for kz = 1:Nzr
        to_layer = layer_mask_r(kz, kx);
        if to_layer == 1; continue; end % skip lens
        x_target = xr(kx);
        z_target = zr(kz);

        for ks = 1:BFC.NSources
            [~, tof] = rt.ray_bending(BFC.XS(ks), BFC.ZS(ks), x_target, z_target, BFC, to_layer);
            LUT_T(kz, kx, ks) = tof;
        end

        for ke = 1:1 %BFC.NElems
            [~, tof] = rt.ray_bending(BFC.XPiezo(ke), BFC.ZPiezo(ke), x_target, z_target, BFC, to_layer);
            LUT_R(kz, kx, ke) = tof;
        end
    end
    rt.util.progressbar_ui(kx, Nxr);
end

% interpolate to original image grid
for ks = 1:BFC.NSources
    LUT_T_interp(:, :, ks) = interp2(Xr, Zr, LUT_T(:, :, ks), X, Z, 'cubic');
end
for ke = 1:BFC.NElems
    LUT_R_interp(:, :, ke) = interp2(Xr, Zr, LUT_R(:, :, ke), X, Z, 'cubic');
end

%%
ks = 1; ke = 1;

figure(99); clf;
subplot(121)
imagesc(BFC.XRecon * 1e3, BFC.ZRecon * 1e3, LUT_T_interp(:, :, ks))
hold on
contour(BFC.XRecon * 1e3, BFC.ZRecon * 1e3, LUT_T_interp(:, :, ks), 50, 'k')
rt.plot.medium(BFC)
daspect([1 1 1])

subplot(122)
imagesc(BFC.XRecon * 1e3, BFC.ZRecon * 1e3, LUT_R_interp(:, :, ke))
hold on
contour(BFC.XRecon * 1e3, BFC.ZRecon * 1e3, LUT_R_interp(:, :, ke), 50, 'k')
daspect([1 1 1])
rt.plot.medium(BFC)

figure(101); clf;
imagesc(BFC.XRecon * 1e3, BFC.ZRecon * 1e3, layer_mask)
hold on
rt.plot.medium(BFC);
daspect([1 1 1])
% scatter(Xr(:)*1e3, Zr(:)*1e3,'kx')
