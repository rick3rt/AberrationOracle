clc
clear all
% close all

% set transducer frequency
P.Fc = 7.5e6;
P.wavespeed_skin = 1600; % m/s
P.wavespeed_brain = 1570; % m/s

% bone properties
bone.thickness = 2e-3; % m
bone.wavespeed = 3000; % m/s
bone.curvature = 0; % m^-1
bone.distance = 2e-3; % m

% define transducer
Trans.Fc = P.Fc;
Trans.lambda = 1540 / P.Fc;
Trans.NElems = 128;
Trans.XPiezo = (0:Trans.NElems-1) .* Trans.lambda - (Trans.NElems-1) / 2 * Trans.lambda;
Trans.ZPiezo = 0 * Trans.XPiezo;
Trans.lens_thickness = 5 * Trans.lambda;
Trans.lens_wavespeed = 1000; % m/s

% from bone, make realization of interfaces
xax = linspace(-Trans.NElems * Trans.lambda, Trans.NElems * Trans.lambda, 100); %
bone.outer_surface = (Trans.lens_thickness + bone.distance) + bone.curvature * xax .^ 2;
bone.inner_surface = (Trans.lens_thickness + bone.distance + bone.thickness) + bone.curvature * xax .^ 2;

% plot transducer and bone
figure(1); clf;
scatter(Trans.XPiezo, Trans.ZPiezo, 'rx');
yline(0);
yline(Trans.lens_thickness)
hold on
plot(xax, bone.outer_surface, 'k-')

plot(xax, bone.inner_surface, 'k-')
%daspect([1 1 1])
set(gca, 'YDir', 'reverse')
axis equal

%% Prepare ray tracer.
% define sound speed per layer
BFC.medium_soundspeeds = [Trans.lens_wavespeed, 
                          P.wavespeed_skin
                          bone.wavespeed
                          P.wavespeed_brain];

% define  interfaces between tissue layers
BFC.medium_interfaces = {Trans.lens_thickness,
                         [bone.curvature, 0, Trans.lens_thickness + bone.distance],
                         [bone.curvature, 0, Trans.lens_thickness + bone.distance + bone.thickness]};

BFC.medium_soundspeeds = BFC.medium_soundspeeds(:).'; % guarantee row vector
BFC.medium_interfaces = BFC.medium_interfaces(:).'; % guarantee row vector

% deremine interface derivatives for normal computation
BFC.medium_interface_derivatives = cellfun(@rt.util.interface_derivative, BFC.medium_interfaces, 'UniformOutput', false);

% define a starting point and a image point, and end point for ray tracing
ke_start = 1; % start at first element
x_start = Trans.XPiezo(ke_start);
z_start = Trans.ZPiezo(ke_start);

x_target = 0e-3; % target point
z_target = 50 * Trans.lambda; % target point

ke_end = Trans.NElems; % end at last element
x_end = Trans.XPiezo(ke_end);
z_end = Trans.ZPiezo(ke_end);


% define reconstruction grid
BFC.XPiezo = Trans.XPiezo;
BFC.ZPiezo = Trans.ZPiezo;
BFC.XRecon = Trans.XPiezo;
BFC.ZRecon = 0:Trans.lambda/2:100*Trans.lambda;
BFC.LensThickness =  Trans.lens_thickness;

% plot the medium
figure(1); clf;
hold on
rt.plot.medium(BFC, x_start, z_start, x_target, z_target, x_end, z_end);
set(gca,'YDir','reverse')

%% Trace rays 
to_layer = 4;
[rays_tx, tof_tx, theta_tx] = rt.ray_bending(x_start, z_start, x_target, z_target, BFC, to_layer);
[rays_rx, tof_rx, theta_rx] = rt.ray_bending(x_end, z_end, x_target, z_target, BFC, to_layer);

tof_tx_c0 = vecnorm([x_target;z_target]-[x_start;z_start])/1540;
tof_rx_c0 = vecnorm([x_end;z_end]-[x_target;z_target])/1540;

%% plot rays


figure(2); clf;
hold on
rt.plot.medium(BFC, x_start, z_start, x_target, z_target, x_end, z_end);
set(gca,'YDir','reverse')
rt.plot.rays(rays_tx);
rt.plot.rays(rays_rx);



figure(3);
bar([tof_tx tof_rx;
    tof_tx_c0 tof_rx_c0 ]*1e6)
ylabel('TOF \mus')
xticklabels({'Ray Tracing', 'Conventional'})
legend('TX','RX')



%% Determine for one source and all receivers

x_start = -0.05;
z_start = -1;

x_target = 0e-3; % target point
z_target = (bone.distance + bone.thickness) + 20 * Trans.lambda; % target point

[rays_tx, tof_tx, theta_tx] = rt.ray_bending(x_start, z_start, x_target, z_target, BFC, to_layer);

tof_rx_all = zeros(1, Trans.NElems);
theta_rx_all = zeros(1, Trans.NElems);
rays_rx_all = cell(1, Trans.NElems);
for ke = 1:Trans.NElems
    x_end = Trans.XPiezo(ke);
    z_end = Trans.ZPiezo(ke);
    [rays_rx, tof_rx, theta_rx] = rt.ray_bending(x_end, z_end, x_target, z_target, BFC, to_layer);
    tof_rx_all(ke) = tof_rx;
    theta_rx_all(ke) = theta_rx;
    rays_rx_all{ke} = rays_rx;
end


ind_valid = find(~isnan(tof_rx_all));
tof_round_trip = tof_tx + tof_rx_all;




% homogenous assumption
c0 = 1540;
tof_round_trip_c0 =  vecnorm([x_target;z_target]-[x_start;z_start])/c0 + ...
    vecnorm([Trans.XPiezo;Trans.ZPiezo]-[x_target;z_target])/c0;

rx_vec = [x_target;z_target]-[Trans.XPiezo;Trans.ZPiezo];
theta_rx_all_c0 = atan2(rx_vec(2,:), rx_vec(1,:)) - pi/2;



figure(2); clf;
hold on
rt.plot.medium(BFC, x_start, z_start, x_target, z_target, x_end, z_end);
set(gca,'YDir','reverse')
rt.plot.rays(rays_tx);
rt.plot.rays(rays_rx_all{ind_valid(1)});
rt.plot.rays(rays_rx_all{ind_valid(end)});
daspect([1 1 1])
xlim(rt.util.minmax(BFC.XRecon)*1e3); ylim(rt.util.minmax(BFC.ZRecon)*1e3)

% rt.plot.rays(rays_rx);

FNumber = 1.5; 
half_opening_angle_rad = atan(1/2/FNumber);

f_number_mask = find(abs(theta_rx_all) < half_opening_angle_rad);


figure(10);clf;
plot(Trans.XPiezo, tof_round_trip_c0)
hold on 
plot(Trans.XPiezo, tof_round_trip)


% get rid of axial shift:
tof_round_trip0 = tof_round_trip - min(tof_round_trip);
tof_round_trip_c0_0 = tof_round_trip_c0 - min(tof_round_trip_c0);


figure(10);clf;
wave_period = 1/P.Fc;


subplot(321)
plot(Trans.XPiezo, tof_round_trip_c0_0)
hold on 
plot(Trans.XPiezo, tof_round_trip0)
xlim(minmax(Trans.XPiezo))

subplot(322)
plot(Trans.XPiezo, tof_round_trip_c0_0 - tof_round_trip0)
yline(wave_period)
hold on 
xlim(minmax(Trans.XPiezo))

subplot(312)
plot(Trans.XPiezo, (tof_round_trip_c0_0 - tof_round_trip0)./wave_period*100)
xline(Trans.XPiezo(f_number_mask([1 end])),'r')
xlim(minmax(Trans.XPiezo))
title('Relative Error wrt Wave Period (%)')
ylabel('Relative error in %')

subplot(313)
plot(Trans.XPiezo, rad2deg(abs(theta_rx_all)))
hold on 
plot(Trans.XPiezo, rad2deg(abs(theta_rx_all_c0)),'r')
legend('Reception Angle Refraction Corrected', 'Reception Angle Conventional DAS')

yline(rad2deg(half_opening_angle_rad),'r', 'HandleVisibility','off')
xline(Trans.XPiezo(f_number_mask([1 end])),'r','HandleVisibility','off')
xlim(minmax(Trans.XPiezo))

