clear
clc

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
P.medium_density = [1.2, 1.02, 2.0, 1.001]; % kg/m^3
P.medium_attenuation = [0.0, 0.54, 6.9, 0.6]; % dB/(MHz*cm)

% relative distances
dist_trans_bone = 10 * P.lambda; % distance from transducer to bone
dist_bone_pix = 10 * P.lambda;
bone_thickness = 0.6e-3; % m
bone_curvature = 5; % 1/m

P.medium_interfaces = ...
    {P.lens_thickness,
 [bone_curvature, 0.02, P.lens_thickness + dist_trans_bone],
 [bone_curvature, 0.035, P.lens_thickness + dist_trans_bone + bone_thickness]};

P.z_pixel = P.lens_thickness + dist_trans_bone + bone_thickness + dist_bone_pix;

P = rt2_derived_parameters(P);

% =============================================================================
% Plot medium
% =============================================================================

figure(1); clf;
rt.plot.medium(P);
yline(P.z_pixel * 1e3)
xline(P.x_pixel * 1e3)
% yline(z_c0*1e3)

%% Compute
data = struct();
fprintf('Computing Max Intensity...\n')

modes = {'AC','NC','LC'};
NModes = numel(modes);
for km  = 1:NModes 
    m = modes{km};
    data.(m) = rt2_compute(P, m);
end

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
for km  = 1:NModes 
    m = modes{km};
    plot(P.x_piezo * 1e3, data.(m).tof_round_trip * 1e6, LOPTS{:},'displayname',m);
end
legend

%% plot medium and rays

rays_tx = data.AC.rays_tx;
rays_rx = data.AC.rays_rx(data.AC.f_number_msk);
f_number_idx_c0 = find(data.NC.f_number_msk);

figure(1); clf;
rt.plot.medium(P);
% yline(P.z_pixel * 1e3)
% xline(P.x_pixel * 1e3)
% yline(z_c0*1e3)

% plot refracted rays
rt.plot.rays(rays_tx);
rt.plot.rays(rays_rx{1});
rt.plot.rays(rays_rx{end});

% plot homogeneous paths
plot([P.x_piezo(f_number_idx_c0(1)) P.x_pixel] * 1e3, [P.z_piezo(f_number_idx_c0(1)) P.z_pixel] * 1e3, 'k--')
plot([P.x_piezo(f_number_idx_c0(end)) P.x_pixel] * 1e3, [P.z_piezo(f_number_idx_c0(end)) P.z_pixel] * 1e3, 'k--')
plot([P.x_source P.x_pixel] * 1e3, [P.z_source P.z_pixel] * 1e3, 'k--')

%% Plot error in time of flight 
ref = data.AC;
cmap = lines(NModes);

figure(31);clf;
subplot(311)
imagesc(P.x_piezo*1e3, rf_data.t_vec*1e6, rf_data.RF)
colormap bone 
hold on ;
for km = 1:NModes 
    m = modes{km};
    error_tof = data.(m).tof_round_trip-ref.tof_round_trip;
    msk = data.(m).f_number_msk;
    plot(P.x_piezo*1e3,(error_tof+rf_data.ttp)*1e6, '--', 'DisplayName',m, 'color', cmap(km,:),'HandleVisibility','off'); 
    plot(P.x_piezo(msk)*1e3,(error_tof(msk)+rf_data.ttp)*1e6,'linewidth',2, 'DisplayName',m, 'color', cmap(km,:)); 
    idx = find(msk);
    xline(P.x_piezo(idx([1 end]))*1e3,'color',cmap(km,:),'HandleVisibility','off');
end
ylim([0 rf_data.ttp*5e6])
legend
ylabel('Time (\mus)')
title('Delayed RF with true time of flight')

wave_period = 1/P.Fc; max_err= 0; min_err = 0;
subplot(312); hold on 
for km = 1:NModes 
    m = modes{km};
    error_tof = data.(m).tof_round_trip-ref.tof_round_trip;
    error_tof_wvl = (error_tof ./ wave_period)*100; % in %
    msk = data.(m).f_number_msk;
    plot(P.x_piezo*1e3,error_tof_wvl, '--', 'DisplayName',m, 'color', cmap(km,:),'HandleVisibility','off'); 
    plot(P.x_piezo(msk)*1e3,error_tof_wvl(msk),'linewidth',2, 'DisplayName',m, 'color', cmap(km,:)); 

    max_err = max(max(error_tof_wvl(msk)), max_err);
    min_err = min(min(error_tof_wvl(msk)), min_err);
end
legend
ylim([min_err max_err]*1.5)
xlim(P.x_piezo([1 end])*1e3)
ylabel('Error ToF (% wave period)')
title('Relative error in time of flight')

subplot(313); hold on 
half_opening_angle_deg = atand(1/2 / P.f_number);
for km = 1:NModes 
    m = modes{km};
    theta_rx = rad2deg(abs(data.(m).theta_rx));
    msk = data.(m).f_number_msk;
    plot(P.x_piezo*1e3,theta_rx, '--', 'DisplayName',m, 'color', cmap(km,:),'HandleVisibility','off'); 
    plot(P.x_piezo(msk)*1e3,theta_rx(msk),'linewidth',2, 'DisplayName',m, 'color', cmap(km,:)); 

    idx = find(msk);
    xline(P.x_piezo(idx([1 end]))*1e3,'color',cmap(km,:),'HandleVisibility','off');
end
legend
xlim(P.x_piezo([1 end])*1e3)
ylim([0 half_opening_angle_deg*2])
yline(half_opening_angle_deg,'DisplayName','F number threshold')
ylabel('Reception angle (deg)')
xlabel('Lateral element position (mm)')
title('Reception angle at element')

%% reconstruct PSF
fprintf('Computing PSFs...\n')
clear data_psf
for k  = 1:NModes 
    m = modes{km};
    data_psf.(m) = rt2_compute_psf(P, rf_data_full, m);
end
fprintf('Computing PSFs Done!\n')

%% Plot PSFs

modes = fieldnames(data_psf);
NModes = numel(modes);

figure(21); clf; hold on
for km = 1:NModes
    m = modes{km};
    plot(data_psf.(m).xv * 1e3, abs(data_psf.(m).IQ_line), 'DisplayName', m)
end
legend

getter_fun = @(field) cellfun(@(m)data_psf.(m).(field), fieldnames(data_psf));
res_x = getter_fun('res_x');
res_z = getter_fun('res_z');
figure(22); clf; hold on
subplot(121)
bar(res_x * 1e3); xticklabels(modes)
ylabel('lateral resolution (mm)')
subplot(122)
bar(res_z * 1e3); xticklabels(modes)
ylabel('axial resolution (mm)')

%%

getter_fun = @(field) cellfun(@(m)data_psf.(m).(field), fieldnames(data_psf), 'UniformOutput', false);
IQ_all = getter_fun('IQ_grid');
IQ_all = cat(3, IQ_all{:});
BMode = flogc(IQ_all);

DR = [-40 0];
figure(23); clf;
for km = 1:NModes
    subplot(1, NModes, km)
    m = modes{km};
    %imagesc(data_psf.(m).xv*1e3,data_psf.(m).zv*1e3, abs(data_psf.(m).IQ_grid))
    imagesc(data_psf.(m).xv * 1e3, data_psf.(m).zv * 1e3, BMode(:, :, km));
    hold on 
    scatter(P.x_pixel*1e3, P.z_pixel*1e3,'rx')
    caxis(DR);
    colorbar
    title(m)
    colormap bone
    daspect([1 1 1])
end

