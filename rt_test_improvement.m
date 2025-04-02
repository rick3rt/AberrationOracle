function [metrics, data] = rt_test_improvement(P, out, BFC, apply_attenuation)

if ~exist('apply_attenuation','var'); apply_attenuation = true; end
fprintf('Applying Att: %i\n',apply_attenuation)

%% Element dictivity

lambda = P.lambda;
W = lambda*0.9; % small kerf
directivity_fun = @(theta) cos(theta).*sinc(pi*W/lambda.*sin(theta));

% theta_test = linspace(-pi/2, pi/2, 128);
% figure(10);clf; 
% plot(rad2deg(theta_test), directivity_fun(theta_test))

%% Total attenuation along path
BFC.medium_attenuation = [1.0, 0.54, 6.9, 0.6]; % dB/(MHz*cm)
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


%% Generate RF
P.Fs = 20 * P.Fc; % Sampling frequency
P.num_cycles = 3;

% Set the impulse response and excitation of the emit aperture
image_pulse = sin(2 * pi * P.Fc * (0:1 / P.Fs:P.num_cycles / P.Fc)); % BP66 %,  sin(2*pi*fc*(0:1/fs:1/fc)); % BP100 %
image_pulse = image_pulse .* hamming(length(image_pulse))';
image_pulse_env = envelope(image_pulse);

[peak_pulse, imax] = max(image_pulse_env);
ttp = imax / P.Fs;

% figure(99);
% plot(image_pulse)
% hold on
% plot(image_pulse_env)


error_rt_tof = out.error_tof_round_trip;
error_rt_tof(error_rt_tof==0) = NaN;
error_rt_tof = fillmissing(error_rt_tof, 'spline');
error_rt_tof = error_rt_tof-min(error_rt_tof);



t_min = 0;
t_max = max(error_rt_tof) * 20; % delay with error
t_max = 20/P.Fc;

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
    RF_pulse = RF_pulse .* directivity_fun( out.ac_theta_rx); % element directivity
end
RF_delayed_nc = pulse_delaying_RF(RF_pulse, error_rt_tof, P.Fs);


cmap_lines = lines(2);
figure(102); clf;
subplot(121)
imagesc(P.x_piezo * 1e3, t_vec * 1e6, RF_delayed_nc)
colormap bone
xline(P.x_piezo(out.ac_f_number_idx([1 end])) * 1e3,'color',cmap_lines(1,:));
xline(P.x_piezo(out.nc_f_number_idx([1 end])) * 1e3,'color',cmap_lines(2,:));

RF_pulse_sum_ac = sum(RF_pulse(:, out.ac_f_number_idx), 2);
RF_pulse_sum_nc = sum(RF_delayed_nc(:, out.nc_f_number_idx), 2);

env_sum_ac = envelope(RF_pulse_sum_ac);
env_sum_nc = envelope(RF_pulse_sum_nc);
peak_ac = max(env_sum_ac);
peak_nc = max(env_sum_nc);
fprintf('Intensity improvement f-masked:     %.2fx\n', peak_ac/peak_nc);

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
metrics.improvement = peak_ac/peak_nc;
metrics.ac_RF_sum = RF_pulse_sum_ac;
metrics.nc_RF_sum = RF_pulse_sum_nc;
% metrics.image_pulse = image_pulse;


%% OPTIONAL SHOW FULL RF

    
    if nargout == 2
    
        RF_full = RF_pulse;
        t_max_full = max(out.ac_tof_round_trip(:)) * 1.5; % full RF
        t_vec_full = t_min:1 / P.Fs:t_max_full;
        RF_full(numel(t_vec_full),1) = 0; % pad with zeros
    
        RF_full = RF_full .* attenuation_rt_mag; % attenuation
        RF_full = RF_full .* TF_rt; % transmission factors
        RF_full = RF_full .* directivity_fun( out.ac_theta_rx); % element directivity
        RF_delayed_full = pulse_delaying_RF(RF_full, out.ac_tof_round_trip-ttp, P.Fs);
        cmap_lines = lines(4);
        cmap_lines = flip(cmap_lines);
        
        figure(101); clf;
        imagesc(P.x_piezo * 1e3, t_vec_full * 1e6, RF_delayed_full)
        colormap bone
        hold on
    
        plot(P.x_piezo * 1e3, out.ac_tof_round_trip * 1e6, '--','color',cmap_lines(1,:), 'LineWidth', 1);
        plot(P.x_piezo(out.ac_f_number_idx) * 1e3, out.ac_tof_round_trip(out.ac_f_number_idx) * 1e6, '-','color',cmap_lines(1,:), 'LineWidth', 2);
        plot(P.x_piezo * 1e3, out.nc_tof_round_trip * 1e6,'--','color',cmap_lines(2,:), 'LineWidth', 1);
        plot(P.x_piezo(out.nc_f_number_idx) * 1e3, out.nc_tof_round_trip(out.nc_f_number_idx) * 1e6,'-','color',cmap_lines(2,:), 'LineWidth', 2);
    
    
    
        data.RF = RF_delayed_full;
        data.tvec = t_vec_full;
        data.ttp = ttp;
    end


end

% RF_pulse_sum_ac = sum(RF_pulse, 2);
% RF_pulse_sum_nc = sum(RF_delayed_nc, 2);
% env_sum_ac = envelope(RF_pulse_sum_ac);
% env_sum_nc = envelope(RF_pulse_sum_nc);
% peak_ac = max(env_sum_ac);
% peak_nc = max(env_sum_nc);
% fprintf('Intensity improvement all elements ( no f number applied): %.2fx\n', peak_ac/peak_nc);

% subplot(211)
% plot(t_vec * 1e6, RF_pulse_sum_ac)
% hold on
% plot(t_vec * 1e6, RF_pulse_sum_nc)
% title('Sum all elements')
