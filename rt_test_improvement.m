function [metrics, data, figs] = rt_test_improvement(P, BFC,out, test_resolution, apply_attenuation)

    if ~exist('test_resolution', 'var'); test_resolution = false; end
    if ~exist('apply_attenuation', 'var'); apply_attenuation = true; end
    figs = {}; 
    
    %% Element dictivity
    
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
    % % plot(P.x_piezo * 1e3, attenuation_rx, '-', 'LineWidth', 2)
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

    % time vector for delaying
    t_min = 0;
    t_max = max(error_rt_tof) * 20; % delay with error
    t_max = 20 / P.Fc;
    t_vec = t_min:1 / P.Fs:t_max;
    N = numel(t_vec);

    RF_pulse = repmat(image_pulse, P.num_elements, 1).';
    RF_pulse(N, 1) = 0; % pad with zeros

    % apply 'real' attenuation etc.
    if apply_attenuation
        RF_pulse = RF_pulse .* attenuation_rt_mag; % attenuation
        RF_pulse = RF_pulse .* TF_rt; % transmission factors
        RF_pulse = RF_pulse .* directivity_fun(out.ac_theta_rx); % element directivity
    end
    % delay RF
    RF_delayed_nc = pulse_delaying_RF(RF_pulse, error_rt_tof, P.Fs);

    % plot RF error 
    cmap_lines = lines(2);
    figure(102); clf;
    subplot(121)
    imagesc(P.x_piezo * 1e3, t_vec * 1e6, RF_delayed_nc)
    colormap bone
    xline(P.x_piezo(out.ac_f_number_idx([1 end])) * 1e3, 'color', cmap_lines(1, :));
    xline(P.x_piezo(out.nc_f_number_idx([1 end])) * 1e3, 'color', cmap_lines(2, :));

    RF_pulse_sum_ac = sum(RF_pulse(:, out.ac_f_number_idx), 2);
    RF_pulse_sum_nc = sum(RF_delayed_nc(:, out.nc_f_number_idx), 2);
    % RF_pulse_sum_nc = sum(RF_delayed_nc(:, out.ac_f_number_idx), 2); % apply AC aperture to NC

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

    if nargout < 2 && ~test_resolution
        return
    end

    RF_full = RF_pulse;
    t_max_full = max(out.ac_tof_round_trip(:)) * 1.5; % full RF
    t_vec_full = t_min:1 / P.Fs:t_max_full;
    RF_full(numel(t_vec_full), 1) = 0; % pad with zeros

    RF_full = RF_full .* attenuation_rt_mag; % attenuation
    RF_full = RF_full .* TF_rt; % transmission factors
    RF_full = RF_full .* directivity_fun(out.ac_theta_rx); % element directivity
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

    % Aberration Corrected lateral resolution
    span_lambda = 5;
    xvec = P.x_pixel + (-span_lambda * P.lambda:P.lambda / 4:span_lambda * P.lambda);
    zvec = P.z_pixel + (-span_lambda * P.lambda:P.lambda / 4:span_lambda * P.lambda);
    Nx = numel(xvec); Nz = numel(zvec);

    IQ_full = hilbert(data.RF);% analytic signal

    IQ_xax_ac = zeros(1, Nx);
    IQ_xax_nc = zeros(1, Nx);
    tof_ac = zeros(Nx, P.num_elements);

    time_remaining_progbar_ui(0, Nx)
    for kx = 1:Nx
        % if 1; kx =  ceil(Nx/2);
        tof_ac(kx,:) = rt_trace_tof(P, BFC, xvec(kx), P.z_pixel);
        IQ_interp = interp1_per_channel(data.tvec, IQ_full, tof_ac(kx,:),'cubic');
        IQ_xax_ac(kx) = sum(IQ_interp, 'omitmissing');
        time_remaining_progbar_ui(kx, Nx)
    end

    % Homogeneous
    span_lambda_x = 10;
    span_lambda_z = 20;
    xvec_nc = P.x_pixel + (-span_lambda_x * P.lambda:P.lambda / 4:span_lambda_x * P.lambda);
    zvec = P.z_pixel + (-span_lambda_z * P.lambda:P.lambda / 4:span_lambda_z * P.lambda);
    Nx = numel(xvec_nc); Nz = numel(zvec);
    [X, Z] = meshgrid(xvec_nc, zvec);

    tau_tx = hypot(X - P.x_source, Z - P.z_source) / P.c0 - P.tx_add_to_nc;
    tau_rx = hypot(X - reshape(P.x_piezo, 1, 1, []), Z - reshape(P.z_piezo, 1, 1, [])) / P.c0;
    theta_rx = atan2(Z - reshape(P.z_piezo, 1, 1, []), X - reshape(P.x_piezo, 1, 1, [])) - pi / 2;
    half_opening_angle_rad = atan(1/2 / P.f_number);

    IQ_nc = zeros(Nz, Nx);
    for kx = 1:Nx
        tau = squeeze(tau_tx(:, kx, :) + tau_rx(:, kx, :));
        f_mask = squeeze(abs(theta_rx(:, kx, :)) < half_opening_angle_rad);
        IQ_interp = interp1_per_channel(data.tvec, IQ_full, tau, 'cubic');
        IQ_line = sum(IQ_interp .* f_mask, 2,'omitmissing');
        IQ_nc(:, kx) = IQ_line;
    end

    figure(990);clf;
    imagesc(P.x_piezo,data.tvec,real(IQ_full))
    hold on 
    kx = ceil(Nx/2); kz = ceil(Nz/2);
    tau = squeeze(tau_tx(kz, kx, :) + tau_rx(kz, kx, :));
    plot(P.x_piezo,tau)



    [v, kz, kx] = maxij(abs(IQ_nc));
    % debug plot
    f = figure(123); clf;
    f.Name = ['psf_homo_' P.name];
    imagesc(xvec_nc * 1e3, zvec * 1e3, flogc(IQ_nc), [-60 0])
    daspect([1 1 1])
    hold on
    scatter(xvec_nc(kx) * 1e3, zvec(kz) * 1e3)
    colormap bone
    title(['psf homo ' P.name])

    figs{end+1} = f;

    % interpolate to same grid as AC
    IQ_xax_nc = IQ_nc(kz, :);
    IQ_xax_nc = interp1(xvec_nc, IQ_xax_nc, xvec);



    %% determine and show resolution.
    [w_ac, y50_ac, x1_ac, x2_ac] = fwhm2(xvec * 1e3, abs(IQ_xax_ac));
    [w_nc, y50_nc, x1_nc, x2_nc] = fwhm2(xvec * 1e3, abs(IQ_xax_nc));

    cmap = lines(2);

    f = figure(69); clf
    f.Name = ['resolution_' P.name];
    figs{end+1} = f;
    plot(xvec * 1e3, abs(IQ_xax_nc))
    hold on
    plot(xvec * 1e3, abs(IQ_xax_ac))
    legend('NC', 'AC')
    plot([x1_nc x2_nc], [y50_nc y50_nc], 'Color', cmap(1, :), 'HandleVisibility', 'off')
    plot([x1_ac x2_ac], [y50_ac y50_ac], 'Color', cmap(2, :), 'HandleVisibility', 'off')
    title(sprintf('%s\\newline Lateral Resolution - NC: %.3f mm - AC: %.3f mm (%.2fx)',P.name, w_nc, w_ac, w_ac/w_nc))
    
    % resolution metrics
    metrics.ac_res_x = w_ac; % in mm
    metrics.nc_res_x = w_nc; % in mm
    metrics.ac_IQ_x = IQ_xax_ac;
    metrics.nc_IQ_x = IQ_xax_nc;
    metrics.ac_IQ_x_peak = max(abs(IQ_xax_ac));
    metrics.nc_IQ_x_peak = max(abs(IQ_xax_nc));

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
