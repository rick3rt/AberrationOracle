function data = rt2_compute_psf(P, rfdata, mode)

    fprintf('=============================\n')
    fprintf('COMPUTING PSF FOR MODE: %s\n', mode);
    fprintf('=============================\n')

    switch mode
        case {'AC', 'LC'} % no correction, homogeneous medium

            if strcmp(mode, 'LC')
                tmp = rt2_compute(P, 'AC');
                tau = tmp.tof_round_trip;
                taum = min(tau);
                tau2 = taum - 2 * P.lens_thickness / P.medium_soundspeeds(1);
                z2 = tau2 / 2 * P.medium_soundspeeds(2);
                P.z_pixel = z2 + P.lens_thickness;
            end

            % Aberration Corrected lateral resolution
            % span_lambda = 6;
            % xv = P.x_pixel + (-span_lambda * P.lambda:P.lambda / 2:span_lambda * P.lambda);
            % zv = P.z_pixel; %  + (-span_lambda * P.lambda:P.lambda / 4:span_lambda * P.lambda);
            % Nx = numel(xv); Nz = numel(zv);
            % [X, Z] = meshgrid(xv, zv);

            % data = rt2_compute_tof_grid(P, mode, Z, X);
            data = rt2_compute_tof_grid_optim(P, mode, rfdata);

            IQ = hilbert(rfdata.RF); % analytic signal
            tvec = rfdata.t_vec;

            xv = data.xv;
            zv = data.zv;
            Nx = numel(xv); Nz = numel(zv);

            IQ_grid = zeros(Nz, Nx);

            for kx = 1:Nx
                tau = reshape(data.tof_round_trip(:, kx, :), Nz, P.num_elements);
                IQ_interp = rt.util.interp1_per_channel(tvec, IQ, tau, 'cubic');
                IQ_grid(:, kx) = sum(IQ_interp, 2, 'omitmissing');
            end

        case 'NC'

            span_lambda_x = 15;
            span_lambda_z = 50;
            xv = P.x_pixel + (-span_lambda_x * P.lambda:P.lambda / 4:span_lambda_x * P.lambda);
            zv = P.z_pixel + (-span_lambda_z * P.lambda:P.lambda / 4:span_lambda_z * P.lambda);
            Nx = numel(xv); Nz = numel(zv);
            [X, Z] = meshgrid(xv, zv);

            data = rt2_compute_tof_grid(P, mode, Z, X);
            IQ = hilbert(rfdata.RF); % analytic signal
            tvec = rfdata.t_vec;

            IQ_grid = zeros(size(X));

            for kx = 1:Nx
                tau = squeeze(data.tof_round_trip(:, kx, :));
                IQ_interp = rt.util.interp1_per_channel(tvec, IQ, tau, 'cubic');
                IQ_grid(:, kx) = sum(IQ_interp, 2, 'omitmissing');
            end
    end

    [peak, kz, kx] = maxij(abs(IQ_grid));

    span_x = 10 * 4; %
    span_z = 10 * 4; % assuming pixel size lambda/4
    ind_x = kx + (-span_x:span_x);
    ind_z = kz + (-span_z:span_z);

    data.IQ_grid = IQ_grid(ind_z, ind_x);
    data.xv = xv(ind_x);
    data.zv = zv(ind_z);
    data.IQ_line = IQ_grid(kz, ind_x);
    data.IQ_line_z = IQ_grid(ind_z, kx);
    data.kz = span_z + 1;
    data.kz = span_x + 1;
    data.peak = abs(peak);

    % determine fwhm
    data.res_x = fwhm(abs(data.IQ_line), data.xv);
    data.res_z = fwhm(abs(data.IQ_line_z), data.zv);

end

function data = rt2_compute_tof_grid_optim(P, mode, rfdata)

    switch mode
        case 'LC' % lens corrected
            recon_to = 1;
            tissue_c = P.c0;
        case 'AC' % full aberration corrected
            recon_to = numel(P.medium_soundspeeds)-1;
            tissue_c = P.medium_soundspeeds(2);
    end

    v_source = [P.x_source; P.z_source];
    v_piezo = [P.x_piezo; P.z_piezo];
    dist_se = vecnorm(v_piezo - v_source);
    tx_delay_lens = dist_se / P.medium_soundspeeds(1); % lens speed of sound
    delay_in_lens = min(tx_delay_lens);

    Pars = struct();
    span = 20 * P.lambda;
    Pars.xv_recon = P.x_pixel + (-span:P.lambda / 4:span);
    Pars.zv_recon = P.z_recon(1):P.lambda / 4:P.z_recon(end);

    Pars.recon_to = recon_to;

    p = P.medium_interfaces{2};
    if ~isstruct(p)
        zp = rt.util.segeval(p, Pars.xv_recon);
        idx_min = find(Pars.zv_recon >= min(zp), 1, 'first');
        idx_max = find(Pars.zv_recon <= max(zp), 1, 'last');
        PeriParabIn = [p idx_min - 1 idx_max - 1];
        PeriPolyOrder = numel(p) - 1;
    else
        PeriParabIn = p;
        PeriPolyOrder = p.pieces;
    end

    % if recon_to == 3
    p = P.medium_interfaces{3};
    if ~isstruct(p)
        zp = rt.util.segeval(p, Pars.xv_recon);
        idx_min = find(Pars.zv_recon >= min(zp), 1, 'first');
        idx_max = find(Pars.zv_recon <= max(zp), 1, 'last');
        EndoParabIn = [p idx_min - 1 idx_max - 1];
        EndoPolyOrder = numel(p) - 1;
    else
        EndoParabIn = p;
        EndoPolyOrder = p.pieces;
    end
    %

    Pars.shape_factor_aniso = 1;
    Pars.lens_thickness = P.lens_thickness;
    Pars.pitch = P.lambda;
    Pars.max_depth_bone = max(P.z_recon);
    Pars.min_thickness_bone = 1e-6;
    Pars.max_thickness_bone = 10e-3;
    Pars.lens_half_opening_angle = atan(1/2 / P.f_number);
    Pars.min_Fnumber = P.f_number;
    Pars.max_error_receive_angle = deg2rad(2.5);
    Pars.Fs = P.Fs;
    Pars.lens_c = P.medium_soundspeeds(1);
    Pars.tissue_c = tissue_c;
    if recon_to == 3
        Pars.brain_c = P.medium_soundspeeds(4);
    else
        Pars.brain_c = P.medium_soundspeeds(3);
    end
    Pars.axial_c = P.medium_soundspeeds(3);
    Pars.radial_c = P.medium_soundspeeds(3);
    Pars.n_cores = 12;
    Pars.NeedTravelTime = 1;
    Pars.SubApertureApodis = 1;
    Pars.PeriPolyOrder = PeriPolyOrder;
    Pars.EndoPolyOrder = EndoPolyOrder;
  
    Pars.X_El = P.x_piezo;
    Pars.Z_El = P.z_piezo;
    Pars.PeriParabIn = PeriParabIn;
    Pars.EndoParabIn = EndoParabIn;
    Pars.XS = P.x_source;
    Pars.ZS = P.z_source;
    Pars.tx_to_use = 1:numel(P.x_source);
    Pars.add_to_delay_firing = delay_in_lens;
    Pars.ReceiveAngleRad = 0;

    if Pars.recon_to == 1
        Pars.tissue_c = P.c0;
    end

    [data3l, ~] = rt.recon.three_layer_recon_wrapper(Pars, rfdata.RF, rfdata.RF);

    perm = @(x) permute(x, [2 3 1]);
    if Pars.recon_to == 1
        tof_rt = perm (data3l.Time_T_Tissue) + perm(data3l.Time_R_Tissue) - delay_in_lens;
        theta_rx = perm(data3l.Angle_R_Tissue);
    elseif Pars.recon_to == 2
        tof_rt = perm (data3l.Time_T_Bone) + perm(data3l.Time_R_Bone) - delay_in_lens;
        theta_rx = perm(data3l.Angle_R_Bone);
    elseif Pars.recon_to == 3
        tof_rt = perm (data3l.Time_T_Marrow) + perm(data3l.Time_R_Marrow) - delay_in_lens;
        theta_rx = perm(data3l.Angle_R_Marrow);
    end

    msk1 = tof_rt < 0;
    tof_rt(msk1) = NaN;
    theta_rx(msk1) = NaN;
    theta_rx(abs(theta_rx) > Pars.lens_half_opening_angle) = NaN;
    f_number_mask = ~isnan(theta_rx);
    tof_rt(~f_number_mask) = NaN;

    data.tof_round_trip = tof_rt;
    data.zv = Pars.zv_recon;
    data.xv = Pars.xv_recon;
    data.f_number_mask = f_number_mask;
    data.theta_rx = theta_rx;

end
