function result = rt2_simulate_multi_source(P)

    P = rt2_derived_parameters(P);

    modes = {'LC', 'AC'};

    NModes = numel(modes);
    %% Compute

    nsources = numel(P.theta_source);

    for ks = 1:nsources
        Ps = P;
        Ps.theta_source = P.theta_source(ks);
        Ps.x_source = P.x_source(ks);
        Ps.z_source = P.z_source(ks);

        [~, rf_data_full] = rt2_generate(Ps);

        % reconstruct PSF
        fprintf('Computing PSFs...\n')
        for km = 1:NModes
            m = modes{km};
            data_psf(ks).(m) = rt2_compute_psf(Ps, rf_data_full, m);
        end

        fprintf('Computing PSFs Done!\n')

    end

    % compound
    for km = 1:NModes
        m = modes{km};
        xv = data_psf(1).(m).xv;
        zv = data_psf(1).(m).zv;

        Nx = numel(xv); Nz = numel(zv);
        IQ_grid = zeros(Nz, Nx, nsources);

        [X, Z] = meshgrid(xv, zv);

        for ks = 1:nsources
            xvl = data_psf(ks).(m).xv;
            zvl = data_psf(ks).(m).zv;
            [Xl, Zl] = meshgrid(xvl, zvl);
            IQl = data_psf(ks).(m).IQ_grid;
            IQint = interp2(Xl, Zl, IQl, X, Z, 'linear', 0);
            IQ_grid(:, :, ks) = IQint;
        end

        IQ_grid = sum(IQ_grid, 3, 'omitnan');
        [peak, kz, kx] = maxij(abs(IQ_grid));

        data_psf_c.(m).IQ_grid = IQ_grid;
        data_psf_c.(m).xv = xv;
        data_psf_c.(m).zv = zv;
        data_psf_c.(m).IQ_line = IQ_grid(kz, :);
        data_psf_c.(m).IQ_line_z = IQ_grid(:, kx);
        data_psf_c.(m).kz = kz;
        data_psf_c.(m).kz = kx;
        data_psf_c.(m).peak = abs(peak);

        % determine fwhm
        data_psf_c.(m).res_x = fwhm(abs(data_psf_c.(m).IQ_line), data_psf_c.(m).xv);
        data_psf_c.(m).res_z = fwhm(abs(data_psf_c.(m).IQ_line_z), data_psf_c.(m).zv);
    end

    %% pack
    result.data_psf = data_psf_c;
    result.data_psf_all = data_psf;

end
