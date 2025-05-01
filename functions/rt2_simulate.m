function result = rt2_simulate(P)

    P = rt2_derived_parameters(P);

    % modes = {'AC', 'NC', 'LC'};
    modes = {'LC', 'AC'};

    %% Compute
    data = struct();
    fprintf('Computing Max Intensity...\n')
    NModes = numel(modes);
    for km = 1:NModes
        m = modes{km};
        data.(m) = rt2_compute(P, m);
    end

    [rf_data, rf_data_full] = rt2_generate(P);
    metrics = rt2_compute_max_intensity(P, data, rf_data);

    %% reconstruct PSF
    fprintf('Computing PSFs...\n')
    clear data_psf
    for km = 1:NModes
        m = modes{km};
        data_psf.(m) = rt2_compute_psf(P, rf_data_full, m);
    end
    fprintf('Computing PSFs Done!\n')

    %% pack
    result.data = data;
    result.metrics = metrics;
    result.rf_data = rf_data;
    result.rf_data_full = rf_data_full;
    result.data_psf = data_psf;

end
