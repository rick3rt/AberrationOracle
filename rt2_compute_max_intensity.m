function metrics = rt2_compute_max_intensity(P, data, rf_data)


    ref = data.AC;
    tof_ref = ref.tof_round_trip;
    tof_ref = tof_ref - min(tof_ref); % make apex 0
    
    modes = fieldnames(data);
    NModes = numel(modes);

    for k = 1:NModes
        m = modes{k};
        tof_mode = data.(m).tof_round_trip - min(data.(m).tof_round_trip); % make apex 0
        
        tof_error = tof_mode - tof_ref;
        
        RF_delayed = rt.util.pulse_delaying_RF(rf_data.RF, tof_error, P.Fs);
        RF_pulse_sum = sum(RF_delayed(:, data.(m).f_number_msk), 2);

        env_sum = envelope(RF_pulse_sum);
        peak = max(env_sum);

        metrics.(m).peak = peak;
        metrics.(m).RF_pulse_sum = RF_pulse_sum;
    end



end