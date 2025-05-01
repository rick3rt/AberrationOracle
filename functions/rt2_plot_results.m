function figs = rt2_plot_results(P, result)
    %%
    unpackStruct(result);

    modes = fieldnames(metrics);
    NModes = numel(modes);
    refIdx = find(strcmp(modes, 'AC'));

    figs = {};

    % ===========================================================
    % getter_fun = @(field) cellfun(@(m)metrics.(m).(field), fieldnames(metrics));
    % f = figure(11); clf;
    % figs{end + 1} = f;
    % bar(getter_fun('peak'))
    % xticklabels(modes)

    rf_plot = rf_data_full;
    % rf_plot = rf_data;
    LOPTS = {'linewidth', 1.5};

    f = figure(12); clf;
    f.Name = 'rf_full';
    figs{end + 1} = f;
    imagesc(P.x_piezo * 1e3, rf_plot.t_vec * 1e6, rf_plot.RF)
    ylabel('Time of flight (\mus)')
    xlabel('Element position (mm)')
    colormap bone
    hold on
    for km = 1:NModes
        m = modes{km};
        plot(P.x_piezo * 1e3, data.(m).tof_round_trip * 1e6, LOPTS{:}, 'displayname', m);
    end
    legend

    %% plot medium and rays

    rays_tx = data.AC.rays_tx;
    rays_rx = data.AC.rays_rx(data.AC.f_number_msk);
    f_number_idx_c0 = find(data.LC.f_number_msk);

    f = figure(13); clf;
    f.Name = 'medium';
    figs{end + 1} = f;
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

    f = figure(14); clf;
    f.Name = 'tof_error';
    f.Position = [713 180 511 714];
    figs{end + 1} = f;
    subplot(311)
    imagesc(P.x_piezo * 1e3, rf_data.t_vec * 1e6, rf_data.RF)
    colormap bone
    hold on;
    for km = 1:NModes
        m = modes{km};
        error_tof = data.(m).tof_round_trip - ref.tof_round_trip;
        msk = data.(m).f_number_msk;
        plot(P.x_piezo * 1e3, (error_tof + rf_data.ttp) * 1e6, '--', 'DisplayName', m, 'color', cmap(km, :), 'HandleVisibility', 'off');
        plot(P.x_piezo(msk) * 1e3, (error_tof(msk) + rf_data.ttp) * 1e6, 'linewidth', 2, 'DisplayName', m, 'color', cmap(km, :));
        idx = find(msk);
        xline(P.x_piezo(idx([1 end])) * 1e3, 'color', cmap(km, :), 'HandleVisibility', 'off');
    end
    ylim([0 rf_data.ttp * 5e6])
    legend
    ylabel('Time (\mus)')
    title('Delayed RF with true time of flight')

    wave_period = 1 / P.Fc; max_err = 0; min_err = 0;
    subplot(312); hold on
    for km = 1:NModes
        m = modes{km};
        error_tof = data.(m).tof_round_trip - ref.tof_round_trip;
        error_tof_wvl = (error_tof ./ wave_period) * 100; % in %
        msk = data.(m).f_number_msk;
        plot(P.x_piezo * 1e3, error_tof_wvl, '--', 'DisplayName', m, 'color', cmap(km, :), 'HandleVisibility', 'off');
        plot(P.x_piezo(msk) * 1e3, error_tof_wvl(msk), 'linewidth', 2, 'DisplayName', m, 'color', cmap(km, :));

        max_err = max(max(error_tof_wvl(msk)), max_err);
        min_err = min(min(error_tof_wvl(msk)), min_err);
    end
    legend
    ylim([min_err max_err] * 1.5)
    xlim(P.x_piezo([1 end]) * 1e3)
    ylabel('Error ToF (% wave period)')
    title('Relative error in time of flight')

    subplot(313); hold on
    half_opening_angle_deg = atand(1/2 / P.f_number);
    for km = 1:NModes
        m = modes{km};
        theta_rx = rad2deg(abs(data.(m).theta_rx));
        msk = data.(m).f_number_msk;
        plot(P.x_piezo * 1e3, theta_rx, '--', 'DisplayName', m, 'color', cmap(km, :), 'HandleVisibility', 'off');
        plot(P.x_piezo(msk) * 1e3, theta_rx(msk), 'linewidth', 2, 'DisplayName', m, 'color', cmap(km, :));

        idx = find(msk);
        xline(P.x_piezo(idx([1 end])) * 1e3, 'color', cmap(km, :), 'HandleVisibility', 'off');
    end
    legend
    xlim(P.x_piezo([1 end]) * 1e3)
    ylim([0 half_opening_angle_deg * 2])
    yline(half_opening_angle_deg, 'DisplayName', 'F number threshold')
    ylabel('Reception angle (deg)')
    xlabel('Lateral element position (mm)')
    title('Reception angle at element')

    %% Plot PSFs

    f = figure(15); clf; hold on
    f.Name = 'lateral_res_PSF';
    f.Position = [100 100 600 400];
    figs{end + 1} = f;
    for km = 1:NModes
        m = modes{km};
        plot(data_psf.(m).xv * 1e3, abs(data_psf.(m).IQ_line), 'DisplayName', m)
    end
    legend

    getter_fun = @(field) cellfun(@(m)data_psf.(m).(field), fieldnames(data_psf));
    res_x = getter_fun('res_x');
    res_z = getter_fun('res_z');

    f = figure(16); clf; hold on
    f.Name = 'lateral_res_PSF_metrics';
    f.Position = [100 100 600 400];
    figs{end + 1} = f;
    subplot(121)
    bar(res_x * 1e3); xticklabels(modes)
    ylabel('lateral resolution (mm)')
    subplot(122)
    bar(res_z * 1e3); xticklabels(modes)
    ylabel('axial resolution (mm)')
    %
    %
    % imp = (res_x-res_x(refIdx))./res_x;
    % mt = modes{refIdx};
    % for k = 1:NModes
    %     m = modes{k};
    %     fprintf('lateral res %s-%s: %.2f\n', m, mt, imp(k));
    % end

    imp = (res_x - res_x(refIdx)) ./ res_x;
    mt = modes{refIdx};
    lbls = {};
    for k = 1:NModes
        m = modes{k};
        lbls{k} = sprintf('lateral res %s-%s: %.2f\n', m, mt, imp(k));
        disp(lbls{k});
    end

    %%

    getter_fun = @(field) cellfun(@(m)data_psf.(m).(field), fieldnames(data_psf), 'UniformOutput', false);
    IQ_all = getter_fun('IQ_grid');
    IQ_all = cat(3, IQ_all{:});
    BMode = flogc(IQ_all);

    DR = [-40 0];
    f = figure(17); clf;
    f.Name = 'full_psf';
    f.Position = [100 100 1000 400];
    figs{end + 1} = f;
    for km = 1:NModes
        subplot(1, NModes, km)
        m = modes{km};
        %imagesc(data_psf.(m).xv*1e3,data_psf.(m).zv*1e3, abs(data_psf.(m).IQ_grid))
        imagesc(data_psf.(m).xv * 1e3, data_psf.(m).zv * 1e3, BMode(:, :, km));
        hold on
        scatter(P.x_pixel * 1e3, P.z_pixel * 1e3, 'rx')
        caxis(DR);
        colorbar
        title([m ' - ' lbls{km}])
        colormap bone
        daspect([1 1 1])
    end

end
