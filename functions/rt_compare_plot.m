function f_all = rt_compare_plot(P, BFC, out)

    f_all = {};
    cmap_lines = lines(2);

    % rays in speed of sound map
    % ===========================================================

    f = figure(2); clf;
    f.Name = sprintf('rays_%s', P.name);
    f.Position = [100 300 400 600];
    hold on
    rt.plot.medium(BFC)
    set(gca, 'YDir', 'reverse')
    rt.plot.rays(out.ac_rays_tx);
    rt.plot.rays(out.ac_rays_rx{out.ac_f_number_idx(1)});
    rt.plot.rays(out.ac_rays_rx{out.ac_f_number_idx(end)});
    daspect([1 1 1])
    xlim(rt.util.minmax(P.x_recon) * 1e3); ylim(rt.util.minmax(P.z_recon) * 1e3)
    title(['  ' P.name])

    % plot homogeneous paths
    plot([P.x_piezo(out.nc_f_number_idx(1)) P.x_pixel] * 1e3, [P.z_piezo(out.nc_f_number_idx(1)) P.z_pixel] * 1e3, 'k--')
    plot([P.x_piezo(out.nc_f_number_idx(end)) P.x_pixel] * 1e3, [P.z_piezo(out.nc_f_number_idx(end)) P.z_pixel] * 1e3, 'k--')
    plot([P.x_source P.x_pixel] * 1e3, [P.z_source P.z_pixel] * 1e3, 'k--')

    f_all{end + 1} = f;

    % plot error relative to wave period
    % ===========================================================
    f = figure(3); clf;
    f.Name = sprintf('error_%s', P.name);
    f.Position = [500 300 400 600];
    subplot(311)
    plot(P.x_piezo * 1e3, out.ac_tof_round_trip * 1e6, 'linewidth', 2, 'DisplayName', 'Aberration Corrected')
    hold on
    plot(P.x_piezo * 1e3, out.nc_tof_round_trip * 1e6, 'linewidth', 2, 'DisplayName', 'No Correction')
    legend
    ylabel('Travel time \mus')
    xlabel('Lateral Element Position (mm)')
    title(['  ' P.name])

    subplot(312)
    plot(P.x_piezo * 1e3, out.error_tof_round_trip_relative * 100, 'k', 'linewidth', 2)
    ylabel({'Relative Error', 'Wave period'})
    xlabel('Lateral Element Position (mm)')

    xline(P.x_piezo(out.ac_f_number_idx([1 end])) * 1e3, '--', 'Color', cmap_lines(1, :), 'linewidth', 1.5, 'HandleVisibility', 'off');
    xline(P.x_piezo(out.nc_f_number_idx([1 end])) * 1e3, '--', 'Color', cmap_lines(2, :), 'linewidth', 1.5, 'HandleVisibility', 'off');

    subplot(313)
    plot(P.x_piezo * 1e3, rad2deg(abs(out.ac_theta_rx)), 'linewidth', 2, 'DisplayName', 'Aberration Corrected')
    hold on
    plot(P.x_piezo * 1e3, rad2deg(abs(out.nc_theta_rx)), 'linewidth', 2, 'DisplayName', 'No Correction')
    ylabel({'Reception angle', 'Wave period'})
    xlabel('Lateral Element Position (mm)')
    % legend

    yline(rad2deg(out.half_opening_angle_rad), 'linewidth', 1.5, 'HandleVisibility', 'off');
    xline(P.x_piezo(out.ac_f_number_idx([1 end])) * 1e3, '--', 'Color', cmap_lines(1, :), 'linewidth', 1.5, 'HandleVisibility', 'off');
    xline(P.x_piezo(out.nc_f_number_idx([1 end])) * 1e3, '--', 'Color', cmap_lines(2, :), 'linewidth', 1.5, 'HandleVisibility', 'off');
    ylim([0 rad2deg(out.half_opening_angle_rad) * 1.55]);

    f_all{end + 1} = f;

end
