function medium(BFC, x_start, z_start, x_target, z_target, x_end, z_end)

    plot(BFC.XPiezo * 1e3, BFC.ZPiezo * 1e3, 'k', 'linewidth', 2);
    plotbox(BFC.XRecon * 1e3, BFC.ZRecon * 1e3, 'k');
    cellfun(@(p) plot(BFC.XRecon * 1e3, rt.util.segeval(p, BFC.XRecon) * 1e3, 'k-'), BFC.medium_interfaces);

    if nargin > 1
        scatter(x_start * 1e3, z_start * 1e3, 'ro', 'filled', 'linewidth', 2);
        scatter(x_target * 1e3, z_target * 1e3, 'ro', 'filled', 'linewidth', 2);
        scatter(x_end * 1e3, z_end * 1e3, 'ro', 'filled', 'linewidth', 2);
    end

end
