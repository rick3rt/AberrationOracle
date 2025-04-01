function medium(BFC, x_start, z_start, x_target, z_target, x_end, z_end)

    % create map with wavespeeds per layer
    [X, Z] = meshgrid(BFC.XRecon, BFC.ZRecon);
    sos_map = ones(size(X)) * BFC.medium_soundspeeds(1);

    for k = 1:numel(BFC.medium_soundspeeds) - 1
        mask = Z > rt.util.segeval(BFC.medium_interfaces{k}, X);
        sos_map(mask) = BFC.medium_soundspeeds(k + 1);
    end

    [~, ~, is] = unique(sos_map(:));
    sos_mask = reshape(is, size(sos_map));

    imagesc(BFC.XRecon * 1e3, BFC.ZRecon * 1e3, sos_mask);
    colormap(my_colormap())
    hold on
    plot(BFC.XPiezo * 1e3, BFC.ZPiezo * 1e3, 'k', 'linewidth', 2);
    rt.plot.plotbox(BFC.XRecon * 1e3, BFC.ZRecon * 1e3, 'k');
    cellfun(@(p) plot(BFC.XRecon * 1e3, rt.util.segeval(p, BFC.XRecon) * 1e3, 'k-'), BFC.medium_interfaces);
    daspect([1 1 1])

    if nargin > 1
        scatter(x_start * 1e3, z_start * 1e3, 'ro', 'filled', 'linewidth', 2);
        scatter(x_target * 1e3, z_target * 1e3, 'ro', 'filled', 'linewidth', 2);
        scatter(x_end * 1e3, z_end * 1e3, 'ro', 'filled', 'linewidth', 2);
    end

    % figure(); imagesc(sos_mask)

end

function cmap = my_colormap()
    cmap = [251, 180, 174;
            179, 205, 227;
            204, 235, 197;
            222, 203, 228] ./ 255;
end
