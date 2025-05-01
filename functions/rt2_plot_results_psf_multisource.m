function figs = rt2_plot_results_psf_multisource(P, result)
    
    unpackStruct(result);
    modes = fieldnames(data_psf);
    refIdx = find(strcmp(modes,'AC'));

    figs = {}; 
    NModes = numel(modes);

    f = figure(15); clf; hold on
    f.Name = 'lateral_res_PSF_multisource';
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
    f.Name = 'lateral_res_PSF_metrics_multisource';
    f.Position = [100 500 600 400];
    figs{end + 1} = f;
    subplot(121)
    bar(res_x * 1e3); xticklabels(modes)
    ylabel('lateral resolution (mm)')
    subplot(122)
    bar(res_z * 1e3); xticklabels(modes)
    ylabel('axial resolution (mm)')

    
    imp = (res_x-res_x(refIdx))./res_x;
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
    f.Name = 'full_psf_multisource';
    f.Position = [800 500 1000 400];
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