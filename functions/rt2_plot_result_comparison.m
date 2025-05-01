function figs = rt2_plot_result_comparison(P, P_all, result)

    % modes = {'AC', 'NC'};
    modes = fieldnames(result.data_psf);
    refIdx = find(strcmp(modes, 'AC'));
    
    psf_metrics = [result.data_psf];
    getter_fun = @(field) arrayfun(@(pm) cellfun(@(m) pm.(m).(field), modes), psf_metrics, 'UniformOutput', 0);
    
    res_x = getter_fun('res_x');
    res_x = cat(1, res_x{:});
    
    
    figs = {};
    
    f = figure(1); clf;
    f.Name = 'Resolution_lateral';
    f.Position = [200 100 600 600];
    figs{end + 1} = f;
    subplot(211)
    barh(res_x * 1e3)
    set(gca, 'YDir', 'reverse')
    yticklabels({P_all.name})
    
    subplot(212)
    

    barh(abs((res_x - res_x(refIdx)) ./ res_x) * 100)
    set(gca, 'YDir', 'reverse')
    yticklabels({P_all.name})
    
    peak_val = getter_fun('peak');
    peak_val = cat(1, peak_val{:});
    
    f = figure(2); clf;
    f.Name = 'Peak_pixel';
    f.Position = [800 100 600 600];
    figs{end + 1} = f;
    subplot(211)
    barh(peak_val)
    set(gca, 'YDir', 'reverse')
    yticklabels({P_all.name})
    
    subplot(212)
    

    barh(abs((peak_val - peak_val(refIdx)) ./ peak_val) * 100)
    set(gca, 'YDir', 'reverse')
    yticklabels({P_all.name})
end