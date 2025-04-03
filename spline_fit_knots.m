
function [sp, xknots ] = spline_fit_knots(xdata, ydata, nknots)

    if numel(nknots) == 1; assert(nknots > 1,'Number of knots must be > 1'); end
    
    % setup knots
    if numel(nknots) == 1
        xknots = linspace(xdata(1), xdata(end), nknots); % The vector of knots that you want it to be
    else
        xknots = nknots;
        nknots = numel(xknots);
    end
        
    % Find out y vector at xknots such that the spline will fit the data
    ndata = length(xdata);
    M = zeros(ndata,nknots);
    for k=1:nknots
        byk = accumarray(k,1,[nknots,1]);
        M(:,k) = spline(xknots,byk,xdata);
    end
    y = M\ydata(:);
    % Check the spline model
    sp = spline(xknots, y);
    
end