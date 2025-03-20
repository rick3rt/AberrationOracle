function h = plotbox(x, y, varargin)
    xm = [min(x), max(x)];
    ym = [min(y), max(y)];
    h = plot(xm([1 end end 1 1]), ym([1 1 end end 1]), varargin{:});
end
