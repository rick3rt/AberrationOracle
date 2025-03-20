function varargout = plot_vec(s, e, varargin)
    % plot_vec = @(s, e, varargin) plot([s(1) e(1)], [s(2) e(2)], varargin{:});
    h = plot([s(1) e(1)], [s(2) e(2)], varargin{:});

    if nargout > 0
        varargout{1} = h;
    end

end
