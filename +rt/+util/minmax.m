function varargout = minmax(x)
    
    minv = min(x(:));
    maxv = max(x(:));

    if nargout <= 1 
        varargout{1} = [minv, maxv];
    elseif nargout == 2
        varargout{1} = minv;
        varargout{2} = maxv;
    end
    
end
