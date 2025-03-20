
function z = segeval(abstractline, x)
    if isstruct(abstractline)
        z = ppval(abstractline, x);
    else
        z = polyval(abstractline(1:end), x);
    end
end
