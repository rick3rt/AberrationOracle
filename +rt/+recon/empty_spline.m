function sp = empty_spline()

    sp = spline([0 1],[0 1]);
    sp.coefs = zeros(sp.pieces, 4);
    sp.order = 4;
    
end