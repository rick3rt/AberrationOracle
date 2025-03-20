

function xroots = poly_roots_real(p)
    xroots = roots(p);
    xroots = xroots(imag(xroots) == 0);
end
