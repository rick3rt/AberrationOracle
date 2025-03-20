function xroot = poly_roots_valid(p, ray)

    import rt.raytheory.poly_roots_real
    
    xroots = poly_roots_real(p);

    if ray.dir(1) > 0
        xroots(xroots < ray.start(1)) = [];
    else
        xroots(xroots > ray.start(1)) = [];
    end

    [~, idx] = min(abs(xroots - ray.start(1)));
    xroot = xroots(idx);

end
