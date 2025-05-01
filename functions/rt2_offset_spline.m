function sp_out = rt2_offset_spline(P, sp, z_offset)
    % utility function to vertically offset a piecewise cubic spline
    % P:        parameters structure
    % sp:       piecewise cubic spline structure
    % z_offset: offset in meters
    % sp_out:   new piecewise cubic spline structure

    xv = P.x_piezo;
    z_sp = rt.util.segeval(sp, xv);
    z_sp0 = z_sp - min(z_sp);
    z_sp_new = z_sp0 + z_offset;

    % fit new spline to offset data
    sp_out = rt.util.spline_fit_knots(xv, z_sp_new, sp.pieces + 1);

end
