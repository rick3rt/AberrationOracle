function sp_out = rt2_offset_spline(P, sp, z_offset)

    xv = P.x_piezo; 

    z_sp = rt.util.segeval(sp, xv);
    z_sp0 = z_sp - min(z_sp);
    z_sp_new = z_sp0 + z_offset;

    sp_out = rt.util.spline_fit_knots(xv, z_sp_new, sp.pieces+1);
    
    
    % z_sp_out = rt.util.segeval(sp_out, xv);   
    % sp
    % sp_out
    % 
    % figure(99);clf;
    % plot(xv, z_sp)
    % hold on 
    % plot(xv, z_sp0)
    % plot(xv, z_sp_out)
    % set(gca,'YDir','reverse')

end