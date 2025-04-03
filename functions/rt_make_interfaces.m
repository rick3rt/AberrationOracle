
function BFC = rt_make_interfaces(P, f_outer, f_inner)

    % functions must be will be evaluated an shifted to zero and then aligned
    
    xv = P.x_piezo;
    z_outer = f_outer(xv);
    z_inner = f_inner(xv);
    
    z_outer = z_outer - min(z_outer); %0 alignment
    z_inner = z_inner - min(z_inner); %0 alignment
    
    
    z_outer = z_outer + P.lens_thickness + P.distance_trans_bone;
    z_inner = z_inner + P.lens_thickness + P.distance_trans_bone + P.bone_thickness;
    
    sp_outer = rt.util.spline_fit_knots(xv, z_outer, 10);
    sp_inner = rt.util.spline_fit_knots(xv, z_inner, 10);
    
    BFC.medium_interfaces = {
        P.lens_thickness,
        sp_outer,
        sp_inner,
    };

end


% BFC.medium_interfaces = {P.lens_thickness,
% [P.bone_curvature, 0, P.lens_thickness + P.distance_trans_bone],
% [P.bone_curvature, 0, P.lens_thickness + P.distance_trans_bone + P.bone_thickness]};