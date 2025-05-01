function [data, Pars] = three_layer_recon_wrapper(Pars, RF_I, RF_Q)

    % if nargin == 2 then RF_I is data struct, so only use plot routines
    if nargin == 2 && isstruct(RF_I)
        data = RF_I;
    end

    % make all vectors in Pars row vectors
    Pars = fixVecPars(Pars);
    Pars.SutureCheck = getFieldIfExist(Pars, 'SutureCheck', 0, 0); % backward compat
    % create Setup
    Setup = [
             Pars.shape_factor_aniso;
             Pars.lens_thickness;
             Pars.pitch;
             Pars.max_depth_bone;
             Pars.min_thickness_bone;
             Pars.max_thickness_bone;
             Pars.lens_half_opening_angle;
             Pars.min_Fnumber;
             Pars.max_error_receive_angle;
             Pars.Fs;
             Pars.lens_c;
             Pars.tissue_c;
             Pars.brain_c;
             Pars.axial_c;
             Pars.radial_c;
             Pars.n_cores;
             Pars.recon_to;
             Pars.NeedTravelTime; % flag = 0, 1, 2
             Pars.SubApertureApodis;
             Pars.PeriPolyOrder;
             Pars.EndoPolyOrder;
             Pars.SutureCheck;
             getFieldIfExist(Pars, 'SutureX1', Pars.SutureCheck, 0);
             getFieldIfExist(Pars, 'SutureX2', Pars.SutureCheck, 0);
             ];

    Setup = Setup(:).'; % make column vector

    % setup imaging grid
    if ~isfield(Pars, 'xv_recon')
        xv_recon = -Pars.img_width / 2:Pars.pixel_size:Pars.img_width / 2;
        zv_recon = 0:Pars.pixel_size:Pars.img_height;
        Pars.xv_recon = xv_recon;
        Pars.zv_recon = zv_recon;
    else
        disp('Reusing reconstruction grid.')
    end

    % select sub aperture if needed
    if isfield(Pars, 'cropAperture') && Pars.cropAperture
        x_ind2 = find(Pars.X_El <= max(Pars.xv_recon), 1, 'last');
        x_ind1 = find(Pars.X_El >= min(Pars.xv_recon), 1, 'first');
        x_ind = x_ind1:x_ind2;
        Pars.x_ind = x_ind;
        RF_I = RF_I(:, x_ind, :);
        RF_Q = RF_Q(:, x_ind, :);
    else
        x_ind = 1:numel(Pars.X_El);
    end

    % init poly vectors if order not right
    if ~isstruct(Pars.PeriParabIn) && all(Pars.PeriParabIn == 0)
        Pars.PeriParabIn = zeros(1, Pars.PeriPolyOrder + 3);
    end
    if ~isstruct(Pars.EndoParabIn) && all(Pars.EndoParabIn == 0)
        Pars.EndoParabIn = zeros(1, Pars.EndoPolyOrder + 3);
    end

    % transducer position
    X_El = Pars.X_El(x_ind);
    Z_El = Pars.Z_El(x_ind);

    % coordinates virtual sources
    XS = Pars.XS;
    ZS = Pars.ZS;

    % reconstruction image grid
    X = Pars.xv_recon;
    Z = Pars.zv_recon;
    [XV, ZV] = meshgrid(X, Z);


    ProcessFunction = @rt.recon.three_layer_recon;

    if Pars.recon_to == 1
        tic
        [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue, ...
             Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue, APixList_R_Tissue] = ...
            ProcessFunction( ...
            Setup, X_El, Z_El, XS, ZS, X, Z, ...
            RF_I, RF_Q, Pars.tx_to_use, Pars.add_to_delay_firing, ...
            Pars.ReceiveAngleRad, Pars.PeriParabIn, Pars.EndoParabIn);
        fprintf('ThreeLayerRecon Tissue time: %f\n', toc);

        mask_tissue = true(size(ZV));

        % Create a B-mode image from I/Q images:
        % sum over transmissions for image with fixed receive angle and F-number
        % if multiple receive angles, then receive angle must be chosen with
        % dimension 3 in Beam_I and Beam_Q
        ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue, 4) .^ 2 + sum(Beam_Q_Tissue, 4) .^ 2));
        % Put in data
        IQ_Tissue = squeeze(complex(sum(Beam_I_Tissue, 4), sum(Beam_Q_Tissue, 4)));
        IQabb = IQ_Tissue;
        IQcorr = IQ_Tissue;
        FullImage = ImageTissue;

    elseif Pars.recon_to == 2
        tic
        [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue, ...
             Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue, APixList_R_Tissue, ...
             ImageBone, Time_T_Bone, Time_R_Bone, Angle_T_Bone, ...
             Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone, ...
             APix_T_Bone, APix_R_Bone, APixList_R_Bone, PeriParab, EndoParab] = ...
            ProcessFunction( ...
            Setup, X_El, Z_El, XS, ZS, X, Z, ...
            RF_I, RF_Q, Pars.tx_to_use, Pars.add_to_delay_firing, ...
            Pars.ReceiveAngleRad, Pars.PeriParabIn, Pars.EndoParabIn);
        fprintf('ThreeLayerRecon Tissue and Bone time: %f\n', toc);

        P_Periosteum = PeriParab;
        P_Endosteum = EndoParab;
        Z_Peri = segeval(P_Periosteum, X);
        Z_Endo = segeval(P_Endosteum, X);

        % create masks
        mask_tissue = ZV < Z_Peri;
        mask_bone = ZV >= Z_Peri; %& ZV < Z_Endo;
        mask_brain = false(size(mask_bone)); % ZV >= Z_Endo;

        % Create a B-mode image from I/Q images:
        % sum over transmissions for image with fixed receive angle and F-number
        % if multiple receive angles, then receive angle must be chosen with
        % dimension 3 in Beam_I and Beam_Q
        ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue, 4) .^ 2 + sum(Beam_Q_Tissue, 4) .^ 2));
        ImageBone_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Bone, 4) .^ 2 + sum(Beam_Q_Bone, 4) .^ 2));

        pixel_size = diff(Pars.zv_recon([1 2]));
        FullImage = ImageTissue;
        FullImage(mask_bone) = ImageBone(mask_bone);
        
        FullImage_FixedReceiveAng = ImageTissue_FixedReceiveAng;
        FullImage_FixedReceiveAng(mask_bone) = ImageBone_FixedReceiveAng(mask_bone);
        
        IQ_Tissue = squeeze(complex(sum(Beam_I_Tissue, 4), sum(Beam_Q_Tissue, 4)));
        IQ_Bone = squeeze(complex(sum(Beam_I_Bone, 4), sum(Beam_Q_Bone, 4)));
        IQcorr = IQ_Tissue;
        IQcorr(mask_bone) = IQ_Bone(mask_bone);
        IQabb = IQ_Tissue;

    elseif Pars.recon_to == 3

        tic
        [ImageTissue, Time_T_Tissue, Time_R_Tissue, Angle_T_Tissue, ...
             Angle_R_Tissue, Beam_I_Tissue, Beam_Q_Tissue, APix_T_Tissue, APix_R_Tissue, APixList_R_Tissue, ...
             ImageBone, Time_T_Bone, Time_R_Bone, Angle_T_Bone, ...
             Angle_R_Bone, Beam_I_Bone, Beam_Q_Bone, APix_T_Bone, APix_R_Bone, APixList_R_Bone, ...
             ImageMarrow, Time_T_Marrow, Time_R_Marrow, Angle_T_Marrow, ...
             Angle_R_Marrow, Beam_I_Marrow, Beam_Q_Marrow, APix_T_Marrow, APix_R_Marrow, APixList_R_Marrow, ...
             PeriParab, EndoParab] = ...
            ProcessFunction( ...
            Setup, X_El, Z_El, XS, ZS, X, Z, ...
            RF_I, RF_Q, Pars.tx_to_use, Pars.add_to_delay_firing, ...
            Pars.ReceiveAngleRad, Pars.PeriParabIn, Pars.EndoParabIn);
        fprintf('ThreeLayerRecon Tissue, Bone and Marrow time: %f\n', toc);

        P_Periosteum = PeriParab;
        P_Endosteum = EndoParab;
        Z_Peri = segeval(P_Periosteum, X);
        Z_Endo = segeval(P_Endosteum, X);

        % create masks
        mask_tissue = ZV < Z_Peri;
        mask_bone = ZV >= Z_Peri & ZV < Z_Endo;
        mask_brain = ZV >= Z_Endo;

        % Create a B-mode image from I/Q images:
        % sum over transmissions for image with fixed receive angle and F-number
        % if multiple receive angles, then receive angle must be chosen with
        % dimension 3 in Beam_I and Beam_Q
        ImageTissue_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Tissue, 4) .^ 2 + sum(Beam_Q_Tissue, 4) .^ 2));
        ImageBone_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Bone, 4) .^ 2 + sum(Beam_Q_Bone, 4) .^ 2));
        ImageMarrow_FixedReceiveAng = squeeze(sqrt(sum(Beam_I_Marrow, 4) .^ 2 + sum(Beam_Q_Marrow, 4) .^ 2));

        FullImage = ImageTissue;
        FullImage(mask_bone) = ImageBone(mask_bone);
        FullImage(mask_brain) = ImageMarrow(mask_brain);

        FullImage_FixedReceiveAng = ImageTissue_FixedReceiveAng;
        FullImage_FixedReceiveAng(mask_bone) = ImageBone_FixedReceiveAng(mask_bone);
        FullImage_FixedReceiveAng(mask_brain) = ImageMarrow_FixedReceiveAng(mask_brain);

        IQ_Tissue = squeeze(complex(sum(Beam_I_Tissue, 4), sum(Beam_Q_Tissue, 4)));
        IQ_Bone = squeeze(complex(sum(Beam_I_Bone, 4), sum(Beam_Q_Bone, 4)));
        IQ_Marrow = squeeze(complex(sum(Beam_I_Marrow, 4), sum(Beam_Q_Marrow, 4)));
        IQcorr = IQ_Tissue;
        IQcorr(mask_bone) = IQ_Bone(mask_bone);
        IQcorr(mask_brain) = IQ_Marrow(mask_brain);
        IQabb = IQ_Tissue;

        IQ_Tissue_angles = squeeze(complex(Beam_I_Tissue, Beam_Q_Tissue));
        IQ_Bone_angles = squeeze(complex(Beam_I_Bone, Beam_Q_Bone));
        IQ_Marrow_angles = squeeze(complex(Beam_I_Marrow, Beam_Q_Marrow));

        mask_bone2 = repmat(mask_bone, 1, 1, size(IQ_Tissue_angles, 3));
        mask_brain2 = repmat(mask_brain, 1, 1, size(IQ_Tissue_angles, 3));

        IQcorr_angles = IQ_Tissue_angles;
        IQcorr_angles(mask_bone2) = IQ_Bone_angles(mask_bone2);
        IQcorr_angles(mask_brain2) = IQ_Marrow_angles(mask_brain2);
        IQabb_angles = IQ_Tissue_angles;

    end

    % set output arguments
    if exist('ImageTissue', 'var'); data.ImageTissue = ImageTissue; end
    if exist('Time_T_Tissue', 'var'); data.Time_T_Tissue = Time_T_Tissue; end
    if exist('Time_R_Tissue', 'var'); data.Time_R_Tissue = Time_R_Tissue; end
    if exist('Angle_T_Tissue', 'var'); data.Angle_T_Tissue = Angle_T_Tissue; end
    if exist('Angle_R_Tissue', 'var'); data.Angle_R_Tissue = Angle_R_Tissue; end
    if exist('Beam_I_Tissue', 'var'); data.Beam_I_Tissue = Beam_I_Tissue; end
    if exist('Beam_Q_Tissue', 'var'); data.Beam_Q_Tissue = Beam_Q_Tissue; end
    if exist('APix_T_Tissue', 'var'); data.APix_T_Tissue = APix_T_Tissue; end
    if exist('APix_R_Tissue', 'var'); data.APix_R_Tissue = APix_R_Tissue; end
    if exist('APixList_R_Tissue', 'var'); data.APixList_R_Tissue = APixList_R_Tissue; end
    if exist('ImageBone', 'var'); data.ImageBone = ImageBone; end
    if exist('Time_T_Bone', 'var'); data.Time_T_Bone = Time_T_Bone; end
    if exist('Time_R_Bone', 'var'); data.Time_R_Bone = Time_R_Bone; end
    if exist('Angle_T_Bone', 'var'); data.Angle_T_Bone = Angle_T_Bone; end
    if exist('Angle_R_Bone', 'var'); data.Angle_R_Bone = Angle_R_Bone; end
    if exist('Beam_I_Bone', 'var'); data.Beam_I_Bone = Beam_I_Bone; end
    if exist('Beam_Q_Bone', 'var'); data.Beam_Q_Bone = Beam_Q_Bone; end
    if exist('APix_T_Bone', 'var'); data.APix_T_Bone = APix_T_Bone; end
    if exist('APix_R_Bone', 'var'); data.APix_R_Bone = APix_R_Bone; end
    if exist('APixList_R_Bone', 'var'); data.APixList_R_Bone = APixList_R_Bone; end
    if exist('ImageMarrow', 'var'); data.ImageMarrow = ImageMarrow; end
    if exist('Time_T_Marrow', 'var'); data.Time_T_Marrow = Time_T_Marrow; end
    if exist('Time_R_Marrow', 'var'); data.Time_R_Marrow = Time_R_Marrow; end
    if exist('Angle_T_Marrow', 'var'); data.Angle_T_Marrow = Angle_T_Marrow; end
    if exist('Angle_R_Marrow', 'var'); data.Angle_R_Marrow = Angle_R_Marrow; end
    if exist('Beam_I_Marrow', 'var'); data.Beam_I_Marrow = Beam_I_Marrow; end
    if exist('Beam_Q_Marrow', 'var'); data.Beam_Q_Marrow = Beam_Q_Marrow; end
    if exist('APix_T_Marrow', 'var'); data.APix_T_Marrow = APix_T_Marrow; end
    if exist('APix_R_Marrow', 'var'); data.APix_R_Marrow = APix_R_Marrow; end
    if exist('APixList_R_Marrow', 'var'); data.APixList_R_Marrow = APixList_R_Marrow; end
    if exist('PeriParab', 'var'); data.PeriParab = PeriParab; end
    if exist('EndoParab', 'var'); data.EndoParab = EndoParab; end
    if exist('Z_Peri', 'var'); data.Z_Peri = Z_Peri; end
    if exist('Z_Endo', 'var'); data.Z_Endo = Z_Endo; end

    if exist('FullImage', 'var'); data.FullImage = FullImage; end
    if exist('FullImage_FixedReceiveAng', 'var'); data.FullImage_FixedReceiveAng = FullImage_FixedReceiveAng; end
    if exist('ImageTissue_FixedReceiveAng', 'var'); data.ImageTissue_FixedReceiveAng = ImageTissue_FixedReceiveAng; end
    if exist('ImageBone_FixedReceiveAng', 'var'); data.ImageBone_FixedReceiveAng = ImageBone_FixedReceiveAng; end
    if exist('ImageMarrow_FixedReceiveAng', 'var'); data.ImageMarrow_FixedReceiveAng = ImageMarrow_FixedReceiveAng; end

    if exist('IQcorr', 'var'); data.IQcorr = IQcorr; end
    if exist('IQabb', 'var'); data.IQabb = IQabb; end
    if exist('IQangles', 'var'); data.IQangles = IQangles; end

    if exist('IQ_Tissue', 'var'); data.IQ_Tissue = IQ_Tissue; end
    if exist('IQ_Bone', 'var'); data.IQ_Bone = IQ_Bone; end
    if exist('IQ_Marrow', 'var'); data.IQ_Marrow = IQ_Marrow; end

    if exist('IQ_Tissue_angles', 'var'); data.IQ_Tissue_angles = IQ_Tissue_angles; end
    if exist('IQ_Bone_angles', 'var'); data.IQ_Bone_angles = IQ_Bone_angles; end
    if exist('IQ_Marrow_angles', 'var'); data.IQ_Marrow_angles = IQ_Marrow_angles; end

    if exist('IQcorr_angles', 'var'); data.IQcorr_angles = IQcorr_angles; end
    if exist('IQabb_angles', 'var'); data.IQabb_angles = IQabb_angles; end

    if exist('mask_tissue', 'var'); data.mask_tissue = mask_tissue; end
    if exist('mask_bone', 'var'); data.mask_bone = mask_bone; end
    if exist('mask_brain', 'var'); data.mask_brain = mask_brain; end

    data.X_El = X_El;
    data.Z_El = Z_El;
    data.XS = XS;
    data.ZS = ZS;
    data.X = X;
    data.Z = Z;

end

function z = segeval(abstractline, x)
    if isstruct(abstractline)
        z = ppval(abstractline, x);
    else
        z = polyval(abstractline(1:end - 2), x);
    end
end

function value = getFieldIfExist(P, field, required, default)
    if ~isfield(P, field) && required
        error('Field %s is required.', field);
    elseif ~isfield(P, field)
        value = default;
    else
        value = P.(field);
    end
end
