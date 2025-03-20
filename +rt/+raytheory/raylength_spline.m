function [lambda,xintersect] = raylength_spline(sp, start, dir, debug)

    % convert to line equation
    a = dir(2) / dir(1);
    b = start(2) - a * start(1);
    p_line = [a b];

    % check which sections to check
    if dir(1) == 0
        % vertical ray, xintersect is at xstart
        lambda = ppval(sp, start(1))-start(2);
        xintersect = start(1);
        return 
    elseif dir(1) < 0 % CHECK SIGN
        if debug; disp('right to left'); end
        step = -1;
        imax = 1;
        [~, imin] = min(abs(sp.breaks - start(1)));
        imin = min(imin, sp.pieces);
    elseif dir(1) > 0
        if debug; disp('left to right'); end
        step = 1;
        imax = sp.pieces;
        [~, imin] = min(abs(sp.breaks - start(1)));
        imin = imin - 1;
        imin = max(1, imin);
    end
    
    if debug; fprintf('[raylength_spline] imin:step:imax - %i:%i:%i\n',  imin,step,imax); end
    
    % k_piece = -1;
    % find section that will intersect
    for k = imin:step:imax
        x1 = sp.breaks(k);
        x2 = sp.breaks(k + 1);

        % otherwise find the root with the intersecting section
        p_line_shift = p_line + [0 p_line(1) * x1]; % shift line by break
        p = intersect_cubic_line(sp.coefs(k, :), p_line_shift);
        
        if (abs(p(1)) > 1e-10)
            potential_roots = cubic_real_roots(p) + x1; % shift roots back by break
        else
            potential_roots = quad_real_roots(p(2:end)) + x1; % shift roots back by break
        end

        % check if roots are within breaks
        root = potential_roots(potential_roots >= x1 & potential_roots <= x2);

        if numel(root) > 0
            [~, idx] = min(abs(root - start(1)));
            root = root(idx);
            break;
        end

    end
    if debug; fprintf('[raylength_spline] intersection piece %i, %i roots found\n',  k, numel(root) ); end


    if numel(root) == 0; lambda = NaN; xintersect = NaN; return ; end
    % and compute the final ray length
    xintersect = root;
    zintersect = ppval(sp, xintersect);
    lambda = vecdist(start(1), start(2), xintersect, zintersect);
    

end

function roots = quad_real_roots(p)

    a=  p(1); b = p(2); c = p(3);
    %  compute discriminant
    D = b * b - 4 * a * c;

    % find the roots
    if (D > 0.0)
        % discriminant is positive, two real roots
        roots(1) = (-b + sqrt(D)) / (2 * a);
        roots(2) = (-b - sqrt(D)) / (2 * a);
    elseif (D == 0)
        % discriminant is zero, one real root
        roots(1) = -b / (2 * a);
    else
        roots = [];
    end

end

function real_roots = cubic_real_roots(p)
    % https://www.particleincell.com/2013/cubic-line-intersection/

    a = p(1); b = p(2); c = p(3); d = p(4);
    A = b / a; B = c / a; C = d / a;
    Q = (3 * B - A ^ 2) / 9;
    R = (9 * A * B - 27 * C - 2 * A ^ 3) / 54;
    D = Q ^ 3 + R ^ 2;
    
    if (D >= 0)
        % fprintf('D>=0  %.2f\n', D)
        S = sign(R + sqrt(D)) * abs(R + sqrt(D)) ^ (1/3);
        T = sign(R - sqrt(D)) * abs(R - sqrt(D)) ^ (1/3);

        t(1) = -A / 3 + (S + T);
        %t(2) = -A / 3 - (S + T) / 2; % real part of complex root
        %t(3) = -A / 3 - (S + T) / 2; % real part of complex root
        %Im = abs(sqrt(3) * (S - T) / 2);
    else
        % fprintf('D<0  %.2f\n', D)
        th = acos(R / sqrt(- (Q ^ 3)));

        t(1) = 2 * sqrt(-Q) * cos(th / 3) - A / 3;
        t(2) = 2 * sqrt(-Q) * cos((th + 2 * pi) / 3) - A / 3;
        t(3) = 2 * sqrt(-Q) * cos((th + 4 * pi) / 3) - A / 3;

    end
    
    real_roots = t;
end

function p = intersect_cubic_line(p_cube, p_line)
    % given the a cubic polynomial of the form
    %   y = a x^3 + b x^2 + c x + d
    % and a line of the form
    %   y = e x + f
    % construct the polynomial
    %   y = a x^3 + b x^2 + (c - e) x + (d - f)
    %
    p = p_cube;
    p(3) = p_cube(3) - p_line(1);
    p(4) = p_cube(4) - p_line(2);
end

function dist = vecdist(x1, y1, x2, y2)
    dist = sqrt((x1 - x2) .^ 2 + (y1 - y2) .^ 2);
end