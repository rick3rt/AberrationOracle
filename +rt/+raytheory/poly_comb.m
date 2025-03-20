

function [p] = poly_comb(p1, p2)
    n1 = numel(p1);
    n2 = numel(p2);

    if (n1 > n2)
        p2 = [zeros(1, n1 - n2) p2];
    elseif (n1 < n2)
        p1 = [zeros(1, n2 - n1) p1];
    end

    p = p2 - p1;
end
