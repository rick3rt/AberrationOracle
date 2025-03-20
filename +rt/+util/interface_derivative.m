function z = interface_derivative(abstractline)
    if isstruct(abstractline)
        z = fnder(abstractline);
    else
        z = polyder(abstractline);
    end
end
