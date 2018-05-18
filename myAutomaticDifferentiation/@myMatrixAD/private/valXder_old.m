function der = valXder(val, der)

    if length(val) == 1
        der = val * der;
    else
        der = val(:,:,ones(size(der,3),1)) .* der;
    end
