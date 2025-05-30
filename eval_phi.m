%% ----------------------------- PHI(xi) ----------------------------------

function phi = eval_phi(GP,i)
    
    % cubic Hermite Polynomials p(t) = (2t^3 - 3t^2 + 1)p0 + (t^3 - 2t^2 + t)m0 + (-2t^3 + 3t^2)p1 + (t^3 - t^2)m1
    % here, t (0,1) is a parametrization of a specific domain interval
    % p0: value at left boundary
    % p1: value at right boundary
    % m0: tangent at left boundary
    % m1: tangent at right boundary    

    gp = GP(1);
    pCoord = GP([2,3]);
    pID = GP([4,5]);

    % cubic Hermite Polynomials p(t) = (2t^3 - 3t^2 + 1)p0 + (t^3 - 2t^2 + t)m0 + (-2t^3 + 3t^2)p1 + (t^3 - t^2)m1

    t = (gp-pCoord(1))/(pCoord(2)-pCoord(1)); % local variable
    
    % evaluate cubic Hermite Polynomials
    if i == pID(1)
        phi  = 2*t^3 - 3*t^2 + 1;
    else
        phi = -2*t^3 + 3*t^2;
    end
end