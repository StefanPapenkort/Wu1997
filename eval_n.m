%% --------------------- n(a_i,a_ii,b_i,b_ii,xi) --------------------------

function [n] = eval_n(a, b, GP)
    % p (1,2): left and right boundary point of the interval
    % a (1,2): scaling factor of phi
    % b (1,2): scaling factor of psi
    % gP: gauss point
    gp = GP(1);
    pCoord = GP([2,3]);
    pID = GP([4,5]);

    % cubic Hermite Polynomials p(t) = (2t^3 - 3t^2 + 1)p0 + (t^3 - 2t^2 + t)m0 + (-2t^3 + 3t^2)p1 + (t^3 - t^2)m1

    t = (gp-pCoord(1))/(pCoord(2)-pCoord(1)); % local variable
    dt_dxi = 1/(pCoord(2)-pCoord(1));

    % evaluate cubic Hermite Polynomials
    phi_i  = 2*t^3 - 3*t^2 + 1;
    phi_ii = -2*t^3 + 3*t^2;
    
    psi_i  = (t^3 - 2*t^2 + t)*1/dt_dxi;
    psi_ii = (t^3 - t^2)*1/dt_dxi;
    
    % determine n_hat
    n = a(pID(1))*phi_i + b(pID(1))*psi_i + a(pID(2))*phi_ii + b(pID(2))*psi_ii;
end
