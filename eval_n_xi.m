%% -------------------- n_xi(a_i,a_ii,b_i,b_ii,xi) ------------------------

function [n_xi] = eval_n_xi(a, b, GP)
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
    phi_i_dt  = 6*t^2 - 6*t;
    phi_ii_dt = -6*t^2 + 6*t;
    
    psi_i_dt  = (3*t^2 - 4*t + 1)*1/dt_dxi;
    psi_ii_dt = (3*t^2 - 2*t)*1/dt_dxi;
    
    % determine n_hat
    n_xi = (a(pID(1))*phi_i_dt + b(pID(1))*psi_i_dt + a(pID(2))*phi_ii_dt + b(pID(2))*psi_ii_dt)*dt_dxi;
end