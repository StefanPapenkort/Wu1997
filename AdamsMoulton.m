%% ------------------ NON-LINEAR SYSTEM OF EQUATIONS ----------------------

function res = AdamsMoulton(x, a_prev, dAdt_prev, b_prev, dBdt_prev, GaussPoints, step, dt)

    numGaussPts = size(GaussPoints,1);
    numPts = numGaussPts/2 + 1;
    % extract variables
    a       = x(1:numPts);
    a_prime = x(numPts+1:2*numPts);
    b       = x(2*numPts+1:3*numPts);
    b_prime = x(3*numPts+1:4*numPts);

    res = zeros(4*numPts,1); % create residual vector (Unknowns: a, a', b, b' in every point xi)

    % ----------------------- calculate residuum --------------------------
    % eqilibrium at Gauss Points
    for j=1:numGaussPts
        % i = floor(j/2)+1;
        RHS = eval_RHS(a,b,GaussPoints(j,:));
        ID1 = GaussPoints(j,4);
        ID2 = GaussPoints(j,5);

        phi_i  = eval_phi(GaussPoints(j,:),ID1);
        phi_ii = eval_phi(GaussPoints(j,:),ID2);
        psi_i  = eval_psi(GaussPoints(j,:),ID1);
        psi_ii = eval_psi(GaussPoints(j,:),ID2);

        res(j) = RHS - a_prime(ID1)*phi_i - b_prime(ID1)*psi_i - a_prime(ID2)*phi_ii - b_prime(ID2)*psi_ii;
    end

    % additional equations
    res(numGaussPts+1) = a_prime(1);
    res(numGaussPts+2) = a_prime(end);

    % Adams-Moulton equations
    for i = 1:numPts
        if step == 1
            % ------- Step 1
            res(2*numPts+2*i-1) = a_prev(i) + dt*a_prime(i) - a(i);
            res(2*numPts+2*i)   = b_prev(i) + dt*b_prime(i) - b(i);

        elseif step == 2
            % ------- Step 2
            res(2*numPts+2*i-1) = a_prev(i) + 0.5*dt*(a_prime(i) + dAdt_prev(i,1)) - a(i);
            res(2*numPts+2*i)   = b_prev(i) + 0.5*dt*(b_prime(i) + dBdt_prev(i,1)) - b(i);

        elseif step == 3
            % ------- Step 3
            res(2*numPts+2*i-1) = a_prev(i) + dt*(5/12*a_prime(i) + 8/12*dAdt_prev(i,1) - 1/12*dAdt_prev(i,2)) - a(i);
            res(2*numPts+2*i)   = b_prev(i) + dt*(5/12*b_prime(i) + 8/12*dBdt_prev(i,1) - 1/12*dBdt_prev(i,2)) - b(i);

        elseif step == 4
            % ------- Step 4
            res(2*numPts+2*i-1) = a_prev(i) + dt*(9/24*a_prime(i) + 19/24*dAdt_prev(i,1) - 5/24*dAdt_prev(i,2) + 1/24*dAdt_prev(i,3)) - a(i);
            res(2*numPts+2*i)   = b_prev(i) + dt*(9/24*b_prime(i) + 19/24*dBdt_prev(i,1) - 5/24*dBdt_prev(i,2) + 1/24*dBdt_prev(i,3)) - b(i);

        else
            % ------- Step >= 5
            res(2*numPts+2*i-1) = a_prev(i) + dt*(251/720*a_prime(i) + 646/720*dAdt_prev(i,1) - 264/720*dAdt_prev(i,2) + 106/720*dAdt_prev(i,3) - 19/720*dAdt_prev(i,4)) - a(i);
            res(2*numPts+2*i)   = b_prev(i) + dt*(251/720*b_prime(i) + 646/720*dBdt_prev(i,1) - 264/720*dBdt_prev(i,2) + 106/720*dBdt_prev(i,3) - 19/720*dBdt_prev(i,4)) - b(i);
        end
    end

end
