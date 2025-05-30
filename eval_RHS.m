%% -------------------- evaluate Right Hand Side ------------------------

function rhs = eval_RHS(a,b,GP)
    h = 1;
    v_t = 2/h; % !!! contraction velocity input !!!
    n = eval_n(a,b,GP);
    n_xi = eval_n_xi(a,b,GP);
    f = eval_f(GP(1));
    g = eval_g(GP(1));
    rhs = v_t*n_xi + f - (f + g)*n;
end