function f = eval_f(xi)
    if xi < 0 || xi > 1
       f = 0;
    else
       f = 43.3*xi;
    end
end