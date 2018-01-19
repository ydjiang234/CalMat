function fit = fit_fun(vector)
    a = vector(1);
    b = vector(2);
    c = vector(3);
    d = vector(4);
    fit = (a-2)^2 + (b-3)^2 + (c-14)^2 + (d-8)^2;
    fit = fit * -1;
end

