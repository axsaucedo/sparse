function regfunc = reg_funcs
    regfunc.abs = @absval;
end

function a = absval(I, R, pimat)
    a = sum(abs(I - R'*pimat));
end