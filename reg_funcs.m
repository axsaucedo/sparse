function regfunc = reg_funcs
    regfunc.abs     = @absval;
    regfunc.squares = @squares;
    regfunc.nccvar  = @nccvar;
    regfunc.cvar    = @cvar;
    regfunc.ridge   = @ridge;
end


function [ pimat_a, value_a ] = absval(I, R, T, n)

    cvx_begin quiet
        variable pimat(n)
        minimize( sum(abs(I - R'*pimat)) )

        subject to
            pimat >= 0
            sum(pimat) == 1
    cvx_end
    
    value_a = sum(abs(I - R'*pimat));
    pimat_a = pimat;
    
end

function [ pimat_q, value_q ] = squares(I,R,T,n)
     
    cvx_begin quiet
        variable pimat(n)

        minimize(sum(norm(I - R'*pimat)))

        subject to
            sum(pimat) <= 1
            pimat >= 0
    cvx_end
    
    value_q = sum(abs(I - R'*pimat));
    pimat_q = pimat;
end

function [ pimat_n, value_n ] = nccvar(I,R,T,n)

    Beta = 0.8;
    k = 1;
    C = 1/sqrt(n) + k*(1-1/sqrt(n))/(10*sqrt(n));
    divcoef = 1 / (T*(1-Beta));
     
    cvx_begin quiet
        variable z_n(T)
        variable Alpha_n
        variable pimat(n)
        minimize( Alpha_n + divcoef * sum(z_n) )

        subject to
            z_n >= 0
            z_n - abs(I - R'*pimat) + Alpha_n >= 0

            pimat >= 0
            sum(pimat) == 1
    cvx_end
    
    value_n = sum(abs(I - R'*pimat));
    pimat_n = pimat;
end

function [ pimat_c, value_c ] = cvar(I,R,T,n)

    Beta = 0.8;
    k = 1;
    C = 1/sqrt(n) + k*(1-1/sqrt(n))/(10*sqrt(n));
    divcoef = 1 / (T*(1-Beta));
     
    cvx_begin quiet
        variable z_c(T)
        variable Alpha_c
        variable pimat(n)
        minimize( Alpha_c + divcoef * sum(z_c) )

        subject to
            z_c >= 0
            z_c - abs(I - R'*pimat) + Alpha_c >= 0

            pimat >= 0
            sum(pimat) == 1
    cvx_end
    
    value_c = sum(abs(I - R'*pimat));
    pimat_c = pimat;
end

function [ pimat_r, value_r ] = ridge(I, R, T, n)

    cvx_begin quiet
        variable pimat(n)
        minimize( sum(norm(I - R'*pimat)) + sum(norm(pimat)) )

        subject to
            pimat >= 0
            sum(pimat) == 1
    cvx_end
    
    value_r = sum(abs(I - R'*pimat));
    pimat_r = pimat;
    
end