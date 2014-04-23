function regfunc = reg_funcs
    regfunc.abs     = @absval;
    regfunc.squares = @squares;
    regfunc.cvar    = @cvar;
    regfunc.ridge   = @ridge;
    regfunc.lasso   = @lasso;
    regfunc.group_lasso         = @group_lasso;
    regfunc.sparse_group_lasso  = @sparse_group_lasso;
end

function [ pimat_a, value_a ] = nccvar(I, R, T, n) end

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

function [ pimat_l, value_l ] = lasso(I, R, T, n)

    cvx_begin quiet
        variable pimat(n)
        minimize( sum(norm(I - R'*pimat)) + sum(norm(pimat)) )

        subject to
            pimat >= 0
            sum(pimat) == 1
    cvx_end
    
    value_l = sum(abs(I - R'*pimat));
    pimat_l = pimat;
    
end

function [ pimat_gl, value_gl ] = group_lasso(I, R, T, n, groups, lambda)

    not_groups = ~groups;
    group_n = size(groups,1);

    cvx_begin quiet
        variable pimat(n,group_n)

        minimize(square_pos(norm(I - sum(R'*pimat,2))) + lambda*sum(sqr_group_sizes'*diag(norms(pimat,2,1))) )

        subject to
            pimat >= 0
            pimat(not_groups') == 0
    cvx_end
    
    % Getting results to add to 1
    pimat=pimat/sum(pimat(:));
    
    value_gl = sum(abs(I - R'*pimat));
    pimat_gl = pimat;
end

function [ pimat_sgl, value_sgl ] = sparse_group_lasso(I, R, T, n, groups, lambda, alpha)

    not_groups = ~groups;
    group_n = size(groups,1);

    cvx_begin quiet
        variable pimat(n,group_n)

        minimize(square_pos(norm(I - sum(R'*pimat,2)))  + (1-alpha)*lambda*sum(sqr_group_sizes'*diag(norms(pimat,2,1))) + alpha*lambda*sum(sum(abs(pimat))) )

        subject to
            pimat >= 0
            pimat(not_groups') == 0
    cvx_end
    
    % Getting results to add to 1
    pimat=pimat/sum(pimat(:));
    
    value_sgl = sum(abs(I - R'*pimat));
    pimat_sgl = pimat;
end