function regfunc = reg_funcs
    regfunc.abs     = @absval;
    regfunc.squares = @squares;
    regfunc.cvar    = @cvar;
    regfunc.ridge   = @ridge;
    regfunc.lasso   = @lasso;
    regfunc.sglasso  = @sparsegrouplasso;
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

    cvx_begin %quiet
        variable pimat(n)
        minimize( sum(norm(I - R'*pimat)) + sum(norm(pimat)) )

        subject to
            pimat >= 0
            sum(pimat) == 1
    cvx_end
    
    value_r = sum(abs(I - R'*pimat));
    pimat_r = pimat;
    
end

function [ pimat_l, value_l ] = lasso(I, R, T, n, l)
    groups = ones(1:n);
    pimat_l = sparsegrouplasso(I,R,groups,0,l);
    
end

function [ pimat_gl, value_gl ] = group_lasso(I, R, T, n, groups, lambda)
    pimat_l = sparsegrouplasso(I,R,groups,l,0);
end

function [ pimat ] = sparsegrouplasso( I, R, groups, l1, l2)
% SPARSEGROUPLASSO - Compute the sparse group lasso minimization
%       I       -   Target Data
%       R       -   Matrix of - Rows:Features, Cols:Observations
%       groups  -   Vector containing the group number each element belongs
%       l1      -   Regularizer which controls group sparsity
%       l2      -   Regularizer which controls feature sparsity

    m = size(unique(groups),1); % Number of groups
    n = size(R,1);              % Number of Assets
    t = size(R,2);              % Time window

    group_idx   = zeros(m,n);   % Each row is a group, columns are logical index
    group_sizes = zeros(m,1);   % Size of each group
    
    % Generating the logical indexes from the 'groups' vector
    for g=1:m
        tmp_logical_idx = groups==g;
        group_idx(g,:) = tmp_logical_idx;
        group_sizes(g) = sum(tmp_logical_idx);
    end
    
    group_idx = logical(group_idx');
    not_group_idx = ~group_idx;             % This is the variable for our constraint
    sqr_group_sizes = sqrt(group_sizes);    % To speed up execution
    
    cvx_begin quiet 
        variable pimat_sgl(n,m)

        % Zero-Constrained Sparse Group Lasso Model
        minimize(sum_square_pos(norm(I - sum(R'*pimat_sgl,2))) + l1*norms(pimat_sgl,2,1)*sqr_group_sizes + l2*sum(abs(pimat_sgl(:))) )

        subject to
            pimat_sgl >= 0 
            % Zero-Constraint
            pimat_sgl(not_group_idx) == 0
    cvx_end
    
    pimat = full(pimat_sgl);
end

