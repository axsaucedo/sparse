function [ pimat ] = SparseGroupLasso( I, R, groups, l1, l2)
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

