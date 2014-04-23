randn('state',23432);
rand('state',3454);

%% MAIN PROGRAM %%
d = datasets;
rf = reg_funcs;

lambda = 5;
l1_groups = 0;
l2_features = 0.1;

% Get Data
[ I, R, totalassets, T ] = d.ftse100();
% [ I, R, totalassets, T ] = d.naive(200,300);


%% MAIN CVX PROGRAM %%


% Calculating n groups with knn + kmeans
n = 5;
[index, groups] = knn(R, n);
combs = [];
size_combs = [];

logical_index = logical(zeros(n,totalassets));
for i=1:5
    logical_index(i,:) = logical(index==i)';
end

comb_idx = [];
group_idx = [];

% Creating indexes for combinations
for i = 1:n
    comb = combnk(1:n, i);    
    
    comb_rn = size(comb,1);
    comb_cn = size(comb,2);
    
    combs = [combs; [comb , zeros(comb_rn,n-comb_cn) ]];
    size_combs = [size_combs; ones(comb_rn,1)*i];
    
    for j=1:comb_rn
        tmp_comb_idx = logical(zeros(1,n));
        tmp_comb_idx(comb(j,:)) = true;
        
        tmp_group_idx = sum(logical_index(comb(j,:),:),1);
        
        comb_idx = [ comb_idx; tmp_comb_idx ];
        group_idx = [ group_idx; tmp_group_idx ];
    end
end


group_idx = logical(group_idx);
group_idx = group_idx(1:n,:);
% n=2;
% group_idx = group_idx(1:2,:);
% group_idx(2,:) = ~group_idx(1,:);

not_group_idx = ~group_idx;
group_sizes = sum(group_idx,2);
sqr_group_sizes = sqrt(group_sizes);



% %% Values for iterating on Group Lasso Lambdas
% lambda = 0;
% limit = 50;
% iterations = 200;
% 
% errors_gl = [];
% zeros_gl = [];
% lambdas_gl = [];
% elapsed_gl = [];
% pimat_sums_gl = [];
% 
% 
% for i=1:iterations
%     i
%     
%     % CVX to find optimal value for Group Lasso
%     tic
%     cvx_begin quiet
%         variable pimat_gl(totalassets,n)
% 
%         minimize(sum_square_pos(norm(I - sum(R'*pimat_gl,2))) + lambda*norms(pimat_gl,2,1)*sqr_group_sizes )
% 
%         subject to
%             pimat_gl >= 0
%             pimat_gl(not_group_idx') == 0
%     cvx_end
%     elapsed=toc
% 
%     pimat_gl        = full(pimat_gl);
%     pimat_norm_gl   = full(pimat_gl/sum(pimat_gl(:)));
% 
% %     sum(pimat_norm_gl,2)
% %     sum(pimat_norm_gl,1)
% %     sum(pimat_gl(:))
% %     sum(abs(I-sum(R'*pimat_norm_gl,2)))
%     
%     
%     error_gl = sum(abs(I - sum(R'*pimat_norm_gl,2)));
%     zero_gl = sum(sum(pimat_norm_gl,2)<0.0001);
%     
%     errors_gl    = [ errors_gl, error_gl ];
%     zeros_gl     = [ zeros_gl, zero_gl ];
%     lambdas_gl    = [ lambdas_gl, lambda ];
%     elapsed_gl   = [ elapsed_gl,  elapsed];
%     pimat_sums_gl = [ pimat_sums_gl;  sum(pimat_norm_gl,1)];
%     
%     [elapsed, zero_gl, lambda]
%     
%     lambda = lambda + limit/iterations;
% end



%% Values for iterating on Sparse Group Lasso Lambdas
% l2_features = 4;
% l1_groups = 0.0;

% l2_features = 0;
% l1_groups = 0;

divcoef=10000;

l2_features = 1;
l1_groups = 0;


iterations = 1;

errors_sgl = [];
zeros_sgl = [];
lambdas_sgl = [];
elapsed_sgl = [];
pimat_sums_sgl = [];
pimats_sgl = [];
zeros_groups_sgl = [];



% CVX to find optimal value for Lasso
for i=1:iterations
    i
    
    % CVX to find optimal value for Group Lasso
    tic
    cvx_begin %quiet
        variable pimat_sgl(totalassets)

        minimize(sum_square_pos(norm(I - R'*pimat_sgl)) + l1_groups*norms(pimat_sgl,2,1)*sqr_group_sizes + l2_features*sum(abs(pimat_sgl(:))) )

        subject to
    %         sum(pimat_sgl(:)) == 1
    %         pimat_sgl >= 0 
            pimat_sgl(not_group_idx') == 0
    cvx_end
    elapsed=toc

    pimat_sgl        = full(pimat_sgl);
    pimat_norm_sgl   = full(pimat_sgl/sum(pimat_sgl(:)));

%     sum(pimat_norm_gl,2)
%     sum(pimat_norm_gl,1)
%     sum(pimat_gl(:))
%     sum(abs(I-sum(R'*pimat_norm_gl,2)))
    
    
    
    zeros_groups = [];
    for i=1:n
        zeros_groups = [ zeros_groups sum(pimat_norm_sgl(group_idx(i,:),i)<0.0001)' ]; 
    end

    error_sgl       = sum(abs(I - sum(R'*pimat_norm_sgl,2)))
    zero_sgl        = sum(sum(pimat_norm_sgl,2)<0.0001)
    
    errors_sgl      = [ errors_sgl, error_sgl ];
    zeros_sgl       = [ zeros_sgl, zero_sgl ];
    elapsed_sgl     = [ elapsed_sgl,  elapsed];
    pimat_sums_sgl  = [ pimat_sums_sgl;  sum(pimat_norm_sgl,1)];
    pimats_sgl      = [ pimats_sgl sum(pimat_norm_sgl,2) ];
    zeros_groups_sgl = [zeros_groups_sgl; zeros_groups];
    
    [elapsed, zero_sgl, lambda]
    
    lambdas_sgl     = [ lambdas_sgl, l2_features ];
%     l2_features = l2_features + 0.05;
    
%     l1_groups = l1_groups + 0.00005;
%     l1_groups = l1_groups + 0.5;

    l2_features = l2_features + 0.5;
%     l1_groups = l2_features/divcoef;
    
%     l1_groups       = l1_groups + 100;
%     l2_features     = l2_features + 100;
    
end

% figure(1);
% h1=plot(lambdas_sgl, pimat_sums_sgl)
% legend('Group 1','Group 2','Group 3','Group 4','Group 5')
% savegraph(h1,'Lambdas','Group Weights','','')
% 
% 
% figure(2);
% h2=plot(lambdas_sgl, pimats_sgl')
% savegraph(h2,'Lambdas','Feature Weights','','')
% 
% 
% figure(3);
% h3=plot(lambdas_sgl,zeros_groups_sgl)
% legend('Group 1','Group 2','Group 3','Group 4','Group 5')
% savegraph(h3,'Lambdas','# Zeros per Group','','')
% 
% 
% figure(4);
% h4=plot(lambdas_sgl,zeros_sgl)
% savegraph(h4,'Lambdas','# Zeros','','')
% 
% 
% t=table([zeros_sgl', errors_sgl']);
% writetable(t,strcat('table',num2str(divcoef),'.csv'));



% pimat_sgl = full(pimat_sgl);
% pimat_norm_sgl=full(pimat_sgl/sum(pimat_sgl(:)));
% sum(pimat_norm_sgl,2)
% sum_of_each_group = sum(pimat_sgl,1)
% original_pimat_result = sum(pimat_sgl(:))
% tracking_error = sum(abs(I-sum(R'*pimat_norm_sgl,2)))
% total_number_of_zeros = sum(sum(pimat_norm_sgl,2)<0.0001)
% zeros_groups = [];
% for i=1:n
%     zeros_groups = [ zeros_groups; sum(group_idx(i,:)) sum(pimat_norm_sgl(group_idx(i,:),i)<0.0001) ]; 
% end
% zeros_groups


% pimat_norm_gl = full(pimat_norm_gl);
% pimat_norm_sgl = full(pimat_norm_sgl);
% 
% p = [ full(sum(pimat_norm_gl,2)) full(sum(pimat_norm_sgl,2)) ]
% t = [ sum(pimat_norm_gl,1) sum(pimat_norm_sgl,1) ]
% e = [ sum(abs(I-sum(R'*pimat_norm_gl,2))) sum(abs(I-sum(R'*pimat_norm_sgl,2))) ]
% z = [ sum(sum(pimat_norm_gl,2)<0.0001) sum(sum(pimat_norm_sgl,2)<0.0001) ]









%% Lasso SPAMS
% 
% 
% % clear all;
% format compact;
% randn('seed',0);
% param.numThreads=-1; % all cores (-1 by default)
% param.verbose=true;   % verbosity, false by default
% param.lambda=0.05; % regularization parameter
% param.it0=10;      % frequency for duality gap computations
% param.max_it=200; % maximum number of iterations
% param.L0=0.1;
% param.tol=1e-3;
% param.intercept=false;
% param.pos=true;
% 
% param.loss='square';
% param.compute_gram=true;


% X=randn(100,200);
% X=X-repmat(mean(X),[size(X,1) 1]);
% X=mexNormalize(X);
% Y=randn(100,1);
% Y=Y-repmat(mean(Y),[size(Y,1) 1]);
% Y=mexNormalize(Y);
% W0=zeros(size(X,2),size(Y,2));
% groups=int32(randi(5,1,size(X,2)));


% X = mexNormalize(R');
% Y = mexNormalize(I);
% W0=zeros(size(X,2),size(Y,2));
% 
% groups=int32(index);
% 
% lambda = 0.0000;
% 
% zeros = [];
% lambdas = [];
% te = [];
% 
% 
% fprintf('\nFISTA + Group Lasso L2 with variable size of groups \n');
% param.regul='group-lasso-l2';
% param2=param;
% param2.groups=groups;
% 
% while lambda < 1 
% param2.lambda=lambda;
% lambda
% [ignore, W_lasso_var, optim_info]= evalc('mexFistaFlat(Y,X,W0,param2);');
% zeros = [zeros sum(W_lasso_var==0)];
% lambdas = [lambdas lambda ];
% te = [te sum(abs(I-R'*W_lasso_var))];
% lambda = lambda + 0.0001;
% end
% 
% fprintf('\nFISTA + Group Lasso L2 with variable size of groups \n');
% param.regul='sparse-group-lasso-l2';
% param2=param;
% param2.groups=groups;  % all the groups are of size 2
% param2.lambda=0.0143295;
% param2.lambda2=0.07;
% % tpm_param=param2
% tic
% [W_sparse_var optim_info]=mexFistaFlat(Y,X,W0,param2);
% t=toc;
% 
% 
% [ W_lasso_var W_sparse_var ]
% [ sum(W_lasso_var==0) sum(W_sparse_var==0) ]
% [ sum(abs(I-R'*W_lasso_var)) sum(abs(I-R'*W_sparse_var)) ]


% sum_groups      = []; % 5 x 99
% zero_groups    = []; % 5 x 99
% 
% for i=1:99
%     curr = pimats_a(:,i)
%     
%     zero_group_row = [];
%     sum_group_row = [];
%     for j=1:5
%         curr_group = curr(group_idx(j,:));
%         zero_group_row = [zero_group_row, sum(curr_group==0) ];
%         sum_group_row = [sum_group_row,  sum(curr_group) ];
%     end
%     
%     zero_groups = [ zero_groups; zero_group_row ];
%     sum_groups = [ sum_groups; sum_group_row ];
% end