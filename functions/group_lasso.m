randn('state',23432);
rand('state',3454);

%% MAIN PROGRAM %%
ds = datasets;
rf = reg_funcs;

lambda = 5;
l1_groups = 0;
l2_features = 0.1;

% Get Data
[ I, R, totalassets, T, index ] = ds.ftse100();            % FTSE with sectors as groups
% [ I, R, totalassets, T, index ] = ds.ftse100(9);             % FTSE with 3 groups created randomly
% [ I, R, totalassets, T, index ] = ds.naive(200,300,9);    % Stochastic with n number of groups
% [ I, R, totalassets, T, index ]=ds.monte_carlo(100,300,.1,.3,9);
n = size(unique(index),1);




T_train     = floor(0.75*T);
T_test      = T - T_train;

R_train     = R(:,1:T_train);
R_test      = R(:,T_train+1:end);

I_train     = I(1:T_train,:);
I_test      = I(T_train+1:end,:);


%% MAIN CVX PROGRAM %%


% Calculating n groups with knn + kmeans
combs = [];
size_combs = [];

logical_index = logical(zeros(n,totalassets));
for i=1:n
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


%% Values for iterating on Sparse Group Lasso Lambdas
% l2_features = 4;
% l1_groups = 0.0;

% l2_features = 0;
% l1_groups = 0;

divcoef=500;

l2_features =  0;
l1_groups = 0;


iterations = 50;

errors_sgl = [];
zeros_sgl = [];
l1_sgl = [];
l2_sgl = [];
elapsed_sgl = [];
pimat_sums_sgl = [];
pimats_sgl = [];
zeros_groups_sgl = [];


ridge_equivalent = [];


% CVX to find optimal value for Lasso
for i=1:iterations
    i
    
    % CVX to find optimal value for Group Lasso
    tic
    cvx_begin %quiet
        variable pimat_sgl(totalassets,n)

        minimize(sum_square_pos(norm(I_train - sum(R_train'*pimat_sgl,2),2)) + l1_groups*norms(pimat_sgl,2,1)*sqr_group_sizes + l2_features*sum(abs(pimat_sgl(:))) )

        subject to
    %         sum(pimat_sgl(:)) == 1
            pimat_sgl >= 0 
            pimat_sgl(not_group_idx') == 0
    cvx_end
    elapsed=toc

    pimat_sgl        = full(pimat_sgl);
    pimat_norm_sgl   = full(pimat_sgl/sum(pimat_sgl(:)));
    
    
    zeros_groups = [];
    for i=1:n
        zeros_groups = [ zeros_groups sum(pimat_norm_sgl(group_idx(i,:),i)<0.0001)' ]; 
    end

    error_sgl       = sum(abs(I_test - sum(R_test'*pimat_norm_sgl,2)))
    zero_sgl        = sum(sum(pimat_norm_sgl,2)<0.0001)
    
    errors_sgl      = [ errors_sgl, error_sgl ];
    zeros_sgl       = [ zeros_sgl, zero_sgl ];
    elapsed_sgl     = [ elapsed_sgl,  elapsed];
    pimat_sums_sgl  = [ pimat_sums_sgl;  sum(pimat_norm_sgl,1)];
    pimats_sgl      = [ pimats_sgl sum(pimat_norm_sgl,2) ];
    zeros_groups_sgl = [zeros_groups_sgl; zeros_groups];
    
    [ignore, ridgeq] = rf.ridge(I,R(sum(pimat_norm_sgl,2)>.0001,:),T,totalassets-zero_sgl);
    ridge_equivalent = [ ridge_equivalent, ridgeq ]
    
    [elapsed, zero_sgl, lambda, error_sgl]
    
    l1_sgl     = [ l1_sgl, l1_groups ];
    l2_sgl     = [ l2_sgl, l2_features ];
%     l2_features = l2_features + 0.05;
    
%     l1_groups = l1_groups + 0.00005;
%     l1_groups = l1_groups + 0.5;

    l2_features = l2_features + .1;
%     l1_groups = l2_features/divcoef;
    
%     l1_groups       = l1_groups + 30;
%     l2_features     = l2_features + 30;
    
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





%% Load all values

% all_err = [];
% all_zeros = [];
% 
% omegas = [500,1000,5000,10000];
% colours = ['b';'r';'c';'m'];
% filename='ftse_sector';
% 
% h1=figure(1)
% for n_iter=1:4
%     load(strcat('/Users/axsauze/IdeaProjects/matlab/SPAMS/savedworkspaces/correctsparsegrouplasso/proportionate/',int2str(omegas(n_iter)),'/',filename,'.mat'));
%     all_err = [all_err; errors_sgl];
%     all_zeros = [all_zeros; zeros_sgl];
%     plot(zeros_sgl',errors_sgl', colours(n_iter));
%     hold on;
% end
% 
% % Loading Lasso
% load('/Users/axsauze/IdeaProjects/matlab/SPAMS/savedworkspaces/correctlasso/ftse_lasso_sector.mat')
% all_err = [all_err; errors_sgl];
% all_zeros = [all_zeros; zeros_sgl];
% plot(zeros_sgl',errors_sgl', 'k');
% hold on;
% 
% % Loading Group Lasso
% load(strcat('/Users/axsauze/IdeaProjects/matlab/SPAMS/savedworkspaces/correctgrouplasso/grouplasso_',filename,'.mat'))
% all_err = [all_err; errors_sgl];
% all_zeros = [all_zeros; zeros_sgl];
% plot(zeros_sgl',errors_sgl', 'b --');
% hold on;
% 
% %  Loading Feature Selection Ridge
% load('/Users/axsauze/IdeaProjects/matlab/SPAMS/savedworkspaces/featureselection/ftse_feature_selection.mat')
% ridge=flip(all_err(:,3));
% ridge=ridge(1:80,:);
% zeros_ridge=linspace(1,80,80);
% plot(zeros_ridge,ridge,'k .');
% % hold on;
% 
% legend('500','1000','5000','10000','0','Infinite', 'Ridge')
% savegraph(h1,'# Zero Weights','Tracking Error')


