% Setting random state to get constant answers
randn('state',23432);
rand('state',3454);


%%%%%%%%%%%%%%%%%%% From File %%%%%%%%%%%%%%%%%%
fI = importdata('FPSEOnlyData');
fR = importdata('FPSESecuritiesData');
% fI = fI(1:100,1);
% fR = fR(1:100,1:50
% fR = fR(:,10:20);
T = size(fI, 1) - 1;
totalassets = size(fR,2);

R = [];
for i = 1:totalassets
    
    returns = [];
    
    for curr = 2:(T+1)
        returns = [returns, fR(curr, i) / fR(curr-1, i)];
    end
    
    R = [R; returns];
end

I = [];
for curr = 2:(T+1)
    I = [I; fI(curr) / fI(curr-1)];
end

% Calculating returns for Index
I = mean(R)';



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Random modelled data %%%%%%%%%%%%%%
totalassets = 200;
T = 1000;
assets = [];
for i=1:totalassets
    volatility = .1 + rand()/5-.1;
    old_price = 100 - rand()*50-25;

    price = [];
    for j=1:T
        rnd = rand(); % generate number, 0 <= x < 1.0
        change_percent = 2 * volatility * rnd;
        if (change_percent > volatility)
            change_percent = change_percent-(2 * volatility);
        end
        change_amount = old_price * change_percent;
        new_price = old_price + change_amount;
        price = [price; old_price/new_price];
    end
    assets = [ assets, price ];
end
R = assets';
I = R'*(ones(totalassets,1)/totalassets);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Tracking error of Portfolio Subset %%%%%%%
%Initialize

Beta = 0.8;
k = 1;
C = 1/sqrt(totalassets) + k*(1-1/sqrt(totalassets))/(10*sqrt(totalassets));
divcoef = 1 / (T*(1-Beta));

delta = 1.0;

cvx_n       = totalassets;
subset_n    = cvx_n;
cvx_R       = R(1:cvx_n,:);
cvx_I       = I;


% % % % % % % % % % % Variables % % % % % % % % % % 
% This array will contain the subset_n of stocks chosen from the index
%   Avaiable : are the stocks that haven't been added to the set
%   Values   : are the values obtained when computing our regression function
%              after adding our selected stock
rf                = reg_funcs;


% % % % % % % ABS % % % % % % % % %
pimats_a          = zeros(cvx_n);
chosen_a          = logical(zeros(1,cvx_n));
chosen_order_a    = [];
chosen_error_a    = [];
elapsed_a         = [];
available_a   = linspace(1,cvx_n,cvx_n);
% % % % % LEAST SQUARES % % % % % %
pimats_q          = zeros(cvx_n);
chosen_q          = logical(zeros(1,cvx_n));
chosen_order_q    = [];
chosen_error_q    = [];
elapsed_q         = [];
available_q   = linspace(1,cvx_n,cvx_n);
% % % % % % % NCCVaR % % % % % % % %
pimats_r          = zeros(cvx_n);
chosen_r          = logical(zeros(1,cvx_n));
chosen_order_r    = [];
chosen_error_r    = [];
elapsed_r         = [];
available_r   = linspace(1,cvx_n,cvx_n);
% % % % % % % CVaR % % % % % % % % %
pimats_c          = zeros(cvx_n);
chosen_c          = logical(zeros(1,cvx_n));
chosen_order_c    = [];
chosen_error_c    = [];
elapsed_c         = [];
available_c   = linspace(1,cvx_n,cvx_n);


for i=1:subset_n
    i
    
    % % % % % % % % % % % % % % % %     
    % % % % % % % ABS
    tic
    [ optimal_available_a, optimal_error_a, pimat_a ] = next_optimal(rf.abs, cvx_I, cvx_R, available_a, chosen_a)
    elapsed=toc
    
    optimal_index_a             = available_a(optimal_available_a);
    chosen_a(optimal_index_a)   = true;
    chosen_order_a              = [chosen_order_a optimal_index_a];
    chosen_error_a              = [chosen_error_a, optimal_error_a];
    elapsed_a                   = [ elapsed_a, elapsed ];
    
    pimats_a(chosen_a,i)= pimat_a;
    
    % Removing the index from the available stock subset
    available_a(optimal_available_a) = [];
    
    
    
    % % % % % % % % % % % % % % % %     
    % % % % % % % LEAST SQUARES
    tic
    [ optimal_available_q, optimal_error_q, pimat_q ] = next_optimal(rf.squares, cvx_I, cvx_R, available_q, chosen_q)
    elapsed=toc
    
    optimal_index_q             = available_q(optimal_available_q);
    chosen_q(optimal_index_q)   = true;
    chosen_order_q              = [chosen_order_q optimal_index_q];
    chosen_error_q              = [chosen_error_q, optimal_error_q];
    elapsed_q                   = [ elapsed_q, elapsed ];
    
    pimats_q(chosen_q,i)= pimat_q;
    
    % Removing the index from the available stock subset
    available_q(optimal_available_q) = [];
    
    
    
    % % % % % % % % % % % % % % % %     
    % % % % % % % Ridge Regression
    tic
    [ optimal_available_r, optimal_error_r, pimat_r ] = next_optimal(rf.ridge, cvx_I, cvx_R, available_r, chosen_r)
    elapsed=toc
    
    optimal_index_r             = available_r(optimal_available_r);
    chosen_r(optimal_index_r)   = true;
    chosen_order_r              = [chosen_order_r optimal_index_r];
    chosen_error_r              = [chosen_error_r, optimal_error_r];
    elapsed_r                   = [ elapsed_r, elapsed ];
    
    pimats_r(chosen_r,i)= pimat_r;
    
    % Removing the index from the available stock subset
    available_r(optimal_available_r) = [];
    
    
    
    % % % % % % % % % % % % % % % %     
    % % % % % % % CVaR
    tic
    [ optimal_available_c, optimal_error_c, pimat_c ] = next_optimal(rf.cvar, cvx_I, cvx_R, available_c, chosen_c)
    elapsed=toc
    
    optimal_index_c             = available_c(optimal_available_c);
    chosen_c(optimal_index_c)   = true;
    chosen_order_c              = [chosen_order_c optimal_index_c];
    chosen_error_c              = [chosen_error_c, optimal_error_c];
    elapsed_c                   = [ elapsed_c, elapsed ];
    
    pimats_c(chosen_c,i)= pimat_c;
    
    % Removing the index from the available stock subset
    available_c(optimal_available_c) = [];
    
end

all_elapsed   = [elapsed_a' elapsed_q' elapsed_r' elapsed_c']
all_err     = [chosen_error_a' chosen_error_q' chosen_error_r' chosen_error_c']
all_order   = [chosen_order_a' chosen_order_q' chosen_order_r' chosen_order_c']

l = linspace(1,134,134);

total_same = zeros(1,134);
for i=1:134
    u=union(union(union(all_order(1:i,1),all_order(1:i,2)),all_order(1:i,3)),all_order(1:i,4))
    total_same(i) = size(u,1)
end
total_same_pct =  l./total_same;

% path='Graphs/feature_selection/';
% 
h1 = figure(1);

idx = linspace(1,100,100);
idx = idx(20:70);
plot(idx, all_err(idx,1));
hold on
plot(idx, all_err(idx,2), 'red')
hold on
plot(idx, all_err(idx,3), 'black')
hold on
plot(idx, all_err(idx,4), 'magenta')
% 
% lgnd = repmat([' '],10)
% lgnd(1,1:3) = 'Abs';
% lgnd(2,1:7) = 'Squares';
% lgnd(3,1:5) = 'Ridge';
% lgnd(4,1:4) = 'CVaR';
% 
% savegraph(h1,'Non-Zero Values','Error',lgnd,fullfile(path,'ftse100_feature_seleciton_zoom'))


