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
%%%%%% Tracking error of Portfolio Subset %%%%%%%
%Initialize

Beta = 0.8;
k = 1;
C = 1/sqrt(totalassets) + k*(1-1/sqrt(totalassets))/(10*sqrt(totalassets));
divcoef = 1 / (T*(1-Beta));

delta = 1.0;

cvx_n       = 10;
subset_n    = cvx_n;
cvx_R       = R(1:cvx_n,:);
cvx_I       = I;


% % % % % % % % % % % Variables % % % % % % % % % % 
% This array will contain the subset_n of stocks chosen from the index
%   Avaiable : are the stocks that haven't been added to the set
%   Values   : are the values obtained when computing our regression function
%              after adding our selected stock
rf                = reg_funcs;

pimats_a          = zeros(cvx_n);
chosen_a          = logical(zeros(1,cvx_n));
chosen_order_a    = [];
chosen_error_a    = [];
elapsed_a         = [];

available_a   = linspace(1,cvx_n,cvx_n);

for i=1:subset_n
    
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
    
end


