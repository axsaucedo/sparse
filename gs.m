randn('state',23432);
rand('state',3454);

fR = importdata('FPSESecuritiesData');
fR = fR(1:270,:);

T = size(fR, 1)-1;
totalassets = size(fR,2);

% Calculating returns for all individual assets
R = [];
for i = 1:totalassets
    
    returns = [];
    
    for curr = 2:(T+1)
        returns = [returns, fR(curr, i) / fR(curr-1, i)];
    end
    
    R = [R; returns];
end

% Calculating returns for Index
I = mean(R)';


%% MAIN CVX PROGRAM %%
Beta = 0.8;
k = 1;
C = 1/sqrt(totalassets) + k*(1-1/sqrt(totalassets))/(10*sqrt(totalassets));
divcoef = 1 / (T*(1-Beta));

% Calculating n groups with knn + kmeans
n = 5;
[index, groups] = knn(R, n);
combs = [];
size_combs = [];

% for i = 1:n
%     
% end

% Creating indexes for combinations
for i = 1:n
    comb = combnk(1:n, i);
    combs = [combs; [comb , zeros(size(comb,1),n-size(comb,2)) ]];
    size_combs = [size_combs; ones(size(comb, 1),1)*i];
end

% Variables to keep track of results
totalcombinations = size(combs,1)

currRs = cell(1, totalcombinations);
pimats_n = cell(1, totalcombinations);
te_n = [];


% Loop for all the possible combinations
for i = 1:totalcombinations
    i
    
    currR = [];
    
    % Generate our dataset for all possible combinations
    for j = 1:size_combs(i)
        currR = [currR; groups{combs(i,j)}];
    end
    
    curr_totalassets = size(currR,1);
    
    C = 1/sqrt(curr) + k*(1-1/sqrt(curr))/(10*sqrt(curr));
    
    % CVX to find optimal value for NCCVAR
    cvx_begin quiet
        variable z_n(T)
        variable Alpha_n
        variable pimat_n(curr_totalassets)
        minimize( Alpha_n + divcoef * sum(z_n) )
        subject to
            z_n >= 0
            z_n - abs(I - transpose(currR)*pimat_n) + Alpha_n >= 0
            
            pimat_n >= 0
            sum(pimat_n) == 1
            
%             % This constraint does not hold
%             norm(pimat_n) <= C

    cvx_end
    
    % Calculating Tracking Error
    currRs{i} = currR;
    pimats_n{i} = pimat_n;
    
    Ret_n = abs(I - transpose(currR) * pimat_n);
    te_n = [te_n; sum(Ret_n)];
        
end