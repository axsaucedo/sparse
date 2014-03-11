function [groups] = knn(k)

randn('state',23432);
rand('state',3454);
    
% Getting Data 

%%%%%%%%%%%%%%%%%%% From File %%%%%%%%%%%%%%%%%%
fI = importdata('../takeda_cvx/FPSEOnlyData');
fR = importdata('../takeda_cvx/FPSESecuritiesData');
% fI = fI(1:100,1);
% fR = fR(1:100,1:50
% fR = fR(:,10:20);
T = size(fI, 1) - 1;
totalassets = size(fR,2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Returns and Proportions %%%%%%%%%%%%

% Generating Portfolio proportion random matrix with rows adding to 1
% pimat = rand(totalassets, 1); 
% pimat = pimat / sum(pimat, 1); % Each value over the total of the sum of
% columns

% Calculating returns for all individual assets
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

l = sqrt(k);

% Creating correlation matrix
corr_mat = corrcoef(R');
corr_mat = 1-corr_mat;
Alower = tril(corr_mat, -1);
Aupper = triu(corr_mat,  1);
G_S = Alower(2:end, 1:end) + Aupper(1:end-1, 1:end);
s = size(G_S,2);
s

sorted = [];
index = [];
% Sorting and keeping track of indexes
for i = 1:s
    
    [std, idx] = sort(G_S(:,i), 1);
   
    sorted = [sorted, std];
    index = [index, idx];
end

% Now we only take the top k neighbors
k_sorted = sorted(1:k,:);
k_index = index(1:k,:);

% We build the Adjacency Matrix A and Degree Matrix D
W_S = zeros(s);
D_S = zeros(s);
for i = 1:s
    
    for j = 1:k
        
        W_S(i, k_index(j,i)) = exp(-k_sorted(j,i)/l);
    end
    
    D_S(i,i) = sum(W_S(i,:));
end

% Computing Laplacian Matrix
size(D_S)
size(W_S)
L_S = D_S - W_S;

% Get eigen values of L
evals_s = eig(L_S);
[ eVecs_s, diagEVals_s ] = eig(L_S);

W = W_S + W_S';
D = diag(sum(W));

% Computing modified Laplacian Matrix
L = D - W;

% Eigen values of Laplacian
evals = eig(L);
[ eVecs, diagEVals ] = eig(L);


three_eVecs = eVecs(:, 1:k);

% Find clusters
group_index = kmeans(three_eVecs, k);


results = cell(1,k); 

for i = 1:s
    curr = group_index(i);
    results{curr} = [ results{curr} ; R(i,:)];
end

groups = results;

% plot(groups);
