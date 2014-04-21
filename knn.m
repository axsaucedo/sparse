function [index, groups] = knn(data, k)

l = sqrt(k);

% Creating correlation matrix
corr_mat = corrcoef(data');
corr_mat = 1-corr_mat;
Alower = tril(corr_mat, -1);
Aupper = triu(corr_mat,  1);
G_S = Alower(2:end, 1:end) + Aupper(1:end-1, 1:end);
s = size(G_S,2);

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
    results{curr} = [ results{curr} ; data(i,:)];
end

groups = results;
index = group_index;
