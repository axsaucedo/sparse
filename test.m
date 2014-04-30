function test(which, iters)
% TEST - This file provides a test for the sparse group
% lasso to show some basic functionality of this model.
% PARAMETERS:
%   omega   -   The Omega that represents the ratio for l1 and l2 (Low
%                   omega would induce group lasso. High omega would induce
%                   simple lasso.
%   l       -   The lambda used for the tests
%   iters   -   The number of iterations on each test (default=3)
%   which   -   an int that specifies which test to run.
%           Values for which:
%               0   -   Show the transition from the Group Lasso model 
%                           to the Simple Lasso model.
%               1   -   Test Group Lasso Behaviour (i.e. l1=0, l2=n)
%               2   -   Test Lasso Behaviour (i.e l1=0, l2=n)
%               3   -   Test Sparse Group Lasso Behaviour (i.e. l1=n, l2=n)
%               4   -   Test the least squares Behaviour (i.e. l1=0, l2=0)



disp(sprintf('Loading FTSE100 stock data\n'));
% The anonymous function ds.ftse100() returns a matrix of FTSE100 
% returns R, the Market Index of these returns I = mean(R,2), the number of
% assets n, the time window T and the groupings contained in the variable
% groups
ds = datasets;
[ I, R, n, T, groups ] = ds.ftse100();

if nargin<2
    iters = 5;
end

limit = 20;

l = limit/iters;

switch which
    case 0      
        disp('Running transition from Group Lasso to Lasso');
        disp('Test 1 - Group Lasso Behaviour');
        test_sparsegrouplasso(I, R, groups, iters, l, 0);
        disp('Test 2 - Sparse Group Lasso Behaviour. Omega=300');
        test_sparsegrouplasso(I, R, groups, iters, l/300, l);
        disp('Test 2 - Sparse Group Lasso Behaviour. Omega=3000');
        test_sparsegrouplasso(I, R, groups, iters, l/3000, l);
        disp('Test 3 - Lasso Behaviour');
        test_sparsegrouplasso(I, R, groups, iters, 0, l);
    case 1
        disp('Test 1 - Group Lasso Behaviour');
        test_sparsegrouplasso(I, R, groups, iters, l, 0);
    case 2
        disp('Test 2 - Sparse Group Lasso Behaviour');
        test_sparsegrouplasso(I, R, groups, iters, l/300, l);
    case 3
        disp('Test 3 - Lasso Behaviour');
        test_sparsegrouplasso(I, R, groups, iters, 0, l);
    case 4
        disp('Test 4 - Least Squares Behaviour');
        test_sparsegrouplasso(I, R, groups, iters, 0, 0);
        
        
        
    otherwise
        error('Please specify a valid test')
end

end


