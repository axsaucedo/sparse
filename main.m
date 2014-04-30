function main(which, iters)
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

% Loading all required files to path for this session
addpath('functions/')
addpath('data/')
addpath('helperfuncs/')
addpath('tests/')


disp(sprintf('Loading FTSE100 stock data\n'));
% The anonymous function ds.ftse100() returns a matrix of FTSE100 
% returns R, the Market Index of these returns I = mean(R,2), the number of
% assets n, the time window T and the groupings contained in the variable
% groups
ds = datasets;
[ I, R, n, T, groups ] = ds.ftse_yahoo();

if nargin<2
    iters = 5;
end

limit = 20;

l = limit/iters;

switch which
    case 0      
        disp('Test 1 - Transition from Group Lasso to Lasso');
        disp('Testing Group Lasso Behaviour');
        test_sparsegrouplasso(I, R, groups, iters, l, 0);
        disp('Testing Sparse Group Lasso Behaviour. Omega=300');
        test_sparsegrouplasso(I, R, groups, iters, l/300, l);
        disp('Testing Sparse Group Lasso Behaviour. Omega=3000');
        test_sparsegrouplasso(I, R, groups, iters, l/3000, l);
        disp('Testing Lasso Behaviour');
        test_sparsegrouplasso(I, R, groups, iters, 0, l);
    case 1
        disp('Test 2 - Feature Selection')
        test_feature_selection(I, R, n, T);
    case 2
        disp('Test 3 - Group Selection')
        test_group_selection(I, R, n, T, groups);
        
        
        
    otherwise
        error('Please specify a valid test')
end

end


