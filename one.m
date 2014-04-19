% 
% Usage:  [V [val_regularizer]]=mexProximalFlat(U,param);
%
% Name: mexProximalFlat
%
% Description: mexProximalFlat computes proximal operators. Depending
%         on the value of param.regul, it computes 
%
%         Given an input matrix U=[u^1,\ldots,u^n], it computes a matrix 
%         V=[v^1,\ldots,v^n] such that
%         if one chooses a regularization functions on vectors, it computes
%         for each column u of U, a column v of V solving
%         if param.regul='l0'
%             argmin 0.5||u-v||_2^2 + lambda||v||_0
%         if param.regul='l1'
%             argmin 0.5||u-v||_2^2 + lambda||v||_1
%         if param.regul='l2'
%             argmin 0.5||u-v||_2^2 + 0.5lambda||v||_2^2
%         if param.regul='elastic-net'
%             argmin 0.5||u-v||_2^2 + lambda||v||_1 + lambda_2||v||_2^2
%         if param.regul='fused-lasso'
%             argmin 0.5||u-v||_2^2 + lambda FL(v) + ...
%                               ...  lambda_2||v||_1 + lambda_3||v||_2^2
%         if param.regul='linf'
%             argmin 0.5||u-v||_2^2 + lambda||v||_inf
%         if param.regul='l1-constraint'
%             argmin 0.5||u-v||_2^2 s.t. ||v||_1 <= lambda
%         if param.regul='l2-not-squared'
%             argmin 0.5||u-v||_2^2 + lambda||v||_2
%         if param.regul='group-lasso-l2'  
%             argmin 0.5||u-v||_2^2 + lambda sum_g ||v_g||_2 
%             where the groups are either defined by param.groups or by param.size_group,
%         if param.regul='group-lasso-linf'
%             argmin 0.5||u-v||_2^2 + lambda sum_g ||v_g||_inf
%         if param.regul='sparse-group-lasso-l2'  
%             argmin 0.5||u-v||_2^2 + lambda sum_g ||v_g||_2 + lambda_2 ||v||_1
%             where the groups are either defined by param.groups or by param.size_group,
%         if param.regul='sparse-group-lasso-linf'
%             argmin 0.5||u-v||_2^2 + lambda sum_g ||v_g||_inf + lambda_2 ||v||_1
%         if param.regul='trace-norm-vec' 
%             argmin 0.5||u-v||_2^2 + lambda ||mat(v)||_* 
%            where mat(v) has param.size_group rows
%
%         if one chooses a regularization function on matrices
%         if param.regul='l1l2',  V= 
%             argmin 0.5||U-V||_F^2 + lambda||V||_{1/2}
%         if param.regul='l1linf',  V= 
%             argmin 0.5||U-V||_F^2 + lambda||V||_{1/inf}
%         if param.regul='l1l2+l1',  V= 
%             argmin 0.5||U-V||_F^2 + lambda||V||_{1/2} + lambda_2||V||_{1/1}
%         if param.regul='l1linf+l1',  V= 
%             argmin 0.5||U-V||_F^2 + lambda||V||_{1/inf} + lambda_2||V||_{1/1}
%         if param.regul='l1linf+row-column',  V= 
%             argmin 0.5||U-V||_F^2 + lambda||V||_{1/inf} + lambda_2||V'||_{1/inf}
%         if param.regul='trace-norm',  V= 
%             argmin 0.5||U-V||_F^2 + lambda||V||_*
%         if param.regul='rank',  V= 
%             argmin 0.5||U-V||_F^2 + lambda rank(V)
%         if param.regul='none',  V= 
%             argmin 0.5||U-V||_F^2 
%         
%         for all these regularizations, it is possible to enforce non-negativity constraints
%         with the option param.pos, and to prevent the last row of U to be regularized, with
%         the option param.intercept
%
% Inputs: U:  double m x n matrix   (input signals)
%               m is the signal size
%         param: struct
%               param.lambda  (regularization parameter)
%               param.regul (choice of regularization, see above)
%               param.lambda2  (optional, regularization parameter)
%               param.lambda3  (optional, regularization parameter)
%               param.verbose (optional, verbosity level, false by default)
%               param.intercept (optional, last row of U is not regularized,
%                 false by default)
%               param.transpose (optional, transpose the matrix in the regularization function)
%               param.size_group (optional, for regularization functions assuming a group
%                 structure). It is a scalar. When param.groups is not specified, it assumes
%                 that the groups are the sets of consecutive elements of size param.size_group
%               param.groups (int32, optional, for regularization functions assuming a group
%                 structure. It is an int32 vector of size m containing the group indices of the
%                 variables (first group is 1).
%               param.pos (optional, adds positivity constraints on the
%                 coefficients, false by default)
%               param.numThreads (optional, number of threads for exploiting
%                 multi-core / multi-cpus. By default, it takes the value -1,
%                 which automatically selects all the available CPUs/cores).
%
% Output: V: double m x n matrix (output coefficients)
%         val_regularizer: double 1 x n vector (value of the regularization
%         term at the optimum).

clear all;
format compact;
randn('seed',0);
param.numThreads=-1; % all cores (-1 by default)
param.verbose=true;   % verbosity, false by default
param.lambda=0.05; % regularization parameter
param.it0=10;      % frequency for duality gap computations
param.max_it=200; % maximum number of iterations
param.L0=0.1;
param.tol=1e-3;
param.intercept=false;
param.pos=true;

param.loss='square';
param.compute_gram=true;


X=randn(100,200);
X=X-repmat(mean(X),[size(X,1) 1]);
X=mexNormalize(X);
Y=randn(100,1);
Y=Y-repmat(mean(Y),[size(Y,1) 1]);
Y=mexNormalize(Y);
W0=zeros(size(X,2),size(Y,2));

groups=int32(randi(5,1,size(X,2)));



% Group Lasso

% fprintf('\nFISTA + Group Lasso L2\n');
% param.regul='group-lasso-l2';
% param2=param;
% param2.size_group=2;  % all the groups are of size 2
% tic
% [W_lasso_fix optim_info]=mexFistaFlat(Y,X,W0,param2);
% t=toc;
% fprintf('mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n',mean(optim_info(1,:)),mean(optim_info(3,:)),t,mean(optim_info(4,:)));


% fprintf('\nFISTA + Group Lasso L2 with variable size of groups \n');
% param.regul='group-lasso-l2';
% param2=param;
% param2.groups=groups;  % all the groups are of size 2
% param2.lambda=10*param2.lambda;
% % tpm_param=param2
% tic
% [W_lasso_var optim_info]=mexFistaFlat(Y,X,W0,param2);
% t=toc;
% fprintf('mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n',mean(optim_info(1,:)),mean(optim_info(3,:)),t,mean(optim_info(4,:)));


% fprintf('\nFISTA + Group Lasso L2 with variable size of groups \n');
% param.regul='sparse-group-lasso-l2';
% param2=param;
% param2.size_group=5;  % all the groups are of size 2
% % param2.lambda=10*param2.lambda;
% param2.lambda=param.lambda/2;
% param2.lambda2=param2.lambda/2;
% % tpm_param=param2
% tic
% [W_sparse_fix optim_info]=mexFistaFlat(Y,X,W0,param2);
% t=toc;
% fprintf('mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n',mean(optim_info(1,:)),mean(optim_info(3,:)),t,mean(optim_info(4,:)));



fprintf('\nFISTA + Group Lasso L2 with variable size of groups \n');
param.regul='sparse-group-lasso-l2';
param2=param;
param2.groups=groups;  % all the groups are of size 2
param2.lambda=param.lambda/2;
param2.lambda2=param.lambda/2;
% tpm_param=param2
tic
[W_sparse_var optim_info]=mexFistaFlat(Y,X,W0,param2);
t=toc;
fprintf('mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n',mean(optim_info(1,:)),mean(optim_info(3,:)),t,mean(optim_info(4,:)));



clear all;

randn('seed',0);
fprintf('test mexLasso\n');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decomposition of a large number of signals
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data are generated
X=randn(100,100000);
X=X./repmat(sqrt(sum(X.^2)),[size(X,1) 1]);
D=randn(100,200);
D=D./repmat(sqrt(sum(D.^2)),[size(D,1) 1]);

% parameter of the optimization procedure are chosen
%param.L=20; % not more than 20 non-zeros coefficients (default: min(size(D,1),size(D,2)))
param.lambda=0.15; % not more than 20 non-zeros coefficients
param.numThreads=-1; % number of processors/cores to use; the default choice is -1
                     % and uses all the cores of the machine
param.mode=2;        % penalized formulation

tic
alpha=mexLasso(X,D,param);
t=toc
fprintf('%f signals processed per second\n',size(X,2)/t);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Regularization path of a single signal 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X=randn(64,1);
D=randn(64,10);
D=D./repmat(sqrt(sum(D.^2)),[size(D,1) 1]);
param.lambda=0;
[alpha path]=mexLasso(X,D,param);
