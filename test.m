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
% randgroups = ceil(linspace(0.1,9,99));
% groups = randgroups(randperm(length(randgroups)))';


% Lambdas (Inputs)
l_l     = 0.01;
l_gl    = 0.7;
l_sgl1  = 0.5;
l_sgl2  = 0.05;

omega1 = 500;
omega2 = 5000;


% Get all anonymousregression functions
rs=reg_funcs;

% Selecting only a subset 100 days to run this test:
T = 100;
R = R(:,1:T);
I = I(1:T);
wsize   = 60;


% Lasso
tes_l       = [];
pimats_l    = [];
zeros_l     = [];
gsums_l     = [];
% Group Lasso
tes_gl       = [];
pimats_gl    = [];
zeros_gl     = [];
gsums_gl     = [];
% Sparse Group Lasso Omega=300
tes_sgl300       = [];
pimats_sgl300    = [];
zeros_sgl300     = [];
gsums_sgl300     = [];
% Sparse Group Lasso Omega=3000
tes_sgl3000       = [];
pimats_sgl3000    = [];
zeros_sgl3000     = [];
gsums_sgl3000     = [];


for i=1:T-wsize
    fprintf('Running #%d of %d...\n',i,T-wsize);
    % Setting up test-train data for window
    R_train = R(:,i:i+wsize-1);
    R_test  = R(:,i+wsize);
    I_train = I(i:i+wsize-1);
    I_test  = I(i+wsize);
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
    % Lasso
    fprintf('\t> Starting Lasso.\n');
    pimat_l     = rs.sglasso(I_train,R_train,groups,0,l_l);
    
    gsum        = sum(pimat_l);
    pimat       = pimat_l/gsum;
    error       = norm(I_test - R_test'*pimat,1);
    nzeros      = sum(pimat<0.0001);
    
    tes_l       = [ tes_l, error ];
    pimats_l    = [ pimats_l, pimat ];
    zeros_l     = [ zeros_l, nzeros ];
    gsums_l     = [ gsums_l, gsum ];
    fprintf('\t\t* Done Lasso. Err: %f, Zeros: %f.\n',error,nzeros);
    
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
%     % Group Lasso
%     fprintf('\t> Starting Group Lasso.\n');
%     pimat_gl     = rs.sglasso(I_train,R_train,groups,l_gl,0);
%     
%     gsum        = sum(pimat_gl);
%     pimat       = pimat_gl/gsum;
%     error       = norm(I_test - R_test'*pimat,1);
%     nzeros      = sum(pimat<0.0005);
%     
%     tes_gl       = [ tes_gl, error ];
%     pimats_gl    = [ pimats_gl, pimat ];
%     zeros_gl     = [ zeros_gl, nzeros ];
%     gsums_gl     = [ gsums_gl, gsum ];
%     fprintf('\t\t* Done Group Lasso. Err: %f, Zeros: %f.\n',error,nzeros);
%     
%     
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
%     % Spares Group Lasso, Omega=300
%     fprintf('\t> Starting Sparse Group Lasso. Omega=500\n');
%     pimat_sgl300     = rs.sglasso(I_train,R_train,groups,l_sgl1/omega1,l_sgl1);
%     
%     gsum        = sum(pimat_sgl300);
%     pimat       = pimat_sgl300/gsum;
%     error       = norm(I_test - R_test'*pimat,1);
%     nzeros      = sum(pimat<0.0005);
%     
%     tes_sgl300       = [ tes_sgl300, error ];
%     pimats_sgl300    = [ pimats_sgl300, pimat ];
%     zeros_sgl300     = [ zeros_sgl300, nzeros ];
%     gsums_sgl300     = [ gsums_sgl300, gsum ];
%     fprintf('\t\t* Done Sparse Group Lasso. Err: %f, Zeros: %f.\n',error,nzeros);
    
    
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
    % Sparse Group Lasso, Omega=3000
    fprintf('\t> Starting Sparse Group Lasso. Omega=3000\n');
    pimat_sgl3000     = rs.sglasso(I_train,R_train,groups,l_sgl2/omega2,l_sgl2);
    
    gsum        = sum(pimat_sgl3000);
    pimat       = pimat_sgl3000/gsum;
    error       = norm(I_test - R_test'*pimat,1);
    nzeros      = sum(pimat<0.0005);
    
    tes_sgl3000       = [ tes_sgl3000, error ];
    pimats_sgl3000    = [ pimats_sgl3000, pimat ];
    zeros_sgl3000     = [ zeros_sgl3000, nzeros ];
    gsums_sgl3000     = [ gsums_sgl3000, gsum ];
    fprintf('\t\t* Done Sparse Group Lasso. Err: %f, Zeros: %f.\n',error,nzeros);
end


plot([tes_l', tes_gl', tes_sgl300', tes_sgl3000'])
sum([zeros_l'<=zeros_sgl3000' tes_sgl3000' <= tes_l'  ])
