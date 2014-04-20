% SFE - Single Feature Equations
% This file contains all the equations with single features including:
% CVX:
% > NCCVAR Minimization
% > CVAR Minimization
% > Least Squares Minimization
% > Lasso Minimization
% SPAMS:
% > Lasso Minimization

clear all;

% Setting random state to get constant answers
randn('state',23432);
rand('state',3454);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Variables that will be used throughout the
% Program can be either set up randomly, obtained
% through the Yahoo Finance API, or read from a 
% file. Please comment/uncomment the following code
% as required to enable/disable these features


%%%%%%%%%%%%%%%% SET UP RANDOMLY %%%%%%%%%%%%%%%%
% roundcoeficient = 99;
% T = 100; % Total number of samples per hist set
% totalassets = 5;
% % Matrix containing historical return for all current assets
% fR = floor(roundcoeficient * rand(T,totalassets)+1);
% % Matrix containing historical returns for index assets
% fI = floor(roundcoeficient * rand(T,totalassets)+1);
% T = T-1;


%%%%%%%%%%%%%%%%%%% Yahoo Finance %%%%%%%%%%%%%%%
% qtimeout = 3;
% y = yahoo;
% % RIGHT NOW WE HAVE AN FPSE 99 BECAUSE 'CCH.L' DOES NOT FETCH!
% % securities = {'AAL.L','ABF.L','ADM.L'}%,'ADN.L','AGK.L','AMEC.L','ANTO.L','ARM.L','AV.L','AZN.L','BA.L','BAB.L','BARC.L','BATS.L','BG.L','BLND.L','BLT.L','BNZL.L','BP.L','BRBY.L','BSY.L','BT-A.L','CCL.L','CNA.L','CPG.L','CPI.L','CRDA.L','CRH.L','DGE.L','EXPN.L','EZJ.L','FRES.L','GFS.L','GKN.L','GLEN.L','GSK.L','HL.L','HMSO.L','HSBA.L','IAG.L','IHG.L','IMI.L','IMT.L','ITRK.L','ITV.L','JMAT.L','KGF.L','LAND.L','LGEN.L','LLOY.L','LSE.L','MGGT.L','MKS.L','MNDI.L','MRO.L','MRW.L','NG.L','NXT.L','OML.L','PFC.L','PRU.L','PSN.L','PSON.L','RB.L','RBS.L','RDSA.L','RDSB.L','REL.L','REX.L','RIO.L','RR.L','RRS.L','RSA.L','RSL.L','SAB.L','SBRY.L','SDR.L','SGE.L','SHP.L','SL.L','SMIN.L','SN.L','SPD.L','SSE.L','STAN.L','SVT.L','TATE.L','TLW.L','TPK.L','TSCO.L','TT.L','ULVR.L','UU.L','VED.L','VOD.L','WEIR.L','WMH.L','WOS.L','WPP.L'};
% securities = {'BBRY','YHOO'}
% index = '^FTSE';
% totalassets = size(securities,2);
% timeformat = 'd';
% fieldname = 'close';
% dtTimeFrom = '11/1/2011';
% dtTimeUntil = '12/1/2012';
% 
% % Fetching Index prices
% disp(index);
% fetched = fetch(y,index,fieldname,dtTimeFrom,dtTimeUntil,timeformat);
% fI = fetched(:,2);
% 
% fR = [];
% for i = 1:totalassets
%     symbol = securities(i);
%     disp(symbol);
%     
%     fetched = fetch(y,symbol,fieldname,dtTimeFrom,dtTimeUntil,timeformat);
%     
%     fR = [fR, fetched(:, 2)];
% end
% 
% T = min(size(fI, 1),size(fR, 1)) - 1; % Minus one because returns will have one less

%%%%%%%%%%%%%%%%%%% From File %%%%%%%%%%%%%%%%%%
fI = importdata('FPSEOnlyData');
fR = importdata('FPSESecuritiesData');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Hard Coded Test Data %%%%%%%%%%%%%%
% T = 12;
% totalassets = 100;
% 
% I = 0.95+(.1).*rand(T,1);
% R = [];
% for i = 1:totalassets
%     R = [R , I + (0.025 - (0.05).*rand(T,1))];
% end
% R = transpose(R)
% R = fliplr(R)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Random modelled data %%%%%%%%%%%%%%
% totalassets = 200;
% T = 1000;
% assets = [];
% for i=1:totalassets
%     volatility = .1 + rand()/5-.1;
%     old_price = 100 - rand()*50-25;
% 
%     price = [];
%     for j=1:T
%         rnd = rand(); % generate number, 0 <= x < 1.0
%         change_percent = 2 * volatility * rnd;
%         if (change_percent > volatility)
%             change_percent = change_percent-(2 * volatility);
%         end
%         change_amount = old_price * change_percent;
%         new_price = old_price + change_amount;
%         price = [price; old_price/new_price];
%     end
%     assets = [ assets, price ];
% end
% R = assets';
% I = R'*(ones(totalassets,1)/totalassets);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% MAIN CVX PROGRAM %%%%%%%%%%%%%%%
Beta = 0.8;
k = 1;
C = 1/sqrt(totalassets) + k*(1-1/sqrt(totalassets))/(10*sqrt(totalassets));
divcoef = 1 / (T*(1-Beta));

%%%%%%%%%%%%%%%%%% Calculations %%%%%%%%%%%%%%%%%
% C2 = 20
% 
% % CVX to find optimal value for NCCVAR
% cvx_begin % quiet
%     variable z_n(T)
%     variable Alpha_n
%     variable pimat_n(totalassets)
%     minimize( Alpha_n + divcoef * sum(z_n) )
%     subject to
%         z_n >= 0
%         z_n - abs(I - transpose(R)*pimat_n) + Alpha_n >= 0
% %         norm(pimat_n) < C
%         nnz(lt(.0000001,pimat_n)) > C2
%         % QUESTION: If the following constraint is not set, pi is given
%         % negative numbers, should this constrain be added? the constraints
%         % are is specified in the paper when tracking variance is
% %         introduced 
%         pimat_n >= 0
%         sum(pimat_n) == 1
% cvx_end
% 
% % CVX to find optimal value for CVAR
% cvx_begin 
%     variable z_c(T)
%     variable Alpha_c
%     variable pimat_c(totalassets)
%     minimize( Alpha_c + divcoef * sum(z_c) )
%     subject to
%         z_c >= 0
% %         transpose(z_c)*z_c <= power(C,2)
%         transpose(R)*pimat_c + Alpha_c + z_c >= 0
%         norm(pimat_c) < C
%         pimat_c >= 0
%         sum(pimat_c) == 1
% cvx_end
% 
% % CVX to find optimal value for Tracking Error Abs
% cvx_begin 
%     variable z_a(T)
%     variable pimat_a(totalassets)
%     minimize( (1/T) * sum(z_a) )
%     subject to
%         z_a >= 0
% %         transpose(z_c)*z_c <= power(C,2)
%         z_a - abs(I - transpose(R)*pimat_n) >= 0
%         norm(pimat_a) < C
%         pimat_a >= 0
%         sum(pimat_a) == 1
% cvx_end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% TRACKING ERROR %%%%%%%%%%%%%%%
% % In this case our distribution for the ratio 
% % invested on each stock is equal, adding to 1
% balancedpimat = ones(totalassets,1) / totalassets; 
% % Formula to Calculate Tracking Error
% totalSum = 0;
% z_tracking_err = [];
% for i = 1:T
%     iRet = I(i);
%     currR = R(:,i);
%     mRet = transpose(currR) * balancedpimat;
%     currErr = iRet - mRet;
%     absSum = abs(currErr);
%     z_tracking_err = [z_tracking_err ; absSum]
% end
% 
% pimats = [pimat_n, pimat_c];
% allz = [z_n, z_c, z_tracking_err];




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Tracking error of Portfolio Subset %%%%%%%
%Initialize

Beta = 0.8;
k = 1;
C = 1/sqrt(totalassets) + k*(1-1/sqrt(totalassets))/(10*sqrt(totalassets));
divcoef = 1 / (T*(1-Beta));

delta = 1.0;

cvx_n = 30;
cvx_R = R(1:cvx_n,:);
cvx_I = I;


% CVX to find optimal value for NCCVAR
cvx_begin quiet
    variable z_n(T)
    variable Alpha_n
    variable pimat_n(cvx_n)
    minimize( Alpha_n + divcoef * sum(z_n) )
    
    subject to
        z_n >= 0
        z_n - abs(cvx_I - cvx_R'*pimat_n) + Alpha_n >= 0

        pimat_n >= 0
        sum(pimat_n) == 1
cvx_end

% CVAR Minimization
cvx_begin quiet
    variable z_c(T)
    variable Alpha_c
    variable pimat_c(cvx_n)
    minimize( Alpha_c + divcoef * sum(z_c) )
    
    subject to
        z_c >= 0
        z_c - (cvx_I - transpose(cvx_R)*pimat_c) + Alpha_c >= 0
        
        pimat_c >= 0
        sum(pimat_c) == 1
        
cvx_end

% CVX to find optimal value for Tracking Error Abs
cvx_begin quiet
    variable z_a(T)
    variable pimat_a(cvx_n)
    minimize( (1/T) * sum(z_a) )
    subject to
        z_a >= 0
%         transpose(z_c)*z_c <= power(C,2)
        z_a - abs(cvx_I - transpose(cvx_R)*pimat_a) >= 0
        pimat_a >= 0
        sum(pimat_a) == 1
        % QUESTION: Should this norm constraint be in CVAR as well?
%         norm(pimat_a) <= C
cvx_end

% CVX to find least squares
cvx_begin quiet
    variable z_q(T)
    variable pimat_q(cvx_n)
    
    minimize(sum(power(cvx_I - cvx_R'*pimat_q, 2)))
    
    subject to
        sum(pimat_q) <= 1
        pimat_q >= 0
cvx_end


% CVX to find optimal value for Lasso
cvx_begin quiet
    variable z_l(T)
    variable pimat_l(totalassets)
    
    minimize(sum(power(I - transpose(R)*pimat_l, 2)) + delta*sum(pimat_l))
    
    subject to
        sum(pimat_l) <= 1
        pimat_l >= 0
cvx_end

pimats_cvx = [pimat_n pimat_c pimat_a ];

% Calculating Tracking Error
cvx_Ret = [ 
            abs(I - cvx_R' * pimat_n)   ...
            abs(I - cvx_R' * pimat_c)   ...
            abs(I - cvx_R' * pimat_a)   ...
            abs(I - cvx_R' * pimat_q)       ...
            abs(I - R' * pimat_l)
          ]

% Tracking Errors
te = sum(cvx_Ret)




% % SPAMS Variables
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


% X=randn(100,200);
% X=X-repmat(mean(X),[size(X,1) 1]);
% X=mexNormalize(X);
% Y=randn(100,1);
% Y=Y-repmat(mean(Y),[size(Y,1) 1]);
% Y=mexNormalize(Y);
% W0=zeros(size(X,2),size(Y,2));


X=R';
X=mexNormalize(X);
Y=I;
Y=mexNormalize(I);
W0=zeros(size(X,2),size(Y,2));


groups = int32(randi(5,1,size(X,2)));
group_size = 5;



% Group Lasso

% fprintf('\nFISTA + Group Lasso L2\n');
% param.regul='group-lasso-l2';
% param2=param;
% param2.lambda=0.5;
% param2.size_group=group_size;  % all the groups are of size 2
% tic
% [full_pimat_glf optim_info]=mexFistaFlat(Y,X,W0,param2);
% t=toc;
% pimat_glf=full_pimat_glf/sum(full_pimat_glf);
% fprintf('mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n',mean(optim_info(1,:)),mean(optim_info(3,:)),t,mean(optim_info(4,:)));
% 
% 
% fprintf('\nFISTA + Group Lasso L2 with variable size of groups \n');
% param.regul='group-lasso-l2';
% param2=param;
% param2.groups=groups;  % all the groups are of size 2
% param2.lambda=10*param2.lambda;
% % tpm_param=param2
% tic
% [full_pimat_glv optim_info]=mexFistaFlat(Y,X,W0,param2);
% t=toc;
% pimat_glv=full_pimat_glv/sum(full_pimat_glv);
% fprintf('mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n',mean(optim_info(1,:)),mean(optim_info(3,:)),t,mean(optim_info(4,:)));
% 
% 
% fprintf('\nFISTA + Group Lasso L2 with variable size of groups \n');
% param.regul='sparse-group-lasso-l2';
% param2=param;
% param2.size_group=group_size;  % all the groups are of size 2
% % param2.lambda=10*param2.lambda;
% param2.lambda=0.01;
% param2.lambda2=0.01;
% % tpm_param=param2
% tic
% [full_pimat_sglf optim_info]=mexFistaFlat(Y,X,W0,param2);
% t=toc;
% pimat_sglf=full_pimat_sglf/sum(full_pimat_sglf);
% fprintf('mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n',mean(optim_info(1,:)),mean(optim_info(3,:)),t,mean(optim_info(4,:)));
% 
% 
% 
% fprintf('\nFISTA + Group Lasso L2 with variable size of groups \n');
% param.regul='sparse-group-lasso-l2';
% param2=param;
% param2.groups=groups;  % all the groups are of size 2
% param2.lambda=param.lambda/2;
% param2.lambda2=param.lambda/2;
% % tpm_param=param2
% tic
% [full_pimat_sglv optim_info]=mexFistaFlat(Y,X,W0,param2);
% t=toc;
% pimat_sglv=full_pimat_sglv/sum(full_pimat_sglv);
% fprintf('mean loss: %f, mean relative duality_gap: %f, time: %f, number of iterations: %f\n',mean(optim_info(1,:)),mean(optim_info(3,:)),t,mean(optim_info(4,:)));
% 
% 
% pimats_g = [pimat_glf, pimat_glv, pimat_sglf, pimat_sglv];



clear param;
% parameter of the optimization procedure are chosen
param.pos=true;
% param.L=50; % not more than 20 non-zeros coefficients (default: min(size(D,1),size(D,2)))
param.lambda=0.03; % not more than 20 non-zeros coefficients
param.numThreads=-1; % number of processors/cores to use; the default choice is -1
                     % and uses all the cores of the machine
param.mode=2;        % penalized formulation

tic
alpha=mexLasso(Y,X,param);
t=toc
fprintf('%f signals processed per second\n',size(X,2)/t);
alpha_full = full(alpha/sum(alpha));
te_g = sum(abs(I-R'*alpha_full));

alpha_full_zeros=sum(alpha_full==0)
pimat_l_zeros=sum(pimat_l<0.001)
% pimats_cvx
te
te_g