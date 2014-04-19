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

cvx_n = 80;
cvx_R = R(1:cvx_n,:);
cvx_I = I;



cvx_begin quiet
    variable z_a(T)
    variable pimat_a(cvx_n)
    minimize( (1/T) * sum(abs(cvx_I - transpose(cvx_R)*pimat_a)) )
    
    subject to
        pimat_a >= 0
        sum(pimat_a) == 1
cvx_end















% % CVX to find optimal value for NCCVAR
% cvx_begin quiet
%     variable z_n(T)
%     variable Alpha_n
%     variable pimat_n(cvx_n)
%     minimize( Alpha_n + divcoef * sum(z_n) )
%     subject to
%         z_n >= 0
%         z_n - abs(cvx_I - cvx_R'*pimat_n) + Alpha_n >= 0
% 
%         % QUESTION: This constrain ensures that less assets have a position of 0?
%         % (Also performs better without this constrain)
% %         norm(pimat_n) <= C
%         pimat_n >= 0
%         sum(pimat_n) == 1
% cvx_end
% 
% cvx_begin quiet
%     variable z_c(T)
%     variable Alpha_c
%     variable pimat_c(cvx_n)
%     minimize( Alpha_c + divcoef * sum(z_c) )
%     subject to
%         z_c >= 0
%         z_c - (cvx_I - transpose(cvx_R)*pimat_c) + Alpha_c >= 0
%         pimat_c >= 0
%         sum(pimat_c) == 1
% 
%         % QUESTION: Should this norm constraint be in CVAR as well?
%         % Note: This performs better without this constraint
% %         norm(pimat_c) <= C
% cvx_end
% 
% % CVX to find optimal value for Tracking Error Abs
% cvx_begin quiet
%     variable z_a(T)
%     variable pimat_a(cvx_n)
%     minimize( (1/T) * sum(z_a) )
%     subject to
%         z_a >= 0
% %         transpose(z_c)*z_c <= power(C,2)
%         z_a - abs(cvx_I - transpose(cvx_R)*pimat_a) >= 0
%         pimat_a >= 0
%         sum(pimat_a) == 1
%         % QUESTION: Should this norm constraint be in CVAR as well?
% %         norm(pimat_a) <= C
% cvx_end
% 
% % CVX to find optimal value for Lasso
% cvx_begin quiet
%     variable z_l(T)
%     variable pimat_l(totalassets)
%     
%     minimize(sum(power(I - transpose(R)*pimat_l, 2)) + delta*sum(pimat_l))
%     
%     subject to
%         sum(pimat_l) <= 1
%         pimat_l >= 0
% cvx_end
% 
% pimats_cvx = [pimat_n pimat_c pimat_a ];
% 
% % Calculating Tracking Error
% cvx_Ret = [ 
%             abs(I - cvx_R' * pimat_n) ...
%             abs(I - cvx_R' * pimat_c) ...
%             abs(I - cvx_R' * pimat_a) ...
%             abs(I - R' * pimat_l)
%           ]
% 
% % Tracking Errors
% te = sum(cvx_Ret)