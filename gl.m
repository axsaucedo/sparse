d = 6;      % dimension of groups
m = 8;      % number of groups
n = 274;   % number of data points
r = m/2;      % number of selected groups

% Getting data
fI = importdata('FPSEOnlyData');
fR = importdata('FPSESecuritiesData');

% Calculating returns for all individual assets
R = [];
for i = 1:size(fR,2)

    returns = [];

    for curr = 2:(n+1)
        returns = [returns, fR(curr, i) / fR(curr-1, i)];
    end

    R = [R; returns];
end
% Calculating returns for Index
I = mean(R)';
fR = fR(:,1:d*m);

R = [];
for i = 1:d*m

    returns = [];

    for curr = 2:(n+1)
        returns = [returns, fR(curr, i) / fR(curr-1, i)];
    end

    R = [R; returns];
end