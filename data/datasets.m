function data = datasets
% DATASETS - contains several anonymous functions that can be
% used to import market or simulated data.
% 
% The anonymous functions available to retreive data are:
%       ftse100 - retreives the 100 stocks from the FTSE100 market index
% In this case, we call the ftse100() function, which returns a matrix of
% returns R, the Market Index of these returns I = mean(R,2), the number of
% assets n, the time window T and the groupings contained in the variable
% groups

    data.ftse_yahoo         = @ftse_yahoo;
    data.ftsegoogle         = @ftse_google;
    data.monte_carlo        = @monte_carlo;
end

function [ I, R, n, T, index, groups ] = ftse_yahoo(varargin)
    if nargin > 1
        error('Wrong Number of Arguments: Optional group number argument');
    end

    fR = importdata('FPSESecuritiesData');
    fR = fR(1:270,:);

    T = size(fR, 1)-1;
    n = size(fR,2);

    % Calculating returns for all individual assets
    R = [];
    for i = 1:n

        returns = [];

        for curr = 2:(T+1)
            returns = [returns, fR(curr, i) / fR(curr-1, i)];
        end

        R = [R; returns];
    end

    % Calculating returns for Index
    I = mean(R)';
    
    % This is the group index based in Stock Sector taken from Google and Yahoo finance
    index = [1, 3, 5, 5, 9, 9, 4, 1, 7, 5, 6, 9, 5, 5, 3, 4, 5, 1, 9, 4, 2, 2, 1, 3, 8, 9, 9, 1, 3, 9, 9, 1, 9, 2, 4, 6, 1, 5, 5, 1, 2, 9, 3, 9, 2, 1, 2, 5, 5, 5, 5, 9, 2, 1, 9, 3, 8, 2, 5, 1, 5, 2, 2, 2, 5, 4, 4, 9, 1, 1, 9, 9, 1, 5, 5, 2, 3, 5, 7, 6, 5, 9, 6, 2, 8, 5, 8, 3, 1, 9, 3, 2, 2, 8, 7, 9, 2, 2, 7]';
    
    if nargin == 1
        [index, groups] = knn(R, varargin{1});
    end
end

%% Load FTSE 100 retreived from Google Finance from disk
function [ c, d, i ] = ftse_google
    load('FPSEGoogleData.mat');
    d = dates_data;
    i = in_data;
    c = close_data;
end

function [I, R, n, T, index, groups] = monte_carlo(numberStocks,rate,volbase,volvar,varargin)

    if nargin > 5
        error('Wrong Number of Arguments');
    end

    sampleRate   = 1/rate; % Assume 250 working days per year.
    % volatility   = 0.10;  % 20% volatility.

    stockPrice   = 100;   % Stock price starts at $100.
    timeToExpiry = 1;     % Length of time.
    dividend     = 0.01;  % 1% annual dividend yield.
    riskFreeRate = 0.005; % 0.5 percent.


    price = stockPrice;
    all_prices = [];
    
    for i=1:numberStocks
        time = 0;
        volatility = volbase + rand()*volvar;

    %     times = [];
        prices = [];

        while time < timeToExpiry
            time = time + sampleRate;
            drift = (riskFreeRate - dividend - volatility*volatility/2)*sampleRate;
            perturbation = volatility*sqrt( sampleRate )*randn();
            price = price*exp(drift + perturbation);
    %         times = [times time];
            prices = [prices, price];
        end

        all_prices = [all_prices; prices];
    end

    R = all_prices(:,2:end)./all_prices(:,1:end-1);
    I = mean(R)';
    n = numberStocks;
    T = rate;

    if nargin == 5
        [index, groups] = knn(R, varargin{1});
    end

end


