function data = datasets
% DATASETS - contains several anonymous functions that can be
% used to import market or simulated data.
% 
% The anonymous functions available to retreive data are:
%       ftse100 - retreives the 100 stocks from the market
% In this case, we call the ftse100() function, which returns a matrix of
% returns R, the Market Index of these returns I = mean(R,2), the number of
% assets n, the time window T and the groupings contained in the variable
% groups

    data.ftse100            = @ftse100;
    data.naive              = @naive;
    data.monte_carlo        = @monte_carlo;
    data.fetch_google       = @fetch_google;
    data.ftse100_google     = @ftse100_google;
end

function [ I, R, n, T, index, groups ] = ftse100(varargin)
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

function [ I, R, n, T, index, groups ] = naive(n, T, varargin)
    if nargin > 3
        error('Wrong # Arguments: Args are n=asset number, T=time window, (Optional) gropup_n=number of groups');
    end

    assets = [];
    for i=1:n
        volatility = .1 + rand()/5-.1;
        old_price = 100 - rand()*50-25;
    
        price = [];
        for j=1:T
            rnd = rand(); % generate number, 0 <= x < 1.0
            change_percent = 2 * volatility * rnd;
            if (change_percent > volatility)
                change_percent = change_percent-(2 * volatility);
            end
            change_amount = old_price * change_percent;
            new_price = old_price + change_amount;
            price = [price; old_price/new_price];
        end
        assets = [ assets, price ];
    end
    R = assets';
    I = R'*(ones(n,1)/n);
    
    if nargin == 3
        [index, groups] = knn(R, varargin{1});
    end
end


function [ close, dates ] = fetch_google()
    
    securities = {'CCH','AAL','ABF','ADM','ADN','AGK','AMEC','ANTO','ARM','AV','AZN','BA','BAB','BARC','BATS','BG','BLND','BLT','BNZL','BP','BRBY','BSY','BT','CCL','CNA','CPG','CPI','CRDA','CRH','DGE','EXPN','EZJ','FRES','GFS','GKN','GLEN','GSK','HL','HMSO','HSBA','IAG','IHG','IMI','IMT','ITRK','ITV','JMAT','KGF','LAND','LGEN','LLOY','LSE','MGGT','MKS','MNDI','MRO','MRW','NG','NXT','OML','PFC','PRU','PSN','PSON','RB','RBS','LON:RDSA','RDSB','REL','REX','RIO','RR','RRS','RSA','RSL','SAB','SBRY','SDR','SGE','SHP','SL','SMIN','SN','SPD','SSE','STAN','SVT','TATE','TLW','TPK','TSCO','TT','ULVR','UU','VED','VOD','WEIR','WMH','WOS','WPP'};
    startDate = 'Dec+01,+2010';
    endDate = 'Apr+22,+2014';


    n = size(securities,2)
    all_stocks = cell(1,n);


    for i=1:n
        i
        disp(securities(i));

        test = 1
        while test
            try
importdata('FPSEGoogleData.mat');                [ ds ] = googleprices(securities{i}, startDate, endDate);
                all_stocks{i} = ds;
                test = 0;
            catch error
                disp('Probably a network error, try again');
                test = 1;
            end
        end
    end


    T = size(all_stocks{2},1);
    close = zeros(n,T);
    in = zeros(n,1);

    for i=1:n
        [i size(all_stocks{i},1)]
        in(i) = size(all_stocks{i},1);

        try
            for j=1:T
                close(i,j) = all_stocks{i}{j,5};
            end
        catch
           disp('out'); 
        end
    end

    dates = cell(n,1);
    for i=1:T
        dates{i} = all_stocks{2}{i,1};
    end

end

function [ c, d, i ] = ftse100_google()
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


