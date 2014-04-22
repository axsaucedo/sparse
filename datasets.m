function data = datasets
    data.ftse100            = @ftse100;
    data.naive              = @naive;
    data.fetch_google       = @fetch_google;
    data.ftse100_google     = @ftse100_google;
end

function [ I, R, n, T ] = ftse100()
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
end

function [ I, R, n, T ] = naive(n, T)

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

