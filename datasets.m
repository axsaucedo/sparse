function data = datasets
    data.yahoo_finance = @yf;
    data.naive = @naive;
end

function [ I, R, n, T ] = yf()
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