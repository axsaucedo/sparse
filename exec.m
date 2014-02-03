%%%%%%%%%%%%%%%%%%% Yahoo Finance %%%%%%%%%%%%%%%
qtimeout = 3;
y = yahoo;
% for later: ,'CCH.L'
securities = {'CCH'};
index = '^FTSE';
totalassets = size(securities,2);
timeformat = 'd';
fieldname = 'close';
dtTimeFrom = '11/1/2011';
dtTimeUntil = '12/1/2012';

% Fetching Index prices
disp(index);
fetched = fetch(y,index,fieldname,dtTimeFrom,dtTimeUntil,timeformat);
fI = fetched(:,2);

T = size(fI, 1) - 1; % Minus one because returns will have one less

fR = [];
for i = 1:totalassets
    symbol = securities(i);
    disp(symbol);
    
    fetched = fetch(y,symbol,fieldname,dtTimeFrom,dtTimeUntil,timeformat);
    
    fR = [fR, fetched(:, 2)];
end








% qtimeout = 3;
% y = yahoo;
% 
% securities = {'TSLA','AAPL','NOK','DDD','SSYS'};
% timeformat = 'd';
% fieldname = 'close';
% dtTimeFrom = '1/1/2011';
% dtTimeUntil = '12/1/2012';
% 
% cols = 1000;
% 
% [m,n] = size(securities);
% pricedata = [];
% returnsdata = [];
% vardata = [];
% 
% for k1 = 1:n
%     symbol = securities{1:k1};
%     
%     fetched = fetch(y,symbol,fieldname,dtTimeFrom,dtTimeUntil,timeformat);
%     s = size(fetched);
%     price = fetched(1:s, 2);
%     
%     returns = zeros(s-1);
% 
%     for curr = 2:s
%         returns(curr-1) = price(curr) / price(curr-1);
%     end
%     
%     volatility = std(returns);
% 
%     valueatrisk = zeros(100);
% 
%     for curr = 1:100
%         valueatrisk(curr) = 2.58 * volatility * sqrt(curr);
%     end
%     
%     pricedata = [pricedata; price];
%     returnsdata = [returnsdata; returns];
%     vardata = [vardata; valueatrisk];
% end
% 
% size(pricedata);
% 
% 
% data = fetch(y,symbol,fieldname,dtTimeFrom,dtTimeUntil,timeformat);
%  
% s = size(data);
% 
% price = data(1:s, 2);
% 
% returns = zeros(s-1);
% 
% for curr = 2:s
%     returns(curr-1) = price(curr) / price(curr-1);
% end
% 
% volatility = std(returns);
% 
% valueatrisk = zeros(100);
% 
% for curr = 1:100
%     valueatrisk(curr) = 2.58 * volatility * sqrt(curr);
% end
% 
% plot(valueatrisk);
% 
% close(y);
% clearvars;

