function [ ds, close, dates ] = googleprices(stockTicker, startDate, endDate)
% PURPOSE: Download the historical prices for a given stock from Google
% Finance and converts it into a MATLAB dataset format.
%---------------------------------------------------
% USAGE: ds = googleprices(stockTicker, startDate, endDate)
% where: stockTicker = Google stock ticker (ExchangeSymbol:SecuritySymbol),
%                      ex. NASDAQ:CSCO for Cisco Stocks.
%        startDate: start date of the prices series. It could be either in
%                   serial matlab form or in Google Date form (mmm+dd,yyyy).
%        endDate: end date of the prices series. It could be either in
%                   serial matlab form or in Google Date form
%                   (mmm+dd,yyyy).
%---------------------------------------------------
% RETURNS: A dataset representing the retrieved prices.
%---------------------------------------------------
% REFERENCES:  a references for the google formats could be found here:
% http://computerprogramming.suite101.com/article.cfm/an_introduction_to_go
% ogle_finance
%---------------------------------------------------

% Version: 1.0
% Written by:
% Display Name: El Moufatich, Fayssal
% Windows: Microsoft Windows NT 5.2.3790 Service Pack 2
% Date: 15-Jun-2010 17:38:18

if isnumeric(startDate)
    startDate = datestr(startDate, 'mmm+dd%2C+yyyy');
end

if ~exist('exportFormat', 'var')
    exportFormat = 'csv';
end

strcat('http://finance.google.com/finance/historical?q=', stockTicker, '&startdate=', startDate, '&enddate=', endDate, '&output=', exportFormat)
% Download the data
fileName = urlwrite(['http://finance.google.com/finance/historical?q=' stockTicker '&startdate=' startDate '&enddate=' endDate '&output=' exportFormat], ['test.' exportFormat]);



% Import the file as a dataset.
ds = dataset('file', fileName, 'delimiter', ',');

n = size(ds,1);
close = zeros(1,n);
dates = cell(1,n);
for i=1:n
    dates{i} = ds{i,1};
    close(i) = ds{i,5};
end

end