function [ f ] = fetch
    f.fetch_google   = @fetch_google;
    f.fetch_yahoo    = @fetch_yahoo;
end

function fetch_yahoo
    qtimeout = 3;
    y = yahoo;
    % RIGHT NOW WE HAVE AN FPSE 99 BECAUSE 'CCH.L' DOES NOT FETCH!
    securities = {'AAL.L','ABF.L','ADM.L','ADN.L','AGK.L','AMEC.L','ANTO.L','ARM.L','AV.L','AZN.L','BA.L','BAB.L','BARC.L','BATS.L','BG.L','BLND.L','BLT.L','BNZL.L','BP.L','BRBY.L','BSY.L','BT-A.L','CCL.L','CNA.L','CPG.L','CPI.L','CRDA.L','CRH.L','DGE.L','EXPN.L','EZJ.L','FRES.L','GFS.L','GKN.L','GLEN.L','GSK.L','HL.L','HMSO.L','HSBA.L','IAG.L','IHG.L','IMI.L','IMT.L','ITRK.L','ITV.L','JMAT.L','KGF.L','LAND.L','LGEN.L','LLOY.L','LSE.L','MGGT.L','MKS.L','MNDI.L','MRO.L','MRW.L','NG.L','NXT.L','OML.L','PFC.L','PRU.L','PSN.L','PSON.L','RB.L','RBS.L','RDSA.L','RDSB.L','REL.L','REX.L','RIO.L','RR.L','RRS.L','RSA.L','RSL.L','SAB.L','SBRY.L','SDR.L','SGE.L','SHP.L','SL.L','SMIN.L','SN.L','SPD.L','SSE.L','STAN.L','SVT.L','TATE.L','TLW.L','TPK.L','TSCO.L','TT.L','ULVR.L','UU.L','VED.L','VOD.L','WEIR.L','WMH.L','WOS.L','WPP.L'};
    % securities = {'BBRY','YHOO'}
    index = '^FTSE';
    totalassets = size(securities,2);
    timeformat = 'd';
    fieldname = 'close';
    % dtTimeFrom = '11/1/2011';
    % dtTimeUntil = '12/1/2012';

    dtTimeFrom = '12/1/2012';
    dtTimeUntil = '12/30/2013';

    % Fetching Index prices
    disp(index);
    fetched = fetch(y,index,fieldname,dtTimeFrom,dtTimeUntil,timeformat);
    fI = fetched(:,2);

    fR = [];
    for i = 1:totalassets
        symbol = securities(i);
        disp(symbol);

        fetched = fetch(y,symbol,fieldname,dtTimeFrom,dtTimeUntil,timeformat);

        fR = [fR, fetched(:, 2)];
    end

    T = min(size(fI, 1),size(fR, 1)) - 1; % Minus one because returns will have one less
end

function fetch_google
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
                [ ds ] = googleprices(securities{i}, startDate, endDate);
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

function [ ds, close, dates ] = googleprices(stockTicker, startDate, endDate)

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