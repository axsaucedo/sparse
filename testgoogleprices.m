


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

