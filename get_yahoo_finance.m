%%%%%%%%%%%%%%%%%%% Yahoo Finance %%%%%%%%%%%%%%%
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
