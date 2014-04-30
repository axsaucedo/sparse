randn('state',23432);
rand('state',3454);

%% MAIN PROGRAM %%

ds = datasets;   % Datasets
rf = test_reg_funcs; % Regression Functions

% Get Data
[ I, R, totalassets, T, index ] = ds.ftse100();            % FTSE with sectors as groups
% [ I, R, totalassets, T, index ] = ds.ftse100(9);             % FTSE with 3 groups created randomly
% [ I, R, totalassets, T, index ]   = ds.monte_carlo(100,300,.1,.3,9);
n = size(unique(index),1);


combs = [];
size_combs = [];
comb_idx = [];
group_idx = [];
te_n = [];
all_size = []
all_elapsed = [];


% Creating indexes for combinations
for i = 1:n
    comb = combnk(1:n, i);    
    
    comb_rn = size(comb,1);
    comb_cn = size(comb,2);
    
    combs = [combs; [comb , zeros(comb_rn,n-comb_cn) ]];
    size_combs = [size_combs; ones(comb_rn,1)*i];
end

totalcombinations = size(combs,1)
currRs = cell(1, totalcombinations);
pimats_n = cell(1, totalcombinations);


% Loop for all the possible combinations
for i = 1:totalcombinations
    i
    
    currR = [];
    
    % Generate our dataset for all possible combinations
    tic
    for j = 1:size_combs(i)
        currR = [currR;  R(index==combs(i,j),:)];
    end
    
    curr_totalassets = size(currR,1);
    
    [ pimat_n, error ] = rf.ridge(I,currR,T,curr_totalassets);
    
    
    elapsed=toc
    all_elapsed = [all_elapsed, elapsed];
    
    
    % Calculating Tracking Error
    currRs{i} = currR;
    pimats_n{i} = pimat_n;
   
    te_n = [te_n; error];
    all_size = [ all_size; curr_totalassets ];
        
end

% compound = [];
% for i=1:size(all_elapsed,1)
%     compound = [ compound; sum(all_elapsed(1:i,:),1)]
% end

[srt_te srt_te_idx] = sort(te_n);
[ totalassets-all_size(srt_te_idx) srt_te ]

[srt_size srt_size_idx] = sort(all_size);
[ totalassets-srt_size te_n(srt_size_idx) ]

 xx = flip(totalassets-srt_size);
 yy = flip(te_n(srt_size_idx));

