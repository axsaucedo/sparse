function test_sparsegrouplasso(I, R, groups, iterations, l1add, l2add)
% TEST_SPARSEGROUPLASSO - Test the Sparse Group Lasso model for x
% iterations and adding l1add and l2add to l1 and l2 respectively on each iteration.
% 
% By default this function will run for 20 iterations with l2+=.5, l1+=l2/5000
    l1=0;
    l2=0;
    
    if nargin<6
        l1add = 0.0005;
    end
    if nargin<5
        l2add = .5;
    end
    if nargin<4
        iterations = 20;
    end
    
    fweights = [];  % Feature Weghts
    gweights = [];  % Group Weights
    
    elapsed = [0];
    
    for i=1:iterations
        tic
        
        disp(sprintf('\t> SparseGroupLasso #%s out of %s iterations.',num2str(i),num2str(iterations)));
        disp(sprintf('\t> L1=%f, L2=%f.',l1,l2));
        
        pimat = SparseGroupLasso(I, R, groups,l1,l2); 
        l2=l2+l2add;
        l1=l1+l1add;
        
        fweights =  [ fweights, sum(pimat,2) ];
        gweights =  [ gweights; sum(pimat,1) ];
        
        elapsed  =  [ elapsed, toc ];
        disp(sprintf('\t\tAverage time=%s',num2str(mean(elapsed))));
    end
   
    figure;
    plot(fweights');
    title(sprintf('l1+=%f, l2+=%f, omega=%f',l1add, l2add, l2add/l1add));
    ylabel('Feature Weights');
    xlabel('Iterations');
%     figure;
%     plot(gweights);
%     ylabel('Group Weights');
%     xlabel('Lambda');
end
