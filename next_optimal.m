function [ optimal_available, optimal_error, pimat_ret ] = next_optimal( rf, I, R, available, chosen )
%NEXT_OPTIMAL Function gives the index, error and weights vector of the
% local one-step optimal value for stock selection from a Market Index  
%   Detailed explanation goes here

clear values T n available_n;

available_n = size(available,2);
n           = (sum(chosen)+1);
T           = size(R,2);
values      = zeros(1,available_n);

for j=1:available_n

    selected = available(j);
    chosen(selected) = true;

    curr_R = R(chosen,:);

%     if ~mod(j,10)  
%         disp([n, j, selected, find(chosen==1)]);
%     end

    % Absolute Value Optimization (Abs)
    cvx_begin quiet
        variable z(T)
        variable pimat(n)
        minimize( rf(I, curr_R, pimat) )

        subject to
            pimat >= 0
            sum(pimat) == 1
    cvx_end

    values(j) = rf(I, curr_R, pimat);
    chosen(selected) = false;
end

[values, sorted_index] = sort(values);

% Return values
optimal_available = sorted_index(1);
optimal_error = values(1);
pimat_ret = pimat; 


end

