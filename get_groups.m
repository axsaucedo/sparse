function [ B ] = get_groups( A, g )
n = size(g,1);
t = size(A,2);

B = zeros(n,t);
for i=1:n
    B(g(i),:) = A(g(i));
end


end

