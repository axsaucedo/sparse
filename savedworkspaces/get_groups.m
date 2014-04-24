function [ B ] = get_groups( a, g )
n = size(g,1);
t = size(a,2);

G = zeros(n,t);
for i=1:n
    G(i,g(i,:)) = a(g(i,:))
end
B=G';

end

