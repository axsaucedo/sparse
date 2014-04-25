function [ v ] = get_norms( A, G )
n = size(G,1);
t = size(A,1);


disp(n);
disp(t);
disp(size(G));
disp(size(A));

v = zeros(1,n);

for i=1:n
    v(i) = norm(A(G(i,:)));
end


end

