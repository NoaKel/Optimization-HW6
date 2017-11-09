function [ sol ] = Frobenius()

[ V, V0 ] = my_example();
[~, ~, dimention] = size(V); 

A = zeros(dimention, dimention);
b = zeros(dimention, 1);

for i=1:dimention
    for j=1:dimention
        A(i,j) = trace(V(:,:,i)' * V(:,:,j));
    end
    b(i,1) = trace(V(:,:,i)' * V0);
end

sol = A \ b;
% 
% linear_comb=0;
% for i=1:dimention
%     linear_comb = x(i)*V(:,:,i);
% end
% 
% f = (norm(linear_comb - V0, 'fro'))^2;
% 
% g = 2 * (A * x - b);

end

