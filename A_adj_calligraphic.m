function [ x ] = A_adj_calligraphic( Y, A )
% A adjoint operator
% Recieves a matrix 
% Returns a vector 

[row, col, dimention] = size(A); 


x = zeros(dimention+1,1);
x(1)=trace(Y);

for i=1:dimention
    x(i+1)=trace((Y'*A(:,:,i)));
end


end

