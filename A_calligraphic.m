function [ Y ] = A_calligraphic( x , V )
% A caligraphic operator 
% Recievs xvector and subspace of matrices
% reutrns matrix 

[row, col, dimention] = size(V); 

Y11 = x(1)*eye(col); % tI
Y22 = x(1)*eye(row); % tI

Y21 = 0; 
for i=1:dimention
    Y21 = Y21 + x(1+i)*V(:,:,i); % sum(xV)
end

Y12=Y21';

Y = [ Y11, Y12;...
      Y21, Y22];
  

    
end

