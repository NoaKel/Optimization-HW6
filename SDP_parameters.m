function [ A, B, c ] = SDP_parameters( V , V0 )
% breaks problem into SDP parameters

[row, col, dimention] = size(V); 

V_trans = permute(V(:,:,:),[2 1 3]); % change dimentions 1 & 2


c = [1; zeros(dimention,1)]; % assume x = [t, x1, x2, .., xn]

B = [ zeros(col) , V0'; ...
      V0, zeros(row)];

A = zeros ( (row+col), (row+col), dimention );

A ( col+1:end , 1:col, :) = V;
A ( 1:col, col+1:end, :) = V_trans; 





end

