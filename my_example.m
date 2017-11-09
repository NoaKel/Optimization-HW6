function [ V, V0 ] = my_example( )
% My subspace and V0

V = zeros(2,3,3); % two first indexes - matrix dimention, last index - subspace dimention
V(:,:,1) = [1,0,0;...
            0,1,0];
V(:,:,2) = [0,0,1;...
            1,1,1];
V(:,:,3) = [0,1,0;...
            0,1,1];
        
% V0= [1,3,2;...
%      2,6,5];

V0= [1.1,3.5,1.5;...
     2.3,6.2,5.04];


end

