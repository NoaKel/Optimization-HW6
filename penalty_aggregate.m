function [ f, g, eig_min ] = penalty_aggregate( x, p, A_calligraphic, A_adj_calligraphic, A_ca, B, c, V, penalty )
% Penalty Aggregate

A = A_calligraphic( x, V )-B;
[S,LAMBDA] = eig(A);
LAMBDA_values = diag(LAMBDA);

[phi,dphi] = penalty( LAMBDA_values , p );

f = c'*x + trace(phi);
g = c+ A_adj_calligraphic( S*dphi*S' , A_ca ); % trace'=A*(S*dphi*S')
eig_min=min(LAMBDA_values);

end
