function [ phi, dphi ] = penalty( vX , p  )


% calculate phi
val = (1/p) * phi_aux(p * vX);

% calculate dphi
vX = vX * p;
dval = zeros(size(vX));
dval(vX <= 0.5) = vX(vX <= 0.5) - 1;
dval(vX >  0.5) = -1 ./ (4 * vX(vX >  0.5));

% turn into diagonal matrix
phi = diag(val);
dphi = diag(dval);

end


function val = phi_aux(vX)

val = zeros(size(vX));

val(vX <= 0.5) = 0.5 * vX(vX <= 0.5).^2 - vX(vX <= 0.5);
val(vX >  0.5) = -0.25 * log(2*vX(vX >  0.5)) - 3/8;

end