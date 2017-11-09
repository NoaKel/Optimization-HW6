% HW6 

clear all; close all;
% spanning subspace and target V0
[ V, V0 ] = my_example();
[row, col, dimention] = size(V); 


%% SDP

% convert to SDP parameters
[ A_ca, B, c ] = SDP_parameters( V , V0 );

% use SDP solver 
[ x_opt_SDP ,iteration, g_convergence, p_convergence, eig_min ] = SDP_solver ( A_ca, B, c, V , @A_calligraphic , @A_adj_calligraphic, @penalty );

% plots
figure;
semilogy(iteration,eig_min);
title('Minimal Eigenvalue as a Function of p');
xlabel('Iterations');
ylabel('Minimal Eigenvalue');

figure;
semilogy(1:1:length(g_convergence),g_convergence);
title('norm(g) as a function of total quasi newton iterations');
xlabel('Qausi Newton Iterations');
ylabel('norm(g)');


%% Forbenius

[ x_opt_F ] = Frobenius();

%% Compare
linear_comb1=0; linear_comb2=0;
for i=1:dimention 
    linear_comb1=linear_comb1+x_opt_SDP(i+1)*V(:,:,i);
    linear_comb2=linear_comb2+x_opt_F(i)*V(:,:,i);
end

norm2_SDP = norm(linear_comb1-V0); 
norm2_F = norm(linear_comb2-V0);
normF_SDP = norm((linear_comb1-V0),'fro');
normF_F = norm((linear_comb2-V0),'fro');

% plot results
f = figure('Position',[400 200 500 100]);
d = [norm2_SDP, normF_SDP ; norm2_F ,  normF_F ];
cnames = {'L2 Norm','Forbenious Norm'};
rnames = {'SDP','Forbenious'};
t = uitable(f,'Data',d,'ColumnName',cnames,'RowName',rnames,'ColumnWidth',{90});
set(t,'Position',[0 0 500 100]);

%% Dual Problem

p = 10^5;
A = A_calligraphic(x_opt_SDP, V) - B;
[S, A_LAMBDA] = eig(A);
A_LAMBDA_values = diag(A_LAMBDA);

[phi,dphi] = penalty( A_LAMBDA_values , p );
Y = -S*dphi*S';
[~, Y_LAMBDA] = eig(Y);
Y_LAMBDA_values = diag(Y_LAMBDA);

%  a) Positive semidefiniteness of Y; 
    a_ans = min(Y_LAMBDA_values);
%  b) Value of the dual objective function <B,Y> (compare it to the optimal value of primal objective function) 
    b_ans = c'*x_opt_SDP-trace(B'*Y);
%  c) KKT equality   \cal{A*}Y=c 
    c_ans = A_adj_calligraphic(Y, A_ca)-c;
%  d) Complementarity slackness <Y, A(x)>  = 0 
    d_ans = trace(Y'*A);

% plot results    
f = figure('Position',[400 200 500 100]);
d = [a_ans; b_ans ; norm(c_ans) ;  d_ans ];
cnames = {'Dual Problem'};
rnames = {'min Lambda','Duality Gap','norm(A*Y-c)','<Y, A(x)>'};
t = uitable(f,'Data',d,'ColumnName',cnames,'RowName',rnames,'ColumnWidth',{90});
set(t,'Position',[0 0 500 100]);