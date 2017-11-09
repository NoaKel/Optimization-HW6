function [ x_opt, iteration, g_convergence, p_convergence, eig_min ] = SDP_solver ( A, B, c, V , A_calligraphic , A_adj_calligraphic, penalty )
% General SDP Solver

% Initializing
x = zeros(length(c),1);
alpha = 2;
p = 1;
p_max = 10^5;
epsilon = 10^-5;
f_prev = 0;
k = 1;
diff = inf; 
max_iter=10^5;
iteration = zeros(max_iter,1);
eig_min = zeros(max_iter,1);
g_convergence = [];
p_convergence = zeros(max_iter,1);

% Iterations
while( diff>=epsilon || k>=max_iter) 
    aggregate_function = @(x) (penalty_aggregate( x, p, A_calligraphic, A_adj_calligraphic, A, B, c, V, penalty ));
    
    [ x , g_norm, iter ] = BFGS( aggregate_function, x );
    if k>1
        iteration(k)=iteration(k-1)+iter;
    else
        iteration(k)=iter;
    end
    g_convergence(1+end:end+length(g_norm))=g_norm;
    
    p_convergence(k)=p;
    p = min(p*alpha,p_max);
    
    
    [f, ~, eig_min(k) ] = penalty_aggregate( x, p, A_calligraphic, A_adj_calligraphic, A, B, c, V, penalty );

    diff = abs(f-f_prev);
    f_prev = f;
    k = k + 1;
end

x_opt = x;
iteration = iteration(1:k-1,1);
eig_min = eig_min(1:k-1,1);
p_convergence = p_convergence(1:k-1,1);
end

