function [ x_min, norm_g, k ] = BFGS( func, x )
%Quasi Newton BFGS method with Amrijo inexact line search

% Initializing
n=length(x);
sigma=0.25;
beta=0.5;
epsilon=10^-5;
a0=1;
B=eye(n); % initial guess
k=0;
[f,g,~]=func(x);

% Quasi Newton BFGS
while(norm(g)>=epsilon)
    x_prev=x;
    g_prev=g;
    % step 1: obtain d
    d=-B*g;
    % step 2: line search
    alpha=armijo_step(a0,x,func,f,g,d,sigma,beta);
    % step 3: calculate new x,f,g
    x=x+alpha*d;
    [f,g,~]=func(x);
    p=x-x_prev;
    q=g-g_prev;
    % step 4: updating B
    d_norm=d/norm(d);
    if ((g_prev'*d_norm)<(g'*d_norm) || k==0)
        s=B*q;
        tau=s'*q;
        meu=p'*q;
        v=(p/meu)-(s/tau);
        B=B+((p*p')/meu) -((s*s')/tau)+(tau*(v*v'));
    end
    k=k+1;
    norm_g(k)=norm(g);
end

x_min=x;

end

