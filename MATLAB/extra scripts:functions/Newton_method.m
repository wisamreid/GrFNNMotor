function [z0,iter,res,his] = Newton_method(fun,dfun,x0,tol,Nmax)

% This function implements the secant method to find the zeros of f(x)
% Inputs:
% ------
%        fun  -> function handle representing f(x)
%        dfun -> derivative of f(x)
%        x0  -> initial guess on the zero
%        tol  -> maximum tolerance accepted on the increment
%               abs(x(i+1)-x(i))
%       Nmax -> maximum number of iterations allowed
%
% Outputs:
% -------
%        z0  -> approximation of the zero of f(x)
%       iter -> number of interation to determine the zero
%       res  -> |f(z0)|
%       his  -> vector collecting all the elements of the sequence 
%               his = [x(1) x(2) x(3) ..... ]

x1     = x0-fun(x0)/dfun(x0);
his(1) = x1;

iter=1;

if abs(x1-x0)<=tol
    z0=x1;
    res=abs(fun(z0));
    return;
end

% USING WHILE LOOP
while abs(x1-x0)>tol && iter<=Nmax
    x0        = x1;
    x1        = x0-fun(x0)/dfun(x0);
    iter      = iter +1;
    his(iter) = x1;
    
end

% % USING FOR LOOP
% for k=3:Nmax
%     x0        = x1;
%     x1        = x0-fun(x0)/q;
%     his(k)    = x1;
% if abs(x0-x1)<=tol    
%     z0=his(k);
%     iter=k;
%     res=abs(fun(z0));
%     break
% end
% end

z0 = his(end);
res=abs(fun(z0));

if iter==Nmax
    warning('Perhaps you need more iterations')
end




