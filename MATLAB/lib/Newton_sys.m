function [z,iter,res] = Newton_sys(fun,Jfun,x0,tol,Nmax)
% This function computes the solution to a nonlinear system
% of dimension n. 

% Input: fun  -> Function handle representing the function
%        Jfun -> Function handle representing the Jacobian
%        x0   -> Initial guess
%        tol  -> tolerange on the 2-norm of the increment
%        Nmax -> maximum number of iteration

% Output  z    -> zero of the nonlinear system
%         iter -> number of iterations to get z
%         res  -> 2-norm of the residual
% Jfun(x0)
e = norm(x0);
iter = 0;

while e >= tol && iter<Nmax
    
x1= x0 - inv(Jfun(x0))*fun(x0);

e   = norm(x1-x0);
res = norm(fun(x1));
iter = iter+1;
x0=x1;
end

if iter==Nmax & e > tol
    warning('Maximum number of iterations reached without achieving e<tol');
end

z=x1;
res = norm(fun(x1));




end

