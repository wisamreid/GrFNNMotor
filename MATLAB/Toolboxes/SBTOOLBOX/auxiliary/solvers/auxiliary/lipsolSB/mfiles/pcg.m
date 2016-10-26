function [x, iter] = pcg(P,U,b,tolcg,flag,x)

% pcg       Solve Ax = b using preconditioned conjugate gradient,
%           where A = P + U*U', P is sparse and U is dense.
%           P is used as the preconditioner.

n = length(b); bnrm = norm(b); iter = 0; 
maxiter = min(250,.5*n);
if (nargin < 6)
   x = zeros(n,1); r = b; rnrm = bnrm;
else
   r = b - P*x;
   if ~isempty(U) r = r - U*(U'*x); end;
   rnrm = norm(r); 
end;
while (rnrm > tolcg) & (iter <= maxiter)
    if (flag == 1) z = linsys(P, r); end;
    if (flag ~= 1) z = linsys([],r); end;
    iter = iter + 1;
    if iter == 1
       p = z; rtz = r'*z;
    else 
       oldrtz = rtz; rtz = r'*z;
       beta = rtz / oldrtz;
       p = z + beta*p;
    end;
    Ap = P*p; if ~isempty(U) Ap = Ap + U*(U'*p); end;
    alpha = rtz / (p'*Ap);
    x = x + alpha*p;
    r = r - alpha*Ap;
    rnrm = norm(r);
end;
if iter == maxiter error('Too many iterations in PCG solver'); end;
