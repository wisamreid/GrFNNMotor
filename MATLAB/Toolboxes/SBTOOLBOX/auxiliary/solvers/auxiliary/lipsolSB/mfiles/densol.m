function x = densol(P,U,b,flag)
% DENSOL      - Solve linear system with dense columns 
%               using either the Sherman-Morrison formula
%               or a preconditioned conjugate gradient method

% Yin Zhang, April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global Sherman_OK Hist
if ~exist('Sherman_OK') Sherman_OK = 1; end;

bnrm = norm(b); x = zeros(size(b));
tol1  = min(1.e-2, 1.e-7*(1+bnrm));
tol2  = 10*tol1; tolcg = tol1;

iter = size(Hist,2);
if iter > 0 
   trerror = Hist(1,iter);
   tolcg = min(trerror, tolcg);
end;

if Sherman_OK
   x = sherman(P,U,b,flag); end;
   resid = norm(b - P*x - U*(U'*x));
   if (resid > tol1) Sherman_OK = 0; end;
   if (resid < tol2) return; end;
end;
if ~Sherman_OK 
   x = pcg(P,U,b,tolcg,flag,x);
end;
