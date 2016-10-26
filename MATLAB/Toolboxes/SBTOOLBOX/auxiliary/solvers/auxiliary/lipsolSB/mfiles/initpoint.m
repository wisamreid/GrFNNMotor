function [x,y,z,s,w] = initpoint(A,b,c,ubounds,bnrm,cnrm)
% INITPOINT   - Specifying the starting point.
%          Usage: [x,y,z,s,w] = initpoint(A,b,c,ubounds,bnrm,cnrm)

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global Ubounds_exist
global Dense_cols_exist

[m,n] = size(A); y = zeros(m,1);
pmin = max(bnrm/100, 100);
dmin = cnrm*.425; dmin = max(dmin, pmin/40);
pmin = min(pmin, 1.e+4); dmin = min(dmin, 1.e+3);

e = ones(n,1); 
[P,U] = getpu(A,e);
if Dense_cols_exist
   x = A'*densol(P,U,full(b),1); 
else
   rho = min(100,bnrm);
   x = A'*linsys(P,full(b-rho*A*e),1) + rho*e; 
end
pmin = max(pmin, -min(x)); x = max(pmin,x); 
z = full((c+dmin).*(c > 0) + dmin*(-dmin < c & c <= 0) - c.*(c <= -dmin));

s = []; w = [];
if (Ubounds_exist)
   s = spones(ubounds).*max(pmin, ubounds-x);
   w = spones(ubounds).*( dmin*(c > 0) + ...
       (dmin - c).*(-dmin < c & c <= 0) - 2*c.*(c <= -dmin) );
end
