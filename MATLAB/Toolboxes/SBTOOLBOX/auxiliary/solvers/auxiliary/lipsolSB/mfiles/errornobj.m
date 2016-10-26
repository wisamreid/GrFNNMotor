function [trerror,rrb,rrc,rru,rdgap,objp,objd] = errornobj...
                  (b,rb,bnrm,c,rc,cnrm,ubounds,ru,unrm,dgap)
% RELERROR    - Calculate the total relative error and objective values.

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global x y z s w
global Ubounds_exist nt

rrb = rb/max(1,bnrm); rrc = rc/max(1,cnrm); rru = 0;
objp = full(c'*x); objd = full(b'*y);
if (Ubounds_exist) 
   rru = ru/(1+unrm); 
   objd = objd - w'*ubounds;
end
rdgap = dgap/nt;
%rdgap = abs(objp-objd)/(1 + min(abs(objp),abs(objd)));
trerror = max([rrb rrc rru rdgap]);
