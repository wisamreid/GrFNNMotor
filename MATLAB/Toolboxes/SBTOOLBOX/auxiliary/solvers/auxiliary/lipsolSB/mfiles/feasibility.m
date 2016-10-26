function [Rb,Rc,Ru,rb,rc,ru] = feasibility(A,b,c,ubounds)
% FEASIBILITY - Evaluate feasibility residual vectors and their norms.
%	Usage:  [Rb,Rc,Ru,rb,rc,ru] = feasibility(A,b,c,ubounds)

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global x y z s w
global Ubounds_exist

Rb = A*x - b; Rc = A'*y + z - c; Ru = [];
if (Ubounds_exist)
   Rc = Rc - w;
   Ru = spones(s).*x + s - ubounds;
end
rb = norm(Rb); rc = norm(Rc); ru = 0;
if Ubounds_exist ru = norm(Ru); end
if any(isnan([rb rc ru])) error(' NaN occured.'); end
