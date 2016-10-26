function [ap,ad] = ratiotest(dx,dz,ds,dw)
% RATIOTEST   - Ratio test.
%          Usage: [ap,ad] = ratiotest(dx,dz,ds,dw)

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global x y z s w
global Ubounds_exist

% ----- ratio test -----
ap = -1/min([dx./x; -0.1]);
ad = -1/min([dz./z; -0.1]);
if (Ubounds_exist)
   as = -1/min([ds(find(s))./nonzeros(s); -0.1]);
   aw = -1/min([dw(find(w))./nonzeros(w); -0.1]);
   ap = min(ap, as); ad = min(ad, aw);
end
