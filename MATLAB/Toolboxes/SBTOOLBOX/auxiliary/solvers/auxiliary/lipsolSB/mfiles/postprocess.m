function [x, objp] = postprocess(x,lbounds)
% POSTPROCESS - Recover the original variables for the solution.
%	Usage: [x, objp] = postprocess(x,lbounds)

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global b_orig c_orig
global Lbounds_non0
global col_scaled colscl
global Fixed_exist ifix infx xfix
global Zrcols_exist izrcol inzcol xzrcol
global Sgtons_exist isolved insolved xsolved

if (col_scaled)   x = x.*colscl;   end;
if (Lbounds_non0) x = x + lbounds; end;

if (Sgtons_exist)
   tmp(insolved) = x;
   tmp(isolved) = xsolved;
   x = sparse(tmp);
end;

if (Zrcols_exist)
   tmp(inzcol) = x;
   tmp(izrcol) = xzrcol;
   x = sparse(tmp);
end;

if (Fixed_exist)
   tmp(infx) = x; 
   tmp(ifix) = xfix;
   x = sparse(tmp);
end;

if (size(x,1) < size(x,2))
   x = x';
end;

objp = full(c_orig'*x);
