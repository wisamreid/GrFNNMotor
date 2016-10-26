function [A,b,c,ubounds] = scaling(A,b,c,ubounds)
% SCALING     - Scale matrix A.
%           Usage: [A,b,c,ubounds] = scaling(A,b,c,ubounds)

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global data_changed
global col_scaled colscl
global Ubounds_exist nub
global DISPLAY_LIPSOL

[m, n] = size(A);
badscl = 1.e-4;

col_scaled = 0;
absnzs = abs(nonzeros(A));
thescl = min(absnzs)/max(absnzs);
clear absnzs

if (thescl < badscl)
   if DISPLAY_LIPSOL >= 1,
       disp('Scaling ...');
   end

% ----- scaling vectors ------
   AA = abs(A);
   colscl = full(sqrt(max(AA)'));
   rowscl = full(sqrt(max(AA')'));
   clear AA;

% ----- column scaling -----
   if (Ubounds_exist) ubounds = ubounds.*colscl; end
   colscl = reciprocal(colscl);
   A = A*sparse(1:n,1:n,colscl);
   c = c.*colscl;
   col_scaled = 1;

% ----- row scaling -----
   rowscl = reciprocal(rowscl);
   A = sparse(1:m,1:m,rowscl)*A;
   b = b.*rowscl;
   bnrm = norm(b);
   if (bnrm > 0) 
      q = median([1 norm(c)/bnrm 1.e+8]);
      if (q > 10) A = q*A; b = q*b; end
   end
   data_changed = 1;

end
