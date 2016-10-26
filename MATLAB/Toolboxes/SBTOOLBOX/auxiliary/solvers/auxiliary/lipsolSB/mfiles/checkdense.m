function checkdense(A);
% CHECKDENSE  - Determine and locate dense columns of matrix A.

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global Dense_cols_exist idense ispars;
global DISPLAY_LIPSOL

[m, n] = size(A);
idense = []; ispars = 1:n;
Dense_cols_exist = 0; nzratio = 1;

if (m >  500) nzratio = 0.20; end;
if (m > 1000) nzratio = 0.10; end;
if (m > 5000) nzratio = 0.05; end;

if (nzratio < 1)
   checking = sum(spones(A))/m <= nzratio;
   if any(checking == 0)
      Dense_cols_exist = 1;
      idense = find(sparse(1-checking));   % Dense  column indices
      ispars = find(checking);             % Sparse column indices
   end
   clear checking;
end

if DISPLAY_LIPSOL>=1,
    disp(sprintf('Dense columns (nnz/m > %g): %i',nzratio,length(idense)));
end
