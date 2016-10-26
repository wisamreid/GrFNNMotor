function Y = reciprocal(X)
% RECIPROCAL  - Invert the nonzero entries of a matrix elementwise.
%             Y = RECIPROCAL(X) has the same sparsity pattern as X
%	      (except possibly for underflow).

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

if issparse(X)
   [m, n]  = size(X);
   [i,j,Y] = find(X);
   Y = sparse(i,j,1./Y,m,n);
else
   Y = 1./X;
end

ceiling = 1.e+16; Y = min(ceiling,Y);
