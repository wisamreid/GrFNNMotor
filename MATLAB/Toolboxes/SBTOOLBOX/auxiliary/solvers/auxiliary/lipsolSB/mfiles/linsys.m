function [sol] = linsys(P,rhs,save_factor)
% LINSYS     - Solve a positive definite system of linear equations.
%        LINSYS(P,rhs) returns a "solution" to linear system Px = rhs,
%        where the matrix P is sparse symmetric and positive definite.
%        The function SYMBFCT must be called at least once before 
%        LINSYS is called.  LINSYS uses three Fortran MEX programs.
%
%        If P = [], LINSYS will try to use, if any, an existing Cholesky 
%        factor.

% Yin Zhang, July, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global XLNZ NNZL PERM INVP LNZ LOOP_LEVEL
global XSUPER XLINDX LINDX SNODE SPLIT TMPSIZ
global LNZ0

% ----- check input rhs -----
if issparse(rhs)
   error('RHS vector must be in full matrix format.');
end;
m = length(PERM);
if (max(size(rhs)) ~= m) | (min(size(rhs)) ~= 1)
   error('No symbolic factor or input sizes mismatch.');
end;

if ~isempty(P)
% ----- check input matrix P -----
   if (size(P,1) ~= m | size(P,2) ~= m | ~issparse(P))
      error('Input matrix must be square and sparse.');
   end;
% ----- Remove diagonal from P ------
   Pdiag = full(diag(P));
   P = P - sparse(1:m,1:m,Pdiag);
% ----- Cholesky factorization ------
   [LNZ] = inpnv(Pdiag,P,INVP,PERM,XLNZ,XSUPER,XLINDX,LINDX,NNZL);
   [LNZ] = blkfct(XLNZ,XSUPER,SNODE,SPLIT,XLINDX,LINDX,LNZ,TMPSIZ,LOOP_LEVEL);
end;
if nargin == 3 if save_factor LNZ0 = LNZ; end; end

% ------ back solve ------
sol = zeros(m,1);
sol(PERM) = blkslv(XLNZ,XSUPER,XLINDX,LINDX,LNZ,rhs(PERM));
