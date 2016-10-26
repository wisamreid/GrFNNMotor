function symbfct(P)
% SYMBFCT     - Symbolic factorization for sparse Cholesky factorization.
%            symbfct(P) does a minimum degree ordering and a symbolic
%            factorization for a symmetric square sparse matrix P, and 
%            passes the results as global variables.
%            It calls two Fortran MEX files.

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global XLNZ NNZL PERM INVP CACHE_SIZE
global XSUPER XLINDX LINDX SNODE SPLIT TMPSIZ
global DISPLAY_LIPSOL

% ----- check input matrix P -----
m = size(P,1);
if size(P,2) ~= m | ~issparse(P) | isempty(P)
   error('Matrix must be square, sparse and nonempty.');
end;
Pdiag = diag(diag(P)); P = P - Pdiag;

% ----- min-degree ordering -----
t = cputime;
if DISPLAY_LIPSOL >= 1,
    disp('  min-degree ordering ... ');
end
[PERM, INVP] = ordmmd(P);
if DISPLAY_LIPSOL >= 1,
    disp(sprintf('Done. CPU seconds: %g\n',cputime-t));
end

% ----- symbalic factorization -----
t = cputime;
if DISPLAY_LIPSOL >= 1,
    disp('  calling symfct.mex* ... ');
end
[XLNZ,NNZL,XSUPER,XLINDX,LINDX,SNODE,SPLIT,TMPSIZ,PERM,INVP] = ...
                         symfct(P,PERM,INVP,CACHE_SIZE);
if DISPLAY_LIPSOL >= 1,
    fprintf('Done. CPU seconds: %g\n',cputime-t);
end