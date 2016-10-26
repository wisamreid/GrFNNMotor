function [sol] = project(rhs)

% Yin Zhang, 1997

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

% ------ back solve ------
sol = zeros(m,1);
sol(PERM) = blkslv(XLNZ,XSUPER,XLINDX,LINDX,LNZ0,rhs(PERM));
