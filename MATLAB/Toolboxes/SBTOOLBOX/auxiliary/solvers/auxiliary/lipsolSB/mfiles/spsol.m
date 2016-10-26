function x = spsol(P,b,flag)
% SPSOL      - Solve linear system without dense columns

% Yin Zhang, April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

if (flag == 1) x = linsys(P, b); end;
if (flag ~= 1) x = linsys([],b); end;
