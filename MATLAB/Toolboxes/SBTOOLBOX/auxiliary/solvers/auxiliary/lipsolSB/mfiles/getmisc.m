function [vmid,xn1,sn1] = getmisc(x,z,s,w);
% GETMISC     - Compute three quantities.
%             Usage: [vmid,xn1,sn1] = getmisc

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County
global Ubounds_exist

xn1 = reciprocal(x); sn1 = [];
if (~Ubounds_exist)
   vmid = reciprocal(z.*xn1);
else
   sn1 = reciprocal(s);
   vmid = reciprocal(z.*xn1 + w.*sn1);
end;
cap = 1.e+11; vmid = full(min(cap,vmid));
