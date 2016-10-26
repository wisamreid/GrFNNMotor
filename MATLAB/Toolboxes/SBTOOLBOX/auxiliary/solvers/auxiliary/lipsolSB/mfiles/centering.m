function mu = centering(dx,dz,ds,dw,ap,ad,dgap,trerror)
% CENTERING   - Computing centering parameter mu.
%             Usage: mu = centering(dx,dz,ds,dw,sn1,wn1,dgap)

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global x y z s w
global Ubounds_exist nt

newdgap = (x + min(1,ap)*dx)'*(z + min(1,ad)*dz);
if (Ubounds_exist) 
   newdgap = newdgap + (s + min(1,ap)*ds)'*(w + min(1,ad)*dw); 
end
sigmak = (newdgap/dgap)^2;

sigmin = 0; sigmax = .208; % Don't ask why.
p = ceil(log10(trerror)); 
if (p < -2 & dgap < 1.e+3) sigmax = 10^(p+1); end;
sigmak = max(sigmin, min(sigmax, sigmak));

mu = sigmak*dgap/nt;
