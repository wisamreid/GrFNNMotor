function [Rxz,Rsw,dgap] = update(ap,ad,dx,dy,dz,ds,dw,trerror,tol)
% UPDATE      - Update the iterates.
%            UPDATE uses backtracking to make iterates stay in an
%            infinity-neighborhood of the central path.

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global x y z s w
global phi0 tau0
global backed nt

tau = tau0;
if ~exist('backed') backed = 0; end; 
if ~backed tau = .9 + 0.1*tau0; end;
k = ceil(log10(trerror)); 
if (k <= -5) tau = max(tau,1-10^k); end;

ap = min(tau*ap,1); ad = min(tau*ad,1); 
xc = x; yc = y; zc = z; sc = s; wc = w;

step = [1 .9975 .95 .90 .75 .50];
for k = 1:length(step)
   x = xc + step(k)*ap*dx;
   y = yc + step(k)*ad*dy;
   z = zc + step(k)*ad*dz;
   s = sc + step(k)*ap*ds;
   w = wc + step(k)*ad*dw;
   [Rxz,Rsw,dgap] = complementy(x,z,s,w);
   phi = nt*full(min([Rxz; Rsw(find(Rsw))]))/dgap;
   if (max(ap,ad) == 1) | (phi >= phi0) break; end;
end;
phi0 = min(phi0, phi); if (k > 1) backed = 1; end;
