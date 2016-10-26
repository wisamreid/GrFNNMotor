function [dx,dy,dz,ds,dw] = ...
         direction(A,P,U,Rb,Rc,Ru,Rxz,Rsw,vmid,xn1,sn1,z,w,mu,bnrm,flag)
% DIRECTION   - Computing search directions

% Yin Zhang, January, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

%global x y z s w
global Ubounds_exist
global Dense_cols_exist

if (mu ~= 0) Rxz = Rxz - mu; end;
Rc0 = Rc;
Rc = Rc - Rxz.*xn1;
if (Ubounds_exist)
   if (mu ~= 0) Rsw = Rsw - mu; end;
   Rc = Rc + (Rsw - Ru.*w).*sn1;
end;
rhs = -(Rb + A*(vmid.*Rc));

if Dense_cols_exist
   dy = densol(P,U,full(rhs),flag);
else
   dy = spsol(P,full(rhs),flag);
end;
dx = vmid.*(A'*dy + Rc);
dz = -(z.*dx + Rxz).*xn1;
ds = []; dw = [];
if (Ubounds_exist)
   ds = -(dx.*spones(w) + Ru);
   dw = -(w.*ds + Rsw).*sn1;
end;

resp = Rb + A*dx;
if norm(resp) > 1.e-1*bnrm
   dx = dx - A'*solaat(full(resp));
end;
