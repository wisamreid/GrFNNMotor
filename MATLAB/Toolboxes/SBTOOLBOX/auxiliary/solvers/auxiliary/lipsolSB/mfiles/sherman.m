function x = sherman(P,U,b,flag)
% CHERMAN    Using the Sherman-Morrison formula to solve
%                 (P + U*U')x = b
%            where P is sparse and U is low-rank.

% Yin Zhang, April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global Mdense

if (flag == 1)
   nd = size(U,2); 
   PinvU = zeros(size(U));
   for i = 1:nd
       PinvU(:,i) = linsys(P,U(:,i)); 
   end;
   Mdense = eye(nd) + U'*PinvU; 
else
   P = [];
end;

x = linsys(P,b); 
tmp = U*(Mdense\(U'*x));
x = x - linsys(P,tmp);
