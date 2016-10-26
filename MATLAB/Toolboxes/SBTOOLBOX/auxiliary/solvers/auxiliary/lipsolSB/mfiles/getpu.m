function [P,U] = getpu(A,vmid)
% GETMISC     - Compute Matrix P and U
%             Usage: [P,U] = getpu(A,vmid)

% Yin Zhang, April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global Dense_cols_exist idense ispars
[m,n] = size(A);

if Dense_cols_exist
   ns = length(ispars); nd = length(idense); 
   Ds = sparse(1:ns,1:ns,vmid(ispars),ns,ns,ns);
   P = A(:,ispars)*Ds*A(:,ispars)';
   Dd = sparse(1:nd,1:nd,sqrt(vmid(idense)),nd,nd,nd);
   U = full(A(:,idense)*Dd);
else
   P = A*sparse(1:n,1:n,vmid,n,n,n)*A'; U = [];
end;
