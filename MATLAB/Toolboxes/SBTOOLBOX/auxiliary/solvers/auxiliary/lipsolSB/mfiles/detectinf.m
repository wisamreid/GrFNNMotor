function message = infeas(tol)
% 

% Yin Zhang, April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global Hist message
msgncov = '<* Not converged *>';
msgpinf = [msgncov ' Primal program infeasible'];
msgdinf = [msgncov ' Dual program infeasible'];

minrp = min(Hist(2,:)+Hist(4,:));
minrd = min(Hist(3,:)); 

if minrp < tol message = msgdinf; return; end;
if minrd < tol message = msgpinf; return; end;

if minrp < 10*tol & minrd > sqrt(tol)
   message = [msgdinf '?']; return;
end;
if minrd < 10*tol & minrp > sqrt(tol)
   message = [msgpinf '?']; return;
end;

iter = size(Hist,2);
if Hist(6,iter) < -1.e+10 & Hist(7,iter) < 1.e+6
   message = [msgdinf '??']; return;
end;
if Hist(7,iter) > 1.e+10 & Hist(6,iter) > -1.e+6
   message = [msgpinf '??']; return;
end;
message = [msgncov ' Both programs infeasible???'];
