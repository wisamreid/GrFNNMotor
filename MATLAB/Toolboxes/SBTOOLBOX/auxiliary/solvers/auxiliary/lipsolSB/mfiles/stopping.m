function [stop, converged] = stopping(tol)
%

% Yin Zhang, April, 1995
% Department of Mathematics and Statistics
% University of Maryland  Baltimore County

global Hist message
stop = 0; converged = 0;

iter = size(Hist,2);
trerror = Hist(1,iter);
if (trerror < tol) 
   stop = 1; converged = 1; 
   message = '<* Converged! *>'; return; 
end;

small = 1.e-5; if (iter <= 2) return; end;
blowup = Hist(1:5,iter) > max(tol,min(Hist(1:5,:)')')/small;
if iter > 5
   blowup = blowup | Hist(1,iter) > 1000*min(Hist(1,:)) | ...
   all(Hist(1,iter-4:iter) > 1.01*(Hist(1,iter-5:iter-1)));
end
if any(blowup)
   stop = 1; converged = 0;
   message = detectinf(tol);
end;

nstall = 5; if (iter <= nstall) return; end;
latest = iter-nstall:iter;
h = Hist(1:5,latest); hmean = mean(h');

for i = 1:5
  stall(i) = all(abs(h(i,:)-hmean(i)) < 1.e-3*hmean(i));
  if i > 1 stall(i) = stall(i) & hmean(i) > small; end;
end
if any(stall)
   if trerror < small
      stop = 1; converged = 1;
      message = 'Likely converged but error > tolerance'; 
   else
      stop = 1; converged = 0;
      message = detectinf(tol);
   end
end;
