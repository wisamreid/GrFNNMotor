function [err,lv,rv] = checkbnds(str,lv,rv)
global variables n_vars

err = 0;
% extract the order number ......
t = find(str == ')');
str = str(t(1)+1:length(str));
str = str(find(str ~= ' ' | str ~= 0));
[substr,lr,j] = parsebnds(str,variables,n_vars);
if(isempty(j))
  err = errmessage('The bounds on this line make no sense');
  return;
end

if(lr(1)*lv(j) == 1 & lr(2)*rv(j) == 1)
   err = errmessage('Both bounds already set'); return; 
elseif(lr(1)*lv(j) == 1)
   err = errmessage('Lower bound already set'); return; 
elseif(lr(2)*rv(j) == 1)
   err = errmessage('Upper bound already set'); return; 
else
   lv(j) = lr(1);
   rv(j) = lr(2);
end
