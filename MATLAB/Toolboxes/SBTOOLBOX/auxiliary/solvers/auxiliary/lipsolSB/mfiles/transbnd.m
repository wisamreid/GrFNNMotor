function [lb,ub,variables,n_vars] = ...
                  transbnd(lb,ub,record,b2e,n_cons,variables,n_vars,BIG)
for (iii = 1:length(b2e))
  b_cont = b2e(iii);
  len = find(record(b_cont+n_cons+5,:) ~= ' ' & ...
             record(b_cont+n_cons+5,:) ~= 0);
  str = record(b_cont+n_cons+5,len);

  t = find(str == ')');
  str = str(t(1)+1:length(str));

  s_gt = find(str == '>');
  s_lt = find(str == '<');
  lcoef = find(str == '{');
  rcoef = find(str ==  '}');
  lencoef = length(lcoef);

  for(i = 1:lencoef)
     coef_number(i) = eval(str(lcoef(i)+1:rcoef(i)-1));
   if(lcoef(i) > 1)
     if(str(lcoef(i)-1) == '-')
       coef_number(i) = -coef_number(i); 
     end
   end
  end


  if(length(s_lt)== 2)
     str2 = str(s_lt(1)+1:s_lt(2)-1);
     j = vsearch(str2, variables, n_vars);
     if(~isempty(j))
       if(lencoef == 0) 
         lb(j) = str2num(str(1:s_lt(1)-1));
         ub(j) = str2num(str(s_lt(2)+1:length(str)));
       elseif(lencoef == 1)
         lb(j) = coef_number(1);
         ub(j) = str2num(str(s_lt(2)+1:length(str)));
       elseif(lencoef == 2)
         lb(j) = coef_number(1);
         ub(j) = coef_number(2);
       end
     end
  elseif(length(s_gt) == 2)
     str2 = str(s_gt(1)+1:s_gt(2)-1);
     j = vsearch(str2, variables, n_vars);
     if(~isempty(j))
       if(lencoef == 0) 
         ub(j) = str2num(str(1:s_gt(1)-1));
         lb(j) = str2num(str(s_gt(2)+1:length(str)));
       elseif(lencoef == 1)
         ub(j) = coef_number(1);
         lb(j) = str2num(str(s_gt(2)+1:length(str)));
       elseif(lencoef == 2)
         ub(j) = coef_number(1);
         lb(j) = coef_number(2);
       end
     end
  elseif(length(s_lt) == 1)
    if(lencoef == 0)
      if((str(1) >= '0' & str(1) <= '9') ...
        | str(1)  == '+' | str(1)  == '-')
        str2 = str(s_lt(1)+1:length(str));
        j = vsearch(str2, variables, n_vars);
        if(~isempty(j)) lb(j) = str2num(str(1:s_lt(1)-1)); end 
      else
        str2 = str(1:s_lt(1)-1);
        j = vsearch(str2, variables, n_vars);
        if(~isempty(j)) ub(j) = str2num(str(s_lt(1)+1:length(str))); end 
      end
    else
      if(lcoef(1)< s_lt(1))
        str2 = str(s_lt(1)+1:length(str));
        j = vsearch(str2, variables, n_vars);
        if(~isempty(j)) lb(j) = coef_number(1); end
      else
        str2 = str(1:s_lt(1)-1);
        j = vsearch(str2, variables, n_vars);
        if(~isempty(j)) ub(j) = coef_number(1); end
      end
    end
  elseif(length(s_gt) == 1)
    if(lencoef == 0)
      if((str(1) >= '0' & str(1) <= '9') ...
        | str(1)  == '+' | str(1)  == '-')
        str2 = str(s_gt(1)+1:length(str));
        j = vsearch(str2, variables, n_vars);
        if(~isempty(j)) ub(j) = str2num(str(1:s_gt(1)-1)); end
      else
        str2 = str(1:s_gt(1)-1);
        j = vsearch(str2, variables, n_vars);
        if(~isempty(j)) lb(j) = str2num(str(s_gt(1)+1:length(str))); end
      end
    else
      if(lcoef(1)< s_gt(1))
        str2 = str(s_gt(1)+1:length(str));
        j = vsearch(str2, variables, n_vars);
        if(~isempty(j)) ub(j) = coef_number(1);end
      else
        str2 = str(1:s_gt(1)-1);
        j = vsearch(str2, variables, n_vars);
        if(~isempty(j)) lb(j) = coef_number(1); end
      end
    end
  end
end
%if(n_vars > length(lb))
%lb(n_vars,1) = 0
%ub(n_vars,1) = 0;
%end

