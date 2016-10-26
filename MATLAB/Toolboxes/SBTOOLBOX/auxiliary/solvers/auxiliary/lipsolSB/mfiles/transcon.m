function [A,b,variables,n_vars] = ...
                      transcon(A,b,record,b2e,cont,variables,n_vars)

for (iii = 1:length(b2e))
  A_cont = b2e(iii);
  A(A_cont,n_vars) = 0;
  b(A_cont) = 0;
  nonzero = find(record(A_cont+4,:) ~= ' ' & ...
             record(A_cont+4,:) ~= 0);
  str = record(A_cont+4,nonzero);


  t = find(str == ')');
  str = str(t(1)+1:length(str));

    if(str(1) ~= '+' & str(1) ~= '-') str = ['+',str] ; end 

    ls = find(str == '{');
    rs = find(str == '}');
    tmpstr = str;
    for(i = 1:length(ls))
       for j = ls(i):rs(i) tmpstr(j) = 0; end
    end

    s = find(tmpstr == '-' | tmpstr == '+' | tmpstr == '=' | ...
             tmpstr == '<' | tmpstr == '>');
    s0 = find(tmpstr == '=' | tmpstr == '>' | tmpstr == '<');
    j0 = find(s == s0);
    s = [s, length(str)+1];
   for i= 2:length(s)
     substr = str(s(i-1):s(i)-1);
     if( i == j0+1) 
       if(length(substr) == 1) substr(1) = '0'; else substr(1) = '+'; end
     end

     coef = find(substr == '{' | substr == '}');
     if(isempty(coef))
       j = min(find(substr >= 'a' & substr <= 'z'));
       if(~isempty(j))
         if(j == 2) 
           substr = [substr(1),'1', substr(2:length(substr))] ;
           j = j+1 ;
         end
       end
       if(isempty(j)) 
         number = str2num(substr(1:length(substr)));
       else 
         num_str = substr(1:j-1);
         num_str = num_str(find(num_str ~= '*'));
         number = str2num(num_str) ;
         string = substr(j:length(substr));
       end
     end
     if(~isempty(coef))
       j = [];
       number1 = str2num([substr(1),'1']);
       number2 = eval(substr(coef(1)+1:coef(2)-1) );
       number = number1*number2;
       if(length(substr) > coef(2))
         j = coef(2)+ 1;
         string = substr(coef(2)+1:length(substr));
       end
     end
     if(isempty(j)) 
       if(i <= j0) number = -number; end
       b(A_cont) = number + b(A_cont);
     else
       if(i > j0) number = -number; end
       jj = vsearch(string, variables, n_vars);
       if(~isempty(jj)) 
         A(A_cont,jj) = number + A(A_cont,jj);
       else
         n_vars = n_vars+1;
         A(A_cont,n_vars) = number;
         variables(n_vars,1:length(string)) = string;
       end
     end
   end
   if(~isempty(find(str == '<')))
     n_vars = n_vars+1;
     string = ['$slack'];
     variables(n_vars, 1:length(string)) = string ;
     A(A_cont,n_vars) = 1;
   elseif(~isempty(find(str == '>')))
     n_vars = n_vars+1;
     string = ['$slack'];
     variables(n_vars, 1:length(string)) = string ;
     A(A_cont,n_vars) = -1;
   end
end    
A = A(1:cont,1:n_vars);
b = b(1:cont);
if size(b,1) < size(b,2) b = b'; end;
variables = variables(1:n_vars,:);
