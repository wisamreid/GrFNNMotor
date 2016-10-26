function [NAME,c,variables,n_vars,BIG] = transobj(record,c,variables,n_vars)
BIG = 1.0e32;

str = record(1,:); str = str(find(str ~= ' '));
spint = find(str == ':');
NAME = deblank(str(1,spint+1:length(str)));

  str = deblank(record(3,find(record(3,:) ~= ' ')));
  if lower(str(1:3)) == 'max' MAX = 1; else MAX = 0; end
  str = str(4:length(str));

  if(str(1) ~= '+' & str(1) ~= '-') str = ['+',str] ; end 
  s = [ find(str == '-' | str == '+') length(str)+1];

  for i= 2:length(s)
     substr = str(s(i-1):s(i)-1);
     sublen = length(substr);

     coef = find(substr == '{' | substr == '}');
     if(isempty(coef))
       j = min(find(substr >= 'a' & substr <= 'z'));
       if(~isempty(j))
         if(j == 2) 
          substr = [substr(1),'1', substr(2:sublen)]; 
          j = j+1 ;
         end
         num_str = substr(1:j-1);
         num_str = num_str(find(num_str ~= '*'));
         number = str2num(num_str);
         string = substr(j:length(substr));
       end
     end
     if(~isempty(coef))
       j = 0;
       if(sublen > coef(2))
         j = coef(2)+ 1;
         number1 = str2num([substr(1),'1']);
         number2 = eval(substr(coef(1)+1:coef(2)-1));
         number = number1*number2; 
         string = substr(coef(2)+1:length(substr));
       end
     end
     if( j ~= 0)
         if MAX number = - number; end
         jj = vsearch(string, variables, n_vars);
         if(~isempty(jj)) 
           c(jj) = number + c(jj);
         else
           n_vars = n_vars+1;
           c(n_vars) = number;
           variables(n_vars,1:length(string)) = string;
         end
     end
  end
if size(c,1) < size(c,2) c = c'; end;
