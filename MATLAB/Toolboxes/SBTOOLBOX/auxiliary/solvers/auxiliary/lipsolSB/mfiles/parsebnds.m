function [substr,lr,j] = parsebnds(str,v,vn_vars)

% do what?
str = str(find(str ~=' ' & str ~= 0 ));
lr = [1 1];
gtorlt = find(str == '>' | str == '<');
s_gt = find(str == '>');

if(length(gtorlt) == 2)
  substr = str(gtorlt(1)+1:gtorlt(2)-1); 
else 
  if((str(1) >= 'A' & str(1) <= 'Z') ...
    |(str(1) >= 'a' & str(1) <= 'z'))  % (x > * or x < *)
    substr = str(1:gtorlt(1)-1);
    if(length(s_gt) == 1) lr(2) = 0; else lr(1) = 0; end
  else                                 % (* > x or * < x)
    substr = str(gtorlt(1)+1:length(str));
    if(length(s_gt) == 1) lr(1) = 0; else lr(2) = 0; end
  end
end
j = vsearch(substr, v, vn_vars);
