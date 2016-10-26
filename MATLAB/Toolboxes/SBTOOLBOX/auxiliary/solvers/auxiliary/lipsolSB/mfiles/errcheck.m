function err=errcheck(str, flag)

global v lv rv  vn_vars  

err = 0;
check_obj    = 1;
check_cons   = 2;
check_bounds = 3;

if (isempty(str)) err = 1; return; end

tmpstr = str(find(str ~= ' ' & str ~= 0));
operator = find(tmpstr == '=' | tmpstr == '>' | tmpstr == '<' );
if (~isempty(operator))
if (operator(1) == 1 | operator(length(operator)) == length(tmpstr))
  err = errmessage('Input incomplete'); return;
end
end

operator = find(tmpstr == '+' | tmpstr == '-' );
if (~isempty(operator))
if ( operator(length(operator)) == length(tmpstr))
  err = errmessage('Input incomplete'); return;
end
end


% coefficient check ......

lbrace = find(str == '{');
rbrace = find(str == '}');
len = length(lbrace);
if ( len ~= length(rbrace) ) 
  errmessage('Braces unmatched'); return;
end
if (len ~= 0)
  tmpstr =['+' tmpstr 'z'];
  tmplbrace = find(tmpstr == '{');
  tmprbrace = find(tmpstr == '}') ;
  patt = [];
  for i = 1:len
    if (~isletter(tmpstr(tmprbrace(i)+1)) &...
         isempty(findstr(tmpstr(tmprbrace(i)+1),'+-<>='))) |...
         isempty(findstr(tmpstr(tmplbrace(i)-1),'+-<>='))
      errstr = 'A coefficient must be a number or an expression in "{}"';
      err = errmessage(errstr);
      return;
    end
    if (rbrace(i)-lbrace(i) <2) 
      errmessage('Empty {}'); return;
    end
    inner = lbrace(i):rbrace(i);
    patt = [patt inner];
    lbracket = find(str(inner) == '(');
    rbracket = find(str(inner) == ')');
    if (length(lbracket) ~= length(rbracket)) 
      err = errmessage(' "(" dose not match ")" ');  return; 
    end
  end
  str(patt) = [];
end

% general checking ......

error = ~( ...
          str == '+' | str == '-' | str == '.' | str == ' ' ...
        | str ==  0  | str == '=' | str == '<' | str == '>' ...
        | (str >= '0' & str <= '9') | isletter(str));
if any(error) 
   errstr = str(find(error));
   err = errmessage(['Illegal charactor(s): ' errstr]); return; end

operator = (find(str == '+' | str == '-' ));
len = length(operator);
if (len > 1)
   error = any(operator(2:len) - operator(1:len-1) == 1);
  if any(error)
    err = errmessage('Successive operators not allowed'); return; 
  end
end

% objective check ......

if (flag == check_obj)
  error = (str == '=' | str == '>' | str == '<'); 
  if any(error) 
    err = errmessage('"=", "<" or ">" not allowed in objective'); 
    return; 
  end

  nonblank = find(str ~= ' ');
  s3 = lower(str(nonblank(1):nonblank(1)+2));
  if any(s3 ~= 'max') & any(s3 ~= 'min')  
    err = errmessage('Max or min must appear first in objective'); 
    return;
  else
     str=str(nonblank(1)+3:length(str));
  end
end
     
% constraints checking ......

if (flag == check_cons)
  operator = find(str == '=' | str == '>' | str == '<');
  if (length(operator) ~= 1)
    err = errmessage('Use one and only one of "=", "<" or ">"'); 
    return;
  end
end

% bounds checking ......

if (flag == check_bounds)
  if (~isempty(find(str == '=')))
    err = errmessage('NO "=" allowed in bounds'); return;
  end
  operator = find(str == '>' | str == '<' );
  len_op = length(find(str == '>' | str == '<' ));
  if (len_op < 1 | len_op > 2)
    err = errmessage('There must be 1 or 2 inequalities in a bound'); 
    return;
  elseif (len_op == 2)
    gt = length(find(str == '>'));
    lt = length(find(str == '<'));
    if (~(gt == 2 | lt == 2))
      err = errmessage('Misused ">" or "<"');
      return;
    end
  end
  
  interval = [0 operator length(str)+1];
  for i = 1:length(interval)-1
    substr = str(interval(i)+1:interval(i+1)-1);
    if (isempty(substr)) substr = '9'; end
    substr = substr(find(substr ~= ' '));
    subletter = find((substr >= 'A' & substr <= 'Z') | ...
                (substr >= 'a' & substr <= 'z') );
    if(~isempty(substr))
      subop = find(substr =='+' | substr =='-' | substr =='*');
    else
      subop = [];
    end
    if (~isempty(subletter))
      if (subletter(1) ~= 1)
        err = errmessage('Coefficients on variables not allowed'); 
        return;
      end
      if (~isempty(subop))
        err = errmessage('Arithmatic operators on variables not allowed'); 
        return;
      end
    end
  end

  letter = find((str >= 'A' & str <= 'Z') | (str >= 'a' & str <= 'z'));
  len_letter = length(letter);
  if (len_letter == 0) err = errmessage('No variable'); return; end

  error = (letter(1) < operator(1) & letter(len_letter) > operator(1)) ...
        | (letter(1) < operator(len_op) & ...
          letter(len_letter) > operator(len_op)) ...
        | (len_op == 2 & letter(1) < operator(1)) ...
        | (len_op == 2 & letter(len_letter) > operator(len_op)); 
  if any(error) 
    err = errmessage('Illegal expression in bounds'); 
    return; 
  end
end

% blank space check ......

nonblank = str ~= ' ';
slen = length(nonblank);
nonblank = [nonblank 1]; 
for i = 1:slen
  if (nonblank(i) == 0 & nonblank(i+1) == 1) nonblank(i) = 1 ; end
end
nonblank = nonblank(1:slen);
str = str(find(nonblank == 1));

str = ['@' str '@'];
blank = find(str == ' ');
if (~isempty(blank))
  b_blank = str(blank -1);
  a_blank = str(blank +1);
  error =  ( (b_blank >= '0' & b_blank <= '9') & ...
              (a_blank >= '0' & a_blank <= '9') ) ...
        | ( isletter(b_blank) & isletter(a_blank) ) ...
        | ( isletter(b_blank) & (a_blank >= '0' & a_blank <= '9') );
if any(error) err = errmessage('Space in wrong place'); return; end 
end
