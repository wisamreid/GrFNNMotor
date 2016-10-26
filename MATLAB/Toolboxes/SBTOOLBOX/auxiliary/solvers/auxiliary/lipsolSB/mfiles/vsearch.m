function j = vsearch(string, strvec, n_vars)

j = [];
for i = 1:n_vars
   s = length(find(strvec(i,:) ~= ' ' & strvec(i,:) ~= 0));
   if strcmp(string,strvec(i,1:s)) j = i; break; end
end
