function aMax = findMaxTimeUnit(asdf)

m = [];
for i=1:N
    m(end+1) = max(asdf{i});
end
aMax = max(m);
