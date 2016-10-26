function [VAL] = extrapolbrentSB(dist,Z,L)

A = L*(L-dist(2)) / (dist(1)*(dist(1)+dist(2)));
B = (L+dist(1))*(dist(2)-L) / (dist(1)*dist(2));
C = L*(L+dist(1)) / (dist(2)*(dist(1)+dist(2)));
VAL = B*Z(:,1) + A*Z(:,2) + C*Z(:,3);

return