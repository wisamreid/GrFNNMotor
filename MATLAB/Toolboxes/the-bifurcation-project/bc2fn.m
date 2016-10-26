function [pl,ql,pr,qr]=bc2fn(xl,ul,xr,ur,t)

% Constants
A = 2;
B = 4; %NOT GIVEN | TO BE CHANGED

% Case 1:- Zero Flux Boundary Conditions (Neumann conditions)
pl= [0;0];
ql=[1;1];
pr =[0;0];
qr =[1;1];

% Case 2:- Fixed Boundary Conditions (Dirichlet conditions)
% pl = [A; B/A];
% ql = [0; 0];
% pr = [A; B/A];
% qr = [0; 0];