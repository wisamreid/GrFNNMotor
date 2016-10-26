function [c,f,s] = pdex4pde(r,t,u,DuDr)

% Diffusion Coefficients
Dy = 1.6*10^(-3);
Dx = 8.0*10^(-3);

% Constants
A = 2;
L = 1;
u1 = u(1);
u2 = u(2);

B = 4; %NOT GIVEN | TO BE CHANGED

% c = [1; 1];
c = [0; 0];
f = [Dx; Dy] .* DuDr;

% Rate equations describing the phenomenon
s1 = A + u1^2*u2 - (B+1)*u1;
s2 = B*u1 - u1^2*u2;

% Linearized equations for the perturbation x and y
% s1 = (B-1)*u1 + A^2*u2;
% s2 = -B*u1 - A^2*u2;

s = [s1; s2];