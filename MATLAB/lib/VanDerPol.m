function xp = VanDerPol(t,x)
% Duffing-VanDerPol equations expressed as a first order system
% (see equations 3-5). Parameter values (mu, f and omega)
% can be chosen to demonstrate various trajectories(see text). 
%
% IN:
% t:time(not used, but necessary for ODEsolver)
% x:Inputvector 
%
% OUT:
% xp:dx/dt
% Current parameters demonstrate chaotic behavior of the oscillator 

mu=0.2; 
f=1.0; 
omega=0.94;

% define equations
xp(1)=x(2);
xp(2)=mu*(1-x(1)^2)*x(2)-x(1)^3+f*cos(x(3));
xp(3)=omega;
% transpose into column vector for ODEsolver
xp = xp';