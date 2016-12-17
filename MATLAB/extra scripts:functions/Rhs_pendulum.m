function F = Rhs_pendulum(y,t)
% This function characterizes the right hand side of the nonlinear 
% prndulum equations written as a first-order system
%
% Input y -> vector collecting the position y(1) and the momentum y(2) of the
%            pendulum
%       l -> length of the pendulum  [m]
%       g -> acceleration of gravity  9.8 m/s^2 

F=zeros(2,1);

F(1) = y(2);
F(2) = -sin(y(1)) -0.01* y(2); 


end

