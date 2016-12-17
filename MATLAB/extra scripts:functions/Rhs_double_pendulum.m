function F = Rhs_double_pendulum(y,t)

% This function represents the right hand side of the double pendulum 
% equations

% y(1) -> theta1
% y(2) -> theta2
% y(3) -> p1
% y(4) -> p2

global m1 m2 l1 l2 g

F=zeros(4,1);  % y(1) the

K1= y(3)*y(4)*sin(y(1)-y(2))/(l1*l2*(m1+m2*sin(y(1)-y(2))^2)); % fine

K2= ( l2^2*m2*y(3)^2+l1^2*(m1+m2)*y(4)^2-l1*l2*m2*y(3)*y(4)*cos(y(1)-y(2)))/ ...
   (2*l1^2*l2^2*(m1+m2*sin(y(1)-y(2))^2))^2*sin(2*(y(1)-y(2))); 

F(1) = (l2*y(3)-l1*y(4)*cos(y(1)-y(2)))/(l1^2*l2*(m1+m2*sin(y(1)-y(2))^2));
F(2) = (l1*(m1+m2)*y(4)-l2*m2*y(3)*cos(y(1)-y(2)))/ ... 
       (l1*l2^2*m2*(m1+m2*sin(y(1)-y(2))^2));
F(3) = -(m2+m1)*g*l1*sin(y(1))-K1+K2;
F(4) = -m2*g*l2*sin(y(2))+K1-K2;

end

