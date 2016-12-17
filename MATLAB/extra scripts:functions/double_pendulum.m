function  double_pendulum()
% This function simulate the dynamics of a duuble pendulum
% by using RK4 method

global m1 m2 l1 l2 g


NSTEPS = 300000;
IOSTEP = 100;
DT     = 1e-4;


m1 = 3.2;   % 1 kg (kilograms)
m2 = 1.3; 
l1 = 0.6; % m (meters)
l2 = 1.3;   % meter

g=9.8; % acceleration of gravity m/s^2

% Initial condition

y0 = [pi/2 3/4*pi 0 0]';

[y,ts] = RK4_method(@ Rhs_double_pendulum,y0,DT,NSTEPS,IOSTEP);



figure(1)
clf
plot(ts,y(1,:),'r')
hold
plot(ts,y(2,:))

figure(2)
for i=1:length(ts)
clf
hold
x1=l1*sin(y(1,i));
y1=-l1*cos(y(1,i));
x2=l1*sin(y(1,i))+l2*sin(y(2,i));
y2=-l1*cos(y(1,i))-l2*cos(y(2,i));

plot(x1,y1,'ko','Markersize',12,'MarkerFaceColor',[0 0 0])
line([0 x1],[0 y1],'Linewidth', 2,'Color',[0 0 0])

plot(x2,y2, ... 
     'ko','Markersize',6,'MarkerFaceColor',[1 0 0])

line([x1 x2 ],[y1 y2],'Linewidth', 2,'Color',[1 0 0])

% Trajectory
plot(l1*sin(y(1,1:i))+l2*sin(y(2,1:i)),-l1*cos(y(1,1:i))-l2*cos(y(2,1:i)),'b');

axis equal 
axis([-3 3 -3 3])
M(i)=getframe;

end

%movie2avi(M,'double_pendulum.avi','fps',50)
%movie(M,1,60);

end

