function [y,ts] = test_solvers_pendulum()


NSTEPS = 100000;
IOSTEP = 100;
DT     = 1e-4;

y0=[ pi/2+pi/4 0]';  % We drop the pendulum from theta0=pi/4 with zero velocity.

%[y,ts] = Euler_forward(@ Rhs_pendulum,y0,DT,NSTEPS,IOSTEP);
%[y,ts] = Heun_method(@ Rhs_pendulum,y0,DT,NSTEPS,IOSTEP);
[y,ts] = RK4_method(@ Rhs_pendulum,y0,DT,NSTEPS,IOSTEP);
%[y,ts] = AB3_method(@ Rhs_pendulum,y0,DT,NSTEPS,IOSTEP);

l = 9.8;  % The pendulum is 9.8 meters (32 feet)

figure(2)
for i=1:length(ts)
clf
hold
plot(l*sin(y(1,i)),-l*cos(y(1,i)),'ko','Markersize',12,'MarkerFaceColor',[0 0 0])
line([0 l*sin(y(1,i))],[0 -l*cos(y(1,i))],'Linewidth', 2,'Color',[0 0 0])

axis equal 
axis([-12 12 -12 12])
M(i)=getframe;

end

movie(M,1,60);

end

