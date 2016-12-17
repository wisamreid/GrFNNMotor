function  phase_plane(f1, f2, ax, ay, bx, by)

% This function plots the phase portrait of an arbitrary 
% 2D nonlinear dynamical system, determines the nullclines 
% and plots the velocity field
%
% Input f1, f2 -> function handles defining the system
%                 f1(x,y,mu) f2(x,y,mu) 
%
%       ax, bx -> these define the rectangle 
%       ay, by      
%                         --------- b=(bx, by)
%                        |         |
%                        |         |
%                         ---------
%                   a=(ax,ay)


%  hereafter we define two different 2D grids
% X1,Y1 represent the initial condition for the 
% trajectories. We also use them to plot the vector field
Nx1=20; 
Ny1=20;
xx1=linspace(ax,bx,Nx1);
yy1=linspace(ay,by,Ny1);
[X1,Y1]=meshgrid(xx1,yy1);

% We use X,Y for computing the nullclines and the trajectories
Nx=200; 
Ny=200;
xx=linspace(ax,bx,Nx);
yy=linspace(ay,by,Ny);
[X,Y]=meshgrid(xx,yy);

fontsize=16;
figure(1)
clf
hold
box on


[C1,h1]=contour(X,Y,f1(X,Y),[0 0],'r','Linewidth',1.3); % nullcline dx/dt=0
[C2,h2]=contour(X,Y,f2(X,Y),[0 0],'k','Linewidth',1.3); % nullcline dy/dt=0
leg=legend('$\dot{x}=0$','$\dot{y}=0$');
set(leg,'Interpreter','Latex')


% Here we plot the trajectories 
p  = haltonset(2,'Skip',1e3,'Leap',1e2);
p  = scramble(p,'RR2');
X0 = net(p,200);
X0(:,1)= ax + (bx-ax)*X0(:,1);
X0(:,2)= ay + (by-ay)*X0(:,2);

streamline(X,Y,f1(X,Y),f2(X,Y),X0(:,1),X0(:,2),[5e-3,1e4])



xlabel('$x$','Interpreter','Latex','Fontsize',fontsize)
ylabel('$y$','Interpreter','Latex','Fontsize',fontsize)
set(gca,'Fontsize',fontsize)
axis tight

% Here we plot the velocity field 
quiver(X1,Y1,f1(X1,Y1),f2(X1,Y1),'k','Linewidth',1.2)











end

