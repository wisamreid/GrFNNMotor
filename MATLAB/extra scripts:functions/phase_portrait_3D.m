function  phase_portrait_3D(f1, f2, f3, ax, ay, az, bx, by, bz)

% This function plots the phase portrait of an arbitrary 
% 3D nonlinear dynamical system, determines the nullclines,
% plots the velocity field and the trajectories
%
% Input f1, f2, f3 -> function handles defining the system
%                     f1(x,y,z)  f2(x,y,z)  f3(x,y,z)
%
%                            _______  (bx, by,bz)
%                          /        /| 
%                         ---------  |
%                        |         | |
%                        |         |/
%                         ---------
%                    (ax,ay,az)


% This grid is for plotting the vector field
Nx1=10; 
Ny1=10;
Nz1=10;
xx1=linspace(ax,bx,Nx1);
yy1=linspace(ay,by,Ny1);
zz1=linspace(az,bz,Nz1);
[X1,Y1,Z1]=meshgrid(xx1,yy1,zz1);

% This grid is computing/plotting nullclines and trajectories
Nx=100; 
Ny=100;
Nz=100;
xx=linspace(ax,bx,Nx);
yy=linspace(ay,by,Ny);
zz=linspace(az,bz,Nz);
[X,Y,Z]=meshgrid(xx,yy,zz);

fontsize=16;
figure(1)
clf
hold
grid

% Here we plot the nullclines (surfaces in 3D)

p1 = patch(isosurface(X, Y, Z, f1(X,Y,Z), 0));
p2 = patch(isosurface(X, Y, Z, f2(X,Y,Z), 0));
p3 = patch(isosurface(X, Y, Z, f3(X,Y,Z), 0));

isonormals(X,Y,Z,f1(X,Y,Z),p1)
set(p1, 'FaceColor', 'red', 'EdgeColor', 'none','FaceAlpha',0.5);
isonormals(X,Y,Z,f2(X,Y,Z),p2)
set(p2, 'FaceColor', 'blue', 'EdgeColor', 'none','FaceAlpha',0.5);
isonormals(X,Y,Z,f3(X,Y,Z),p3)
set(p3, 'FaceColor', 'green', 'EdgeColor', 'none','FaceAlpha',0.5);

leg=legend('$\dot{x}=0$','$\dot{y}=0$','$\dot{z}=0$');
set(leg,'Interpreter','Latex')

% Initial condition for the trajectories 
p  = haltonset(3,'Skip',1e3,'Leap',1e2);
p  = scramble(p,'RR2');
X0 = net(p,300);
X0(:,1)= ax + (bx-ax)*X0(:,1);
X0(:,2)= ay + (by-ay)*X0(:,2);
X0(:,3)= az + (bz-az)*X0(:,3);

streamline(X,Y,Z,f1(X,Y,Z),f2(X,Y,Z),f3(X,Y,Z),X0(:,1),X0(:,2),X0(:,3),[5e-3,1e4])


xlabel('$x$','Interpreter','Latex','Fontsize',fontsize)
ylabel('$y$','Interpreter','Latex','Fontsize',fontsize)
zlabel('$z$','Interpreter','Latex','Fontsize',fontsize)
view(-30,30)
set(gca,'Fontsize',fontsize)
axis tight

% Here we plot the velocity field 
quiver3(X1,Y1,Z1,f1(X1,Y1,Z1),f2(X1,Y1,Z1),f3(X1,Y1,Z1),'k','Linewidth',1)











end

