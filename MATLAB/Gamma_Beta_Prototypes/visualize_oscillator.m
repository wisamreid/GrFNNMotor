% This script shows a limit cycle oscillator in the beta band (20Hz)
% by Irán Román
% implementation of the theory in Kim and Large 2015

clear all
close all
clc

% differential equation for the amplitude of a hopf bifurcation
F = @(t,r,alpha,beta1,beta2,epsilon) r*alpha + beta1*r^3 + (epsilon*beta2*r^5)/(1-epsilon*r^2);

% parameters of the limit cycle oscillator
r0 = 0.9;
alpha = 1;
beta1 = -1;
beta2 = -10;
epsilon = 1;
f = 1;

% parameters for time
fs = 1000;
dur = 1; % in seconds
T = 1/fs;
time = 0:T:dur;

% integrate the ode to obtain the amplitude
r = ode4(F,time,r0,alpha,beta1,beta2,epsilon);
drdt = diff(r);

% visualize the amplitude of the oscillator over time
figure(1)
subplot(2,1,1)
plot(time,r)
ylabel('$A.U.$','Interpreter','Latex')
xlabel('$time(s)$','Interpreter','Latex')
grid on
subplot(2,1,2)
stem((r(1:end-1)),drdt)
xlabel('$r$','Interpreter','Latex')
ylabel('$\frac{dr}{dt}$','Interpreter','Latex')
grid on

% to calculate z we need omega
omega = ones(size(r))*2*pi*f;

% calculate z
z = r.*exp(1i.*omega.*time');
dzdt = diff(z);

% plot z
figure('units','normalized','position',[1 0.5 1 0.5])
subplot(1,3,1)
plot(time,real(z))
axis square
ylabel('$A.U.$','Interpreter','Latex')
xlabel('$time(s)$','Interpreter','Latex')
subplot(1,3,2)
plot(z)
axis square
xlabel('real')
ylabel('imaginary')
subplot(1,3,3)
% plot the phase portrait
plot((abs(z(1:end-1))),abs((dzdt)))
axis square
xlabel('$z$','Interpreter','Latex')
ylabel('$\frac{dz}{dt}$','Interpreter','Latex')

%%
% clear all
% close all
% clc
% % UNCOMMENT to try the same with a saddle node
% Fx = @(t,x,r) r - x^2;
% Fy = @(t,y) -y;
% 
% % bifurcation parameter
% r = 0.1;
% 
% for x0 = -1:0.1:1
%     for y0 = -1:0.1:1
%         
%         % parameters for time
%         fs = 100;
%         dur = 10; % in seconds
%         T = 1/fs;
%         time = 0:T:dur;
%         
%         % integrate the ode to obtain the amplitude
%         x = ode4(Fx,time,x0,r);
%         y = ode4(Fy,time,y0);
%         
%         
%         % plot the phase portrait
%         plot(x,y)
%         hold on
%         grid on       
%         axis square
%         axis([-1 1 -1 1])
%         xlabel('$x$','Interpreter','Latex')
%         ylabel('$y$','Interpreter','Latex')
%         
%     end
% end