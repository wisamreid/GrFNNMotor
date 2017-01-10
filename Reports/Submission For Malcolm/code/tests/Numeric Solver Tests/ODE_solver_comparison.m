% Author: Wisam Reid 
%
% Comparison of ODE solvers on a system of equations.
% 1. Euler
% 2. Runge-Kutta
% 3. ODE23
% 4. ODE45
% 5. ODE113
% 
%% CLEAN AND CLEAR

clc;
clear;
close all;

%% vanderpol

% yprime = @(t,y) vanderpoldemo(t,y,1); 
yprime = @(t,y) vanderpoldemo(t,y,0.4);

%% Intital Condition and Variables

tspan = [0 20];
y0 = [2; 0]; % intital condition
% time step
h = 0.1; 
 
%% Plots

for i=1:length(y0)
    
    figure(i)
    clf; 
    hold on 
    grid on
    xlabel('t')
    ylabel('y(t)')
    title('Comparison of ODE solvers on a system of equations')

    [teuler,yeuler] = euler(yprime, tspan, y0, h);
    plot(teuler,yeuler(i,:),'r.-','markersize',10);

    [trunge,yrunge] = runge_kutta_4order(yprime, tspan, y0, h);
    plot(trunge,yrunge(i,:),'g.-','markersize',10);
    
    [tode23,yode23] = ode23(yprime, tspan, y0);
    plot(tode23,yode23(:,i),'c.-','markersize',10);
    
    [tode45,yode45] = ode45(yprime, tspan, y0);
    plot(tode45,yode45(:,i),'b.-','markersize',10);
 
    [tode113,yode113] = ode113(yprime, tspan, y0);
    plot(tode113,yode113(:,i),'c.-','markersize',10);
    
    legend('Euler','Runge-Kutta','ODE45', 'ODE23', 'ODE113')
end

%%
