% Replicating Figure 1: Kim and Large 2015
% by Wisam Reid
% 
% Abandoning this script for now

%% CLEAN AND CLEAR

clc;
clear;
close all;

%% equations to solve

% two oscillators with recurrent coupling
drdt = @(t,r,alpha1,beta11,beta12,epsilon1, ... % osc 1 params
             alpha2,beta21,beta22,epsilon2, ... % osc 2 params
             alpha3,beta31,beta32,epsilon3, ... % osc 3 params
             alpha4,beta41,beta42,epsilon4) ... % osc 4 params 
    ... % osc 1
    [alpha1*r(1) + beta11*r(1)^3 + (epsilon1*beta12*r(1)^5)/(1-epsilon1*r(1)^2); ...
    ... % osc 2
     alpha2*r(2) + beta21*r(2)^3 + (epsilon2*beta22*r(2)^5)/(1-epsilon2*r(2)^2);
    ... % osc 3
     alpha3*r(3) + beta31*r(3)^3 + (epsilon3*beta32*r(3)^5)/(1-epsilon3*r(3)^2);
    ... % osc 4
     alpha4*r(4) + beta41*r(4)^3 + (epsilon4*beta42*r(4)^5)/(1-epsilon4*r(4)^2)];
 
%% Oscillator parameters

% oscillator 1 parameters: a critical Hopf regime
%   a = 0, b1 < 0, b2 = 0
r01 = 0.2;
alpha1 = 0;
beta11 = -1.0;
beta12 = 0;
epsilon1 = 0;

% oscillator 2 parameters: a supercritical Hopf regime 
%   a > 0, b1 < 0, b2 = 0
r02 = 0.2;
alpha2 = 1.0;
beta21 = -1.0;
beta22 = 0;
epsilon2 = 0;

% oscillator 3 parameters: a super critical double limit cycle regime
%   a < 0, b1 > 0, b2 < 0, localmax > 0
r03 = 1;
alpha3 = -1;
beta31 = 1;
beta32 = -1;
epsilon3 = 1;

% oscillator 4 parameters: a subcritical double limit cycle regime: 
%   a < 0, b1 > 0, b2 < 0, local max < 0
r04 = 1;
alpha4 = -1;
beta41 = 1;
beta42 = -10;
epsilon4 = 1;

%% time parameters

fs = 5000;
dur = 10; % in seconds
T = 1/fs;
time = 0:T:dur;

%% Integrate the ode to obtain the amplitude

[t,r] = ode45(@(t,r) ...
        drdt(t,r,alpha1,beta11,beta12,epsilon1,  ... % osc 1 params
                 alpha2,beta21,beta22,epsilon2,  ... % osc 2 params
                 alpha3,beta31,beta32,epsilon3,  ... % osc 3 params
                 alpha4,beta41,beta42,epsilon4), ... % osc 4 params
                 time,[r01,r02,r03,r04]);

%% grab the derivative

drdt = diff(r);

%% Plot Responses

figure(1)
% critical Hopf
plot(r(1:end-1,1),drdt,'-k','LineWidth',2)
hold on
% super critical Hopf
plot(r(1:end-1,2),drdt,'-r','LineWidth',2)
hold on
% super critical double limit cycle
plot(r(1:end-1,3),drdt,'-g','LineWidth',2)
hold on
% double limit cycle (subcritical)
plot(r(1:end-1,4),drdt,'-b','LineWidth',2)
hold on
line(xlim,[0 0],'Color','c','LineStyle','--') 
title('Autonomous behavior in different parameter regimes')
xlabel('$r$','Interpreter','Latex')
ylabel('$\frac{dr}{dt}$','Interpreter','Latex')
legend('Critical Hopf','Super Critical Hopf','SC Double Limit Cycle', ...
    'Sub-Critical DLC','Location', 'South')