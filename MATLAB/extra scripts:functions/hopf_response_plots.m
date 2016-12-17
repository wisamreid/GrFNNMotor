% Hopf Responses      
% 
% Author: Wisam Reid 
%
%% CLEAN AND CLEAR

clc;
clear;
close all;

%% equations to solve

% % two oscillators with recurrent coupling
% dzdt = @(t,z,w1,alpha1,beta11,beta12,epsilon1, ... % osc 1 params
%              w2,alpha2,beta21,beta22,epsilon2, ... % osc 2 params
%              w3,alpha3,beta31,beta32,epsilon3, ... % osc 3 params
%              w4,alpha4,beta41,beta42,epsilon4, ... % osc 4 params
%              F,step) ... % input 
%     ... % osc 1
%     [z(1)*(alpha1 + 1i*w1 + beta11*abs(z(1))^2 + ...     
%      (epsilon1*beta12*abs(z(1))^4)/(1-epsilon1*abs(z(1))^2)) + F*step(t); ...
%     ... % osc 2
%      z(2)*(alpha2 + 1i*w2 + beta21*abs(z(2))^2 + ...
%      (epsilon2*beta22*abs(z(2))^4)/(1-epsilon2*abs(z(2))^2)) + F*step(t); ...
%     ... % osc 3
%      z(3)*(alpha3 + 1i*w3 + beta31*abs(z(3))^2 + ...
%      (epsilon3*beta32*abs(z(3))^4)/(1-epsilon3*abs(z(3))^2)) + F*step(t); ...
%     ... % osc 4
%      z(4)*(alpha4 + 1i*w4 + beta41*abs(z(4))^2 + ...
%      (epsilon4*beta42*abs(z(4))^4)/(1-epsilon4*abs(z(4))^2)) + F*step(t)];

%% Oscillator parameters

% a critical Hopf regime: 
%   a = 0, b1 < 0, b2 = 0
% a supercritical Hopf regime: 
%   a > 0, b1 < 0, b2 = 0
% a super critical double limit cycle regime:
%   a < 0, b1 > 0, b2 < 0, localmax > 0
% a subcritical double limit cycle regime: 
%   a < 0, b1 > 0, b2 < 0, local max < 0

% % oscillator 1 parameters
% % linear response
% z01 = 0.0001;
% alpha1 = -1;
% beta11 = 0;
% beta12 = 0;
% epsilon1 = 0;
% f1 = 1;

% oscillator 2 parameters
% critical Hopf
z02 = 0.0001;
alpha2 = 0;
beta21 = -10;
beta22 = 0;
epsilon2 = 1;
f2 = 1;
% 
% % oscillator 3 parameters
% % limit cycle
% z03 = 0.0001;
% alpha3 = -5;
% beta31 = 10;
% beta32 = 0;
% epsilon3 = 1;
% f3 = 1;
% 
% % oscillator 4 parameters
% % double limit cycle (subcritical)
% z04 = 0.0001;
% alpha4 = 0.01;
% beta41 = -5;
% beta42 = 0;
% epsilon4 = 1;
% f4 = 1;
%  
% oscillator tuning
% w1 = 2*pi*f1; % osc 1
% w2 = 2*pi*f2; % osc 2
% w3 = 2*pi*f3; % osc 1
% w4 = 2*pi*f4; % osc 2


% % oscillator parameters
% % linear response
% alpha = -1;
% beta1 = 0;
% beta2 = 0;
% epsilon = 0;
% f = 1;

% [stim_strengths, stead_state_amp] = plotSteadyState(f1,alpha1,beta11,beta12,epsilon1,0.1);
[stim_strengths, stead_state_amp] = plotSteadyState(f2,alpha2,beta21,beta22,epsilon2,0.1);

% [stim_strengths, stead_state_amp] = plotSteadyState(f,alpha,beta1,beta2,epsilon,0.1);


figure(1)
% plot result
plot(stim_strengths,stead_state_amp,'-k','LineWidth',2)
axis([0 1 0 1])


%% time parameters
% 
% fs = 1000;
% dur = 10; % in seconds
% T = 1/fs;
% time = 0:T:dur;

%% parameters for input

% f0 = 0.25;
% omega0 = 2*pi*f0;
% F = 1;

% ramp input 
% ramp = @(t) max(0,t);
% step = @(t) heaviside(t);
% test = @(t) exp(1i*omega0*t);
% plot(time,ramp(time))
 

%% Integrate the ode to obtain the amplitude

% [t,z] = ode45(@(t,z) ...
%         dzdt(t,z,w1,alpha1,beta11,beta12,epsilon1, ... % osc 1 params
%                  w2,alpha2,beta21,beta22,epsilon2, ... % osc 2 params
%                  w3,alpha3,beta31,beta32,epsilon3, ... % osc 3 params
%                  w4,alpha4,beta41,beta42,epsilon4, ... % osc 4 params
%                  F,step),time,[z01,z02,z03,z04]);

%% grab the amplitudes

% r = abs(z);
% 
% %% Plot Responses
% 
% figure(1)
% % linear response
% plot(t,r(:,1),'-k','LineWidth',2)
% hold on
% % % critical Hopf
% % plot(t,r(:,2),'-r','LineWidth',2)
% % hold on
% % % limit cycle
% % plot(t,r(:,3),'-g','LineWidth',2)
% % hold on
% % % double limit cycle (subcritical)
% % plot(t,r(:,4),'-b','LineWidth',2)
% hold on
% % plot(t,ramp(t),'--m')
% title('Intrinsic Dynamics')
% xlabel('Stimulus Strength')
% ylabel('Oscillator Amplitude, r')
% legend('Linear','Critical Hopf','Limit Cycle','Double Limit Cycle')

