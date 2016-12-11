% coupled oscillators
% 
% Author: Wisam Reid 
%
%% CLEAN AND CLEAR

clc;
clear;
close all;

%% equations to solve

% two oscillators with recurrent coupling
dzdt = @(t,z,w1,alpha1,beta11,beta12,epsilon1,C12,  ... % osc 1 params
             w2,alpha2,beta21,beta22,epsilon2,C21,  ... % osc 2 params
             F,input) ... % input 
    ... % osc 1
    [z(1)*(alpha1 + 1i*w1 + beta11*abs(z(1))^2 + ...     
    (epsilon1*beta12*abs(z(1))^4)/(1-epsilon1*abs(z(1))^2)) + ...
    C21*z(2) + ... % coupling from osc 2 
    F*input(t); ...   % input
    ... % osc 2
    z(2)*(alpha2 + 1i*w2 + beta21*abs(z(2))^2 + ...
    (epsilon2*beta22*abs(z(2))^4)/(1-epsilon2*abs(z(2))^2)) + ...
    C12*z(1)]; % coupling from osc 1

%% Oscillator parameters

% a critical Hopf regime: 
%   a = 0, b1 < 0, b2 = 0
% a supercritical Hopf regime: 
%   a > 0, b1 < 0, b2 = 0
% a super critical double limit cycle regime:
%   a < 0, b1 > 0, b2 < 0, localmax > 0
% a subcritical double limit cycle regime: 
%   a < 0, b1 > 0, b2 < 0, local max < 0

% coupling parameters 
C12 = 1;
C21 = -10;

% oscillator 1 parameters
z01 = 0.0001;
alpha1 = -5;
beta11 = 0.1;
beta12 = 0;
epsilon1 = 1;
f1 = 79.5;

% oscillator 2 parameters
z02 = 0.0001;
alpha2 = 0.01;
beta21 = -5;
beta22 = 0;
epsilon2 = 1;
f2 = 20;

% oscillator tuning
w1 = 2*pi*f1; % osc 1
w2 = 2*pi*f2; % osc 2

%% time parameters

fs = 1000;
dur = 30; % in seconds
T = 1/fs;
time = 0:T:dur;

%% parameters for input

f0 = 0.25;
omega0 = 2*pi*f0;
F = 1;
D = 0 : 1/f0 : dur;

% input 
% input = @(t) 0;
% input = @(t) exp(1i*omega0*t);
input = @(t) 0.5*square(omega0*t,2) + 0.5;  
% input = @(t) pulstran(t,D,'tripuls',0.1,-1);

%% integrate the ode to obtain the amplitude
[t,z] = ode45(@(t,z) dzdt(t,z,w1,alpha1,beta11,beta12,epsilon1,C12, ...
                              w2,alpha2,beta21,beta22,epsilon2,C21, ...
                              F,input),time,[z01,z02]);

%% Plot z
% figure(1)
% plot(z)

%% Plot Oscillators

% extract information
r = abs(z);
dr2dt = diff(r(:,1));
dr1dt = diff(r(:,2));

figure(2)
subplot(4,1,1)
plot(t,r)
title('Magnitude of Oscillators')
ylabel('Amplitude')
xlabel('time(s)')
legend(['Osc 1: ', num2str(f1),' Hz'],['Osc 2: ', num2str(f2), ' Hz'])
subplot(4,1,2)
plot(t,real(z))
title('Real part of Oscillators')
legend('Osc 1','Osc 2')
ylabel('Amplitude')
xlabel('time(s)')
subplot(4,1,3)
plot(t,input(t))
title('Input')
% legend('Beta','Gamma')
ylabel('Amplitude')
xlabel('time(s)')
subplot(4,1,4)
spectrogram((z(:,1)+z(:,2)),kaiser(256,5),220,20000,fs,'yaxis')
title('Spectrum')
ylim([0 100])
xlabel('time(s)')
% title(['$\w_1 = $', num2str(w1),             ...
%        '$~\w_2 = $', num2str(w2),            ...
%        '$~\alpha_1 = $', num2str(alpha1),    ...
%        '$~\alpha_2 = $', num2str(alpha2),    ...
%        '$~\beta_{11} = $', num2str(beta11),  ...
%        '$~\beta_{21} = $', num2str(beta21),  ...
%        '$~\beta_{12} = $', num2str(beta12),  ...
%        '$~\beta_{22} = $', num2str(beta22)], ...
%         'Interpreter','latex')
%% Extra 

figure(3)
plot3(t,real(z),angle(z))
xlabel('time')
ylabel('$\mathcal{R}$','Interpreter','latex')
zlabel('\textit{Im}','Interpreter','latex')

% figure(4)
% plot(r(1:end-1,2),dr1dt)
% hold on
% plot(r(1:end-1,1),dr2dt)

%% Save params 
% lopes

% % coupling parameters 
% C12 = 1;
% C21 = -1;
% 
% % oscillator 1 parameters
% z01 = 0.0001;
% alpha1 = 0;
% beta11 = -0.1;
% beta12 = 0;
% epsilon1 = 1;
% f1 = 40;
% 
% % oscillator 2 parameters
% z02 = 0.0001;
% alpha2 = 0;
% beta21 = -0.1;
% beta22 = 0;
% epsilon2 = 1;
% f2 = 40;
