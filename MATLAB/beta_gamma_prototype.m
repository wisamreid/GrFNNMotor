% This script shows a limit cycle oscillator in the beta band (20Hz)

clear all
close all
clc
% Include the GrFNN toolbox
addpath(genpath('~/Documents/GrFNNToolbox/'))

% the output of beta drives a gamma oscillator
dzdt = @(t,z,alpha,beta1,beta2,epsilon,F,omega0)  ...        
    [z(1)*(100*alpha + 1i*2*pi*40 + 200*beta1*abs(z(1))^2 + ...     
    (epsilon*beta2*abs(z(1))^4)/(1-epsilon*abs(z(1))^2)) + 3.5*(z(2)*exp(pi)); ...
    (z(2)*(alpha + 1i*2*pi*20 + beta1*abs(z(2))^2 + ... % beta layer
    (epsilon*beta2*abs(z(2))^4)/(1-epsilon*abs(z(2))^2)) + ...   
    F*exp(1i*omega0*t))]; % input to beta

% parameters of the limit cycle oscillator
z0 = 0.3;
alpha = 1;
beta1 = -15;
beta2 = -1;
epsilon = 1;
F = 0.5;
f0 = 20.5;
omega0 = 2*pi*f0;

% parameters for time
fs = 1000;
dur = 10; % in seconds
T = 1/fs;
time = 0:T:dur;

% integrate the ode to obtain the amplitude
[t,z] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,F,omega0),time,[z0,z0]);

% extract information
r = abs(z);
dr2dt = diff(r(:,1));
dr1dt = diff(r(:,2));

figure(2)
plot(z)
% plot(r(1:end-1,2),dr1dt)
% hold on
% plot(r(1:end-1,1),dr2dt)

figure(3)
subplot(3,1,1)
plot(time,r)
title('Magnitude of Oscillators')
ylabel('Amplitude')
xlabel('time(s)')
legend('Gamma','Beta')
subplot(3,1,2)
plot(time,real(z))
title('Real part of Oscillators')
legend('Gamma','Beta')
ylabel('Amplitude')
xlabel('time(s)')
subplot(3,1,3)
spectrogram((z(:,1)),kaiser(256,5),220,20000,fs,'yaxis')
title('Spectrum')
ylim([0 100])
xlabel('time(s)')