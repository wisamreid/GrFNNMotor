% This script shows a limit cycle oscillator in the beta band (20Hz)

clear all
close all
clc
% Include the GrFNN toolbox
addpath(genpath('~/Documents/GrFNNToolbox/'))

% differential equation for the amplitude of a hopf bifurcation
% %%% no input %%%
% dzdt = @(t,z,alpha,beta1,beta2,epsilon) z*(alpha + 1i*2*pi + ...
%     beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
% %%% with input %%%
dzdt = @(t,z,alpha,beta1,beta2,epsilon,F,omega0) z*(alpha + 1i*2*pi*20 + ...
    beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)) + ...
    F*exp(1i*omega0*t);

% parameters of the limit cycle oscillator
z0 = 0.3;
alpha = 1;
beta1 = -1;
beta2 = -1;
epsilon = 1;
F = 0.7;
f0 = 20.5;
omega0 = 2*pi*f0;

% parameters for time
fs = 1000;
dur = 20; % in seconds
T = 1/fs;
time = 0:T:dur;

% integrate the ode to obtain the amplitude
% %%% no input %%%
% z = ode4(dzdt,time,z0,alpha,beta1,beta2,epsilon);
% %%% with input %%%
z = ode4(dzdt,time,z0,alpha,beta1,beta2,epsilon,F,omega0);

% extract inormation
r = abs(z);
drdt = diff(r);
dzdt = diff(z);

% visualize the amplitude of the oscillator over time
figure(1)
subplot(2,1,1)
plot(time,r)
ylabel('$A.U.$','Interpreter','Latex')
xlabel('$time(s)$','Interpreter','Latex')
grid on
subplot(2,1,2)
plot((r(1:end-1)),drdt)
xlabel('$r$','Interpreter','Latex')
ylabel('$\frac{dr}{dt}$','Interpreter','Latex')
grid on

% plot z
figure('units','normalized','position',[1 0.5 1 0.5])
subplot(1,3,1)
plot(time,real(z))
hold on
plot(time,imag(z))
title('real part of z')
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
plot((z(1:end-1)),(dzdt))
axis square
xlabel('$z$','Interpreter','Latex')
ylabel('$\frac{dz}{dt}$','Interpreter','Latex')

% uncomment for a quick and dirty fft of z
figure(3)
freqs = linspace(-fs/2,fs/2,2^(nextpow2(length(z))));
plot(freqs,fftshift((fft(z,2^(nextpow2(length(z)))))))
axis tight