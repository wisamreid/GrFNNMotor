%%
clear all;close all;clc;

% parameters of the network
freqs = [0.5,1,2,4,8,16,32,64];
nfreqs = length(freqs);

% parameters of the oscillator
r0 = 0.01;
alpha = 1;
beta1 = -1;
beta2 = -1;
epsilon = 0.11;
psi0 = 0.5*pi; % initial phase of the oscillator

% parameters for the input
F = 1;
f = 1;
w0 = 2*pi*f;
    
% parameters for time
fs = 1000;
end_t = 50; % in seconds
T = 1/fs;
time = 0:T:end_t;

for i=1:nfreqs
    
    f_scale = freqs(i);
    w = 2*pi*f_scale;
    omega = w - w0;
    
    % solve r_dot
    full_sys = @(t,z,alpha,beta1,beta2,epsilon,omega,F,f_scale)  ...
        [z(1)*alpha + beta1*z(1)^3 + (epsilon*beta2*z(1)^5)/(1-epsilon*z(1)^2) + F*cos(z(2));
        omega/f_scale - (F/z(1))*sin(z(2))];
    
    % integrate the ode to obtain the amplitude
    [time,z] = ode45(@(t,z) full_sys(t,z,alpha,beta1,beta2,epsilon,omega,F,f_scale),time,[r0,psi0]);
       
    figure(1);
    subplot(1,nfreqs,i)
    plot(time,z(:,1))
    title(sprintf('%d Hz',f_scale))
    ylim([0 1.3])
    grid on
    xlabel('Time (s)')
    ylabel('Magnitude')    
    figure(2)    
    polar(pi,1.3)
    hold on
    polar(z(:,2),z(:,1))
    hold on
    grid on
    pause
    
end