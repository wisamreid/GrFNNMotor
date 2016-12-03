%%
clear all;close all;clc;

% parameters of the network
freqs = [0.5, 1, 2, 4, 8];
nfreqs = length(freqs);

% parameters of the oscillator
r0 = 0.1;
alpha = 1;
beta1 = -1;
beta2 = -1;
epsilon = 1;
w = 2*pi;

% parameters for the input
F = 1;
f = 1;
w0 = 2*pi*f;
omega = w - w0;
    
psi0 = 0.5*pi;

% parameters for time
fs = 1000;
end_t = 10; % in seconds
T = 1/fs;
time = 0:T:end_t;

for i=1:nfreqs
    
    f_scale = freqs(i);
    
    % solve r_dot
    full_sys = @(t,z,alpha,beta1,beta2,epsilon,w,omega,F,w0,f_scale)  ...
        [z(1)*alpha + beta1*z(1)^3 + (epsilon*beta2*z(1)^5)/(1-epsilon*z(1)^2) + F*cos(z(2));
        (omega/f_scale - (F/z(1))*sin(z(2)))];
    
    % integrate the ode to obtain the amplitude
    [time,z] = ode45(@(t,z) full_sys(t,z,alpha,beta1,beta2,epsilon,w,omega,F,w0,f_scale),time,[r0,psi0]);
    
    figure(1)
    subplot(3,nfreqs,i)
    plot(time,(1/f_scale)*z(:,1))
    title(sprintf('%d Hz',f_scale))
    ylim([0 1])
    grid on
    xlabel('Time (s)')
    ylabel('Magnitude')
    subplot(3,nfreqs,nfreqs/1+i)
    plot(time,(1/f_scale)*z(:,2))
    grid on
    xlabel('Time (s)')
    ylabel('Relative Phase')
    
%     figure(2)
%     polar(z(:,2),z(:,1))
%     hold on
%     
%     figure(3)
%     plot3(time,real(z(:,1).*exp(1i*z(:,2))),imag(z(:,1).*exp(1i*z(:,2))))
%     grid on
    
end
