close all
clear all
clc

% This is a linear model of oscillators in the beta band, gating activity
% in the gamma band.

fs = 100;
dur = 10;
time = 0:1/fs:dur;
f = 1; % in Hz
omega = 2*pi*f;

r = 0.00000001;
alpha = -3;
% beta_1 = 1;
beta_1 = 0;
beta_2 = -0.002;
% epsilon = 0.0000001;
epsilon = 0;

r_dot = ones(size(time))*r;
phi_dot = ones(size(time))*omega;
z_dot = zeros(size(time));
% this will be a cannonical linear oscllator

% evaluate z_dot with a linear oscillator
for t=1:length(time)-1
    
    r_dot(t+1) = alpha.*r_dot(t) + beta_1.*r_dot(t).^3 + ...
        (epsilon.*beta_2.*r_dot(t).^5)./(1 - epsilon.*r_dot(t).^2);
    
%     r_dot(t+1)
    
    z_dot(t+1) = r_dot(t+1).*exp(1i.*phi_dot(t)) ...
        + r_dot(t).*exp(1i.*phi_dot(t))*1i.*phi_dot(t+1);
    
        z_dot(t+1)
    
    if isnan(z_dot(t+1))
        t
        time(t)
        error('A Nan has been evaluated');
    end
    
    plot(r_dot(t),r_dot(t+1),'o')
    hold on
    grid on
    axis([0 r -0.3 0.3])
    
    if rem(time(t),1) == 0.1 && time(t) > 0
        pause
        
    end
    
    %
    %     if rem(time(t),1) == 0 && time(t) > 0
    %         subplot(1,2,1)
    %         plot(time(1:end-1),real(z_dot(2:end)))
    %         grid on
    %         subplot(1,2,2)
    %         plot(z_dot(2:t))
    %         grid on
    %         pause;
    %         clf;
    %     end
end

subplot(1,2,1)
plot(time(1:end-1),abs(z_dot(2:end)))
grid on
subplot(1,2,2)
plot(z_dot(2:end))
grid on