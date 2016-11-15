% implementation of a double-limit cycle oscillator

clear all;close all;clc;

for z0 = [0.45 0.46 0.99]
    % parameters of the gamma oscillator
    % z0 = 0.1;
    alpha = -100;
    beta1 = 500;
    beta2 = -50;
    epsilon = 1;
    f1 = 100;
    
    % parameters for time
    fs = 100000;
    dur = 0.1; % in seconds
    T = 1/fs;
    time = 0:T:dur;
    
    ntime = length(time);
    
    % initialize empty arrays
    t_out = zeros(size(time));
    z_out_1 = zeros(size(time));
    
    % initialization values
    t_out(1) = time(1);
    z_out_1(1) = z0;
    
    input = zeros(size(time));
    
    %     for i=1:ntime-1
    %
    %         %%% first layer %%%
    %
    %         k1 = z_out_1(i)*(alpha + 1i*2*pi*f1 + beta1*abs(z_out_1(i))^2 + ...
    %             ((epsilon*beta2*abs(z_out_1(i)))^4)/(1-epsilon*abs(z_out_1(i))^2)) + input(i);
    %
    %         z1 = z_out_1(i)+(T/2)*k1;
    %
    %         k2 = z1*(alpha + 1i*2*pi*f1 + beta1*abs(z1)^2 + ...
    %             ((epsilon*beta2*abs(z1))^4)/(1-epsilon*abs(z1)^2)) + (input(i) + input(i+1))/2;
    %
    %         z2 = z_out_1(i)+(T/2)*k2;
    %
    %         k3 = z2*(alpha + 1i*2*pi*f1 + beta1*abs(z2)^2 + ...
    %             ((epsilon*beta2*abs(z2))^4)/(1-epsilon*abs(z2)^2)) + (input(i) + input(i+1))/2;
    %
    %         z3 = z_out_1(i)+T*k3;
    %
    %         k4 = z3*(alpha + 1i*2*pi*f1 + beta1*abs(z3)^2 + ...
    %             ((epsilon*beta2*abs(z3))^4)/(1-epsilon*abs(z3)^2)) + (input(i));
    %
    %         z_out_1(i+1) = z_out_1(i) + (1/6)*T*(k1+(2*k2)+(2*k3)+k4);
    %
    %         t_out(i+1) = t_out(i) + T;
    %
    %     end
    %
    % the output of beta drives a gamma oscillator
    dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
        z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); % input to beta
    
    % integrate the ode to obtain the amplitude
    [t_out_2,z_out_2] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);
    
    figure(1)
    subplot(2,2,1)
    %     plot((z_out_1))
    hold on
    plot(z_out_2)
    axis square
    title('Oscillator 1')
    grid on
    subplot(2,2,2)
    %     plot((abs(z_out_1(1:end-1))),diff(abs(z_out_1)));
    axis square
    hold on
    plot((abs(z_out_2(1:end-1))),diff(abs(z_out_2)));
    title('Oscillator 1')
    grid on
    subplot(2,2,3:4)
    %     plot(time,abs(z_out_1))
    grid on
    hold on
    plot(time,abs(z_out_2))
    title('Magnitude')
    
    figure(2)
    plot3(time,real(z_out_2),imag(z_out_2))
    hold on
    grid on
    
end

%%
% Here the same double limit cycle oscillator is forced with external input

% implementation of a double-limit cycle oscillator

clear all;close all;clc;

% parameters of the gamma oscillator
z0 = 0.2;
alpha = -100;
beta1 = 500;
beta2 = -50;
epsilon = 1;
f1 = 10;

% parameters for time
fs = 100000;
dur = 0.1; % in seconds
T = 1/fs;
time = 0:T:dur;

ntime = length(time);

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); % input to beta

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);

figure(1)
subplot(2,2,1)
plot((z_out_1))
hold on
axis square
title('Oscillator 1')
grid on
subplot(2,2,2)
plot((abs(z_out_1(1:end-1))),diff(abs(z_out_1)));
axis square
hold on
title('Oscillator 1')
grid on
subplot(2,2,3:4)
plot(time,abs(z_out_1))
grid on
hold on
title('Magnitude')

figure(2)
plot3(time,real(z_out_1),imag(z_out_1))
hold on
grid on

n_impulses = 4;

for i=1:n_impulses
    
    % new oscillator parameter
    z0 = z_out_1(end) + 0.40;
    
    % parameters for time
    fs = 100000;
    dur = 0.031; % in seconds
    T = 1/fs;
    time = t_out_1(end):T:t_out_1(end)+dur;
    
    ntime = length(time);
    
    dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
        z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); % input to beta
    
    % integrate the ode to obtain the amplitude
    [t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);
    
    figure(1)
    subplot(2,2,1)
    plot((z_out_1))
    hold on
    axis square
    title('Oscillator 1')
    grid on
    subplot(2,2,2)
    plot((abs(z_out_1(1:end-1))),diff(abs(z_out_1)));
    axis square
    hold on
    title('Oscillator 1')
    grid on
    subplot(2,2,3:4)
    plot(time,abs(z_out_1))
    grid on
    hold on
    title('Magnitude')
    
    figure(2)
    plot3(time,real(z_out_1),imag(z_out_1))
    hold on
    grid on
    
    
end

% new oscillator parameter
z0 = z_out_1(end);

% parameters for time
fs = 100000;
dur = 0.1; % in seconds
T = 1/fs;
time = t_out_1(end):T:t_out_1(end)+dur;

ntime = length(time);

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); % input to beta

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);

figure(1)
subplot(2,2,1)
plot((z_out_1))
hold on
axis square
title('Oscillator 1')
grid on
subplot(2,2,2)
plot((abs(z_out_1(1:end-1))),diff(abs(z_out_1)));
axis square
hold on
title('Oscillator 1')
grid on
subplot(2,2,3:4)
plot(time,abs(z_out_1))
grid on
hold on
title('Magnitude')

figure(2)
plot3(time,real(z_out_1),imag(z_out_1))
hold on
grid on