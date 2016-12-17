%%

clear all;close all;clc;

animate = 1; % true or false

% parameters of the gamma oscillator
z0 = 0.1;
alpha = 1;
beta1 = -100;
beta2 = 0;
epsilon = 0;
f1 = 1;

% parameters for the input
F = 0.2;
f0 = 1.55;
w0 = 2*pi*f0;

% parameters for time
fs = 1000;
end_t = 100; % in seconds
T = 1/fs;
time = 0:T:end_t;

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1,F,w0)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + ...
    (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)) + F*exp(1i*w0*t);

% integrate the ode to obtain the amplitude
[time,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1,F,w0),time,z0);

z_all = [];
time_all = [];
z_all = [z_all z_out_1'];
time_all = [time_all time'];

figure(1)
subplot(3,2,1:2)
plot(time_all,abs(z_all))
grid on
xlabel('Time (s)')
xlabel('Magnitude')

if animate == 1
    
    for i=1:length(time_all)
        
        if mod(i,fs/20) == 0
            subplot(3,2,3:6)
            plot3(time_all(1:i),real(z_all(1:i)),imag(z_all(1:i)),'-')
            grid on
            xlabel('Time (s)')
            ylabel('Real part')
            zlabel('Imaginary part')
            pause(0.1)
            
        end
    end
    
else
    
    figure(2)
    plot3(time_all,real(z_all),imag(z_all))
    hold on
    grid on
    xlabel('Time (s)')
    ylabel('Real part')
    zlabel('Imaginary part')
    
end

%%
clear all;close all;clc;

% parameters of the oscillator
r0 = 0.1;
alpha = 0;
beta1 = -100;
beta2 = 0;
epsilon = 0;
f = 1;
w = f*2*pi;

% parameters for the input
F = 0.2;
f0 = 0.9;
w0 = 2*pi*f0;

omega = w - w0;
psi0 = 0;

% parameters for time
fs = 1000;
end_t = 100; % in seconds
T = 1/fs;
time = 0:T:end_t;

% solve r_dot
full_sys = @(t,z,alpha,beta1,beta2,epsilon,w,omega,F,w0)  ...
    [z(1)*alpha + beta1*z(1)^3 + (epsilon*beta2*z(1)^5)/(1-epsilon*z(1)^2) + F*cos(z(2));
    omega - (F/z(1))*sin(z(2))];   

% integrate the ode to obtain the amplitude
[time,r] = ode45(@(t,z) full_sys(t,z,alpha,beta1,beta2,epsilon,w,omega,F,w0),time,[r0,psi0]);

figure(1)
subplot(2,1,1)
plot(time,r(:,1))
grid on
xlabel('Time (s)')
ylabel('Magnitude')
subplot(2,1,2)
plot(time,r(:,2))
grid on
xlabel('Time (s)')
ylabel('Relative Phase')

figure(2)
polar(r(:,2),r(:,1))
hold on

figure(3)
plot3(time,r(:,1).*real(exp(1i*r(:,2))),r(:,1).*imag(exp(1i*r(:,2))))
% plot3(time(fs*5:end),r(fs*5:end,1).*real(exp(1i*r(fs*5:end,2))),r(fs*5:end,1).*imag(exp(1i*r(fs*5:end,2))))
grid on

%%
clear all;close all;clc;

animate = 1; % true or false

% parameters of the gamma oscillator
r0 = 0.05;
alpha = 1;
beta1 = -100;
beta2 = 0;
epsilon = 0;
f = 1;
w = f*2*pi;

% parameters for the input
F = 0.2;
f0 = 0.3;
w0 = 2*pi*f0;

omega = w - w0;
psi0 = 0.9*pi;

% parameters for time
fs = 1000;
end_t = 10; % in seconds
T = 1/fs;
time = 0:T:end_t;

% solve r_dot
full_sys = @(t,z,alpha,beta1,beta2,epsilon,w,omega,F,w0)  ...
    [z(1)*alpha + beta1*z(1)^3 + (epsilon*beta2*z(1)^5)/(1-epsilon*z(1)^2) + F*cos(z(2));
    omega - (F/z(1))*sin(z(2))];   

% integrate the ode to obtain the amplitude
[time,r] = ode45(@(t,z) full_sys(t,z,alpha,beta1,beta2,epsilon,w,omega,F,w0),time,[r0,psi0]);

figure(1)
subplot(2,1,1)
plot(time,r(:,1))
grid on
xlabel('Time (s)')
ylabel('Magnitude')
subplot(2,1,2)
plot(time,r(:,2))
grid on
xlabel('Time (s)')
ylabel('Relative Phase')

figure(2)
polar(r(:,2),r(:,1))
hold on

figure(3)
plot3(time,real(r(:,1).*exp(1i*r(:,2))),imag(r(:,1).*exp(1i*r(:,2))))
grid on

%%
clear all;close all;clc;

% parameters of the gamma oscillator
r0 = 0.9;
alpha = -1;
beta1 = 4;
beta2 = -1;
epsilon = 1;
f = 1;
w = f*2*pi;

% parameters for the input
F = 0.3;
f0 = 0.8;
w0 = 2*pi*f0;

omega = w - w0;
psi0 = 0.9*pi;

% parameters for time
fs = 1000;
end_t = 20; % in seconds
T = 1/fs;
time = 0:T:end_t;

% solve r_dot
full_sys = @(t,z,alpha,beta1,beta2,epsilon,w,omega,F,w0)  ...
    [z(1)*alpha + beta1*z(1)^3 + (epsilon*beta2*z(1)^5)/(1-epsilon*z(1)^2) + F*cos(z(2));
    omega - (F/z(1))*sin(z(2))];   

% integrate the ode to obtain the amplitude
[time,r] = ode45(@(t,z) full_sys(t,z,alpha,beta1,beta2,epsilon,w,omega,F,w0),time,[r0,psi0]);

figure(1)
subplot(2,1,1)
plot(time,r(:,1))
grid on
xlabel('Time (s)')
ylabel('Magnitude')
subplot(2,1,2)
plot(time,r(:,2))
grid on
xlabel('Time (s)')
ylabel('Relative Phase')

figure(2)
polar(r(:,2),r(:,1))
hold on

figure(3)
plot3(time,real(r(:,1).*exp(1i*r(:,2))),imag(r(:,1).*exp(1i*r(:,2))))
grid on