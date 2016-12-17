% implementation of the Runge Kutta method

clear all;close all;clc;

% parameters of the limit cycle oscillator
z0 = 0.1;
alpha = 1;
beta1 = -0.0001;
beta2 = 2.1;
epsilon = 1;

% parameters for time
fs = 1000;
dur = 30; % in seconds
T = 1/fs;
time = 0:T:dur;

ntime = length(time);

% initialize empty arrays
t_out = zeros(size(time));
z_out = zeros(size(time));

t_out(1) = time(1);
z_out(1) = z0;

% add an input
input = zeros(size(time));
input(fs*10:fs*10+10) = 100;

for i=1:ntime-1
        
        k1 = z_out(i)*(alpha + 1i*2*pi + beta1*abs(z_out(i))^2 + ((epsilon*beta2*abs(z_out(i)))^4)/(1-epsilon*abs(z_out(i))^2)) + input(i);                
        
        z1 = z_out(i)+(T/2)*k1;
        
        k2 = z1*(alpha + 1i*2*pi + beta1*abs(z1)^2 + ((epsilon*beta2*abs(z1))^4)/(1-epsilon*abs(z1)^2)) + (input(i) + input(i+1))/2;
        
        z2 = z_out(i)+(T/2)*k2;
        
        k3 = z2*(alpha + 1i*2*pi + beta1*abs(z2)^2 + ((epsilon*beta2*abs(z2))^4)/(1-epsilon*abs(z2)^2)) + (input(i) + input(i+1))/2;
                
        z3 = z_out(i)+T*k3;
        
        k4 = z3*(alpha + 1i*2*pi + beta1*abs(z3)^2 + ((epsilon*beta2*abs(z3))^4)/(1-epsilon*abs(z3)^2)) + (input(i));
                
        z_out(i+1) = z_out(i) + (1/6)*T*(k1+(2*k2)+(2*k3)+k4);
        
        t_out(i+1) = t_out(i) + T;        
    
end

figure(1)
axis square
plot((z_out));
grid on
figure(2)
subplot(2,1,1)
plot(time,abs(z_out))
grid on
subplot(2,1,2)
plot(time,input)
grid on