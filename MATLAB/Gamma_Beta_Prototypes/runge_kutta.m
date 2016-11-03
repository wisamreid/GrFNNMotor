% implementation of the Runge Kutta method

clear all;close all;clc;

% indicate which variables will be used symbolically
% syms z t

% parameters of the limit cycle oscillator
z0 = 0.1;
alpha = 1;
beta1 = -1;
beta2 = 0;
epsilon = 1;

% parameters for time
fs = 1000;
dur = 30; % in seconds
T = 1/fs;
time = 0:T:dur;

ntime = length(time);

h=time(ntime)/(ntime-1);

% initialize empty arrays
t_out = zeros(size(time));
z_out = zeros(size(time));

t_out(1) = time(1);
z_out(1) = z0;

% declare function
fun = @(t,z) (z*(alpha + 1i*2*pi + beta1*abs(z)^2 + ...     
    ((epsilon*beta2*abs(z))^4)/(1-epsilon*abs(z)^2)));

% add an input
input = zeros(size(time));
input(fs*10:fs*10+10) = 0.2;
input(fs*10+1000:fs*10+1010) = 0.2;
input(fs*10+2000:fs*10+2010) = 0.2;
input(fs*10+3000:fs*10+3010) = 0.2;

for i=1:ntime-1    
    
    k1 = fun(t_out(i),z_out(i));
    
    k2 = fun(t_out(i)+0.5*h,z_out(i)+0.5*h*k1);
    
    k3 = fun(t_out(i)+0.5*h,z_out(i)+0.5*h*k2);
    
    k4 = fun(t_out(i)+h,z_out(i)+h*k3);
    
    z_out(i+1) = z_out(i) + (1/6)*h*(k1+(2*k2)+(2*k3)+k4) + input(i);

    t_out(i+1) = t_out(i) + h;
    
end

figure(1)
plot(z_out)
figure(2)
plot(abs(z_out(1:end)));
grid on