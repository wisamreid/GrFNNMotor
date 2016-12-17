% implementation of the Runge Kutta method

clear all;close all;clc;

% parameters of the gamma oscillator
z01 = 0.0001;
alpha = -100;
beta1 = 400;
beta2 = 1;
epsilon1 = 20;

% parameters of the modulating oscillator
z02 = 0.0001;
alpha2 = -0.5;
beta12 = 1;
beta22 = 10;
epsilon2 = 1;

% parameters of the beta oscillator
z03 = 0.5;
alpha3 = 1;
beta13 = -4;
beta23 = 0;

% parameters for time
fs = 1000;
dur = 10; % in seconds
T = 1/fs;
time = 0:T:dur;

ntime = length(time);

% initialize empty arrays
t_out = zeros(size(time));
z_out_1 = zeros(size(time));
z_out_2 = zeros(size(time));
z_out_3 = zeros(size(time));

t_out(1) = time(1);
z_out_1(1) = z01;
z_out_2(1) = z02;
z_out_3(1) = z03;

% coupling parameter
c12 = 10;
c23 = 0;

% add an input
input = zeros(size(time));
ntaps = 4;
for i=1:ntaps
    input(fs*2*i:fs*2*i+25) = 50*hann(25+1);
end

for i=1:ntime-1
    
    %%% first layer %%%
    
    k1 = z_out_1(i)*(alpha + 1i*2*pi*80 + beta1*abs(z_out_1(i))^2 + ...
        ((epsilon1*beta2*abs(z_out_1(i)))^4)/(1-epsilon1*abs(z_out_1(i))^2)) + input(i);
    
    z1 = z_out_1(i)+(T/2)*k1;
    
    k2 = z1*(alpha + 1i*2*pi*80 + beta1*abs(z1)^2 + ...
        ((epsilon1*beta2*abs(z1))^4)/(1-epsilon1*abs(z1)^2)) + (input(i) + input(i+1))/2;
    
    z2 = z_out_1(i)+(T/2)*k2;
    
    k3 = z2*(alpha + 1i*2*pi*80 + beta1*abs(z2)^2 + ...
        ((epsilon1*beta2*abs(z2))^4)/(1-epsilon1*abs(z2)^2)) + (input(i) + input(i+1))/2;
    
    z3 = z_out_1(i)+T*k3;
    
    k4 = z3*(alpha + 1i*2*pi*80 + beta1*abs(z3)^2 + ...
        ((epsilon1*beta2*abs(z3))^4)/(1-epsilon1*abs(z3)^2)) + (input(i));
    
    z_out_1(i+1) = z_out_1(i) + (1/6)*T*(k1+(2*k2)+(2*k3)+k4);
    
    %%% second layer %%%
    
    k1 = z_out_2(i)*(alpha2 + 1i*2*pi*0.5 + beta12*abs(z_out_2(i))^2 + ...
        ((epsilon2*beta22*abs(z_out_2(i)))^4)/(1-epsilon2*abs(z_out_2(i))^2)) + c12*z_out_1(i);
    
    z1 = z_out_2(i)+(T/2)*k1;
    
    k2 = z1*(alpha2 + 1i*2*pi*20 + beta12*abs(z1)^2 + ...
        ((epsilon2*beta22*abs(z1))^4)/(1-epsilon2*abs(z1)^2)) + c12*(z_out_1(i) + z_out_1(i+1))/2;
    
    z2 = z_out_2(i)+(T/2)*k2;
    
    k3 = z2*(alpha2 + 1i*2*pi*0.5 + beta12*abs(z2)^2 + ...
        ((epsilon2*beta22*abs(z2))^4)/(1-epsilon2*abs(z2)^2)) + c12*(z_out_1(i) + z_out_1(i+1))/2;
    
    z3 = z_out_2(i)+T*k3;
    
    k4 = z3*(alpha2 + 1i*2*pi*0.5 + beta12*abs(z3)^2 + ...
        ((epsilon2*beta22*abs(z3))^4)/(1-epsilon2*abs(z3)^2)) + c12*(z_out_1(i));
    
    z_out_2(i+1) = z_out_2(i) + (1/6)*T*(k1+(2*k2)+(2*k3)+k4);
    
    %%% third layer %%%
    
    k1 = z_out_3(i)*(alpha3 + 1i*2*pi*20 + beta13*abs(z_out_3(i))^2 + ...
        ((epsilon1*beta23*abs(z_out_3(i)))^4)/(1-epsilon1*abs(z_out_3(i))^2)) + c23*z_out_2(i);
    
    z1 = z_out_3(i)+(T/2)*k1;
    
    k2 = z1*(alpha3 + 1i*2*pi*20 + beta13*abs(z1)^2 + ...
        ((epsilon1*beta23*abs(z1))^4)/(1-epsilon1*abs(z1)^2)) + c23*(z_out_2(i) + z_out_2(i+1))/2;
    
    z2 = z_out_3(i)+(T/2)*k2;
    
    k3 = z2*(alpha3 + 1i*2*pi*20 + beta13*abs(z2)^2 + ...
        ((epsilon1*beta23*abs(z2))^4)/(1-epsilon1*abs(z2)^2)) + c23*(z_out_2(i) + z_out_2(i+1))/2;
    
    z3 = z_out_3(i)+T*k3;
    
    k4 = z3*(alpha3 + 1i*2*pi*20 + beta13*abs(z3)^2 + ...
        ((epsilon1*beta23*abs(z3))^4)/(1-epsilon1*abs(z3)^2)) + c23*(z_out_2(i));
    
    z_out_3(i+1) = z_out_3(i) + (1/6)*T*(k1+(2*k2)+(2*k3)+k4);
    
    t_out(i+1) = t_out(i) + T;
    
end

figure(1)
axis square
subplot(1,3,1)
plot((z_out_1));
title('Oscillator 1')
grid on
subplot(1,3,2)
plot((z_out_2));
grid on
title('Oscillator 2')
subplot(1,3,3)
plot((z_out_3));
grid on
title('Oscillator 3')

figure(2)
subplot(3,1,1)
plot(time,abs(z_out_1))
grid on
title('Magnitude')
subplot(3,1,2)
plot(time,input)
grid on
title('Input')
subplot(3,1,3)
spectrogram(z_out_1,kaiser(256,5),220,20000,fs,'yaxis')
title('Spectrum')
ylim([0 100])
xlabel('time(s)')

figure(3)
subplot(3,1,1)
plot(time,abs(z_out_2))
grid on
title('Magnitude')
subplot(3,1,2)
plot(time,input)
grid on
title('Input')
subplot(3,1,3)
spectrogram(z_out_2,kaiser(256,5),220,20000,fs,'yaxis')
title('Spectrum')
ylim([0 100])
xlabel('time(s)')

figure(4)
subplot(3,1,1)
plot(time,abs(z_out_3))
grid on
title('Magnitude')
subplot(3,1,2)
plot(time,input)
grid on
title('Input')
subplot(3,1,3)
spectrogram(z_out_3,kaiser(256,5),220,20000,fs,'yaxis')
title('Spectrum')
ylim([0 100])
xlabel('time(s)')