% dynamical systems model of beta gamma coupling in motor cortex
% by Iran Roman

%%
% start with a beta oscillator at a limit cycle. This is the intrinsic
% activity of beta band without stimulation.

clear all; close all; clc

dzdt = @(t,z,alpha,beta1,beta2,epsilon)  ...        
    z(1)*(alpha + 1i*2*pi + beta1*abs(z(1))^2 + ...     
    (epsilon*beta2*abs(z(1))^4)/(1-epsilon*abs(z(1))^2));

% parameters of the limit cycle oscillator
z0 = 0.99;
alpha = 1;
beta1 = -1;
beta2 = 0;
epsilon = 1;

% parameters for time
fs = 1000;
dur = 10; % in seconds
T = 1/fs;
time = 0:T:dur;

% integrate the ode to obtain the amplitude
[time,z] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon),time,z0);

% extract information
figure(1)
plot(z)

%%
% A gamma oscillator will only become transiently active when stimulated.
% For this purpose, we use a subcritical double limit cycle that becomes a
% supercritical doule limit cycle upon stimulation.

clear all; close all; clc

dzdt = @(t,z,alpha,beta1,beta2,epsilon,F,omega0)  ...        
    z(1)*(alpha + 1i*2*pi*79.5 + beta1*abs(z(1))^2 + ...     
    (epsilon*beta2*abs(z(1))^4)/(1-epsilon*abs(z(1))^2)) + ...
    F*exp(1i*omega0*t);

% parameters of the limit cycle oscillator
z0 = 0.5;
alpha = -0.6;
beta1 = 2.5;
beta2 = -1;
epsilon = 1;

% parameters for time
fs = 1000;
dur = 20; % in seconds
T = 1/fs;
time = 0:T:dur;

% parameters for input
f0 = 80;
omega0 = 2*pi*f0;
F = 0.6;

% integrate the ode to obtain the amplitude
[time,z] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,F,omega0),time,z0);

r = abs(z);
Psi = angle(z) - omega0*time;

figure(1)
plot(r,Psi)
xlabel('rho')
ylabel('Psi')
figure(2)
plot(Psi)
xlabel('Time')
ylabel('Relative Phase')
figure(3)
plot(time,abs(z))
xlabel('Time')
ylabel('Magnitude')
figure(4)
plot(z)

%%
clear all; close all; clc

% the output of gamma modulates a betaoscillator
dzdt = @(t,z,alpha,beta1,beta2,epsilon,F,omega0)  ...        
    [z(1)*(alpha + 1i*2*pi*20 + 15*beta1*abs(z(1))^2 + ...     
    (epsilon*beta2*abs(z(1))^4)/(1-epsilon*abs(z(1))^2)) + 600*(z(2)); ...
    (z(2)*(alpha + 1i*2*pi*79.5 + beta1*abs(z(2))^2 + ... % beta layer
    (epsilon*beta2*abs(z(2))^4)/(1-epsilon*abs(z(2))^2)) + ...   
    F*exp(1i*omega0*t))]; % input to beta

% parameters of the limit cycle oscillator
z0 = 0.3;
alpha = -0.1;
beta1 = 2.5;
beta2 = -10;
epsilon = 1;

% parameters for time
fs = 1000;
dur = 30; % in seconds
T = 1/fs;
time = 0:T:dur;

% parameters for input
f0 = 80;
omega0 = 2*pi*f0;
F = 0.6;

% integrate the ode to obtain the amplitude
[t,z] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,F,omega0),time,[z0,z0]);

% extract information
r = abs(z);
dr2dt = diff(r(:,1));
dr1dt = diff(r(:,2));

figure(1)
plot(z)
% plot(r(1:end-1,2),dr1dt)
% hold on
% plot(r(1:end-1,1),dr2dt)

figure(2)
subplot(3,1,1)
plot(time,r)
title('Magnitude of Oscillators')
ylabel('Amplitude')
xlabel('time(s)')
legend('Beta','Gamma')
subplot(3,1,2)
plot(time,real(z))
title('Real part of Oscillators')
legend('Beta','Gamma')
ylabel('Amplitude')
xlabel('time(s)')
subplot(3,1,3)
spectrogram((z(:,1)+z(:,2)),kaiser(256,5),220,20000,fs,'yaxis')
title('Spectrum')
ylim([0 100])
xlabel('time(s)')

% use 80 Hz instead of 40Hz. We just don't know too much about 40 Hz gamma. 
% Talamocortical loop may be the one giving the 40Hz response (Healthy). In
% Parkinson's this loop is not there. 
% 
% You see beta everywhere. Even in the cerebelum and Olivo-cerebelar loop. 
% Next variation: Without sound, a top down (intention to move) cue is
% given to the system, and the system produces movement. Have beta dip
% prior to single movement, and high gamma burst before the dip of beta. 
%
% After having this, we can modulate the parameters for the initial beta,
% with high amplitude beta not letting gamma act. Jenkinson's and Brown hypothesis:
% there is dopaminergic input coming to STN, low level of input, and enhances
% beta. This is the parkinsonism state. Have the model have a parameter:
% High Gamma amplitude, beta level, beta level is modulated by the
% available Dopamine

%%

% notes from Takako
% Gamma does not spontaneously oscillate. Beta does. 
% read: New insights into the relationship between dopamine, beta oscillations and motor function.
% Beta is always oscillating. Gamma gives exogenous or endogenous input to
% affect beta. Think about high gamma. DBS people didn't find beta
% desynchronization when giving low gamma. High gamma will give the
% desynchronization of beta. If you stimulate 20Hz, people are affected, with 40Hz gamma 
% nothing happens, with high gamma 80Hz there is improvement. Dopamine
% depleted brain gives constant beta. 

% Gamma must start with sound. 

% Granger causality. Berhnard. Issues with phase coherence. Bayes might
% adress causality clarification. DCM

%%

% 11/02/2016
% Due to status, these people can say things they think without very clear evicende
% Cautious Note: Not Always possible what they say. 
% Compare to Terry, who is a pioneer, but doesn't show the flag too much.
% Only cited when thinking about predictive coding in the physiological
% system that they are trying to infer. Friston has a visionary position in
% the community, and thus one has to be carefull with their observations. 

% Hippocampal, rat work, gamma band, Rudolfo Llinas work. Very stablished
% animal physiology lab. Everyone ends up coming to MEG. 