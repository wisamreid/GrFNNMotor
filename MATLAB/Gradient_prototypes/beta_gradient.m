%%
clear all;close all;clc;

% parameters for the 1st layer in the network
freqs1 = (0.1:0.1:10)';
nfreqs1 = length(freqs1);
z10 = 0.00001*ones(nfreqs1,1); % initial complex state of the oscillators
alpha1 = 0*ones(nfreqs1,1);
beta11 = -500*ones(nfreqs1,1);
beta12 = 0*ones(nfreqs1,1);
epsilon1 = 1*ones(nfreqs1,1);

% parameters for the 2nd layer in the network
freqs2 = 20*ones(nfreqs1,1);
nfreqs2 = length(freqs2);
z20 = 0.5*ones(nfreqs2,1); % initial complex state of the oscillators
alpha2 = (5:0.1:14.9)'; %*ones(nfreqs2,1);
beta21 = -alpha2 - 40;
beta22 = 0*ones(nfreqs2,1);
epsilon2 = 1*ones(nfreqs2,1);

c = 600; % coupling strength between the two layers

% parameters for time
fs = 1000;
begin_t = 3; % in seconds
T = 1/fs;
time = 0:T:begin_t-T;

z_all = [];
t_all = [];

%%% Integrate the beginning of the function
% note this initial segment has no input
% write the ODE for z dot
full_sys = @(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c)  ...
    [z(1:nfreqs2).*(alpha2 + 1i*2*pi*freqs2 + beta21.*abs(z(1:nfreqs2)).^2 + ...
    (epsilon2.*beta22.*abs(z(1:nfreqs2)).^4)./(1-epsilon2.*abs(z(1:nfreqs2)).^2)) + ...
    c*(z((nfreqs2+1:nfreqs2+nfreqs1)));
    z(nfreqs2+1:nfreqs2+nfreqs1).*(alpha1 + 1i*2*pi*freqs1 + beta11.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2 + ...
    (epsilon1.*beta12.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^4)./(1-epsilon1.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2))];
% integrate the ode to obtain the amplitude
[time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c),time,[z20' z10']);

z_all = [z_all z_i'];
t_all = [t_all time'];

%%% Period of Stimulation 1: at 1.6Hz
% parameters for time during impulse stimulation
num_pulses = 5;
F = 0.1; % forcing amplitude
f = 0.5; % period of stimulation in Hz
stim_t = 1/f; % in seconds
time = t_all(end)+T:T:t_all(end)+stim_t;

for istim = 1:num_pulses
    
    z20 = z_all(1:nfreqs2,end);
    z10 = z_all(nfreqs2+1:nfreqs2+nfreqs1,end) + (F./(1-sqrt(epsilon1).*F)).*(1./(1-sqrt(epsilon1).*conj(z_all(nfreqs2+1:nfreqs2+nfreqs1,end))));
    
    % write the ODE for z dot during stimulation
    full_sys = @(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c)  ...
        [z(1:nfreqs2).*(alpha2 + 1i*2*pi*freqs2 + beta21.*abs(z(1:nfreqs2)).^2 + ...
        (epsilon2.*beta22.*abs(z(1:nfreqs2)).^4)./(1-epsilon2.*abs(z(1:nfreqs2)).^2)) + ...
        c*(z((nfreqs2+1:nfreqs2+nfreqs1)));
        z(nfreqs2+1:nfreqs2+nfreqs1).*(alpha1 + 1i*2*pi*freqs1 + beta11.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2 + ...
        (epsilon1.*beta12.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^4)./(1-epsilon1.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2))];
    % integrate the ode to obtain the amplitude
    [time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c),time,[z20' z10']);
    
    z_all = [z_all z_i'];
    t_all = [t_all time'];
    
    % update time vector
    time = t_all(end)+T:T:t_all(end)+stim_t;
    
end

%%% let the system relax after stimulation
% parameters for time after stimulation
relax_t = 3; % in seconds
time = t_all(end)+T:T:t_all(end)+relax_t;

z20 = z_all(1:nfreqs2,end);
z10 = z_all(nfreqs2+1:nfreqs2+nfreqs1,end);

% write the ODE for z dot during relaxation
full_sys = @(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c)  ...
    [z(1:nfreqs2).*(alpha2 + 1i*2*pi*freqs2 + beta21.*abs(z(1:nfreqs2)).^2 + (epsilon2.*beta22.*abs(z(1:nfreqs2)).^4)./(1-epsilon2.*abs(z(1:nfreqs2)).^2)) + c*(z((nfreqs2+1:nfreqs2+nfreqs1)));
    z(nfreqs2+1:nfreqs2+nfreqs1).*(alpha1 + 1i*2*pi*freqs1 + beta11.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2 + (epsilon1.*beta12.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^4)./(1-epsilon2.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2))];
% integrate the ode to obtain the amplitude
[time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c),time,[z20' z10']);

z_all = [z_all z_i'];
t_all = [t_all time'];

%%% Period of Stimulation 2: at 2.2Hz
% parameters for time during impulse stimulation
f = 1.6; % period of stimulation in Hz
stim_t = 1/f; % in seconds
time = t_all(end)+T:T:t_all(end)+stim_t;

for istim = 1:num_pulses
    
    z20 = z_all(1:nfreqs2,end);
    z10 = z_all(nfreqs2+1:nfreqs2+nfreqs1,end) + (F./(1-sqrt(epsilon1).*F)).*(1./(1-sqrt(epsilon1).*conj(z_all(nfreqs2+1:nfreqs2+nfreqs1)')));
    
    % write the ODE for z dot during stimulation
    full_sys = @(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c)  ...
        [z(1:nfreqs2).*(alpha2 + 1i*2*pi*freqs2 + beta21.*abs(z(1:nfreqs2)).^2 + (epsilon2.*beta22.*abs(z(1:nfreqs2)).^4)./(1-epsilon2.*abs(z(1:nfreqs2)).^2)) + c*(z((nfreqs2+1:nfreqs2+nfreqs1)));
        z(nfreqs2+1:nfreqs2+nfreqs1).*(alpha1 + 1i*2*pi*freqs1 + beta11.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2 + (epsilon1.*beta12.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^4)./(1-epsilon2.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2))];
    % integrate the ode to obtain the amplitude
    [time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c),time,[z20' z10']);
    
    z_all = [z_all z_i'];
    t_all = [t_all time'];
    
    % update time vector
    time = t_all(end)+T:T:t_all(end)+stim_t;
    
end

%%% let the system relax after stimulation
% parameters for time after stimulation
end_t = 5; % in seconds
time = t_all(end)+T:T:t_all(end)+end_t;

z20 = z_all(1:nfreqs2,end);
z10 = z_all(nfreqs2+1:nfreqs2+nfreqs1,end);

% write the ODE for z dot during stimulation
full_sys = @(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c)  ...
    [z(1:nfreqs2).*(alpha2 + 1i*2*pi*freqs2 + beta21.*abs(z(1:nfreqs2)).^2 + (epsilon2.*beta22.*abs(z(1:nfreqs2)).^4)./(1-epsilon2.*abs(z(1:nfreqs2)).^2)) + c*(z((nfreqs2+1:nfreqs2+nfreqs1)));
    z(nfreqs2+1:nfreqs2+nfreqs1).*(alpha1 + 1i*2*pi*freqs1 + beta11.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2 + (epsilon1.*beta12.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^4)./(1-epsilon2.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2))];
% integrate the ode to obtain the amplitude
[time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c),time,[z20' z10']);

z_all = [z_all z_i'];
t_all = [t_all time'];

figure(1)
subplot(2,1,1)
plot(t_all,20*log10(abs(sum(z_all(1:nfreqs2,:)))))
grid on
subplot(2,1,2)
plot(t_all,abs(z_all(nfreqs2+1:nfreqs2+nfreqs1,:)))
grid on

figure(2)
spectrogram(sum(conj(z_all(1:nfreqs2,:)))/max(abs(sum(conj(z_all(1:nfreqs2,:))))),hann(128),100,2^nextpow2(20000),fs,'yaxis')
title('Spectrum of the Beta Oscillator (2nd layer)')
ylim([0 100])
caxis([-70 -10])

figure(3)
polar(angle(sum(z_all(1:nfreqs2,:))),abs(sum(z_all(1:nfreqs2,:))))