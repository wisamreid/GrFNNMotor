%%
clear all;close all;clc;

% parameters for the 1st layer in the network
freqs1 = logspace(-0.1,0.6,30)';
nfreqs1 = length(freqs1);
z10 = 0.00001*ones(nfreqs1,1); % initial complex state of the oscillators
alpha1 = 0*ones(nfreqs1,1);
beta11 = -500*ones(nfreqs1,1);
beta12 = 0*ones(nfreqs1,1);
epsilon1 = 1*ones(nfreqs1,1);

% parameters for the 2nd layer in the network
freqs2 = 20;
nfreqs2 = length(freqs2);
z20 = 0.5*ones(nfreqs2,1); % initial complex state of the oscillators
alpha2 = 10*ones(nfreqs2,1);
beta21 = -15*ones(nfreqs2,1);
beta22 = 0*ones(nfreqs2,1);
epsilon2 = 1*ones(nfreqs2,1);

c = 50; % coupling strength between the two layers

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
    c*sum(z((nfreqs2+1:nfreqs2+nfreqs1)));
    z(nfreqs2+1:nfreqs2+nfreqs1).*(alpha1 + 1i*2*pi*freqs1 + beta11.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2 + ...
    (epsilon1.*beta12.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^4)./(1-epsilon1.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2))];
% integrate the ode to obtain the amplitude
[time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c),time,[z20' z10']);

z_all = [z_all z_i'];
t_all = [t_all time'];

%%% Period of Stimulation 1: at 1.6Hz
% parameters for time during impulse stimulation
num_pulses = 10;
F = 0.1; % forcing amplitude
f = 1; % period of stimulation in Hz
stim_t = 1/f; % in seconds
time = t_all(end)+T:T:t_all(end)+stim_t;

for istim = 1:num_pulses
    
    z20 = z_all(1:nfreqs2,end);
    z10 = z_all(nfreqs2+1:nfreqs2+nfreqs1,end) + (F./(1-sqrt(epsilon1).*F)).*(1./(1-sqrt(epsilon1).*conj(z_all(nfreqs2+1:nfreqs2+nfreqs1,end))));
    
    % write the ODE for z dot during stimulation
    full_sys = @(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c)  ...
        [z(1:nfreqs2).*(alpha2 + 1i*2*pi*freqs2 + beta21.*abs(z(1:nfreqs2)).^2 + (epsilon2.*beta22.*abs(z(1:nfreqs2)).^4)./(1-epsilon2.*abs(z(1:nfreqs2)).^2)) + c*sum(z((nfreqs2+1:nfreqs2+nfreqs1)));
        z(nfreqs2+1:nfreqs2+nfreqs1).*(alpha1 + 1i*2*pi*freqs1 + beta11.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2 + (epsilon1.*beta12.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^4)./(1-epsilon1.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2))];
    % integrate the ode to obtain the amplitude
    [time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c),time,[z20' z10']);
    
    z_all = [z_all z_i'];
    t_all = [t_all time'];
    
    % update time vector
    time = t_all(end)+T:T:t_all(end)+stim_t;
    
end

%%% let the system relax after stimulation
% parameters for time after stimulation
relax_t = 5; % in seconds
time = t_all(end)+T:T:t_all(end)+relax_t;

z20 = z_all(1:nfreqs2,end);
z10 = z_all(nfreqs2+1:nfreqs2+nfreqs1,end);

% write the ODE for z dot during relaxation
full_sys = @(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c)  ...
    [z(1:nfreqs2).*(alpha2 + 1i*2*pi*freqs2 + beta21.*abs(z(1:nfreqs2)).^2 + (epsilon2.*beta22.*abs(z(1:nfreqs2)).^4)./(1-epsilon2.*abs(z(1:nfreqs2)).^2)) + c*sum(z((nfreqs2+1:nfreqs2+nfreqs1)));
    z(nfreqs2+1:nfreqs2+nfreqs1).*(alpha1 + 1i*2*pi*freqs1 + beta11.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2 + (epsilon1.*beta12.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^4)./(1-epsilon2.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2))];
% integrate the ode to obtain the amplitude
[time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c),time,[z20' z10']);

z_all = [z_all z_i'];
t_all = [t_all time'];

%%% Period of Stimulation 2: at 2.2Hz
% parameters for time during impulse stimulation
f = 2.2; % period of stimulation in Hz
stim_t = 1/f; % in seconds
time = t_all(end)+T:T:t_all(end)+stim_t;

for istim = 1:num_pulses
    
    z20 = z_all(1:nfreqs2,end);
    z10 = z_all(nfreqs2+1:nfreqs2+nfreqs1,end) + (F./(1-sqrt(epsilon1).*F)).*(1./(1-sqrt(epsilon1).*conj(z_all(nfreqs2+1:nfreqs2+nfreqs1)')));
    
    % write the ODE for z dot during stimulation
    full_sys = @(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c)  ...
        [z(1:nfreqs2).*(alpha2 + 1i*2*pi*freqs2 + beta21.*abs(z(1:nfreqs2)).^2 + (epsilon2.*beta22.*abs(z(1:nfreqs2)).^4)./(1-epsilon2.*abs(z(1:nfreqs2)).^2)) + c*sum(z((nfreqs2+1:nfreqs2+nfreqs1)));
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
end_t = 3; % in seconds
time = t_all(end)+T:T:t_all(end)+end_t;

z20 = z_all(1:nfreqs2,end);
z10 = z_all(nfreqs2+1:nfreqs2+nfreqs1,end);

% write the ODE for z dot during stimulation
full_sys = @(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c)  ...
    [z(1:nfreqs2).*(alpha2 + 1i*2*pi*freqs2 + beta21.*abs(z(1:nfreqs2)).^2 + (epsilon2.*beta22.*abs(z(1:nfreqs2)).^4)./(1-epsilon2.*abs(z(1:nfreqs2)).^2)) + c*sum(z((nfreqs2+1:nfreqs2+nfreqs1)));
    z(nfreqs2+1:nfreqs2+nfreqs1).*(alpha1 + 1i*2*pi*freqs1 + beta11.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2 + (epsilon1.*beta12.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^4)./(1-epsilon2.*abs(z(nfreqs2+1:nfreqs2+nfreqs1)).^2))];
% integrate the ode to obtain the amplitude
[time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,alpha2,beta21,beta22,epsilon2,freqs1,freqs2,c),time,[z20' z10']);

z_all = [z_all z_i'];
t_all = [t_all time'];

figure(1)
subplot(2,1,1)
plot(t_all,20*log10(abs(z_all(1,:))))
grid on
subplot(2,1,2)
plot(t_all,sum((z_all(2:end,:))))
grid on

figure(2)
spectrogram(conj(z_all(1,:)),hann(128),100,2^nextpow2(20000),fs,'yaxis')
title('Spectrum of the Beta Oscillator (2nd layer)')
ylim([0 100])
caxis([-70 -10])

figure(3)
polar(angle(z_all(1,:)),abs(z_all(1,:)))
%%
% old version starts beyond this line:

for ifreq = 1:nfreqs1
    max_abs = 0.22;
    figure(1)
    subplot(1,nfreqs1,ifreq)
    plot(t_all,abs(z_all(ifreq,:)))
    title(sprintf('%.2f Hz',freqs1(ifreq)))
    axis([t_all(1) t_all(end) 0 max_abs])
    grid on
    xlabel('Time (s)')
    ylabel('Magnitude')
    figure(2)
    subplot(1,nfreqs1,ifreq)
    polar(pi,max_abs)
    hold on
    polar(angle(z_all(ifreq,:)),abs(z_all(ifreq,:)))
    hold on
    grid on
end

%%
clear all;clf(1);clf(2);clc;

% parameters of the network
freqs1 = logspace(0,0.6,20);
nfreqs1 = length(freqs1);


z_all_freqs = [];

for ifreq=1:nfreqs1
    
    
    % parameters of the oscillator
    z10 = 0.001; % initial complex state of the oscillator
    alpha1 = 0;
    beta11 = -100;
    beta12 = 0;
    epsilon1 = 1;
    
    % parameters for time
    fs = 1000;
    begin_t = 2; % in seconds
    T = 1/fs;
    time = 0:T:begin_t-T;
    
    z_all = [];
    t_all = [];
    
    %%% Integrate the beginning of the function
    % note this initial segment has no input
    % write the ODE for z dot
    full_sys = @(t,z,alpha,beta1,beta2,epsilon,freqs)  ...
        z*(alpha + 1i*2*pi*freqs(ifreq) + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
    % integrate the ode to obtain the amplitude
    [time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,freqs1),time,z10);
    z_all = [z_all z_i'];
    t_all = [t_all time'];
    
    %%% Period of Stimulation 1: at 3Hz
    % parameters for time during impulse stimulation
    num_pulses = 12;
    F = 0.1; % forcing amplitude
    f = 1.6; % period of stimulation in Hz
    stim_t = 1/f; % in seconds
    time = t_all(end)+T:T:t_all(end)+stim_t;
    
    for istim = 1:num_pulses
        
        z10 = z_all(end) + (F/(1-sqrt(epsilon1)*F))*(1/(1-sqrt(epsilon1)*conj(z_all(end))));
        
        % write the ODE for z dot during stimulation
        full_sys = @(t,z,alpha,beta1,beta2,epsilon,freqs)  ...
            z*(alpha + 1i*2*pi*freqs(ifreq) + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
        % integrate the ode to obtain the amplitude
        [time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,freqs1),time,z10);
        
        z_all = [z_all z_i'];
        t_all = [t_all time'];
        
        % update time vector
        time = t_all(end)+T:T:t_all(end)+stim_t;
        
    end
    
    %%% let the system relax after stimulation
    % parameters for time after stimulation
    relax_t = 5; % in seconds
    time = t_all(end)+T:T:t_all(end)+relax_t;
    
    z10 = z_all(end);
    % write the ODE for z dot during stimulation
    full_sys = @(t,z,alpha,beta1,beta2,epsilon,freqs)  ...
        z*(alpha + 1i*2*pi*freqs(ifreq) + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
    % integrate the ode to obtain the amplitude
    [time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,freqs1),time,z10);
    
    z_all = [z_all z_i'];
    t_all = [t_all time'];
    
    %%% Period of Stimulation 2: at 4Hz
    % parameters for time during impulse stimulation
    f = 2.2; % period of stimulation in Hz
    stim_t = 1/f; % in seconds
    time = t_all(end)+T:T:t_all(end)+stim_t;
    
    for istim = 1:num_pulses
        
        z10 = z_all(end) + (F/(1-sqrt(epsilon1)*F))*(1/(1-sqrt(epsilon1)*conj(z_all(end))));
        
        % write the ODE for z dot during stimulation
        full_sys = @(t,z,alpha,beta1,beta2,epsilon,freqs)  ...
            z*(alpha + 1i*2*pi*freqs(ifreq) + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
        % integrate the ode to obtain the amplitude
        [time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,freqs1),time,z10);
        
        z_all = [z_all z_i'];
        t_all = [t_all time'];
        
        % update time vector
        time = t_all(end)+T:T:t_all(end)+stim_t;
        
    end
    
    %%% let the system relax after stimulation
    % parameters for time after stimulation
    end_t = 3; % in seconds
    time = t_all(end)+T:T:t_all(end)+end_t;
    
    z10 = z_all(end);
    % write the ODE for z dot during stimulation
    full_sys = @(t,z,alpha,beta1,beta2,epsilon,freqs)  ...
        z*(alpha + 1i*2*pi*freqs(ifreq) + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
    % integrate the ode to obtain the amplitude
    [time,z_i] = ode45(@(t,z) full_sys(t,z,alpha1,beta11,beta12,epsilon1,freqs1),time,z10);
    
    z_all = [z_all z_i'];
    t_all = [t_all time'];
    
    max_abs = 0.22;
    figure(1)
    subplot(1,nfreqs1,ifreq)
    plot(t_all,abs(z_all))
    title(sprintf('%.2f Hz',freqs1(ifreq)))
    axis([t_all(1) t_all(end) 0 max_abs])
    grid on
    xlabel('Time (s)')
    ylabel('Magnitude')
    figure(2)
    subplot(1,nfreqs1,ifreq)
    polar(pi,max_abs)
    hold on
    polar(angle(z_all),abs(z_all))
    hold on
    grid on
    
    z_all_freqs(ifreq,:) = z_all;
    
end

figure(3)
subplot(2,1,1)
plot3(t_all,real(sum(z_all_freqs)),imag(sum(z_all_freqs)))
grid on
subplot(2,1,2)
plot(t_all,abs(sum(z_all_freqs)))