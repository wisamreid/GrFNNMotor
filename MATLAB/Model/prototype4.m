%%
% a double limit cycle oscillator driving a gama oscillator

clear all;close all;clc;

animate = 1; % true or false

% play with the number of pulses (10 won't be enough to learn)
% play with making alpha more negative so that the learning quality of the
% oscillator goes away (make z0 0.9 to see the contour of the phase portrait).

% parameters of the double-limit cycle oscillator
z01 = 0.;
alpha1 = -2;
beta11 = 10;
beta21 = -1;
epsilon1 = 1;
f11 = 1;

% coupling strength from double-limit cycle to the gamma oscillator
c = 1000;

% parameters of the gamma oscillator
z02 = 0;
alpha2 = -15;
beta12 = 0;
beta22 = 0;
epsilon2 = 0;
f21 = 80;

% parameters for time
fs = 1000;
start_t = 1.5; % in seconds
end_t = 3; % in seconds
T = 1/fs;
time = 0:T:start_t;

% input parameters
npulses = 10;
pulse_amp = 0.2945; % in Fujioka et al. 2012 the tempos were 2.5Hz, 1.7Hz, and 1.3Hz
% the corresponding pulse amplitudes are: 2.5Hz : 0.153, 1.7Hz : 0.205, 1.3Hz : 0.247
pul_freq = 1;
pul_step = 1/pul_freq; % in seconds

% the output of beta drives a gamma oscillator
dzdt = @(t,z,alpha1,beta11,beta21,epsilon1,f11,c,alpha2,beta12,beta22,epsilon2,f21)  ...
    [z(1)*(alpha2 + 1i*2*pi*f21 + beta12*abs(z(1))^2 + ...
    (epsilon2*beta22*abs(z(1))^4)/(1-epsilon2*abs(z(1))^2)) + c*z(2); ...
    z(2)*(alpha1 + 1i*2*pi*f11 + beta11*abs(z(2))^2 + ...
    (epsilon1*beta21*abs(z(2))^4)/(1-epsilon1*abs(z(2))^2))];

% integrate the ode to obtain the amplitude
[t_out,z_out] = ode45(@(t,z) dzdt(t,z,alpha1,beta11,beta21,epsilon1,f11,c,...
    alpha2,beta12,beta22,epsilon2,f21),time,[z02 z01]);

z_all_1 = [];
z_all_2 = [];
time_all = [];
z_all_1 = [z_all_1 z_out(:,2)'];
z_all_2 = [z_all_2 z_out(:,1)'];
time_all = [time_all t_out'];


% stimulation
for i=1:npulses
    
    % new oscillator parameter
    z01 = z_all_1(end) + pulse_amp;
    z02 = z_all_2(end);
    
    % parameters for time
    time = time_all(end):T:time_all(end)+pul_step;
    
    % the output of beta drives a gamma oscillator
    dzdt = @(t,z,alpha1,beta11,beta21,epsilon1,f11,c,alpha2,beta12,beta22,epsilon2,f21)  ...
        [z(1)*(alpha2 + 1i*2*pi*f21 + beta12*abs(z(1))^2 + ...
        (epsilon2*beta22*abs(z(1))^4)/(1-epsilon2*abs(z(1))^2)) + c*z(2); ...
        z(2)*(alpha1 + 1i*2*pi*f11 + beta11*abs(z(2))^2 + ...
        (epsilon1*beta21*abs(z(2))^4)/(1-epsilon1*abs(z(2))^2))];
    
    % integrate the ode to obtain the amplitude
    [t_out,z_out] = ode45(@(t,z) dzdt(t,z,alpha1,beta11,beta21,epsilon1,f11,c,...
        alpha2,beta12,beta22,epsilon2,f21),time,[z02 z01]);
    
    z_all_1 = [z_all_1 z_out(:,2)'];
    z_all_2 = [z_all_2 z_out(:,1)'];
    time_all = [time_all t_out'];
    
end

% new oscillator parameter
z01 = z_all_1(end);
z02 = z_all_2(end);

% parameters for time
time = time_all(end):T:time_all(end)+end_t;

% the output of beta drives a gamma oscillator
dzdt = @(t,z,alpha1,beta11,beta21,epsilon1,f11,c,alpha2,beta12,beta22,epsilon2,f21)  ...
    [z(1)*(alpha2 + 1i*2*pi*f21 + beta12*abs(z(1))^2 + ...
    (epsilon2*beta22*abs(z(1))^4)/(1-epsilon2*abs(z(1))^2)) + c*z(2); ...
    z(2)*(alpha1 + 1i*2*pi*f11 + beta11*abs(z(2))^2 + ...
    (epsilon1*beta21*abs(z(2))^4)/(1-epsilon1*abs(z(2))^2))];

% integrate the ode to obtain the amplitude
[t_out,z_out] = ode45(@(t,z) dzdt(t,z,alpha1,beta11,beta21,epsilon1,f11,c,...
    alpha2,beta12,beta22,epsilon2,f21),time,[z02 z01]);

z_all_1 = [z_all_1 z_out(:,2)'];
z_all_2 = [z_all_2 z_out(:,1)'];
time_all = [time_all t_out'];% figure('position', [1, 1, 500, 1000])


if animate
    figure('position', [1, 1, 5000, 1000])
    for i=1:length(time_all)
        
        if mod(i,fs/20) == 0
            subplot(2,2,1:2)
            plot3(time_all(1:i),real(z_all_2(1:i)),imag(z_all_2(1:i)),'-k')
            grid on
            hold on
            xlabel('Time (s)')
            ylabel('Real part')
            zlabel('Imaginary part')
            subplot(2,2,3:4)
            plot3(time_all(1:i),real(z_all_1(1:i)),imag(z_all_1(1:i)),'-b')
            grid on
            hold on
            xlabel('Time (s)')
            ylabel('Real part')
            zlabel('Imaginary part')
            pause
        end
    end
end

figure(2)
plot(linspace(-(fs/2),(fs/2)-1,length(z_all_1)),fftshift(abs(fft(z_all_2))))
axis tight
grid on