%%
% Linear Oscillator

clear all;close all;clc;

animate = 1; % true or false

% parameters of the gamma oscillator
z0 = 0.8;
alpha = -10;
beta1 = 0;
beta2 = 0;
epsilon = 0;
f1 = 10;

% parameters for time
fs = 1000;
start_t = 1; % in seconds
end_t = 1; % in seconds
T = 1/fs;
time = 0:T:start_t;

% input parameter
npulses = 3;
pulse_amp = 0.8;
pul_step = 1; % in seconds

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);

z_all = [];
time_all = [];
z_all = [z_all z_out_1'];
time_all = [time_all t_out_1'];

for i=1:npulses
    
    % new oscillator parameter
    z0 = z_out_1(end) + pulse_amp;
    
    % parameters for time
    time = t_out_1(end):T:t_out_1(end)+pul_step;
    
    dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
        z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
    
    % integrate the ode to obtain the amplitude
    [t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);
    
    z_all = [z_all z_out_1'];
    time_all = [time_all t_out_1'];
    
end

% new oscillator parameter
z0 = z_out_1(end);

% parameters for time
time = t_out_1(end):T:t_out_1(end)+end_t;

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);

z_all = [z_all z_out_1'];
time_all = [time_all t_out_1'];

figure('position', [1, 1, 500, 1000])
subplot(4,2,1)
plot((z_all))
axis square
xlabel('Real Part')
ylabel('Imaginary Part')
title('Complex Plane')
grid on
subplot(4,2,2)
plot((abs(z_all(1:start_t*fs-1))),diff(abs(z_all(1:start_t*fs))));
axis square
xlabel('Magintude')
ylabel('Maginitude Derivative')
title('Phase Portrait')
grid on
subplot(4,2,3:4)
plot(time_all,abs(z_all))
grid on
hold on
xlabel('Time (s)')
ylabel('Amplitude')
title('Magnitude')

if animate == 1
    
    for i=1:length(time_all)
        
        if mod(i,fs/20) == 0
            subplot(4,2,5:8)
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
% Hopf Oscillator

clear all;close all;clc;

animate = 1; % true or false

% parameters of the gamma oscillator
z0 = 0.8;
alpha = 0;
beta1 = -100;
beta2 = 0;
epsilon = 0;
f1 = 10;

% parameters for time
fs = 1000;
start_t = 1; % in seconds
end_t = 5; % in seconds
T = 1/fs;
time = 0:T:start_t;

% input parameter
npulses = 3;
pulse_amp = 0.8;
pul_step = 1; % in seconds

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);

z_all = [];
time_all = [];
z_all = [z_all z_out_1'];
time_all = [time_all t_out_1'];

for i=1:npulses
    
    % new oscillator parameter
    z0 = z_out_1(end) + pulse_amp;
    
    % parameters for time
    time = t_out_1(end):T:t_out_1(end)+pul_step;
    
    dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
        z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
    
    % integrate the ode to obtain the amplitude
    [t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);
    
    z_all = [z_all z_out_1'];
    time_all = [time_all t_out_1'];
    
end

% new oscillator parameter
z0 = z_out_1(end);

% parameters for time
time = t_out_1(end):T:t_out_1(end)+end_t;

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);

z_all = [z_all z_out_1'];
time_all = [time_all t_out_1'];

figure('position', [1, 1, 500, 1000])
subplot(4,2,1)
plot((z_all))
axis square
xlabel('Real Part')
ylabel('Imaginary Part')
title('Complex Plane')
grid on
subplot(4,2,2)
plot((abs(z_all(1:start_t*fs-1))),diff(abs(z_all(1:start_t*fs))));
axis square
xlabel('Magintude')
ylabel('Maginitude Derivative')
title('Phase Portrait')
grid on
subplot(4,2,3:4)
plot(time_all,abs(z_all))
grid on
hold on
xlabel('Time (s)')
ylabel('Amplitude')
title('Magnitude')

if animate == 1
    
    for i=1:length(time_all)
        
        if mod(i,fs/20) == 0
            subplot(4,2,5:8)
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
% Hopf Oscillator

clear all;close all;clc;

animate = 1; % true or false

% parameters of the gamma oscillator
z0 = 0.8;
alpha = 5;
beta1 = -10;
beta2 = 0;
epsilon = 0;
f1 = 10;

% parameters for time
fs = 1000;
start_t = 1; % in seconds
end_t = 1; % in seconds
T = 1/fs;
time = 0:T:start_t;

% input parameter
npulses = 3;
pulse_amp = -0.8;
pul_step = 1; % in seconds

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); 

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);

z_all = [];
time_all = [];
z_all = [z_all z_out_1'];
time_all = [time_all t_out_1'];

for i=1:npulses
    
    % new oscillator parameter
    z0 = z_out_1(end) + pulse_amp;
    
    % parameters for time
    time = t_out_1(end):T:t_out_1(end)+pul_step;
    
    dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
        z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
    
    % integrate the ode to obtain the amplitude
    [t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);
    
    z_all = [z_all z_out_1'];
    time_all = [time_all t_out_1'];
    
end

% new oscillator parameter
z0 = z_out_1(end);

% parameters for time
time = t_out_1(end):T:t_out_1(end)+end_t;

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);

z_all = [z_all z_out_1'];
time_all = [time_all t_out_1'];

figure('position', [1, 1, 500, 1000])
subplot(4,2,1)
plot((z_all))
axis square
xlabel('Real Part')
ylabel('Imaginary Part')
title('Complex Plane')
grid on
subplot(4,2,2)
plot((abs(z_all(1:start_t*fs-1))),diff(abs(z_all(1:start_t*fs))));
axis square
xlabel('Magintude')
ylabel('Maginitude Derivative')
title('Phase Portrait')
grid on
subplot(4,2,3:4)
plot(time_all,abs(z_all))
grid on
hold on
xlabel('Time (s)')
ylabel('Amplitude')
title('Magnitude')

if animate == 1
    
    for i=1:length(time_all)
        
        if mod(i,fs/20) == 0
            subplot(4,2,5:8)
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
% Hopf Oscillator

clear all;close all;clc;

animate = 1; % true or false

% play with the number of pulses (9 won't be enough to learn)
% play with making alpha more negative so that the learning quality of the
% oscillator goes away (make z0 0.9 to see the contour of the phase portrait).

% parameters of the gamma oscillator
z0 = 0.;
alpha = -2;
beta1 = 10;
beta2 = -1;
epsilon = 1;
f1 = 1;

% parameters for time
fs = 1000;
start_t = 3; % in seconds
end_t = 3; % in seconds
T = 1/fs;
time = 0:T:start_t;

% input parameters
npulses = 5;
pulse_amp = 0.2945; % 1Hz : 0.2945, 2Hz : 0.18, 3Hz : 0.135, 4Hz : 0.11, 5Hz : 0.10
% in Fujioka et al. 2012 the tempos were 2.5Hz, 1.7Hz, and 1.3Hz
% the corresponding pulse amplitudes are: 2.5Hz : 0.153, 1.7Hz : 0.205, 1.3Hz : 0.247
pul_freq = 1;
pul_step = 1/pul_freq; % in seconds

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); 

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);

z_all = [];
time_all = [];
z_all = [z_all z_out_1'];
time_all = [time_all t_out_1'];

for i=1:npulses
    
    % new oscillator parameter
    z0 = z_out_1(end) + pulse_amp;
    
    % parameters for time
    time = t_out_1(end):T:t_out_1(end)+pul_step;
    
    dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
        z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); 
    
    % integrate the ode to obtain the amplitude
    [t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);
    
    z_all = [z_all z_out_1'];
    time_all = [time_all t_out_1'];
    
end

% new oscillator parameter
z0 = z_out_1(end);

% parameters for time
time = t_out_1(end):T:t_out_1(end)+end_t;

dzdt = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); 

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f1),time,z0);

z_all = [z_all z_out_1'];
time_all = [time_all t_out_1'];

figure('position', [1, 1, 500, 1000])
subplot(4,2,1)
plot((z_all))
axis square
xlabel('Real Part')
ylabel('Imaginary Part')
title('Complex Plane')
grid on
subplot(4,2,2)
plot((abs(z_all(1:start_t*fs-1))),diff(abs(z_all(1:start_t*fs))));
axis square
grid on
xlabel('Magintude')
ylabel('Maginitude Derivative')
title('Phase Portrait')
hold on
for i=0:npulses-1
    % ugly but works
    pause
    plot(abs(z_all((start_t*fs+i*pul_step*fs)+i+2:(start_t*fs+(i+1)*pul_step*fs)-1)),diff(abs(z_all((start_t*fs+i*pul_step*fs)+i+2:(start_t*fs+(i+1)*pul_step*fs)))),'r');
end
plot(abs(z_all((start_t*fs+(i+1)*pul_step*fs)+i+1+2:end-1)),diff(abs(z_all((start_t*fs+(i+1)*pul_step*fs)+i+1+2:end))),'r');
grid on
subplot(4,2,3:4)
plot(time_all,abs(z_all))
grid on
hold on
xlabel('Time (s)')
ylabel('Amplitude')
title('Magnitude')

if animate == 1
    
    for i=1:length(time_all)
        
        if mod(i,fs/20) == 0
            subplot(4,2,5:8)
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
% bank of three oscillators learning a pulse

clear all;close all;clc;

animate = 1; % true or false

% play with the number of pulses (10 won't be enough to learn)
% play with making alpha more negative so that the learning quality of the
% oscillator goes away (make z0 0.9 to see the contour of the phase portrait).

% parameters of the oscillators
z01 = 0.0;
z02 = 0.0;
z03 = 0.0;
alpha = -2;
beta1 = 10;
beta2 = -1;
epsilon = 1;
f1 = 1.3;
f2 = 1.7;
f3 = 2.5;

% parameters for time
fs = 1000;
start_t = 1.5; % in seconds
end_t = 5; % in seconds
T = 1/fs;
time = 0:T:start_t;

% input parameters
npulses = 10;
pulse_amp = 0.205; % 1Hz : 0.2945, 2Hz : 0.18, 3Hz : 0.135, 4Hz : 0.11, 5Hz : 0.10
% in Fujioka et al. 2012 the tempos were 2.5Hz, 1.7Hz, and 1.3Hz
% the corresponding pulse amplitudes are: 2.5Hz : 0.153, 1.7Hz : 0.205, 1.3Hz : 0.247
pul_freq = 1.7;
pul_step = 1/pul_freq; % in seconds

% bank of oscillators
dzdt1 = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); 
dzdt2 = @(t,z,alpha,beta1,beta2,epsilon,f2)  ...
    z*(alpha + 1i*2*pi*f2 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); 
dzdt3 = @(t,z,alpha,beta1,beta2,epsilon,f3)  ...
    z*(alpha + 1i*2*pi*f3 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); 

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt1(t,z,alpha,beta1,beta2,epsilon,f1),time,z01);
[~,z_out_2] = ode45(@(t,z) dzdt2(t,z,alpha,beta1,beta2,epsilon,f2),time,z02);
[~,z_out_3] = ode45(@(t,z) dzdt3(t,z,alpha,beta1,beta2,epsilon,f3),time,z03);

time_all = [];
time_all = [time_all t_out_1'];
z_all_1 = [];
z_all_1 = [z_all_1 z_out_1'];
z_all_2 = [];
z_all_2 = [z_all_2 z_out_2'];
z_all_3 = [];
z_all_3 = [z_all_3 z_out_3'];

for i=1:npulses
    
    % new oscillator parameter
    z01 = z_out_1(end) + pulse_amp;
    z02 = z_out_2(end) + pulse_amp;
    z03 = z_out_3(end) + pulse_amp;
    
    % parameters for time
    time = t_out_1(end):T:t_out_1(end)+pul_step;
    
    dzdt1 = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
        z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
    dzdt2 = @(t,z,alpha,beta1,beta2,epsilon,f2)  ...
        z*(alpha + 1i*2*pi*f2 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
    dzdt3 = @(t,z,alpha,beta1,beta2,epsilon,f3)  ...
        z*(alpha + 1i*2*pi*f3 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
    
    % integrate the ode to obtain the amplitude
    [t_out_1,z_out_1] = ode45(@(t,z) dzdt1(t,z,alpha,beta1,beta2,epsilon,f1),time,z01);
    [~,z_out_2] = ode45(@(t,z) dzdt2(t,z,alpha,beta1,beta2,epsilon,f2),time,z02);
    [~,z_out_3] = ode45(@(t,z) dzdt3(t,z,alpha,beta1,beta2,epsilon,f3),time,z03);
    
    time_all = [time_all t_out_1'];
    z_all_1 = [z_all_1 z_out_1'];
    z_all_2 = [z_all_2 z_out_2'];
    z_all_3 = [z_all_3 z_out_3'];
    
end

% new oscillator parameter
z01 = z_out_1(end);
z02 = z_out_2(end);
z03 = z_out_3(end);

% parameters for time
time = t_out_1(end):T:t_out_1(end)+end_t;

dzdt1 = @(t,z,alpha,beta1,beta2,epsilon,f1)  ...
    z*(alpha + 1i*2*pi*f1 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
dzdt2 = @(t,z,alpha,beta1,beta2,epsilon,f2)  ...
    z*(alpha + 1i*2*pi*f2 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));
dzdt3 = @(t,z,alpha,beta1,beta2,epsilon,f3)  ...
    z*(alpha + 1i*2*pi*f3 + beta1*abs(z)^2 + (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2));

% integrate the ode to obtain the amplitude
[t_out_1,z_out_1] = ode45(@(t,z) dzdt1(t,z,alpha,beta1,beta2,epsilon,f1),time,z01);
[~,z_out_2] = ode45(@(t,z) dzdt2(t,z,alpha,beta1,beta2,epsilon,f2),time,z02);
[~,z_out_3] = ode45(@(t,z) dzdt3(t,z,alpha,beta1,beta2,epsilon,f3),time,z03);

time_all = [time_all t_out_1'];
z_all_1 = [z_all_1 z_out_1'];
z_all_2 = [z_all_2 z_out_2'];
z_all_3 = [z_all_3 z_out_3'];

figure('position', [1, 1, 5000, 1000])
subplot(4,6,1)
plot((z_all_1))
axis square
xlabel('Real Part')
ylabel('Imaginary Part')
title('Complex Plane')
grid on
subplot(4,6,3)
plot((z_all_2))
axis square
xlabel('Real Part')
ylabel('Imaginary Part')
title('Complex Plane')
grid on
subplot(4,6,5)
plot((z_all_3))
axis square
xlabel('Real Part')
ylabel('Imaginary Part')
title('Complex Plane')
grid on
subplot(4,6,2)
plot((abs(z_all_1(1:start_t*fs-1))),diff(abs(z_all_1(1:start_t*fs))));
axis square
grid on
xlabel('Magintude')
ylabel('Maginitude Derivative')
title('Phase Portrait')
hold on
subplot(4,6,4)
plot((abs(z_all_2(1:start_t*fs-1))),diff(abs(z_all_2(1:start_t*fs))));
axis square
grid on
xlabel('Magintude')
ylabel('Maginitude Derivative')
title('Phase Portrait')
hold on
subplot(4,6,6)
plot((abs(z_all_3(1:start_t*fs-1))),diff(abs(z_all_3(1:start_t*fs))));
axis square
grid on
xlabel('Magintude')
ylabel('Maginitude Derivative')
title('Phase Portrait')
hold on
for i=0:npulses-1
    % ugly but works
    pause
    subplot(4,6,2)
    plot(abs(z_all_1((start_t*fs+i*pul_step*fs)+i+2:(start_t*fs+(i+1)*pul_step*fs)-1)),diff(abs(z_all_1((start_t*fs+i*pul_step*fs)+i+2:(start_t*fs+(i+1)*pul_step*fs)))),'r');
    subplot(4,6,4)
    plot(abs(z_all_2((start_t*fs+i*pul_step*fs)+i+2:(start_t*fs+(i+1)*pul_step*fs)-1)),diff(abs(z_all_2((start_t*fs+i*pul_step*fs)+i+2:(start_t*fs+(i+1)*pul_step*fs)))),'r');
    subplot(4,6,6)
    plot(abs(z_all_3((start_t*fs+i*pul_step*fs)+i+2:(start_t*fs+(i+1)*pul_step*fs)-1)),diff(abs(z_all_3((start_t*fs+i*pul_step*fs)+i+2:(start_t*fs+(i+1)*pul_step*fs)))),'r');
end
subplot(4,6,2)
plot(abs(z_all_1((start_t*fs+(i+1)*pul_step*fs)+i+1+2:end-1)),diff(abs(z_all_1((start_t*fs+(i+1)*pul_step*fs)+i+1+2:end))),'r');
subplot(4,6,4)
plot(abs(z_all_2((start_t*fs+(i+1)*pul_step*fs)+i+1+2:end-1)),diff(abs(z_all_2((start_t*fs+(i+1)*pul_step*fs)+i+1+2:end))),'r');
subplot(4,6,6)
plot(abs(z_all_3((start_t*fs+(i+1)*pul_step*fs)+i+1+2:end-1)),diff(abs(z_all_3((start_t*fs+(i+1)*pul_step*fs)+i+1+2:end))),'r');
subplot(4,6,7:8)
plot(time_all,abs(z_all_1))
grid on
hold on
xlabel('Time (s)')
ylabel('Amplitude')
title('Magnitude')
subplot(4,6,9:10)
plot(time_all,abs(z_all_2))
grid on
hold on
xlabel('Time (s)')
ylabel('Amplitude')
title('Magnitude')
subplot(4,6,11:12)
plot(time_all,abs(z_all_3))
grid on
hold on
xlabel('Time (s)')
ylabel('Amplitude')
title('Magnitude')

if animate == 1
    
    for i=1:length(time_all)
        
        if mod(i,fs/20) == 0
            subplot(4,6,[13,14,19,20])
            plot3(time_all(1:i),real(z_all_1(1:i)),imag(z_all_1(1:i)),'-')
            grid on
            xlabel('Time (s)')
            ylabel('Real part')
            zlabel('Imaginary part')
            subplot(4,6,[15,16,21,22])
            plot3(time_all(1:i),real(z_all_2(1:i)),imag(z_all_2(1:i)),'-')
            grid on
            xlabel('Time (s)')
            ylabel('Real part')
            zlabel('Imaginary part')
            subplot(4,6,[17,18,23,24])
            plot3(time_all(1:i),real(z_all_3(1:i)),imag(z_all_3(1:i)),'-')
            grid on
            xlabel('Time (s)')
            ylabel('Real part')
            zlabel('Imaginary part')
            pause(0.05)
            
        end
    end
    
end