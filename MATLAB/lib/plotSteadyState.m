function [stim_strengths, stead_state_amp] = plotSteadyState(f,alpha,beta1,beta2,epsilon,step_size,fs)
% plotSteadyState: Calculate steady state amplitude 
%                  vs stimulus strength
% Arguments: 
% 
%               f:  Oscillator Frequency
%           alpha: (Hopf Oscillator Parameter)
%           beta1: (Hopf Oscillator Parameter)
%           beta2: (Hopf Oscillator Parameter)
%         epsilon: (Hopf Oscillator Parameter)
%       step_size: (optional, default 10/T, fs = 1000 Hz) 
%                  Step function amplitude incriments
%              fs: (optional, default 1000 Hz) Sample Rate
% Returns:
% 
%  stim_strengths: A vector of stimulus strengths
% stead_state_amp: A vector of steady state gains

max_stim_amp = 1; % initial step amplitude
step = @(t) heaviside(t); % step function
z0 = 1e-5; % initial condition
% threshold = 1e-3; % threshold for the norm of the derivative
w = 2*pi*f; % radians

% plotting flag
plotting = 1;

switch nargin
    case 5
        fs = 1000;
        T = 1/fs;
        steps = 10*T;
        % Displays 
        disp(['Running plotSteadyState with ',num2str(nargin),' parameters:'])
        fprintf('\n')
        disp(['The Amplitude Step is: ',num2str(steps)])
        disp(['The Sample Rate is: ',num2str(fs)])
        disp(['The Oscillator Frequency is: ',num2str(f), ' Hz'])
        disp(['alpha is: ',num2str(alpha)])
        disp(['beta1 is: ',num2str(beta1)])
        disp(['beta2 is: ',num2str(beta2)])
        disp(['epsilon is: ',num2str(epsilon)])
        fprintf('\n')
        
    case 6
        fs = 1000;
        T = 1/fs;
        steps = step_size;
        % Displays 
        disp(['Running plotSteadyState with ',num2str(nargin),' parameters:'])
        fprintf('\n')
        disp(['The Amplitude Step is: ',num2str(steps)])
        disp(['The Sample Rate is: ',num2str(fs)])
        disp(['The Oscillator Frequency is: ',num2str(f), ' Hz'])
        disp(['alpha is: ',num2str(alpha)])
        disp(['beta1 is: ',num2str(beta1)])
        disp(['beta2 is: ',num2str(beta2)])
        disp(['epsilon is: ',num2str(epsilon)])
        fprintf('\n')
    case 7
        T = 1/fs;
        steps = step_size;
        % Displays 
        disp(['Running plotSteadyState with ',num2str(nargin),' parameters:'])
        fprintf('\n')
        disp(['The Amplitude Step is: ',num2str(steps)])
        disp(['The Sample Rate is: ',num2str(fs)])
        disp(['The Oscillator Frequency is: ',num2str(f), ' Hz'])
        disp(['alpha is: ',num2str(alpha)])
        disp(['beta1 is: ',num2str(beta1)])
        disp(['beta2 is: ',num2str(beta2)])
        disp(['epsilon is: ',num2str(epsilon)])
        fprintf('\n')
    otherwise
        disp('You fucked up! Give me the right number of arguments')
end

% Hopf equation
dzdt = @(t,z,w1,alpha1,beta11,beta12,epsilon1,F,step) ...                  % params 
            z(1)*(alpha1 + 1i*w1 + ...                                     % linear
            beta11*abs(z(1))^2 + ...                                       % quadratc
            (epsilon1*beta12*abs(z(1))^4)/(1-epsilon1*abs(z(1))^2)) + ...  % quintic
            F*step(t); % step       

% initialize convergence parameters
prev_z = 1e5;
curr_z = 0;
opts = odeset('Events',@isConvergent);


r = [];
stim_strengths = 0.1:steps:max_stim_amp;
for F = stim_strengths
    
%     disp(['The Current Stimulus Strength is: ',num2str(F)])

    [t,z] = ode45(@(t,z) ...
            dzdt(t,z,w,alpha,beta1,beta2,epsilon, ... % oscillator params
                     F,step),[0 Inf],z0,opts);
                 
    disp(['r is: ',num2str(abs(z(end)))])
    disp(['Step Amplitude is: ',num2str(F)])
    fprintf('\n')
    
    r = [r abs(z(end))];

end

if plotting
    figure
    % plot the last response
    plot(t,abs(z),'-k','LineWidth',2)
    axis tight
end

stead_state_amp = r;

function [value,isterminal,direction] = isConvergent(t,z)
    threshold = 1e-3;
    curr_z = z;
%     abs(curr_z-prev_z)
    value = abs(curr_z-prev_z)-threshold; % Don't try to detect exact zero for asymptotics
    isterminal = 1;
    direction = -1;
    prev_z = curr_z;
end

end

