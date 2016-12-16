function [plot_obj] = plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon, F, f_osc, f_input, figure_number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% plotAmplitudeVectorFeild: 
%         Plots the amplitude vector feild for a given
%         oscillator
%           
%   Note: * This function returns the plot as an object, allowing for 
%           manipulation of the figure outside of the function
%         * This function does not close the figure allowing for plotting
%           and parameter analysis over for loops
%         * This function automatically simulates the oscillator amplitude 
%           behavior using all necessary initial conditions
% 
% Author: Wisam Reid
%  Email: wisam@ccrma.stanford.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Function Definition:
% 
%   plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon, F, f_osc, f_input, figure_number)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Arguments: 
% 
%           alpha:  (float) Hopf Oscillator Parameter
%           beta1:  (float) Hopf Oscillator Parameter
%           beta2:  (float) Hopf Oscillator Parameter
%         epsilon:  (float) Hopf Oscillator Parameter
%               F:  [Optional] (float) 0 means autonomous 
%           f_osc:  [Optional] (float) oscillator frequency
%         f_input:  [Optional] (float) input frequency
%   figure_number:  [Optional] (int) number to be used for the plot
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Returns:  [plot_obj]
%   
%        plot_obj:  Figure (figure_number) with properties
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Example Calls:  
% 
% An autonomous critical Hopf oscillator (See Kim & Large 2015 Fig 1):
% 
%       fig = plotAmplitudeVectorFeild(0, -100, 0, 0)
% 
% A forced supercritical Hopf oscillator (See Kim & Large 2015 Fig 4):
%       
%       fig = plotAmplitudeVectorFeild(1, -100, 0, 0, 0.02, 1, 0.5)
% 
% A forced supercritical DLC oscillator using plot number 3:
%   (See Kim & Large 2015 Fig 7)
%       
%       fig3 = plotAmplitudeVectorFeild(-1, 4, -1, 1, 0.3, 1, 0.8, 3)
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Handle Arguments

if nargin < 4 || nargin > 8  
    disp('Error running plotAmplitudeVectorFeild')
    disp('Please provide the right number of arguments')
    disp('Function Call: ')
    disp('plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon)')
    disp('or')
    disp('plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon, F, f_osc, f_input)')
    fprintf('\n')
    disp('type "help plotAmplitudeVectorFeild" for more details')
    fprintf('\n')
    return
end

if nargin == 4
    f_osc = 1;
    f_input = 0;
    F = 0;
    figure_number = 1;
end

if nargin == 7
    figure_number = 1;
end

% make sure the input is properly zeroed out for autonomous osc
if F == 0.0 && f_input ~= 0.0
    display('Both F and f_input must be zero for autonomous oscillators')
    display('F and f_input have been set to zero')
    f_input = 0;
    F = 0;
elseif  F ~= 0.0 && f_input == 0.0
    display('Both F and f_input must be zero for autonomous oscillators')
    display('F and f_input have been set to zero')
    f_input = 0;
    F = 0;
end 

% Check input arguments
if F < 0
  error('F > 0')
end

if beta2 == 0 && epsilon ~= 0
    disp('set epsilon = 0 because beta2 = 0')
    epsilon = 0; 
end

if epsilon == 0 && beta2 ~= 0
    disp('set beta2 = 0 because epsilon = 0')
    beta2 = 0;
end

%% Get Param Regime and FPs

% Regime:
%   critical Hopf = 1
%   supercritical Hopf = 2
%   supercritical DLC = 3
%   subcritical DLC = 4
%   unknown = 0

% Stability Type:
%   stable node = 1 
%   stable spiral = 2
%   unstable node = 3
%   unstable spiral = 4
%   a saddle point = 5 
%   unknown = 0

[r_star, ~, ~, regime] = getFP(f_osc, f_input, alpha, beta1, beta2, epsilon, F, 0);

%% grab all required initial conditions

switch regime
    % Critical Hopf
    case regime == 1 
        assert(length(r_star) == 1,'Critical Hopf has more than one FP')
        z0 = [1.0]; % initial conditions
    % Supercritical Hopf
    case regime == 2 % Supercritical Hopf
        assert(length(r_star) == 1,'Critical Hopf has more than one FP')
        z0 = [eps,1.0]; % initial conditions
    % Supercritical DLC
    case regime == 3 
        assert(length(r_star) == 3,'Supercritical DLC does not have 3 FPs')
        z0 = [r_star(2)-eps,r_star(2)+eps, 1.0]; % initial conditions
    % Subcritical DLC
    case regime == 4 
        assert(length(r_star) == 1,'Subcritical DLC has more than one FP')
         z0 = [1.0]; % initial conditions
    % Unknown
    case regime == 0 
        display('The oscillators parameter regime is unknown plot was not generated')
        return
end

%% Plot 

% parameters for time
fs = 1000;
dur = 10; % in seconds
T = 1/fs;
time = 0:T:dur;

for initial_conditon  = z0
    

    dzdt = @(t,z,alpha,beta1,beta2,epsilon,f_osc)  ...
        z*(alpha + 1i*2*pi*f_osc + beta1*abs(z)^2 + ...
        (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2)); 

    % integrate the ode to obtain the amplitude
    [t_out,z_out] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f_osc),time,z0);

    plot_obj = figure(figure_number);
    plot((abs(z_out(1:end-1))),diff(abs(z_out)));

    % grid on
end

%% Define r_dot

% % r values
% r = 0:0.01:1;

% r_dot = @(t) alpha*r + beta1*r.^3 + ((epsilon*beta2*r.^5)/(1 - epsilon*r.^2));

% % parameters for time
% fs = 1000;
% end_t = 10; % in seconds
% T = 1/fs;
% time = 0:T:end_t;
% 
% % Initial Condition
% r0 = 0.001;
% 
% % solve r_dot
% r_dot = @(t,r,alpha,beta1,beta2,epsilon)  ...
%     [alpha*r(1) + beta1*r(1)^3 + ((epsilon*beta2*r(1)^5)/(1-epsilon*r(1)^2))];
% 
% % integrate the ode to obtain the amplitude
% [time,r] = ode45(@(t,r) r_dot(t,r,alpha,beta1,beta2,epsilon),time,[r0]);

%% Plot

% regime = getOscParamRegime(alpha, beta1, beta2, epsilon); % Check for a known regime


% figure(1);
% % subplot(1,nfreqs,i)
% plot(time,r)
% title('Yo bro')
% ylim([0 1.3])
% grid on
% xlabel('Time (s)')
% ylabel('Magnitude') 

% figure(2)    
% polar(pi,1.3) 
% hold on
% polar(z(:,2),z(:,1))
% hold on
% grid on
% % pause

% figure(1)
% if regime == 1
%     plot(r,r_dot)
%     hold on 
%     plot([0 max(r)],[0 0],'k-')
%     ylim([min(r_dot),0.2])
% elseif regime == 2
%     plot(r,r_dot)
%     hold on 
%     plot([0 max(r)],[0 0],'k-')
%     ylim([min(r_dot),max(0.2,max(r_dot)+ 0.1)])
% elseif regime == 3
%     plot(r,r_dot)
%     hold on 
%     plot([0 max(r)],[0 0],'k-')
%     ylim([min(-0.1,min(r_dot)),max(0.2,max(r_dot)+ 0.1)])
% elseif regime == 4
%     plot(r,r_dot)
%     hold on 
%     plot([0 max(r)],[0 0],'k-')
%     ylim([min(-0.1,min(r_dot)),max(0.2,max(r_dot)+ 0.1)])
% else
%     display('Warning: The oscillator parameter regime could not be classified')
% end

end

