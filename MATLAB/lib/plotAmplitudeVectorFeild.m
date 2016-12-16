function [plot_obj, FPs] = plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon, F, f_osc, f_input, figure_number, plot_color)
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
%   plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon, F, f_osc, f_input, figure_number, plot_color)
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
%      plot_color:  [Optional] specifiy color for plotting
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Returns:  [plot_obj, FPs]
%   
%        plot_obj:  Figure (figure_number) with properties
%             FPs:  (float) Array of fixed points (r*)
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

if nargin < 4 || nargin > 9  
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
    f_input = 0.0;
    F = 0.0;
    figure_number = 1;
    plot_color = 0;
end

if nargin == 5 || nargin == 6
    display('Please enter the right number of parameters')
    error('Type "help plotAmplitudeVectorFeild" for more information')
end

if nargin == 7
    figure_number = 1;
    plot_color = [];
end

if nargin == 8
    plot_color = [];
end

% make sure the input is properly zeroed out for autonomous osc
if F && ~f_input
    display('Both F and f_input must be zero for autonomous oscillators')
    display('F and f_input have been set to zero')
    f_input = 0;
    F = 0;
elseif  f_input && ~F
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

numFPs = length(r_star);
FPs = r_star;

%% grab all required initial conditions

% Critical Hopf
if regime == 1 
    theta = pi;
    if ~F
        assert(numFPs == 1,'Critical Hopf has more than one FP')
    end
    z0 = [exp(1i*theta)]; % initial conditions
% Supercritical Hopf
elseif regime == 2 % Supercritical Hopf
    if ~F
        assert(numFPs == 2,'Supercritical Hopf does not have 2 FP')
    end
    z0 = [eps,r_star(2)+0.3]; % initial conditions
% Supercritical DLC
elseif regime == 3
    if ~F
        assert(numFPs == 3,'Supercritical DLC does not have 3 FPs')
    end
    z0 = [r_star(2)-0.001,r_star(2)+0.001, 0.99999999]; % initial conditions
% Subcritical DLC
elseif regime == 4 
    if ~F
        assert(numFPs == 1,'Subcritical DLC has more than one FP')
    end
    z0 = [0.99999999]; % initial conditions
% Unknown
elseif regime == 0 
    display('The oscillators parameter regime is unknown, plot was not generated')
    return
end

%% Plot 

% parameters for time
fs = 1000;
dur = 30; % in seconds
T = 1/fs;
time = 0:T:dur;

ic_num = 1;
for ic = z0 % loop through initial conditions

    dzdt = @(t,z,alpha,beta1,beta2,epsilon,f_osc,f_input)  ...
        z*(alpha + 1i*2*pi*f_osc + beta1*abs(z)^2 + ...
        (epsilon*beta2*abs(z)^4)/(1-epsilon*abs(z)^2) + F*exp(1i*2*pi*f_input*t)); 

    % integrate the ode to obtain the amplitude
    [~,z_out] = ode45(@(t,z) dzdt(t,z,alpha,beta1,beta2,epsilon,f_osc,f_input),time,[ic]);
    
    if ic_num == 1
        % store color for first condition
        
        r1 = (abs(z_out(1:end-1)));
        r_dot1 = diff(abs(z_out));
        
        if ~isempty(plot_color)
            % plot
            plot_obj = figure(figure_number);
            plot(r1,r_dot1,'linewidth',2,'Color', plot_color);
            hold on
        else    
            % plot
            plot_obj = figure(figure_number);
            plot(r1,r_dot1,'linewidth',2);
            hold on
        end
    else
        % In order to keep the colors consistent with the legend we 
        % will remove the rest of the initial conditions from the legend
        
        r = (abs(z_out(1:end-1)));
        r_dot = diff(abs(z_out));
    
         if ~isempty(plot_color)
            % plot
            plot_obj = figure(figure_number);
            h_rest = plot(r,r_dot,'linewidth',2,'Color', plot_color);
            hold on
        else    
            % plot
            plot_obj = figure(figure_number);
            h_rest = plot(r,r_dot,'linewidth',2);
            hold on
         end
         % remove the line from the legend 
         hasbehavior(h_rest,'legend',false);
         
    end
    
    ic_num = ic_num + 1;
    
%     figure(2)
%     polar((angle(z_out)),abs(abs(z_out)))
%     hold on
end
